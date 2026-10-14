# Minimal TRES example (runs in seconds; no patient data needed).
# Part 1 fits TRES on synthetic data with a known risk hand-off.
# Part 2 reads the frozen PDAC model and prints its per-interval log hazard ratios and relative risk for an example profile.
# Run from this folder:  Rscript minimal_example.R
suppressPackageStartupMessages(library(survival))
source(file.path("..", "R", "TRES.R"))

# ---------------------------------------------------------------- Part 1: synthetic data
set.seed(2027); n <- 800
d <- data.frame(cohort = sample(c("A", "B", "C"), n, TRUE), T = rnorm(n), N = rnorm(n), PCOSP = rnorm(n), P1 = rnorm(n), P2 = rnorm(n))
h <- tres_transition(12, 4); base <- c(0.05, 0.04, 0.03, 0.022); brk <- c(0, 6, 12, 24, 36)
death <- vapply(seq_len(n), function(i) { for (j in 1:4) {
  e <- 0.2 * d$T[i] + (0.6 - 0.5 * h[j]) * d$PCOSP[i] + (0.1 + 0.3 * h[j]) * d$N[i]       # tumour score fades, N rises
  w <- rexp(1, base[j] * exp(e)); if (w <= brk[j + 1] - brk[j]) return(brk[j] + w) }; Inf }, 0)
cen <- 8 + rexp(n, 1 / 24); d$time <- pmax(pmin(death, cen, 36), 0.01); d$event <- as.integer(death <= pmin(cen, 36))
fit <- tres_fit(d, static = "T", core = c("N", "PCOSP"), sources = list(SRC = c("P1", "P2")))
print(fit)                                   # gate decisions (P1/P2 carry no signal and should not be kept)
print(tres_gate_report(fit)$log_hr)          # per-interval log HR
print(round(tres_predict(fit, d[1:3, ], months = c(12, 24, 36)), 3))   # survival curves for three patients

# ---------------------------------------------------------------- Part 2: frozen PDAC model (coefficients only)
F <- readRDS(file.path("..", "frozen", "TRES_frozen_PDAC_model.rds"))$TN$model
cat(sprintf("\nFrozen PDAC TRES | cell programs kept: %s | temporal gate p = %.2g | kappa/tau set %d | theta %g\n",
            if (length(F$keep)) paste(F$keep, collapse = "+") else "none", F$gate_p, F$transition, F$theta))
L <- t(sapply(F$vars, function(v) F$b[[v]] + F$h * F$d[[v]])); colnames(L) <- c("0-6m", "6-12m", "12-24m", "24-36m")
print(round(L, 3))
# Example profile in model units (T stage 3; N1; PCOSP +1 SD within its cohort). Relative risk vs T0/N0/PCOSP 0, by interval.
# Absolute survival needs the cohort reference batch and the development baseline hazard, which are not in this compact object
# (see docs/TRES_model_card.md, "Deployment limitations").
x <- c(N = 1, PCOSP = 1); eta <- F$bS[["T"]] * 3 + colSums(L * x[rownames(L)])
print(round(exp(eta), 2))
