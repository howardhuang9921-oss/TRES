# Base-R unit tests for TRES (run by R CMD check, or Rscript tests/run_tests.R)
suppressPackageStartupMessages({ library(survival); if (requireNamespace("TRES", quietly = TRUE)) library(TRES) else source(file.path("..", "R", "TRES.R")) })
sim <- function(n, dyn, seed) { set.seed(seed)
  d <- data.frame(cohort = sample(c("A", "B", "C"), n, TRUE), T = rnorm(n), N = rnorm(n), PCOSP = rnorm(n), P1 = rnorm(n), P2 = rnorm(n))
  h <- tres_transition(15, 4); base <- c(0.05, 0.04, 0.03, 0.022); brk <- c(0, 6, 12, 24, 36)
  death <- vapply(seq_len(n), function(i) { for (j in 1:4) { e <- 0.2 * d$T[i] + (0.5 + dyn * (0.1 - 0.5) * h[j]) * d$PCOSP[i] + (0.3 + dyn * 0.3 * h[j]) * d$N[i]
    w <- rexp(1, base[j] * exp(e)); if (w <= brk[j + 1] - brk[j]) return(brk[j] + w) }; Inf }, 0)
  cen <- 8 + rexp(n, 1 / 24); d$time <- pmax(pmin(death, cen, 36), 0.01); d$event <- as.integer(death <= pmin(cen, 36)); d }
# 1. transition is monotone from ~0 to 1
h <- tres_transition(12, 4); stopifnot(all(diff(h) > 0), h[1] >= 0, h[length(h)] <= 1)
# 2. static truth: gate closes, exact fallback (all dynamic coefficients zero)
m0 <- tres_fit(sim(800, 0, 1), static = "T", core = c("N", "PCOSP"), sources = list(SRC = c("P1", "P2")))
stopifnot(!m0$gate$activated, all(m0$d == 0))
# 3. strong hand-off: gate opens and PCOSP effect declines while N effect rises
m1 <- tres_fit(sim(1500, 1, 2), static = "T", core = c("N", "PCOSP"), sources = list(SRC = c("P1", "P2")))
r <- tres_gate_report(m1); stopifnot(m1$gate$activated, r$log_hr["PCOSP", 1] > r$log_hr["PCOSP", 4], r$log_hr["N", 4] > r$log_hr["N", 1])
# 4. predictions: valid, monotone survival for a single new patient (no retraining needed)
nd <- data.frame(T = 0.5, N = 1, PCOSP = 1, P1 = 0, P2 = 0); S <- tres_predict(m1, nd, months = c(6, 12, 24, 36))
stopifnot(all(S > 0 & S <= 1), all(diff(S[1, ]) <= 1e-12))
# 5. single-sample program scores are invariant to adding other samples (deployability)
X <- matrix(rnorm(20 * 200), 20, dimnames = list(NULL, paste0("g", 1:200))); pr <- list(A = paste0("g", 1:10), B = paste0("g", 11:20))
stopifnot(isTRUE(all.equal(tres_score_programs(X[1:5, ], pr)[1, ], tres_score_programs(X, pr)[1, ])))
cat("all TRES tests passed\n")
