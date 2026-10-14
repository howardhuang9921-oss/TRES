# TRES: Time-Resolved Effect Selection (evidence-gated time-varying survival modelling).
# Same algorithm as the frozen analysis engine used in the study; this file is a cleaned, documented, user-facing
# re-implementation for single datasets and new patients.

#' Single-sample program scores (platform-agnostic)
#' @param expr samples x genes matrix (log-scale expression; any platform)
#' @param programs named list of gene vectors (e.g. two programs per cell source)
#' @param universe optional common gene universe; ranks are computed within it
#' @param centre_by optional batch / cohort labels; scores are z-scored within each batch
#' @return samples x programs matrix
tres_score_programs <- function(expr, programs, universe = NULL, centre_by = NULL) {
  if (!is.null(universe)) expr <- expr[, intersect(colnames(expr), universe), drop = FALSE]
  one <- function(x, gs) { gs <- intersect(gs, names(x)); if (length(gs) < 4) return(NA_real_); mean(rank(x, ties.method = "average")[gs] / length(x)) }
  S <- sapply(programs, function(gs) apply(expr, 1, one, gs = gs))
  S <- matrix(S, nrow = nrow(expr), dimnames = list(rownames(expr), names(programs)))
  if (!is.null(centre_by)) for (b in unique(centre_by)) { i <- centre_by == b
    S[i, ] <- apply(S[i, , drop = FALSE], 2, function(v) (v - mean(v)) / stats::sd(v)) }
  S
}

#' Shared smooth transition h at interval mid-points (0 -> 1 between 0 and horizon)
tres_transition <- function(kappa, tau, mid = c(3, 9, 18, 30), horizon = 36) {
  raw <- stats::plogis((c(0, mid, horizon) - kappa) / tau); J <- length(mid)
  (raw[2:(J + 1)] - raw[1]) / (raw[J + 2] - raw[1])
}

.split <- function(d, vars, static, h, brk) {
  s <- survival::survSplit(Surv(time, event) ~ ., data = d[, c("time", "event", "cohort", static, vars)], cut = brk[2:(length(brk) - 1)], episode = "ep")
  for (v in vars) s[[paste0(v, "_D")]] <- s[[v]] * h[s$ep]; s
}
.fit_dyn <- function(d, vars, static, h, theta, brk) {
  s <- .split(d, vars, static, h, brk); Dv <- paste0(vars, "_D")
  rhs <- if (theta > 0) paste0("ridge(", paste(Dv, collapse = ","), ", theta = ", theta, ", scale = FALSE)") else paste(Dv, collapse = "+")
  f <- survival::coxph(stats::as.formula(paste("Surv(tstart, time, event) ~", paste(c(static, vars), collapse = "+"), "+", rhs, "+ strata(cohort)")), data = s, ties = "breslow")
  cf <- stats::coef(f); p0 <- length(static)
  list(bS = cf[static], b = cf[vars], d = stats::setNames(cf[(p0 + length(vars) + 1):(p0 + 2 * length(vars))], vars), h = h, vars = vars, static = static, ll = f$loglik[2])
}
.static <- function(d, vars, static) survival::coxph(stats::as.formula(paste("Surv(time, event) ~", paste(c(static, vars), collapse = "+"), "+ strata(cohort)")), data = d, ties = "breslow")
.eta <- function(m, d) { x <- as.matrix(d[, m$vars, drop = FALSE]); s0 <- if (length(m$static)) drop(as.matrix(d[, m$static, drop = FALSE]) %*% m$bS) else 0
  matrix(sapply(seq_along(m$h), function(j) s0 + drop(x %*% (m$b + m$h[j] * m$d))), nrow = nrow(d)) }
.base_rates <- function(eta, d, brk) { J <- length(brk) - 1
  E <- sapply(seq_len(J), function(j) pmax(0, pmin(d$time, brk[j + 1]) - brk[j])); Dm <- sapply(seq_len(J), function(j) as.integer(d$event == 1 & d$time > brk[j] & d$time <= brk[j + 1]))
  (colSums(Dm) + 1e-6) / (colSums(E * exp(pmin(pmax(eta, -25), 25))) + 1e-6) }
.surv <- function(eta, base, brk, months) matrix(sapply(months, function(t) { w <- pmax(0, pmin(t, brk[-1]) - brk[-length(brk)]); exp(-rowSums(sweep(exp(pmin(pmax(eta, -25), 25)), 2, base * w, "*"))) }), nrow = nrow(eta))
.ibs <- function(S, d, months) { ck <- survival::survfit(survival::Surv(d$time, 1 - d$event) ~ 1)
  G <- function(x) vapply(x, function(u) { i <- which(ck$time < u); if (!length(i)) 1 else max(ck$surv[max(i)], 1e-4) }, 0)
  br <- sapply(seq_along(months), function(k) { h <- months[k]; alive <- d$time > h; w <- ifelse(d$time <= h & d$event == 1, 1 / G(pmin(d$time, h)), ifelse(alive, 1 / G(h), 0)); mean(w * (alive - S[, k])^2) })
  sum(diff(months) * (utils::head(br, -1) + utils::tail(br, -1)) / 2) / (max(months) - min(months)) }

#' Fit TRES
#' @param data data.frame with columns time (months), event (0/1), cohort, and all variables below (standardised)
#' @param static static clinical covariates (e.g. "T")
#' @param core dynamic-eligible clinical + tumour-intrinsic covariates (e.g. c("N", "PCOSP"))
#' @param sources named list of program pairs, e.g. list(CAF = c("CAF1","CAF2"), ...)
#' @param alpha family-wise level for each gate
#' @param transitions data.frame(kappa, tau) of candidate transitions
#' @param thetas ridge grid for the dynamic part
#' @param brk interval breaks in months (horizon = last break)
#' @return object of class "tres"
tres_fit <- function(data, static = character(0), core, sources = list(), alpha = 0.05,
                     transitions = data.frame(kappa = c(12, 15, 18), tau = c(4, 6, 8)), thetas = c(0, 1, 4, 16), brk = c(0, 6, 12, 24, 36)) {
  stopifnot(all(c("time", "event", "cohort") %in% names(data)))
  data$time <- pmin(data$time, max(brk)); data$event <- as.integer(data$event == 1 & data$time < max(brk) + 1e-9)
  mid <- (utils::head(brk, -1) + utils::tail(brk, -1)) / 2; months <- seq(6, max(brk))
  base <- .static(data, core, static)
  src_p <- vapply(names(sources), function(s) { f <- .static(data, c(core, sources[[s]]), static)
    min(1, length(sources) * stats::pchisq(2 * (f$loglik[2] - base$loglik[2]), df = length(sources[[s]]), lower.tail = FALSE)) }, 0)
  keep <- names(sources)[src_p < alpha]; vars <- c(core, unlist(sources[keep], use.names = FALSE)); st <- .static(data, vars, static)
  hs <- lapply(seq_len(nrow(transitions)), function(g) tres_transition(transitions$kappa[g], transitions$tau[g], mid, max(brk)))
  lrt <- vapply(hs, function(h) 2 * (.fit_dyn(data, vars, static, h, 0, brk)$ll - st$loglik[2]), 0)
  gate_p <- min(1, length(hs) * stats::pchisq(max(lrt), df = length(vars), lower.tail = FALSE)); best <- which.max(lrt)
  if (gate_p >= alpha) { cf <- stats::coef(st)
    m <- list(bS = cf[static], b = cf[vars], d = stats::setNames(rep(0, length(vars)), vars), h = hs[[best]], vars = vars, static = static); theta <- NA
  } else { h <- hs[[best]]; cs <- unique(data$cohort)
    folds <- if (length(cs) >= 3) data$cohort else { set.seed(1); sample(rep(1:3, length.out = nrow(data))) }
    ibs <- vapply(thetas, function(th) { w <- sapply(unique(folds), function(q) { a <- data[folds != q, ]; v <- data[folds == q, ]; mm <- .fit_dyn(a, vars, static, h, th, brk)
      c(.ibs(.surv(.eta(mm, v), .base_rates(.eta(mm, a), a, brk), brk, months), v, months), nrow(v)) }); sum(w[1, ] * w[2, ]) / sum(w[2, ]) }, 0)
    theta <- thetas[max(which(ibs <= min(ibs) + 1e-4))]; m <- .fit_dyn(data, vars, static, h, theta, brk) }
  m$base <- .base_rates(.eta(m, data), data, brk); m$brk <- brk
  structure(c(m, list(gate = list(source_p = src_p, sources_kept = keep, temporal_p = gate_p, activated = gate_p < alpha,
                                  transition = transitions[best, ], theta = theta, alpha = alpha))), class = "tres")
}

#' Predict survival probabilities for new patients (single-sample: no retraining)
tres_predict <- function(object, newdata, months = c(12, 24, 36)) {
  S <- .surv(.eta(object, newdata), object$base, object$brk, months); colnames(S) <- paste0("S", months); S
}

#' Log-hazard-ratio trajectory per interval and gate decisions
tres_gate_report <- function(object) {
  L <- t(sapply(object$vars, function(v) object$b[[v]] + object$h * object$d[[v]]))
  colnames(L) <- paste0(utils::head(object$brk, -1), "-", utils::tail(object$brk, -1), "m")
  list(sources_p = object$gate$source_p, sources_kept = object$gate$sources_kept, temporal_p = object$gate$temporal_p,
       time_varying = object$gate$activated, transition = object$gate$transition, theta = object$gate$theta, log_hr = round(L, 3))
}

print.tres <- function(x, ...) { g <- x$gate
  cat(sprintf("TRES | sources kept: %s | temporal gate p = %.3g -> %s\n", if (length(g$sources_kept)) paste(g$sources_kept, collapse = "+") else "none",
              g$temporal_p, if (g$activated) sprintf("time-varying (kappa %g, tau %g, theta %g)", g$transition$kappa, g$transition$tau, g$theta) else "static (exact fallback)"))
  invisible(x) }
