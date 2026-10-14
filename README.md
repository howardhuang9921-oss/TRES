# TASA — evidence-gated time-varying survival modelling

**TASA** (Time-Adaptive Survival Analysis) is a clinically anchored Cox model in which two null-calibrated evidence gates decide

1. **which single-cell-defined source programs enter the model** (source gate: likelihood-ratio test per source, Bonferroni);
2. **whether covariate effects change over follow-up time** (temporal gate: likelihood-ratio test across candidate transitions, Bonferroni).

If the evidence is insufficient, the model reverts to the static Cox model (all time-change coefficients are zero). Time variation is represented by one shared smooth transition `h(t; kappa, tau)`, with ridge shrinkage tuned by leave-one-cohort-out integrated Brier score.

## Why
Time-resolved survival models (for example piecewise per-bin cell coefficients) report early-to-late "risk reversals" even when none exist. In our simulations:
- under a static truth, interval-wise Cox models reported false early-to-late reversals in **80%** of simulations versus **1.7% for TASA**; across 1,800 simulations of varying size, TASA's false-positive rate was **2.0–3.3%** (other methods 4.0–86%);
- power was **≥ 90%** for a moderate hand-off at n ≥ 600.

## Usage
```r
source("R/TASA.R")                       # or library(TASA) after installation
P <- tasa_score_programs(expr, programs, universe = common_genes, centre_by = cohort)   # single-sample scores
fit <- tasa_fit(train, static = "T", core = c("N", "PCOSP"),
                sources = list(CAF = c("CAF1","CAF2"), MYELOID = c("MY1","MY2"), LYMPHOID = c("LY1","LY2")))
print(fit)                               # gate decisions
tasa_gate_report(fit)$log_hr            # log-HR trajectory per interval
tasa_predict(fit, new_patient, months = c(12, 24, 36))   # one new patient, no retraining
```
Input data need the columns `time` (months), `event` (0/1), `cohort`, and standardised covariates.

## Properties
- At most `2·|V| + |static|` parameters (≤ 18 in the PDAC application). Fits in milliseconds to seconds on a laptop CPU.
- **Single-sample scoring, batch-referenced standardisation.** Raw program scores are within-sample ranks and need no retraining. The validated PDAC pipeline, however, standardised scores within cohort, so a new site currently needs a local reference batch; stand-alone one-patient absolute-risk use is not yet validated.
- Tests: `Rscript tests/run_tests.R` checks:
  - the shape of the transition;
  - reversion to the static model under no effect;
  - detection of a hand-off;
  - valid survival curves for a single new patient;
  - invariance of single-sample scores.

## Validation (pre-registered, SHA-256-frozen analyses)
- Simulation: SPEC_60, SPEC_62, SPEC_65.
- PDAC development (4 cohorts, 794 patients): SPEC_61/62.
- 12 TCGA cancers + METABRIC: SPEC_64.
- PDAC external validation (10 cohorts, 616 patients): SPEC_63.
- Head-to-head against full scTREND on its own melanoma data: SPEC_66.
