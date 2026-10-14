# TRES — Time-Resolved Effect Selection

TRES is a prognostic model for postoperative overall survival (0–36 months) in resected pancreatic ductal adenocarcinoma (PDAC).
It decides its own structure with two statistical tests:

1. **Layer 1 (source test):** do single-cell-defined cell programs (CAF, myeloid, lymphoid) add information beyond
   pathological T, N and the PCOSP tumour RNA score? Likelihood-ratio test per cell source, Bonferroni-corrected.
2. **Layer 2 (time test):** do the retained effects really change over follow-up time? If yes, all time-varying effects share one
   smooth transition curve h(t; κ, τ) with ridge shrinkage; if not, the model reverts to a static, cohort-stratified Cox model.

Information that adds nothing is not added; effects that do not change are not forced to change.

Author: Hung-Ying Huang. Supervisor and corresponding author: Tzong-Yi Lee (National Yang Ming Chiao Tung University).
This repository accompanies the 2027 Taiwan International Science Fair research report.

## Repository contents
| Path | Content |
|---|---|
| `R/TRES.R` | The algorithm: `tres_score_programs`, `tres_fit`, `tres_predict`, `tres_gate_report`, `tres_transition` |
| `tests/run_tests.R` | 5 unit tests (transition shape, reversion to the static model under no effect, detection of a hand-off, valid survival curves for one new patient, invariance of single-sample scores) |
| `examples/minimal_example.R` | Runs in seconds on synthetic data, then prints the frozen PDAC model; expected output in `minimal_example_output.txt` |
| `frozen/TRES_frozen_PDAC_model.rds` | Frozen final PDAC model (coefficients, test results, transition and ridge settings) |
| `docs/` | Model card, reproducibility audit (17/17 checks passed), SHA-256 manifest of the frozen analysis files, patient flow, cohort data dictionary |

## Quick start
```bash
cd tests && Rscript run_tests.R          # expect: "all TRES tests passed"
cd ../examples && Rscript minimal_example.R
```
Requires R ≥ 4.2 and the `survival` package (developed with R 4.6.0, survival 3.8.6).

```r
source("R/TRES.R")
fit <- tres_fit(train, static = "T", core = c("N", "PCOSP"),
                sources = list(CAF = c("CAF1","CAF2"), MYELOID = c("MY1","MY2"), LYMPHOID = c("LY1","LY2")))
print(fit)                                  # test results
tres_gate_report(fit)$log_hr                # log hazard ratio per interval
tres_predict(fit, new_patient, months = c(12, 24, 36))
```
Input data need the columns `time` (months), `event` (0/1), `cohort`, and standardised covariates.

## Main results (from the report)
- Simulation with no true time change: false early–late reversals in 80.0% of datasets for the scTREND risk layer vs 1.7% for TRES.
- PDAC development (4 cohorts, 794 patients, leave-one-cohort-out): time-dependent AUC 0.704 and IBS 0.1979, best among nine models.
- Frozen external validation (10 cohorts, 616 patients): time-dependent AUC 0.673 vs 0.575 for AJCC staging; comparable to the static
  T+N+PCOSP model.

## Data availability
No patient-level data are included. All cohorts are public: E-MTAB-6134 and E-MTAB-6830 (ArrayExpress); PCSI, ICGC array, KIRBY, OUH
and Collisson (MetaGxPancreas); TCGA-PAAD (GDC/Xena); GSE21501, GSE57495, GSE62452, GSE224564 and GSE263733 (GEO); CPTAC-PDA and
CPTAC-GDC (PDC/GDC/cBioPortal). `docs/TRES_file_manifest.csv` lists the SHA-256 hashes of the frozen analysis plans, scripts and model
objects used in the study (internal file names keep the development name TASA).

## Limitations
- Research model only; it is not a treatment-response model and not a cleared medical device.
- External scoring used within-cohort standardisation, so a new site currently needs a local reference batch. The frozen object stores
  coefficients and test results but not the reference distribution or development baseline hazard, so absolute survival for a single
  new patient is not yet supported by this package.
- External calibration was not recalibrated (observed/expected 1.15–1.23); prospective validation is required.

## License
MIT (see `LICENSE`).
