# TRES frozen model card

(TRES was called TASA during development; internal file names in `TRES_file_manifest.csv` keep that working name.)

## Intended use
Research prediction of 0–36 month postoperative overall survival in resected PDAC. It is not a treatment-response model or a cleared clinical device.

## Frozen analysis population
- Development: 815 available patients (508 deaths); complete T/N analysis: 794 patients (492 deaths).
- External validation: 616 patients across 10 cohorts after six unusable outcome records were removed from the 622-patient readiness pool.

## Final hazard model
For patient i in cohort h and interval j:

`lambda_ihj(t) = lambda_0hj(t) * exp(beta_T*T_i + sum_v z_iv * (b_v + h_j*d_v))`,

where `v` is N or PCOSP, intervals are 0–6, 6–12, 12–24 and 24–36 months, and the cohort-specific baseline hazard is stratified by development cohort.
- Static T coefficient: 0.095712
- Transition weights h: 0.050440, 0.287754, 0.810596, 0.991039
- N log-HR by interval: 0.085471, 0.149803, 0.291538, 0.340453
- PCOSP log-HR by interval: 0.632457, 0.528649, 0.299942, 0.221011
- Source-gate adjusted p-values (CAF, myeloid, lymphoid): 0.682075, 1, 0.364511
- Temporal gate p-value: 3.77598e-05; activated: TRUE; ridge theta: 16
- No single-cell-derived program was retained in the final PDAC survival layer.

## Preprocessing and leakage control
- PCOSP was retrained inside development folds; held-out patients were scored only by models that did not see them.
- Molecular scores were standardised within cohort, followed by development-derived model fitting.
- External outcomes were not used to select programs, time functions, cutoffs or coefficients.

## Primary performance readouts
- Development LOCO: iAUC 0.7040, IBS 0.1979, Uno C 0.6427.
- External: iAUC 0.6731, IBS 0.1983, Uno C 0.6288.

## Deployment limitations
- The validated external pipeline used cohort-wise molecular standardisation; it is therefore a cross-cohort relative-risk model, not yet a fully validated one-patient absolute-risk tool.
- The compact `frozen/TRES_frozen_PDAC_model.rds` stores coefficients and gate decisions but not a self-contained reference distribution and baseline-hazard deployment bundle. These components are reproducible from the frozen code/data but should be serialized before prospective deployment.
- External absolute calibration was not recalibrated and remained imperfect; prospective validation is required.
