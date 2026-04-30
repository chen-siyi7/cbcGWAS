# cbcGWAS

R code accompanying:

> **Conditional-prevalence correction for liability-scale conversion of
> covariate-adjusted binary GWAS.**
> Siyi Chen, LSU Health Sciences Center New Orleans (manuscript under review).

This repository reproduces every numerical result in the main text and
supplement, and provides the real-data pipeline applying the correction to
published Alzheimer's disease and BMI/T2D GWAS summary statistics.

---

## What the correction does

A standard liability-scale conversion divides a log-odds GWAS coefficient by

```
lambda_marg = phi(t) / [p (1 - p)]
```

where `p` is the population disease prevalence and `t = Phi^{-1}(1 - p)`. This
factor is correct for *unadjusted* binary GWAS. When the GWAS coefficient is
estimated *conditional* on heritable covariates, the correct factor is

```
lambda_cond  >=  lambda_marg
```

a quantity that depends additionally on `rho^2`, the liability variance
explained by the adjusted covariates. The ratio

```
kappa_score = lambda_cond / lambda_marg  >=  1
```

is the bias factor of the standard conversion. This repository computes
`lambda_cond` via Gauss-Hermite quadrature and applies it consistently across
heritability conversion, polygenic-score scaling, genetic correlation, and
univariable / multivariable Mendelian randomization.

For BMI-adjusted T2D at `(p, rho^2) = (0.10, 0.30)` the correction is
sizeable: `kappa_score ≈ 1.148`, a ~13 % bias under the marginal default.
For APOE-adjusted Alzheimer's disease the correction is small:
`kappa_score ≈ 1.01` to `1.04` across plausible APOE liability shares.

---

## Repository layout

```
cbcGWAS/
├── R/
│   ├── sims/
│   │   ├── cbcGWAS_replication.R     # all manuscript / supplement tables
│   │   ├── downstream_simulations.R   # heritability + PGS sims (§3.6, §3.7)
│   │   └── run_high_R.R               # high-precision R = 800 reruns
│   └── realdata/
│       ├── helpers.R                  # column dictionary, preflight, sumstat normalizer
│       ├── ad_ldsc.R                  # LDSC heritability + κ-correction (§4.2.1)
│       └── bmi_t2d_mr.R               # IVW MR + κ-correction (§4.2.2)
├── tests/
│   └── synthetic_test.R               # no-dependency sanity check (recommended first run)
├── examples/
│   ├── 01_basic_factors.R             # minimal κ-correction example
│   └── 02_table_kappa.R               # reproduce manuscript Table 1
├── data/                              # (empty; user supplies GWAS sumstats)
├── DEPENDENCIES.md                    # full package list with versions
├── LICENSE                            # MIT
└── README.md                          # this file
```

All scripts assume you `cd` into the repository root and run them from there
(e.g. `Rscript R/sims/cbcGWAS_replication.R`). Internal `source()` calls use
paths relative to the repo root.

---

## Quick start

### 1. Verify the math on synthetic data (10 seconds, base R only)

```bash
git clone https://github.com/sychen-lsuhsc/cbcGWAS
cd cbcGWAS
Rscript tests/synthetic_test.R
```

You should see four `PASS ... TRUE` lines and `=== ALL TESTS PASS: TRUE ===`.
This confirms the κ-correction logic on simulated ground truth before any
external dependencies are touched.

### 2. Reproduce a manuscript table

```r
# from R, with working directory = repo root
source("R/sims/cbcGWAS_replication.R")
table_kappa()        # Table 1: kappa_score across (p, rho^2)
table_no_collider()  # Table 2: simple-design simulations
```

`run_all()` reproduces every table at the seeds used in the manuscript;
expect 30–60 minutes on a recent laptop.

### 3. Reproduce the §3.6–3.7 downstream simulations

```bash
Rscript R/sims/downstream_simulations.R 2>&1 | tee downstream_sim.log
```

Wall time: ~5 s for heritability, ~3–4 min for the polygenic-score table.

### 4. Real-data pipelines (require external GWAS files; see below)

```bash
# Edit the CONFIG block at the top of the script to point at your local files.
# Run with phase = "PREFLIGHT" first, then "RUN".
Rscript R/realdata/ad_ldsc.R
Rscript R/realdata/bmi_t2d_mr.R
```

---

## Dependencies

**Simulations only** (`R/sims/`, `tests/`): base R ≥ 4.0, plus

* `statmod` (Gauss-Hermite nodes for `lambda_cond`)
* `mvtnorm` (correlated sampling in some simulations)

```r
install.packages(c("statmod", "mvtnorm"))
```

**Real-data pipelines** also need:

* `data.table` — sumstat I/O
* `GenomicSEM` — LDSC munging and ldsc h² estimation
* `TwoSampleMR`, `ieugwasr`, `genetics.binaRies` — MR pipeline (local clumping)
* PLINK v1.9 binary — invoked by `ieugwasr::ld_clump_local`
* 1000 Genomes EUR PLINK reference for clumping (~1 GB)
* HapMap3 reference SNP list `w_hm3.snplist` and `eur_w_ld_chr/` LDSC reference
  (~5 GB unpacked)

Full installation instructions are in `DEPENDENCIES.md`.

The pipelines do **not** redistribute any GWAS summary statistics; the user
downloads them from the original sources (Bellenguez et al. 2022 via
GWAS Catalog GCST90027158; Yengo et al. 2018 from GIANT; Mahajan et al.
2018 from DIAGRAM).

---

## What each script reproduces

| Manuscript element | Script(s) | Function |
|---|---|---|
| Closed-form factors (Table 1) | `R/sims/cbcGWAS_replication.R` | `table_kappa()` |
| Simple-design simulations (Table 2) | `R/sims/cbcGWAS_replication.R` | `table_no_collider()` |
| Collider-aware simulations (Table 3) | `R/sims/cbcGWAS_replication.R` | `table_collider()` |
| MVMR simulations (Table 4) | `R/sims/cbcGWAS_replication.R` | `table_mvmr()` |
| Heritability simulation (Table 5) | `R/sims/downstream_simulations.R` | `table_heritability()` |
| Polygenic-score simulation (Table 6) | `R/sims/downstream_simulations.R` | `table_pgs()` |
| Wald-test invariance (Table S1) | `R/sims/cbcGWAS_replication.R` | `table_wald()` |
| Sensitivity to ρ² misspecification (Table S2) | `R/sims/cbcGWAS_replication.R` | `table_sensitivity()` |
| Prentice–Pyke case-control (Table S3) | `R/sims/cbcGWAS_replication.R` | `table_cc()` |
| Bias-slope estimator recovery (Table S5) | `R/sims/cbcGWAS_replication.R` | `table_slope_estimators()` |
| Non-IVW MVMR estimators (Table S11) | `R/sims/cbcGWAS_replication.R` | `table_non_ivw_mvmr()` |
| AD heritability (§4.2.1) | `R/realdata/ad_ldsc.R` | `Rscript R/realdata/ad_ldsc.R` |
| BMI→T2D MR (§4.2.2) | `R/realdata/bmi_t2d_mr.R` | `Rscript R/realdata/bmi_t2d_mr.R` |

Tables S3, S5, S6 can be re-run at the higher precision used for some
sensitivity analyses via `Rscript R/sims/run_high_R.R` (R = 800 replicates,
~8 minutes).

---

## Reproducibility

Every simulation function fixes a deterministic seed via the `SEED_BASE`
constant at the top of `R/sims/cbcGWAS_replication.R` (and re-seeds inside
`downstream_simulations.R`). Re-running with the same R version and package
versions (see `DEPENDENCIES.md`) reproduces the manuscript numbers
bit-for-bit.

The only stochastic component left is in `R/realdata/bmi_t2d_mr.R`, where
the LD-clumping output depends on the exact PLINK version and 1000 Genomes
reference; we report the harmonized SNP count (K = 505) and the realised
estimates explicitly so a re-run can be checked for consistency rather than
identity.

---

## Citation

If you use this code, please cite the manuscript (currently under review;
this README will be updated with the journal reference upon acceptance).
The bibliography in the manuscript also cites the upstream methods this
work composes with — most importantly Aschard et al. 2015, Lee et al. 2011,
Mahmoud et al. 2022 (Slope-Hunter), Dudbridge et al. 2019 (index-event),
Wang et al. 2024 (MVMR-cML-bias-correction), and the LDSC, GCTA/mtCOJO,
and TwoSampleMR software stacks.

---

## License

MIT. See `LICENSE`.

The author thanks colleagues at the LSU Health Sciences Center New Orleans
School of Public Health, and acknowledges the open-source software developers
whose tools (LDSC, GCTA/mtCOJO, Slope-Hunter, MVMR-cML, MVMR-cML-bias-correction,
MendelianRandomization, TwoSampleMR, REGENIE, SAIGE) made this work possible.

---

## Issues and contact

Please file issues at <https://github.com/sychen-lsuhsc/cbcGWAS/issues>.
For correspondence regarding the manuscript, contact Siyi Chen at LSU
Health Sciences Center New Orleans (email in the manuscript).
