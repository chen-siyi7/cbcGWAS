# `R/realdata/`

Real-data pipeline scripts that apply the conditional-prevalence correction
to public GWAS summary statistics. These cover §4.2 of the manuscript.

| File           | Manuscript section   | What it does                                                                                  |
| -------------- | -------------------- | --------------------------------------------------------------------------------------------- |
| `helpers.R`    | (utility)            | column-name dictionary, `prepare_sumstats()` normalizer, `preflight()`, `sanity_check_sumstats()` |
| `ad_ldsc.R`    | §4.2.1               | LDSC liability-scale h² + κ-correction on APOE-adjusted Alzheimer's disease GWAS              |
| `bmi_t2d_mr.R` | §4.2.2               | IVW Mendelian randomization with κ-correction for BMI → BMI-adjusted T2D                       |

## Preflight is mandatory

Both pipelines have a `CONFIG` block at the top with a `phase` field. Set
`phase = "PREFLIGHT"` first; the script will:

1. Verify every R package is installed and at a version that has the API it expects
2. Verify every reference file exists and is the expected size
3. Open each sumstat file and probe its column names against the dictionary in `helpers.R`
4. Print a diagnostic summary

If any check fails the script stops with a clear error message before
committing to the slow LDSC munging or PLINK clumping steps. Only switch
`phase` to `"RUN"` once preflight is clean.

## Wall-clock budget

| Step                                        | Wall time         |
| ------------------------------------------- | ----------------- |
| Preflight (either pipeline)                 | <30 s             |
| `ad_ldsc.R` munge + LDSC h²                 | 5–10 min          |
| `bmi_t2d_mr.R` clump + harmonize + IVW      | 15–30 min         |

## Output files

Both pipelines write a tidy `.rds` of the realised numbers into the working
directory at the end of the run, plus a human-readable text summary printed
to the console (which you should redirect to a log file via `tee`):

```bash
Rscript R/realdata/ad_ldsc.R     2>&1 | tee realdata_ad.log
Rscript R/realdata/bmi_t2d_mr.R  2>&1 | tee realdata_mr.log
```

## What's intentionally not here

* **Multivariable MR application** — discussed in §5 of the manuscript as
  future work; clean MVMR-bias-correction requires a covariate-adjusted
  binary exposure GWAS that is not yet publicly available at the scale
  needed for instrument selection.
* **AD → CAD MR** — early drafts considered this; we settled on the
  AD-heritability application in §4.2.1 plus the BMI → T2D MR in §4.2.2 as
  the cleanest pair of demonstrations.

## Sanity check before any of this

Before you touch the real-data pipelines at all, run

```bash
Rscript tests/synthetic_test.R
```

from the repo root. This is base-R-only and confirms the κ-correction
math on simulated ground truth. If the synthetic test fails, the real-data
pipeline will too.
