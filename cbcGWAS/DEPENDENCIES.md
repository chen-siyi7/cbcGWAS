# Dependencies

This file lists every external dependency for `cbcGWAS`, separated
into a small "simulations only" set and a larger "real-data pipelines" set
needed for §4.2 of the manuscript.

The versions listed below are those used to produce the manuscript numbers.
Other versions will most likely work, but exact bit-for-bit reproduction
is only guaranteed at these versions.

---

## R itself

R ≥ 4.0 is required. The manuscript was produced under

```
R version 4.3.3 (2024-02-29) -- "Angel Food Cake"
```

Anything in the 4.x series should work. R 3.x is not supported.

---

## Simulations and synthetic test (`R/sims/`, `tests/`)

These scripts are **base R + 2 packages**:

| Package   | Version | Purpose                                    |
| --------- | ------- | ------------------------------------------ |
| `statmod` | 1.5.0   | `gauss.quad.prob()` for `lambda_cond_gauss` |
| `mvtnorm` | 1.2.4   | correlated MVN sampling in some sims        |

Install:

```r
install.packages(c("statmod", "mvtnorm"))
```

Both are pure-R, lightweight, and on CRAN. `statmod` provides Gauss–Hermite
quadrature nodes used inside `lambda_cond_gauss(p, rho2)`; `mvtnorm` is used
for correlated-instrument simulations in `table_mvmr()` and similar.

---

## Real-data pipelines (`R/realdata/`)

These add a substantially larger stack and pull in genetics-specific tools.

### R packages

| Package             | Version  | Source                                    | Purpose                                |
| ------------------- | -------- | ----------------------------------------- | -------------------------------------- |
| `data.table`        | ≥ 1.14   | CRAN                                      | fast TSV/CSV I/O for sumstat files     |
| `remotes`           | ≥ 2.4    | CRAN                                      | install GitHub-only packages           |
| `GenomicSEM`        | ≥ 0.0.5  | GitHub: `GenomicSEM/GenomicSEM`           | LDSC munging + h² estimation           |
| `TwoSampleMR`       | ≥ 0.6.0  | GitHub: `MRCIEU/TwoSampleMR`              | IVW MR, sensitivity analyses, harmonization |
| `ieugwasr`          | ≥ 1.0    | GitHub: `MRCIEU/ieugwasr`                 | local PLINK clumping wrapper           |
| `genetics.binaRies` | ≥ 0.0.6  | GitHub: `MRCIEU/genetics.binaRies`        | bundled PLINK binary                   |

```r
install.packages(c("data.table", "remotes"))
remotes::install_github("GenomicSEM/GenomicSEM")
remotes::install_github("MRCIEU/TwoSampleMR")
remotes::install_github("MRCIEU/genetics.binaRies")
remotes::install_github("MRCIEU/ieugwasr")
```

`GenomicSEM` itself depends on `gdata`, `Matrix`, `lavaan`, `qgraph`,
which install automatically. Total install time: 5–10 minutes on a typical
laptop.

### External binaries

| Tool       | Version | Source                                                  |
| ---------- | ------- | ------------------------------------------------------- |
| PLINK 1.9  | b7.2    | <https://www.cog-genomics.org/plink/1.9/> (or via `genetics.binaRies`) |

PLINK is required for `ieugwasr::ld_clump_local()` in the MR pipeline.
The `genetics.binaRies` package ships a working binary; you can also point
the pipeline at a system-installed PLINK by editing `CONFIG$plink_bin` in
`R/realdata/bmi_t2d_mr.R`.

### Reference data

| File / directory     | Size    | Source                                                                                  | Purpose                       |
| -------------------- | ------- | --------------------------------------------------------------------------------------- | ----------------------------- |
| `eur_w_ld_chr/`      | ~5 GB   | LDSC reference, EUR, 23 chromosome files (mirror via Zenodo: <https://zenodo.org/records/8182036>) | LDSC h² estimation            |
| `w_hm3.snplist`      | ~3 MB   | same source                                                                              | HapMap3 munging reference     |
| 1000G EUR PLINK fileset | ~1 GB | <http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz>                                          | local clumping LD reference   |

The 1000G LD reference is the same one shipped by the MR-Base group and is
what `ieugwasr::ld_clump_local()` expects. The LDSC reference files are
maintained by the original LDSC authors (Bulik-Sullivan et al. 2015); an
up-to-date Zenodo mirror is given above because the Broad-hosted directory
listing has changed format.

### GWAS summary statistics

The manuscript uses three public GWAS releases. None are redistributed
in this repository; download them yourself from the original sources:

| Trait                        | Source / accession                                               |
| ---------------------------- | ---------------------------------------------------------------- |
| Alzheimer's disease          | Bellenguez et al. 2022, GWAS Catalog GCST90027158                |
| BMI (continuous exposure)    | Yengo et al. 2018, GIANT Consortium downloads                    |
| BMI-adjusted T2D (outcome)   | Mahajan et al. 2018, DIAGRAM Consortium downloads                |

After download, edit the `CONFIG` block at the top of each real-data
script to point at the local file paths. Run with `CONFIG$phase = "PREFLIGHT"`
first to verify column names and paths; then switch to `"RUN"`.

---

## Operating system

All development and testing was done on macOS 14 (Apple Silicon and Intel)
and Linux x86_64 (Ubuntu 22.04, Ubuntu 24.04). Windows users should be fine
for the simulation-only pieces; the real-data pipelines have not been
tested on Windows but should work under WSL.

---

## Quick check that your setup is good

From the repo root:

```bash
Rscript tests/synthetic_test.R
```

Should print `=== ALL TESTS PASS: TRUE ===` in about 10 seconds. If this
works, the simulation half of the codebase is fully functional. The
real-data half can be sanity-checked separately by running the relevant
script with `CONFIG$phase = "PREFLIGHT"`, which inspects all paths and
column names without committing to any long downloads or LDSC / MR runs.
