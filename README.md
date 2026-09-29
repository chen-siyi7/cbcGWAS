# cbcGWAS

R functions for converting covariate-adjusted binary GWAS associations to a
liability scale, with simulation code from the revised manuscript.

The package separates two steps: estimating and subtracting a collider-bias
contribution, and converting the remaining association to liability units.
It provides the marginal factor, the conditional-prevalence approximation,
and the working-logistic factor. The latter solves the null logistic score
equations and accounts for the sample case fraction.

## Installation

Install the core dependencies, then install the source archive:

```r
install.packages(c("statmod", "mvtnorm"))
install.packages("cbcGWAS_0.1.0.tar.gz", repos = NULL, type = "source")
library(cbcGWAS)
```

Alternatively, extract the GitHub ZIP and install the `cbcGWAS` directory:

```sh
R CMD INSTALL cbcGWAS
```

R 4.1.0 or later is required. No compilation is needed. A GitHub repository
address can be added to `DESCRIPTION` after the repository has been created.

## Conversion factors

```r
library(cbcGWAS)

# Population prevalence and covariate-explained liability variance
p <- 0.10
rho2 <- 0.30

lambda_marginal(p)
lambda_approx(p, rho2)
lambda_working(p, rho2)$lambda

# Case-control GWAS with half cases
factor <- lambda_working(p, rho2, pi = 0.50)$lambda

# Example coefficient after subtracting a separately estimated bias term
adjusted_log_odds <- 0.12
estimated_bias_term <- 0.02
liability_effect <- (adjusted_log_odds - estimated_bias_term) / factor
```

These calculations use a Gaussian liability model and a first-order expansion
in small genetic effects. `p` is population prevalence; `pi` is the case fraction
in the GWAS sample. The conditional-prevalence approximation and the
working-logistic factor are distinct quantities. Bias identification depends
on the chosen estimator's assumptions.

The functions and their return values are documented in R:

```r
?lambda_working
?bias_cwls
?run_simulations
help(package = "cbcGWAS")
```

| Task | Functions |
| --- | --- |
| Conversion | `lambda_marginal()`, `lambda_approx()`, `lambda_working()`, `kappa_cond()` |
| Individual-level associations | `logistic_gwas()`, `linear_gwas()` |
| Population calibration | `probit_rho2()`, `fitted_variance_factor()` |
| Bias-slope estimation | `bias_cwls()`, `bias_slopehunter()` |
| Mendelian randomization | `mr_ivw()`, `mr_raps()`, `mvmr_ivw()` |
| Simulations | `run_simulations()`, `simulate_collider()`, `simulate_end_to_end()`, `summarize_end_to_end()` |

The original names `lambda_marg()`, `lambda_cond_gauss()`, `kappa_score()`,
`bias_dudbridge()`, and `ivw_mr()` remain available for existing analysis code.
`kappa_score()` refers to the conditional-prevalence approximation.

## Simulations

Start with a small run:

```r
small <- run_simulations("smoke")
names(small)
small$population_factors

# Save generated tables, settings, and session information
run_simulations("smoke", output_dir = "smoke_results")
```

The default smoke run produces nine tables. It uses three replicates and
smaller individual-level samples to check that the code runs. Use the full
settings to assess the statistical results.

| Mode | Experiment | Default replicates |
| --- | --- | --- |
| `factors` | Population and case-control factors | Deterministic |
| `population` | No-collider model | 500 |
| `casecontrol` | Four case fractions | 200 |
| `collider` | Unit-variance collider model | 300 |
| `mr` | Univariable and multivariable MR | 100 and 500 |
| `multicovariate` | One, three, and five covariates | 300 |
| `bias` | Illustrative bias-slope regressions | 200 |
| `end-to-end` | Estimated CWLS and Slope-Hunter correction and MR | 200 per design |
| `all` | All experiments | As above |

```r
population <- run_simulations("population", output_dir = "results/population")
collider <- run_simulations("collider", replicates = 100)
full <- run_simulations("end-to-end", output_dir = "results/end_to_end", cores = 1)
```

`replicates` overrides the default count. Full end-to-end runs use 50,000 GWAS
individuals and independent samples for calibration and other associations.
They can take hours on one core. Use `cores = 1` on Windows; additional workers
on other systems require more memory. Fixed seeds make each run reproducible
within the same software environment. The simulation entry points restore the
caller's random-number state. Files are written only when `output_dir` is set.

The `bias` mode evaluates residual-trimmed and precision-weighted regressions.
The published CWLS and Slope-Hunter methods are evaluated in `end-to-end` mode.
The disease labels in the MR experiments refer to simulated traits. Output
columns containing `score` retain the original conditional-prevalence meaning;
the legacy end-to-end method label `Dudbridge` denotes CWLS.

A command-line runner is included in `inst/scripts/run_simulations.R`:

```sh
Rscript cbcGWAS/inst/scripts/run_simulations.R smoke default smoke_results
```

After installation, locate it with
`system.file("scripts", "run_simulations.R", package = "cbcGWAS")`.

## Optional estimators

CWLS, Slope-Hunter, and MR-RAPS wrappers require their respective packages.
The full end-to-end experiment needs all three. Install them separately:

```r
install.packages(c("remotes", "mr.raps"))
remotes::install_github("DudbridgeLab/indexevent")
remotes::install_github("Osmahmoud/SlopeHunter")
```

The revision used R 4.5.1, statmod 1.5.2, mvtnorm 1.4-2, indexevent 0.2.0,
SlopeHunter 1.1.0, and mr.raps 0.4.3. Core conversion and smoke simulations
do not require the optional estimators.

## Development

The source includes help pages, examples, input checks, and tests for numerical
identities, regression estimates, calibration, and simulation reproducibility.

```sh
R CMD build cbcGWAS
R CMD check --no-manual cbcGWAS_0.1.0.tar.gz
```

Install `testthat` and the suggested estimator packages before a complete
check. Documentation is generated with `roxygen2::roxygenise("cbcGWAS")`.
Loading the package does not run simulations, install dependencies, or write
files. The archive contains source code, documentation, and tests, with no
manuscript files or saved simulation results.

## License

A redistribution license has not yet been selected by the author. See `LICENSE`.
