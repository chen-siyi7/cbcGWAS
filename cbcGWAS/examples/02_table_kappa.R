# =============================================================================
#  examples/02_table_kappa.R
#
#  Reproduce Table 1 of the manuscript -- kappa_score across a grid of
#  (p_Y, rho^2). This is the core summary table showing that the bias of
#  the standard marginal conversion grows with both prevalence and the
#  covariate's liability share.
#
#  Run from the repo root:
#    Rscript examples/02_table_kappa.R
#
#  Wall time: ~1 second. Base R + statmod only (loaded by the source below).
# =============================================================================

source("R/sims/cbcGWAS_replication.R")

# table_kappa() prints the full grid and returns it invisibly
tab <- table_kappa()

# Save a CSV next to this script for easy import into Excel / pandas
out_path <- "examples/table1_kappa.csv"
write.csv(tab, out_path)
cat(sprintf("\nSaved CSV: %s\n", out_path))

# -----------------------------------------------------------------------------
# Reading the table:
#
#   - Each cell is kappa_score = lambda_cond / lambda_marg, the multiplicative
#     bias of the standard liability conversion when the GWAS adjusts for a
#     covariate that explains rho^2 of the liability variance.
#   - kappa = 1 when rho^2 = 0 (no adjustment needed).
#   - kappa grows with both p_Y and rho^2.
#   - For the practical T2D / BMI cell (p_Y = 0.10, rho^2 = 0.30): kappa ~ 1.148,
#     i.e. ~13% upward bias if the marginal factor is used naively.
# -----------------------------------------------------------------------------
