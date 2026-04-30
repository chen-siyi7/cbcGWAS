# =============================================================================
#  examples/01_basic_factors.R
#
#  Minimal, ~30-line example showing how to use the three core functions:
#    lambda_marg(p_Y)          -- standard liability conversion factor
#    lambda_cond_gauss(p, rho) -- conditional-prevalence factor
#    kappa_score(p, rho)       -- ratio of the two (the bias factor)
#
#  Run from the repo root:
#    Rscript examples/01_basic_factors.R
#
#  No external packages required for these three functions.
# =============================================================================

source("R/sims/cbcGWAS_replication.R")

# -----------------------------------------------------------------------------
# Setting: a covariate-adjusted binary GWAS
#   p_Y  = population disease prevalence            (e.g. 0.10 for T2D)
#   rho2 = liability variance explained by the
#          adjusted covariate                       (e.g. 0.30 for BMI's
#                                                   contribution to T2D)
# -----------------------------------------------------------------------------

p_Y  <- 0.10
rho2 <- 0.30

cat(sprintf("Inputs: p_Y = %.2f, rho^2 = %.2f\n\n", p_Y, rho2))

L_marg <- lambda_marg(p_Y)
L_cond <- lambda_cond_gauss(p_Y, rho2)
kappa  <- kappa_score(p_Y, rho2)

cat(sprintf("lambda_marg          = %.6f\n", L_marg))
cat(sprintf("lambda_cond (gauss)  = %.6f\n", L_cond))
cat(sprintf("kappa_score          = %.6f   (= lambda_cond / lambda_marg)\n\n",
            kappa))

# -----------------------------------------------------------------------------
# Worked example: converting a covariate-adjusted log-OR coefficient back to
# the liability scale.
#
# Suppose a SNP's BMI-adjusted T2D log-OR coefficient is beta_logOR = 0.15,
# with standard error 0.020. Then:
# -----------------------------------------------------------------------------

beta_logOR <- 0.15
se_logOR   <- 0.020

beta_liab_marg  <- beta_logOR / L_marg
beta_liab_score <- beta_logOR / L_cond

cat("Worked example -- BMI-adjusted T2D SNP:\n")
cat(sprintf("  beta (log-OR scale)         = %.4f (SE %.4f)\n",
            beta_logOR, se_logOR))
cat(sprintf("  beta (marginal-converted)   = %.4f  <-- biased upward by kappa = %.3f\n",
            beta_liab_marg, kappa))
cat(sprintf("  beta (score-corrected)      = %.4f  <-- recovers the liability-scale truth\n",
            beta_liab_score))

# -----------------------------------------------------------------------------
# Note: the bias is multiplicative, not additive. For heritability (a sum of
# squared liability-scale effects) the bias is kappa^2; for IVW MR (a ratio
# of liability-scale slopes) the bias is kappa. See examples/02_table_kappa.R
# for the full grid of factors used in the manuscript.
# -----------------------------------------------------------------------------
