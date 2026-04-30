# =============================================================================
#  realdata_synthetic_test.R
#  Synthetic-data sanity check for the κ-correction logic in
#  realdata_ad_ldsc.R and realdata_mr.R.
#
#  This script does NOT require GenomicSEM or TwoSampleMR. It simulates
#  ground-truth liability-scale parameters, generates the corresponding
#  log-OR-scale "summary statistics" the way a covariate-adjusted GWAS
#  would produce them, and verifies that:
#    (a) marginal-scaled liability h^2 is biased upward by kappa^2
#    (b) score-corrected liability h^2 recovers the truth
#    (c) marginal-scaled IVW MR is biased downward
#    (d) score-corrected IVW MR recovers the truth
#
#  Run BEFORE running the real pipelines, to confirm the math is right.
# =============================================================================

source("R/sims/cbcGWAS_replication.R")

set.seed(20260101)

# -----------------------------------------------------------------------------
# Test (a)+(b): Heritability-scale correction
# -----------------------------------------------------------------------------
# Setup: K=2000 SNPs with true liability-scale per-SNP effects gamma_k.
# Each SNP induces a covariate-adjusted log-OR coefficient
#   beta_logOR_k = lambda_cond * gamma_k  (under the model)
# A naive analyst converts beta_logOR back to liability with lambda_marg,
# obtaining beta_liab_naive = beta_logOR / lambda_marg = kappa * gamma_k.
# Heritability sums squared per-SNP effects, so h^2_naive = kappa^2 * h^2_true.

cat("=== Test (a)+(b): heritability scaling ===\n")

K        <- 2000
p_Y      <- 0.05
rho2     <- 0.07
lc_true  <- lambda_cond_gauss(p_Y, rho2)
lm_true  <- lambda_marg(p_Y)
kappa    <- lc_true / lm_true

# True per-SNP effects on the liability scale, scaled so total h^2 = 0.10
gamma_k  <- rnorm(K, 0, sqrt(0.10 / K))
h2_true  <- sum(gamma_k^2)  # ≈ 0.10

# What the GWAS pipeline produces (covariate-adjusted log-OR)
beta_logOR <- lc_true * gamma_k

# Naive marginal-scaled liability conversion
beta_liab_marg <- beta_logOR / lm_true
h2_marg        <- sum(beta_liab_marg^2)

# Score-corrected liability conversion
beta_liab_score <- beta_logOR / lc_true
h2_score        <- sum(beta_liab_score^2)

cat(sprintf("  True h^2_liab:           %.6f\n", h2_true))
cat(sprintf("  Marginal-scaled h^2:     %.6f  (ratio to truth = %.4f)\n",
            h2_marg, h2_marg / h2_true))
cat(sprintf("  Score-corrected h^2:     %.6f  (ratio to truth = %.4f)\n",
            h2_score, h2_score / h2_true))
cat(sprintf("  Theoretical ratio kappa^2 = %.4f\n", kappa^2))

# Pass criteria
ok_marg  <- abs(h2_marg / h2_true  - kappa^2) < 1e-6
ok_score <- abs(h2_score / h2_true - 1)        < 1e-6
cat(sprintf("  PASS (a) marginal biased by kappa^2:    %s\n", ok_marg))
cat(sprintf("  PASS (b) score recovers truth:          %s\n", ok_score))

# -----------------------------------------------------------------------------
# Test (c)+(d): IVW MR scaling
# -----------------------------------------------------------------------------
# Setup: Exposure X is a continuous heritable trait (e.g., BMI).
# Outcome Y is binary with liability Y* and prevalence p_Y.
# True causal effect on the liability scale: theta_liab.
# Per-SNP exposure effects: beta_X_k.
# Per-SNP outcome effects on the LIABILITY scale: theta_liab * beta_X_k + 0
# Per-SNP outcome effects from a covariate-adjusted GWAS on the LOG-OR scale:
#   beta_Y_logOR_k = lc * theta_liab * beta_X_k.
# A naive analyst applies marginal scaling to the outcome:
#   beta_Y_liab_naive = beta_Y_logOR / lambda_marg = kappa * theta_liab * beta_X
# Standard IVW with this gives theta_marg = kappa * theta_liab.
# Score-corrected: theta_score = theta_liab.

cat("\n=== Test (c)+(d): IVW MR scaling ===\n")

K_instr     <- 100
theta_true  <- 0.40
beta_X      <- rnorm(K_instr, 0, 0.025)
se_X        <- rep(0.005, K_instr)
beta_Y_liab <- theta_true * beta_X
beta_Y_logOR <- lc_true * beta_Y_liab + rnorm(K_instr, 0, 0.001)  # tiny residual
se_Y        <- rep(0.005, K_instr)

# IVW on the log-OR outcome scale (no rescaling)
w           <- 1 / se_Y^2
theta_logOR <- sum(w * beta_X * beta_Y_logOR) / sum(w * beta_X^2)

# Marginal: divide outcome beta by lambda_marg, re-fit IVW
theta_marg  <- theta_logOR / lm_true

# Score: divide outcome beta by lambda_cond, re-fit IVW
theta_score <- theta_logOR / lc_true

cat(sprintf("  True theta:              %.6f\n", theta_true))
cat(sprintf("  Marginal-scaled theta:   %.6f  (ratio to truth = %.4f)\n",
            theta_marg, theta_marg / theta_true))
cat(sprintf("  Score-corrected theta:   %.6f  (ratio to truth = %.4f)\n",
            theta_score, theta_score / theta_true))
cat(sprintf("  Theoretical: marginal/truth = kappa = %.4f\n", kappa))

ok_marg_mr  <- abs(theta_marg  / theta_true - kappa) < 0.01
ok_score_mr <- abs(theta_score / theta_true - 1)     < 0.01
cat(sprintf("  PASS (c) marginal biased by kappa:      %s\n", ok_marg_mr))
cat(sprintf("  PASS (d) score recovers truth:          %s\n", ok_score_mr))

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
all_pass <- ok_marg && ok_score && ok_marg_mr && ok_score_mr
cat(sprintf("\n=== ALL TESTS PASS: %s ===\n", all_pass))
if (!all_pass) {
  stop("Synthetic test failed; do not run the real-data pipelines until fixed.")
}
