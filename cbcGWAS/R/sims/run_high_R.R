# =============================================================================
#  run_high_R.R
#  Re-run Supplementary Tables S3, S5, S6 at R = 800 replicates
# =============================================================================
#
#  Usage (from the repo root):
#    Rscript R/sims/run_high_R.R 2>&1 | tee high_R.log
#
#  Estimated wall time: ~8 min on a typical laptop.
#  Outputs:
#    high_R_results.rds      -- saved data frames
#    high_R.log              -- full console output
#
#  All seeds are set deterministically via SEED_BASE inside each table
#  function, so the run is fully reproducible.
# =============================================================================

cat("=== loading replication code ===\n")
source("R/sims/cbcGWAS_replication.R")

R_HIGH <- 800L

# -----------------------------------------------------------------------------
#  S3 -- Prentice-Pyke case-control invariance
# -----------------------------------------------------------------------------
cat("\n=== Table S3: case-control invariance at R =", R_HIGH, "===\n")
t0 <- Sys.time()
s3 <- table_cc(R = R_HIGH, n = 20000, p_Y = 0.10, rho2 = 0.30, gamma = 0.05)
cat("Wall:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
print(s3, row.names = FALSE, digits = 6)

# -----------------------------------------------------------------------------
#  S5 -- Bias-slope robustness to SNP-covariate effect magnitude
# -----------------------------------------------------------------------------
cat("\n=== Table S5: sigma_beta robustness at R =", R_HIGH, "===\n")
t0 <- Sys.time()
s5 <- table_sigma_beta(R = R_HIGH, K = 200, n = 50000, b_lin = 0.50)
cat("Wall:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
print(s5, row.names = FALSE, digits = 6)

# -----------------------------------------------------------------------------
#  S6 -- Multi-covariate adjustment
# -----------------------------------------------------------------------------
cat("\n=== Table S6: multi-covariate adjustment at R =", R_HIGH, "===\n")
t0 <- Sys.time()
s6 <- table_multi_covar(R = R_HIGH, n = 20000, p_Y = 0.10, gamma = 0.05)
cat("Wall:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
print(s6, row.names = FALSE, digits = 6)

# -----------------------------------------------------------------------------
#  Save raw outputs
# -----------------------------------------------------------------------------
saveRDS(list(s3 = s3, s5 = s5, s6 = s6, R = R_HIGH),
        file = "high_R_results.rds")
cat("\n=== saved high_R_results.rds ===\n")
cat("=== send the .rds (or paste the printed output) back to update tables ===\n")
