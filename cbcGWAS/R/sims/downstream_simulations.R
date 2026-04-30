# =============================================================================
#  downstream_simulations.R
#
#  Standalone simulations for Sections 3.6 and 3.7 of the manuscript:
#    Section 3.6 -- Liability-scale heritability under marginal vs score-based
#                   correction (table_heritability)
#    Section 3.7 -- Polygenic-score scaling under marginal vs score-based
#                   correction (table_pgs)
#
#  This file is self-contained: it does not require cbcGWAS_replication.R.
#  It re-defines the small set of helpers (lambda_marg, lambda_cond_gauss)
#  needed for these two simulations, with a fixed SEED_BASE for reproducibility.
#
#  USAGE
#  -----
#    source("R/sims/downstream_simulations.R")
#    h2  <- table_heritability(R = 200)        # ~5 s
#    pgs <- table_pgs(R = 100)                 # ~3-4 min
#    saveRDS(list(h2 = h2, pgs = pgs), "downstream_sim_results.rds")
#
#  Or run as a script:
#    Rscript R/sims/downstream_simulations.R 2>&1 | tee downstream_sim.log
#
#  Send the printed tables (or the .rds) back to update Sections 3.6 and 3.7
#  of the manuscript with actual numbers.
# =============================================================================

# ----------------------------------------------------------------------------- 
# Reproducibility seed -- matches cbcGWAS_replication.R
# ----------------------------------------------------------------------------- 
SEED_BASE <- 20260101L

# ----------------------------------------------------------------------------- 
# Helpers (re-exported here so this file is self-contained)
# ----------------------------------------------------------------------------- 

# Marginal liability-to-logit factor:
#   lambda_marg = phi(t_Y) / [p_Y * (1 - p_Y)]
lambda_marg <- function(p_Y) {
  t_Y <- qnorm(1 - p_Y)
  dnorm(t_Y) / (p_Y * (1 - p_Y))
}

# Conditional liability-to-logit factor under Gaussian M(H), via 200-point
# Gauss-Hermite quadrature.
#
#   M ~ N(0, rho2);  tau^2 = 1 - rho2
#   Numerator:    E_M[ phi((M - t_Y) / tau) / tau ]
#   Denominator:  E_M[ Phi((M - t_Y)/tau) * (1 - Phi(...)) ]
#   lambda_cond = Numerator / Denominator
lambda_cond_gauss <- function(p_Y, rho2, n_quad = 200) {
  t_Y     <- qnorm(1 - p_Y)
  tau     <- sqrt(1 - rho2)
  sigma_M <- sqrt(rho2)

  # Probabilist's Hermite nodes (= physicist's nodes / sqrt(2))
  if (requireNamespace("statmod", quietly = TRUE)) {
    gh <- statmod::gauss.quad.prob(n_quad, dist = "normal")
    nodes   <- gh$nodes
    weights <- gh$weights
  } else {
    # Fallback: trapezoidal integration on a fine grid; less accurate but
    # adequate to ~1e-6 over the (p_Y, rho2) ranges of interest.
    M_grid <- seq(-8, 8, length.out = 4000) * sigma_M
    fM     <- dnorm(M_grid, 0, sigma_M)
    weights <- fM * (M_grid[2] - M_grid[1])
    nodes   <- M_grid / sigma_M
  }

  M       <- sigma_M * nodes
  z       <- (M - t_Y) / tau
  pH      <- pnorm(z)
  dpH     <- dnorm(z) / tau

  num <- sum(weights * dpH)
  den <- sum(weights * pH * (1 - pH))
  num / den
}

# ============================================================================= 
#  Section 3.6: Liability-scale per-SNP heritability under marginal vs
#               score-based correction
# ============================================================================= 
#
#  Verifies that summing squared per-SNP liability-scale effects from a
#  covariate-adjusted GWAS gives the right h^2_liab only under score-based
#  rescaling. Marginal scaling biases h^2_liab upward by kappa^2.
#
#  Setup: K SNPs with true gamma_k ~ N(0, h2_true / K). Per-SNP covariate-
#  adjusted log-OR is lc * gamma_k + epsilon_k, where epsilon_k has variance
#  1 / (n * p_Y * (1 - p_Y)). Convert each beta_logOR to liability via
#  (a) marginal: beta_liab = beta_logOR / lambda_marg, and
#  (b) score:    beta_liab = beta_logOR / lambda_cond.
#  Sum squared per-SNP effects to get h2_marg and h2_score, then subtract the
#  per-SNP sampling-noise floor K * (se / lambda)^2 from each.

table_heritability <- function(R = 200, K = 5000, n = 100000,
                               h2_true = 0.20, p_Y = 0.10,
                               rho2_grid = c(0.10, 0.20, 0.30, 0.40)) {
  out <- data.frame()
  for (rho2 in rho2_grid) {
    set.seed(SEED_BASE + 1100L + round(rho2 * 100))
    lc    <- lambda_cond_gauss(p_Y, rho2)
    lm    <- lambda_marg(p_Y)
    kappa <- lc / lm
    se_C  <- 1 / sqrt(n * p_Y * (1 - p_Y))

    h2_marg_vals  <- numeric(R)
    h2_score_vals <- numeric(R)
    for (r in 1:R) {
      gamma_k    <- rnorm(K, 0, sqrt(h2_true / K))
      beta_logOR <- lc * gamma_k + rnorm(K, 0, se_C)
      h2_marg_vals[r]  <- sum((beta_logOR / lm)^2) - K * (se_C / lm)^2
      h2_score_vals[r] <- sum((beta_logOR / lc)^2) - K * (se_C / lc)^2
    }
    out <- rbind(out, data.frame(
      rho2          = rho2,
      kappa_sq      = kappa^2,
      h2_true       = h2_true,
      h2_marg_mean  = mean(h2_marg_vals),
      h2_marg_mcse  = sd(h2_marg_vals) / sqrt(R),
      h2_score_mean = mean(h2_score_vals),
      h2_score_mcse = sd(h2_score_vals) / sqrt(R),
      ratio_marg_to_truth  = mean(h2_marg_vals)  / h2_true,
      ratio_score_to_truth = mean(h2_score_vals) / h2_true
    ))
  }
  cat("\n=== Section 3.6: Liability-scale heritability under marginal vs score correction ===\n")
  cat("(true h^2 =", h2_true, ", p_Y =", p_Y, ", K =", K,
      ", n =", n, ", R =", R, ")\n")
  print(out, row.names = FALSE, digits = 4)
  invisible(out)
}

# ============================================================================= 
#  Section 3.7: Polygenic-score scaling under marginal vs score-based
#               correction
# ============================================================================= 
#
#  Verifies that a polygenic score built from covariate-adjusted GWAS effects
#  regresses on the true liability with a slope of 1 only under score-based
#  rescaling. Under marginal scaling, the regression slope is biased by kappa
#  relative to the score-based version.
#
#  Setup: K SNPs in a discovery GWAS produce per-SNP covariate-adjusted log-OR
#  coefficients beta_logOR_k. In a target sample of size n_target, draw
#  genotypes G_target_k ~ N(0,1), construct the true liability
#  Y_star = sum_k gamma_k * G_target_k + residual, build
#  PGS_marg  = sum_k (beta_logOR_k / lm) * G_target_k and PGS_score similarly,
#  regress Y_star on each PGS, and report the slope and R^2.

table_pgs <- function(R = 100, K = 2000, n_disc = 100000, n_target = 20000,
                      h2_true = 0.20, p_Y = 0.10,
                      rho2_grid = c(0.10, 0.20, 0.30, 0.40)) {
  out <- data.frame()
  for (rho2 in rho2_grid) {
    set.seed(SEED_BASE + 1200L + round(rho2 * 100))
    lc    <- lambda_cond_gauss(p_Y, rho2)
    lm    <- lambda_marg(p_Y)
    kappa <- lc / lm
    se_C  <- 1 / sqrt(n_disc * p_Y * (1 - p_Y))

    slope_marg_vals  <- numeric(R)
    slope_score_vals <- numeric(R)
    r2_marg_vals     <- numeric(R)
    r2_score_vals    <- numeric(R)
    for (r in 1:R) {
      # Discovery: per-SNP covariate-adjusted log-OR coefficients
      gamma_k    <- rnorm(K, 0, sqrt(h2_true / K))
      beta_logOR <- lc * gamma_k + rnorm(K, 0, se_C)

      # Target: build PGS, regress true liability on each variant
      G_target <- matrix(rnorm(n_target * K), n_target, K)
      Y_star   <- as.vector(G_target %*% gamma_k) +
                    rnorm(n_target, 0, sqrt(1 - h2_true))

      pgs_marg  <- as.vector(G_target %*% (beta_logOR / lm))
      pgs_score <- as.vector(G_target %*% (beta_logOR / lc))

      fit_marg  <- lm(Y_star ~ pgs_marg)
      fit_score <- lm(Y_star ~ pgs_score)
      slope_marg_vals[r]  <- coef(fit_marg)[2]
      slope_score_vals[r] <- coef(fit_score)[2]
      r2_marg_vals[r]     <- summary(fit_marg)$r.squared
      r2_score_vals[r]    <- summary(fit_score)$r.squared
    }
    out <- rbind(out, data.frame(
      rho2              = rho2,
      kappa             = kappa,
      slope_marg_mean   = mean(slope_marg_vals),
      slope_score_mean  = mean(slope_score_vals),
      ratio_marg_score  = mean(slope_marg_vals) / mean(slope_score_vals),
      r2_marg_mean      = mean(r2_marg_vals),
      r2_score_mean     = mean(r2_score_vals)
    ))
  }
  cat("\n=== Section 3.7: Polygenic-score scaling under marginal vs score correction ===\n")
  cat("(true h^2 =", h2_true, ", p_Y =", p_Y, ", K =", K,
      ", n_disc =", n_disc, ", n_target =", n_target, ", R =", R, ")\n")
  print(out, row.names = FALSE, digits = 4)
  invisible(out)
}

# ============================================================================= 
#  Auto-run when called via Rscript (not when source()'d)
# ============================================================================= 
#
#  The check is: was this file launched directly via `Rscript file.R`, as
#  opposed to being source()'d from another R session? We detect the former
#  by examining commandArgs(trailingOnly = FALSE) for "--file=" entries.

.is_rscript_launch <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  any(grepl("^--file=", args))
}

if (.is_rscript_launch()) {
  cat("=== Running Section 3.6 (heritability) ===\n")
  t0 <- Sys.time()
  h2  <- table_heritability(R = 200)
  cat("Wall:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
      "s\n")

  cat("\n=== Running Section 3.7 (PGS) ===\n")
  t0 <- Sys.time()
  pgs <- table_pgs(R = 100)
  cat("Wall:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
      "s\n")

  saveRDS(list(h2 = h2, pgs = pgs, timestamp = Sys.time()),
          file = "downstream_sim_results.rds")
  cat("\n=== saved downstream_sim_results.rds ===\n")
  cat("=== send the .rds (or paste the printed tables) back to update ",
      "Sections 3.6 and 3.7 of the manuscript ===\n")
} else {
  cat("Loaded downstream_simulations.R\n")
  cat("Run table_heritability(R = 200)  for Section 3.6 (~5 s)\n")
  cat("Run table_pgs(R = 100)           for Section 3.7 (~3-4 min)\n")
}
