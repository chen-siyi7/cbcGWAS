#  =============================================================================
#  cbcGWAS: Replication code for all simulations
#  =============================================================================
#
#  This script replicates every numerical result in the main manuscript and
#  supplementary materials. Each function corresponds to one table.
#
#  Usage (from repo root):
#    source("R/sims/cbcGWAS_replication.R")
#    run_all()                        # all tables, ~30-60 min on a laptop
#    run_quick()                      # closed-form tables only, < 1 min
#    table_kappa()                    # specific table
#
#  Dependencies: base R, stats, statmod (for Gauss-Hermite nodes), mvtnorm
#  =============================================================================

# ----------------------------------------------------------------------------- 
# Setup and dependencies
# ----------------------------------------------------------------------------- 

# Defensive: remove any stale function definitions from the global environment
# that we are about to redefine. Without this, a half-loaded earlier source()
# call can leave incomplete versions in scope.
.cbcGWAS_funcs <- c(
  "lambda_marg", "lambda_cond_gauss", "kappa_score", "lambda_cond_discrete",
  "lambda_cond_exact", "lambda_cond_exact_logit",
  "sim_one_no_collider", "sim_one_collider",
  "sim_t2d_cad_one", "run_t2d_cad", "sim_mvmr_one",
  "ivw_mr", "mvmr_ivw", ".print_rounded",
  ".mvmr_cml_one", ".mvmr_median_one", ".mvmr_egger_one",
  "diagnose_mvmr_packages",
  "table_kappa", "table_ad_apoe", "table_bias_coverage", "table_rho_sens",
  "table_rho_sens_n20k", "table_t2d_cad", "table_theta_sens",
  "table_residual_decomp", "table_mvmr", "table_mvmr_strength",
  "table_mvmr_non_ivw", "table_pl_gap", "table_nongaussian",
  "table_cc", "table_bias_slope_recovery", "table_sigma_beta",
  "table_multi_covar", "table_heritability", "table_pgs",
  "run_quick", "run_all"
)
rm(list = intersect(.cbcGWAS_funcs, ls(envir = .GlobalEnv)),
      envir = .GlobalEnv)
rm(.cbcGWAS_funcs)

if (!requireNamespace("statmod",  quietly = TRUE)) install.packages("statmod")
if (!requireNamespace("mvtnorm",  quietly = TRUE)) install.packages("mvtnorm")
if (!requireNamespace("MASS",     quietly = TRUE)) install.packages("MASS")

library(statmod)
library(mvtnorm)
library(MASS)

# Optional packages for Supplementary Table S9 (MVMR-cML, MR-Median, MR-Egger):
#   install.packages("MendelianRandomization")
#   # MVMRcML is GitHub-only:
#   # install.packages("remotes")
#   # remotes::install_github("ZhaotongL/MVMR-cML")
# If absent, table_mvmr_non_ivw() reports IVW only with a notice.

# Reproducibility seed bank: each table uses a derived seed for independence.
SEED_BASE <- 20260101L

# Helper: round numeric columns of a data frame; print without rounding factors/strings.
.print_rounded <- function(df, digits = 3) {
  if (is.data.frame(df)) {
    num_cols <- sapply(df, is.numeric)
    df[, num_cols] <- round(df[, num_cols], digits)
    print(df)
  } else {
    print(round(df, digits))
  }
}

# ============================================================================= 
# Core: closed-form lambda_cond and kappa_score
# ============================================================================= 

#' Marginal liability-to-logit factor: phi(t_Y) / [p_Y * (1 - p_Y)]
lambda_marg <- function(p_Y) {
  t_Y <- qnorm(1 - p_Y)
  dnorm(t_Y) / (p_Y * (1 - p_Y))
}

#' Conditional liability-to-logit factor under Gaussian-M assumption.
#' Computed by 200-point Gauss-Hermite quadrature.
lambda_cond_gauss <- function(p_Y, rho2, n_quad = 200) {
  if (rho2 <= 0) return(lambda_marg(p_Y))
  t_Y <- qnorm(1 - p_Y)
  tau <- sqrt(1 - rho2)
  rho <- sqrt(rho2)
  gh  <- statmod::gauss.quad(n_quad, kind = "hermite")
  M   <- sqrt(2) * rho * gh$nodes
  p_H <- pnorm((M - t_Y) / tau)
  D   <- sum(gh$weights * p_H * (1 - p_H)) / sqrt(pi)
  dnorm(t_Y) / D
}

kappa_score <- function(p_Y, rho2, n_quad = 200) {
  lambda_cond_gauss(p_Y, rho2, n_quad) / lambda_marg(p_Y)
}

#' Discrete-M branch for binary or categorical adjustment covariates.
#' Used in supplement equation (S1.21).
lambda_cond_discrete <- function(p_Y, rho2, p_H = 0.5) {
  t_Y     <- qnorm(1 - p_Y)
  tau     <- sqrt(1 - rho2)
  delta_H <- sqrt(rho2 / (p_H * (1 - p_H)))
  m0  <- -delta_H * p_H        # H = 0, centered
  m1  <-  delta_H * (1 - p_H)  # H = 1, centered
  p0  <- pnorm((m0 - t_Y) / tau)
  p1  <- pnorm((m1 - t_Y) / tau)
  num <- (1 - p_H) * dnorm((m0 - t_Y) / tau) / tau +
              p_H  * dnorm((m1 - t_Y) / tau) / tau
  den <- (1 - p_H) * p0 * (1 - p0) +
              p_H  * p1 * (1 - p1)
  num / den
}

# ============================================================================= 
# Table 5 / S5: kappa_score across (p_Y, rho^2)
# ============================================================================= 

table_kappa <- function() {
  p_Y_grid  <- c(0.50, 0.20, 0.10, 0.05, 0.01)
  rho2_grid <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
  out <- outer(rho2_grid, p_Y_grid,
               Vectorize(function(rho2, p_Y) kappa_score(p_Y, rho2)))
  rownames(out) <- paste0("rho2=", rho2_grid)
  colnames(out) <- paste0("p_Y=",  p_Y_grid)
  cat("\n=== Main Table 5: kappa_score ===\n")
  .print_rounded(out, 3)
  invisible(out)
}

# =============================================================================
# Main Table 6: APOE-adjusted Alzheimer's disease application
# =============================================================================
# Closed-form kappa_score across the published range of APOE's contribution
# to AD liability variance, rho^2_APOE in [0.04, 0.13] (Liu et al. 2025
# PLOS Genet, citing studies that directly genotyped epsilon-2 / epsilon-4),
# at AD population prevalence p_Y = 0.05 (the standard value used in Wightman
# 2021, de la Fuente / Grotzinger 2022, and Liu 2025 review). Marginal
# liability-to-logit factor lambda_marg is constant across rows because p_Y
# is fixed.

table_ad_apoe <- function() {
  p_Y       <- 0.05
  rho2_grid <- c(0.04, 0.07, 0.10, 0.13)
  l_marg    <- lambda_marg(p_Y)
  l_cond    <- sapply(rho2_grid, function(r) lambda_cond_gauss(p_Y, r))
  kappa     <- l_cond / l_marg
  out <- data.frame(
    rho2_APOE        = rho2_grid,
    lambda_marg      = round(l_marg, 3),
    lambda_cond      = round(l_cond, 3),
    kappa_score      = round(kappa,  3),
    correction_pct   = round(100 * (kappa - 1), 1)
  )
  cat("\n=== Main Table 6: APOE-adjusted AD, p_Y = 0.05 ===\n")
  print(out, row.names = FALSE)
  invisible(out)
}

# ============================================================================= 
# Helpers: data-generating model under explicit collider structure
# ============================================================================= 

#' Simulate one replicate of the no-collider model and return mean kappa estimates.
#' Population logistic-regression slope of Y ~ G + H, no collider bias.
sim_one_no_collider <- function(n, p_Y, rho2, gamma = 0.05) {
  H <- rnorm(n, 0, sqrt(rho2))
  G <- rnorm(n)              # standardized SNP
  # Liability: rho^2 from H, gamma*G direct, residual to make var = 1
  res_var <- max(1 - rho2 - gamma^2, 1e-8)
  Y_star  <- H + gamma * G + rnorm(n, 0, sqrt(res_var))
  t_Y     <- qnorm(1 - p_Y)
  Y       <- as.numeric(Y_star > t_Y)
  fit     <- glm(Y ~ G + H, family = binomial())
  beta_C  <- coef(fit)["G"]
  se_C    <- summary(fit)$coefficients["G", "Std. Error"]
  list(beta_C = beta_C, se_C = se_C)
}

#' Simulate explicit-collider model with unmeasured confounder U.
sim_one_collider <- function(n, p_Y, rho2, gamma = 0.05,
                                       alpha_U = sqrt(0.30),
                                       alpha_UH = 0.55,
                                       beta_GH  = 0.10,
                                       sigma_zeta2 = 0.45) {
  U     <- rnorm(n)
  G     <- rnorm(n)
  zeta  <- rnorm(n, 0, sqrt(sigma_zeta2))
  H     <- beta_GH * G + alpha_UH * U + zeta
  # Standardize H to var roughly 1 in population for cleanest scaling
  H     <- H / sd(H)
  # Liability: scale so total var = 1 with variance contribution rho^2 from H
  delta_H <- sqrt(rho2)
  res_var <- max(1 - rho2 - gamma^2 - alpha_U^2, 1e-6)
  Y_star  <- gamma * G + delta_H * H + alpha_U * U + rnorm(n, 0, sqrt(res_var))
  t_Y     <- qnorm(1 - p_Y)
  Y       <- as.numeric(Y_star > t_Y)
  fit     <- glm(Y ~ G + H, family = binomial())
  beta_C  <- coef(fit)["G"]
  se_C    <- summary(fit)$coefficients["G", "Std. Error"]
  list(beta_C = beta_C, se_C = se_C)
}

# ============================================================================= 
# Table 6: per-SNP relative bias and coverage (no-collider model)
# ============================================================================= 

table_bias_coverage <- function(R = 500, n = 20000, gamma = 0.05) {
  cells <- list(
    c(0.10, 0.05), c(0.10, 0.10), c(0.20, 0.10),
    c(0.05, 0.30), c(0.10, 0.30), c(0.20, 0.30),
    c(0.50, 0.30), c(0.50, 0.50)
  )
  out <- data.frame(p_Y = numeric(), rho2 = numeric(),
                                      kappa = numeric(),
                                      marg_bias = numeric(), marg_mcse = numeric(),
                                      score_bias = numeric(), score_mcse = numeric(),
                                      marg_cov = numeric(),  score_cov = numeric())
  for (k in seq_along(cells)) {
    set.seed(SEED_BASE + 100L + k)
    p_Y  <- cells[[k]][1]
    rho2 <- cells[[k]][2]
    kap  <- kappa_score(p_Y, rho2)
    lm   <- lambda_marg(p_Y)
    lc   <- lm * kap

    rel_marg  <- numeric(R)
    rel_score <- numeric(R)
    cov_marg  <- logical(R)
    cov_score <- logical(R)
    for (r in 1:R) {
      rep <- sim_one_no_collider(n, p_Y, rho2, gamma)
      beta_marg  <- rep$beta_C / lm
      beta_score <- rep$beta_C / lc
      se_marg    <- rep$se_C   / lm
      se_score   <- rep$se_C   / lc
      rel_marg[r]  <- beta_marg  / gamma - 1
      rel_score[r] <- beta_score / gamma - 1
      cov_marg[r]  <- abs(beta_marg  - gamma) < 1.96 * se_marg
      cov_score[r] <- abs(beta_score - gamma) < 1.96 * se_score
    }
    out[k, ] <- list(
      p_Y, rho2, kap,
      mean(rel_marg),  sd(rel_marg)  / sqrt(R),
      mean(rel_score), sd(rel_score) / sqrt(R),
      mean(cov_marg)  * 100,
      mean(cov_score) * 100
    )
  }
  cat("\n=== Main Table 6: per-SNP relative bias and coverage ===\n")
  .print_rounded(out, 3)
  invisible(out)
}

# ============================================================================= 
# Table 7: sensitivity to misspecified rho^2
# ============================================================================= 

table_rho_sens <- function(R = 300, n = 5000, gamma = 0.05, p_Y = 0.10,
                                                      header = NULL) {
  if (is.null(header)) {
    header <- sprintf("=== Main Table 7: sensitivity to rho^2 misspecification (n=%d) ===",
                                          as.integer(n))
  }
  true_grid    <- c(0.20, 0.30)
  assumed_grid <- c(0.05, 0.10, 0.20, 0.30, 0.40)

  out <- data.frame(true_rho2 = numeric(), assumed_rho2 = numeric(),
                                      rel_bias = numeric(), coverage = numeric())
  for (true_r in true_grid) {
    for (asm_r in assumed_grid) {
      set.seed(SEED_BASE + 200L + 100L * which(true_grid == true_r) +
                       which(assumed_grid == asm_r))
      lm    <- lambda_marg(p_Y)
      lc_a  <- lambda_cond_gauss(p_Y, asm_r)

      rel <- numeric(R); cov <- logical(R)
      for (r in 1:R) {
        rep         <- sim_one_no_collider(n, p_Y, true_r, gamma)
        beta_score  <- rep$beta_C / lc_a
        se_score    <- rep$se_C   / lc_a
        rel[r]      <- beta_score / gamma - 1
        cov[r]      <- abs(beta_score - gamma) < 1.96 * se_score
      }
      out <- rbind(out, data.frame(true_rho2 = true_r,
                                                              assumed_rho2 = asm_r,
                                                              rel_bias = mean(rel),
                                                              coverage = mean(cov) * 100))
    }
  }
  cat("\n", header, "\n", sep = "")
  .print_rounded(out, 3)
  invisible(out)
}

#' Larger-n version (Supp Table S8)
table_rho_sens_n20k <- function(R = 500, n = 20000, gamma = 0.05, p_Y = 0.10) {
  table_rho_sens(R = R, n = n, gamma = gamma, p_Y = p_Y,
                              header = sprintf("=== Supp Table S8: sensitivity at n=%d ===",
                                                                  as.integer(n)))
}

# ============================================================================= 
# IVW MR helpers (univariable and multivariable)
# ============================================================================= 

ivw_mr <- function(beta_X, se_X, beta_Y, se_Y) {
  w     <- 1 / se_Y^2
  num   <- sum(w * beta_X * beta_Y)
  den   <- sum(w * beta_X^2)
  theta <- num / den
  se    <- sqrt(1 / den)
  list(theta = theta, se = se)
}

mvmr_ivw <- function(B_X, B_Y, se_Y) {
  W   <- diag(1 / se_Y^2)
  inv <- solve(t(B_X) %*% W %*% B_X)
  list(theta = as.numeric(inv %*% t(B_X) %*% W %*% B_Y),
              cov   = inv)
}

# ============================================================================= 
# Helpers: simulate one IVW MR replicate for the T2D-CAD plausibility check
# ============================================================================= 

sim_t2d_cad_one <- function(K = 120, n_T2D = 80000, n_BMI = 80000,
                                                      n_CAD = 100000,
                                                      p_T2D = 0.10, p_CAD = 0.06,
                                                      rho2 = 0.30,
                                                      theta = 0.40,
                                                      b_lin = 0.50,
                                                      gamma_sd = 0.025,
                                                      beta_GH_sd = 0.02,
                                                      known_b = FALSE) {
  # True per-instrument liability-scale T2D effects
  gamma_k     <- rnorm(K, 0, gamma_sd)
  beta_GH_k   <- rnorm(K, 0, beta_GH_sd)

  # T2D liability: gamma_k + b_lin * beta_GH_k (collider contamination on liability scale)
  lc_T2D      <- lambda_cond_gauss(p_T2D, rho2)
  lm_T2D      <- lambda_marg(p_T2D)
  lc_CAD      <- lambda_marg(p_CAD)         # CAD unadjusted

  # Adjusted log-OR coefficient: lambda_cond * (gamma + b * beta_GH)
  beta_C_T2D  <- lc_T2D * (gamma_k + b_lin * beta_GH_k)
  se_C_T2D    <- 1 / sqrt(n_T2D * p_T2D * (1 - p_T2D))
  beta_C_T2D  <- beta_C_T2D + rnorm(K, 0, se_C_T2D)

  # SNP-BMI association with its own SE
  beta_GH_obs <- beta_GH_k + rnorm(K, 0, 1 / sqrt(n_BMI))
  se_GH       <- rep(1 / sqrt(n_BMI), K)

  # CAD outcome: log-OR scale, theta * gamma_k for liability-on-liability
  beta_GY_CAD <- lc_CAD * theta * gamma_k + rnorm(K, 0, 1 / sqrt(n_CAD * p_CAD * (1 - p_CAD)))
  se_GY_CAD   <- rep(1 / sqrt(n_CAD * p_CAD * (1 - p_CAD)), K)
  # Convert outcome to liability scale
  beta_Y_liab <- beta_GY_CAD / lc_CAD
  se_Y_liab   <- se_GY_CAD   / lc_CAD

  # Bias slope: known true b_logOR vs estimated via ordinary regression
  b_logOR_true <- lc_T2D * b_lin
  if (known_b) {
    b_hat <- b_logOR_true
  } else {
    b_hat <- coef(lm(beta_C_T2D ~ beta_GH_obs))[2]   # ordinary regression fit
  }

  A           <- beta_C_T2D - b_hat * beta_GH_obs

  # Three exposure rescalings
  beta_X_marg    <- A / lm_T2D
  se_X_marg      <- se_C_T2D / lm_T2D
  beta_X_score   <- A / lc_T2D
  se_X_score     <- se_C_T2D / lc_T2D
  # Uncorrected log-OR: exposure stays on log-OR scale (no rescaling), outcome
  # converted to liability via marginal CAD scaling. This is the "naive" pipeline
  # that an inattentive user would follow: rescale outcome but forget to rescale
  # exposure. The IVW slope under this pipeline targets theta/lambda_cond_T2D
  # at first order, severely underestimating theta.
  beta_X_uncorr  <- A                # log-OR scale, no rescaling
  se_X_uncorr    <- se_C_T2D

  # IVW MR with each rescaling
  marg   <- ivw_mr(beta_X_marg,   se_X_marg,   beta_Y_liab, se_Y_liab)
  score  <- ivw_mr(beta_X_score,  se_X_score,  beta_Y_liab, se_Y_liab)
  uncorr <- ivw_mr(beta_X_uncorr, se_X_uncorr, beta_Y_liab, se_Y_liab)

  list(marg = marg, score = score, uncorr = uncorr)
}

run_t2d_cad <- function(R = 100, rho2_grid = c(0.20, 0.30, 0.40),
                                          theta = 0.40, K = 120,
                                          n_T2D = 80000, known_b = FALSE) {
  out <- data.frame()
  for (rho2 in rho2_grid) {
    set.seed(SEED_BASE + 600L + round(rho2 * 100))
    kap <- kappa_score(0.10, rho2)
    em <- es <- eu <- numeric(R)
    sm <- ss <- su <- numeric(R)
    for (r in 1:R) {
      rep <- sim_t2d_cad_one(K = K, n_T2D = n_T2D, p_T2D = 0.10,
                                                  p_CAD = 0.06, rho2 = rho2,
                                                  theta = theta, known_b = known_b)
      em[r] <- rep$marg$theta
      es[r] <- rep$score$theta
      eu[r] <- rep$uncorr$theta
      sm[r] <- rep$marg$se;   ss[r] <- rep$score$se;  su[r] <- rep$uncorr$se
    }
    cov_m <- mean(abs(em - theta) < 1.96 * sm) * 100
    cov_s <- mean(abs(es - theta) < 1.96 * ss) * 100
    cov_u <- mean(abs(eu - theta) < 1.96 * su) * 100
    out <- rbind(out,
      data.frame(rho2 = rho2, kappa = kap, method = "Marginal-scaled",
                              theta_hat = mean(em), SE = mean(sm),
                              rel_bias = mean(em) / theta - 1, coverage = cov_m),
      data.frame(rho2 = rho2, kappa = kap, method = "Score-based",
                              theta_hat = mean(es), SE = mean(ss),
                              rel_bias = mean(es) / theta - 1, coverage = cov_s),
      data.frame(rho2 = rho2, kappa = kap, method = "Uncorrected log-OR",
                              theta_hat = mean(eu), SE = mean(su),
                              rel_bias = mean(eu) / theta - 1, coverage = cov_u)
    )
  }
  out
}

# ============================================================================= 
# Table 8: BMI-T2D plausibility check (theta = 0.40)
# ============================================================================= 

table_t2d_cad <- function(R = 100) {
  out <- run_t2d_cad(R = R, theta = 0.40)
  cat("\n=== Main Table 8: T2D-CAD plausibility check, theta=0.40 ===\n")
  .print_rounded(out, 3)
  invisible(out)
}

# ============================================================================= 
# Supp Table S7: theta = 0.20 sensitivity
# ============================================================================= 

table_theta_sens <- function(R = 100) {
  out <- run_t2d_cad(R = R, theta = 0.20)
  cat("\n=== Supp Table S7: theta=0.20 sensitivity ===\n")
  .print_rounded(out, 3)
  invisible(out)
}

# ============================================================================= 
# Supp Table S13: residual decomposition (known b vs estimated b)
# ============================================================================= 

table_residual_decomp <- function(R = 200) {
  cat("\n=== Supp Table S13: residual decomposition ===\n")
  configs <- list(
    list(K = 120, n_T2D = 80000),
    list(K = 300, n_T2D = 200000)
  )
  out <- data.frame()
  for (cfg in configs) {
    for (kb in c(TRUE, FALSE)) {
      sub <- run_t2d_cad(R = R, rho2_grid = 0.30, theta = 0.40,
                                          K = cfg$K, n_T2D = cfg$n_T2D, known_b = kb)
      score_row <- sub[sub$method == "Score-based", ]
      out <- rbind(out, data.frame(
        K = cfg$K, n_T2D = cfg$n_T2D,
        b_lin_source = ifelse(kb, "Known true b", "Estimated"),
        theta_hat = score_row$theta_hat,
        rel_bias  = score_row$rel_bias
      ))
    }
  }
  .print_rounded(out, 3)
  invisible(out)
}
#  =============================================================================
#  cbcGWAS: Replication code, Part 2
#  MVMR simulations, supplement tables, and master driver
#  =============================================================================

# ============================================================================= 
# Helpers: simulate a two-exposure MVMR replicate
# ============================================================================= 

#' One MVMR replicate. Returns naive, marginal-scaled, and score-based estimates
#' of the joint causal effects under three rescalings.
sim_mvmr_one <- function(K = 30, n_X = 8000,
                                                  p1 = 0.10, rho2_1 = 0.30,
                                                  p2 = 0.10, rho2_2 = 0.30,
                                                  theta = c(0.4, -0.3),
                                                  gamma_sd = 0.06) {       # match manuscript Sec 3.5
  # Per-instrument liability-scale exposure effects
  gamma_1 <- rnorm(K, 0, gamma_sd)
  gamma_2 <- rnorm(K, 0, gamma_sd)

  # Conversion factors per exposure
  lc1 <- lambda_cond_gauss(p1, rho2_1)
  lc2 <- lambda_cond_gauss(p2, rho2_2)
  lm1 <- lambda_marg(p1)
  lm2 <- lambda_marg(p2)

  # Adjusted log-OR coefficients (collider bias already absorbed; we use A_ij directly)
  se_X1 <- 1 / sqrt(n_X * p1 * (1 - p1))
  se_X2 <- 1 / sqrt(n_X * p2 * (1 - p2))
  A1    <- lc1 * gamma_1 + rnorm(K, 0, se_X1)
  A2    <- lc2 * gamma_2 + rnorm(K, 0, se_X2)

  # Outcome on liability scale: continuous
  beta_Y <- theta[1] * gamma_1 + theta[2] * gamma_2 + rnorm(K, 0, 0.005)
  se_Y   <- rep(0.005, K)

  # Three rescalings of the per-SNP exposure vectors
  X_naive <- cbind(A1,         A2)
  X_marg  <- cbind(A1 / lm1,    A2 / lm2)
  X_score <- cbind(A1 / lc1,    A2 / lc2)

  # Per-SNP exposure SEs (rescaled correspondingly)
  SE_naive <- cbind(rep(se_X1, K),       rep(se_X2, K))
  SE_marg  <- cbind(rep(se_X1 / lm1, K),  rep(se_X2 / lm2, K))
  SE_score <- cbind(rep(se_X1 / lc1, K),  rep(se_X2 / lc2, K))

  list(
    naive = mvmr_ivw(X_naive, beta_Y, se_Y),
    marg  = mvmr_ivw(X_marg,  beta_Y, se_Y),
    score = mvmr_ivw(X_score, beta_Y, se_Y),
    # Raw data: needed by external MVMR estimators (MVMRcML, MR-Egger, etc.)
    raw   = list(
      X_naive = X_naive, SE_naive = SE_naive,
      X_marg  = X_marg,  SE_marg  = SE_marg,
      X_score = X_score, SE_score = SE_score,
      beta_Y  = beta_Y,  se_Y    = se_Y,
      n_X     = n_X
    )
  )
}

# ============================================================================= 
# Table 9: MVMR with two binary exposures
# ============================================================================= 

table_mvmr <- function(R = 500, K = 30, n_X = 8000, theta = c(0.4, -0.3)) {
  scenarios <- list(
    list(p1 = 0.10, r1 = 0.30, p2 = 0.10, r2 = 0.30, label = "S1: equal kappa"),
    list(p1 = 0.10, r1 = 0.30, p2 = 0.05, r2 = 0.10, label = "S2: unequal kappa")
  )
  out <- data.frame()
  for (s in scenarios) {
    set.seed(SEED_BASE + 800L + which(sapply(scenarios, function(x) x$label) == s$label))
    k1 <- kappa_score(s$p1, s$r1); k2 <- kappa_score(s$p2, s$r2)
    naive <- marg <- score <- matrix(NA_real_, nrow = R, ncol = 2)
    for (r in 1:R) {
      rep <- sim_mvmr_one(K = K, n_X = n_X,
                                              p1 = s$p1, rho2_1 = s$r1,
                                              p2 = s$p2, rho2_2 = s$r2,
                                              theta = theta)
      naive[r, ] <- rep$naive$theta
      marg[r,  ] <- rep$marg$theta
      score[r, ] <- rep$score$theta
    }
    out <- rbind(out,
      data.frame(scenario = s$label, kappa1 = k1, kappa2 = k2,
                              method = "Naive (log-OR)",
                              theta1 = mean(naive[, 1]), theta1_mcse = sd(naive[, 1]) / sqrt(R),
                              theta2 = mean(naive[, 2]), theta2_mcse = sd(naive[, 2]) / sqrt(R)),
      data.frame(scenario = s$label, kappa1 = k1, kappa2 = k2,
                              method = "Marginal",
                              theta1 = mean(marg[, 1]),  theta1_mcse = sd(marg[, 1])  / sqrt(R),
                              theta2 = mean(marg[, 2]),  theta2_mcse = sd(marg[, 2])  / sqrt(R)),
      data.frame(scenario = s$label, kappa1 = k1, kappa2 = k2,
                              method = "Score-based",
                              theta1 = mean(score[, 1]), theta1_mcse = sd(score[, 1]) / sqrt(R),
                              theta2 = mean(score[, 2]), theta2_mcse = sd(score[, 2]) / sqrt(R))
    )
  }
  cat("\n=== Main Table 9: MVMR with two binary exposures ===\n")
  .print_rounded(out, 3)
  invisible(out)
}

# ============================================================================= 
# Supp Table S10: MVMR weak-instrument shrinkage at increased K, n
# ============================================================================= 

table_mvmr_strength <- function() {
  configs <- list(
    list(K = 30,  n_X =  8000),
    list(K = 60,  n_X = 20000),
    list(K = 200, n_X = 50000)
  )
  out <- data.frame()
  for (cfg in configs) {
    set.seed(SEED_BASE + 900L + cfg$K)
    R <- 500
    score_t1 <- score_t2 <- numeric(R)
    marg_t1  <- marg_t2  <- numeric(R)
    for (r in 1:R) {
      rep <- sim_mvmr_one(K = cfg$K, n_X = cfg$n_X,
                                              p1 = 0.10, rho2_1 = 0.30,
                                              p2 = 0.10, rho2_2 = 0.30)
      score_t1[r] <- rep$score$theta[1]; score_t2[r] <- rep$score$theta[2]
      marg_t1[r]  <- rep$marg$theta[1];  marg_t2[r]  <- rep$marg$theta[2]
    }
    out <- rbind(out, data.frame(
      K = cfg$K, n_X = cfg$n_X,
      marg_theta1  = mean(marg_t1),  marg_theta2  = mean(marg_t2),
      score_theta1 = mean(score_t1), score_theta2 = mean(score_t2),
      score_resid_pct1 = (mean(score_t1) / 0.4 - 1) * 100,
      score_resid_pct2 = (mean(score_t2) / -0.3 - 1) * 100
    ))
  }
  cat("\n=== Supp Table S10: MVMR weak-instrument shrinkage ===\n")
  .print_rounded(out, 3)
  invisible(out)
}

# ============================================================================= 
# Supp Table S2: probit-logit gap and non-Gaussian-M sensitivity
# ============================================================================= 

#' Solve population score equation exactly via deterministic quadrature.
#' Returns the exact lambda_cond at given p_Y, rho^2, and a custom M distribution.
#' For Gaussian M, this should reproduce the closed form.
lambda_cond_exact <- function(p_Y, rho2, M_density,
                                                          M_grid_lo = -8, M_grid_hi = 8, n_grid = 4000) {
  t_Y <- qnorm(1 - p_Y)
  tau <- sqrt(1 - rho2)
  M   <- seq(M_grid_lo * sqrt(rho2), M_grid_hi * sqrt(rho2), length.out = n_grid)
  dM  <- diff(M)[1]
  fM  <- M_density(M)
  fM  <- fM / sum(fM * dM)        # normalize density
  pH  <- pnorm((M - t_Y) / tau)
  num <- sum(fM * dnorm((M - t_Y) / tau) / tau) * dM
  den <- sum(fM * pH * (1 - pH)) * dM
  num / den
}

table_nongaussian <- function() {
  cells <- list(c(0.10, 0.30), c(0.10, 0.40), c(0.05, 0.40), c(0.05, 0.50))
  rows  <- list()

  for (cell in cells) {
    p_Y <- cell[1]; rho2 <- cell[2]
    closed <- lambda_cond_gauss(p_Y, rho2)

    # Skew-normal: use moments to scale
    sn <- function(x, alpha) {
      delta <- alpha / sqrt(1 + alpha^2)
      omega <- sqrt(rho2 / (1 - 2 * delta^2 / pi))
      xi    <- -omega * delta * sqrt(2 / pi)
      2 * dnorm((x - xi) / omega) * pnorm(alpha * (x - xi) / omega) / omega
    }
    sn05  <- lambda_cond_exact(p_Y, rho2,
                                                            function(x) sn(x, alpha = 4))
    sn10  <- lambda_cond_exact(p_Y, rho2,
                                                            function(x) sn(x, alpha = 8))

    # Zero-inflated Gaussian
    zi <- function(x, p0) {
      sd_used <- sqrt(rho2 / (1 - p0))
      (1 - p0) * dnorm(x, 0, sd_used) + p0 * (abs(x) < 1e-3) / 2e-3
    }
    zi10 <- lambda_cond_exact(p_Y, rho2, function(x) zi(x, 0.10))
    zi30 <- lambda_cond_exact(p_Y, rho2, function(x) zi(x, 0.30))

    # Binary
    b50  <- lambda_cond_discrete(p_Y, rho2, 0.50)
    b30  <- lambda_cond_discrete(p_Y, rho2, 0.30)
    b10  <- lambda_cond_discrete(p_Y, rho2, 0.10)

    rows[[paste(cell, collapse = ",")]] <- c(
      gauss      = 0,
      sn_skew05  = abs(closed - sn05) / sn05 * 100,
      sn_skew10  = abs(closed - sn10) / sn10 * 100,
      zi10       = abs(closed - zi10) / zi10 * 100,
      zi30       = abs(closed - zi30) / zi30 * 100,
      bin50      = abs(closed - b50)  / b50  * 100,
      bin30      = abs(closed - b30)  / b30  * 100,
      bin10      = abs(closed - b10)  / b10  * 100
    )
  }
  out <- do.call(cbind, rows)
  cat("\n=== Supp Table S2: closed-form vs exact, percent discrepancy ===\n")
  .print_rounded(out, 2)
  invisible(out)
}

# ============================================================================= 
# Supp Table S1: probit-logit gap (closed-form vs exact working-model variance)
# ============================================================================= 

#' Compute lambda_cond using the exact logistic working-model variance,
#' i.e. solving the population score equations for (alpha^o, delta^o) at zeroth order
#' and using sigma(...) (1-sigma(...)) instead of p(1-p).
lambda_cond_exact_logit <- function(p_Y, rho2, n_grid = 2000) {
  t_Y <- qnorm(1 - p_Y)
  tau <- sqrt(1 - rho2)
  M   <- seq(-6, 6, length.out = n_grid) * sqrt(rho2)
  dM  <- diff(M)[1]
  fM  <- dnorm(M, 0, sqrt(rho2))
  pH  <- pnorm((M - t_Y) / tau)

  # Find (alpha, delta) minimizing weighted log-loss between probit and logistic.
  # At zeroth order, the working logistic fits expit(alpha + delta * M) to pH.
  obj <- function(par) {
    a <- par[1]; d <- par[2]
    q <- plogis(a + d * M)
    -sum(fM * (pH * log(q + 1e-12) + (1 - pH) * log(1 - q + 1e-12))) * dM
  }
  init <- c(qlogis(p_Y), 1)
  fit  <- optim(init, obj, method = "Nelder-Mead", control = list(reltol = 1e-10))
  a_o  <- fit$par[1]; d_o  <- fit$par[2]
  q    <- plogis(a_o + d_o * M)
  num  <- sum(fM * dnorm((M - t_Y) / tau) / tau) * dM
  den  <- sum(fM * q * (1 - q)) * dM
  num / den
}

# =============================================================================
# Main Table 5: Liability-scale per-SNP heritability under marginal vs
# score-based correction
# =============================================================================
# Verifies that summing squared per-SNP liability-scale effects from a
# covariate-adjusted GWAS gives the right h^2_liab only when score-based
# rescaling is applied. Marginal scaling biases h^2_liab upward by kappa^2.
#
# Setup: K SNPs with true gamma_k ~ N(0, h2_true / K). Per-SNP covariate-
# adjusted log-OR is lc * gamma_k + epsilon_k, where epsilon_k has variance
# 1 / (n * p_Y * (1 - p_Y)). Convert each beta_logOR to liability via
# (a) marginal: beta_liab = beta_logOR / lambda_marg, and
# (b) score:    beta_liab = beta_logOR / lambda_cond.
# Sum squared per-SNP effects to get h2_marg and h2_score.

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
      gamma_k     <- rnorm(K, 0, sqrt(h2_true / K))
      beta_logOR  <- lc * gamma_k + rnorm(K, 0, se_C)
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
  cat("\n=== Main Table 5: Liability-scale heritability under marginal vs score correction ===\n")
  cat("(true h^2 =", h2_true, ", p_Y =", p_Y, ", K =", K, ", n =", n, ", R =", R, ")\n")
  print(out, row.names = FALSE, digits = 4)
  invisible(out)
}

# =============================================================================
# Main Table 6: Polygenic-score scaling under marginal vs score-based
# correction
# =============================================================================
# Verifies that a polygenic score built from covariate-adjusted GWAS effects
# regresses on the true liability with a slope of 1 only under score-based
# rescaling. Under marginal scaling, the regression slope is biased by
# kappa relative to the score-based version.
#
# Setup: K SNPs in a discovery GWAS produce per-SNP covariate-adjusted log-OR
# coefficients beta_logOR_k (as in table_heritability above). In a target
# sample of size n_target, draw genotypes G_target_k ~ N(0,1), construct
# the true liability Y_star = sum_k gamma_k * G_target_k + residual, build
# PGS_marg = sum_k (beta_logOR_k / lm) * G_target_k and PGS_score similarly,
# regress Y_star on each PGS, report the slope and R^2.

table_pgs <- function(R = 200, K = 2000, n_disc = 100000, n_target = 50000,
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
      # Discovery: estimate per-SNP covariate-adjusted log-OR coefficients
      gamma_k     <- rnorm(K, 0, sqrt(h2_true / K))
      beta_logOR  <- lc * gamma_k + rnorm(K, 0, se_C)

      # Target sample: build PGS, regress true liability on each variant
      G_target    <- matrix(rnorm(n_target * K), n_target, K)
      Y_star      <- as.vector(G_target %*% gamma_k) +
                       rnorm(n_target, 0, sqrt(1 - h2_true))

      pgs_marg    <- as.vector(G_target %*% (beta_logOR / lm))
      pgs_score   <- as.vector(G_target %*% (beta_logOR / lc))

      fit_marg    <- lm(Y_star ~ pgs_marg)
      fit_score   <- lm(Y_star ~ pgs_score)
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
  cat("\n=== Main Table 6: Polygenic-score scaling under marginal vs score correction ===\n")
  cat("(true h^2 =", h2_true, ", p_Y =", p_Y, ", K =", K,
      ", n_disc =", n_disc, ", n_target =", n_target, ", R =", R, ")\n")
  print(out, row.names = FALSE, digits = 4)
  invisible(out)
}


table_pl_gap <- function() {
  p_Y_grid  <- c(0.50, 0.20, 0.10, 0.05)
  rho2_grid <- c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
  out <- matrix(NA_real_, length(rho2_grid), length(p_Y_grid))
  for (i in seq_along(rho2_grid)) {
    for (j in seq_along(p_Y_grid)) {
      closed <- lambda_cond_gauss(p_Y_grid[j], rho2_grid[i])
      exact  <- lambda_cond_exact_logit(p_Y_grid[j], rho2_grid[i])
      out[i, j] <- abs(closed - exact) / exact * 100
    }
  }
  rownames(out) <- paste0("rho2=", rho2_grid)
  colnames(out) <- paste0("p_Y=",  p_Y_grid)
  cat("\n=== Supp Table S1: probit-logit gap (percent) ===\n")
  .print_rounded(out, 2)
  invisible(out)
}

# ============================================================================= 
# Supp Table S3: Prentice-Pyke invariance under case-control sampling
# ============================================================================= 

table_cc <- function(R = 800, n = 20000, p_Y = 0.10, rho2 = 0.30, gamma = 0.05) {
  pi_grid <- c(0.05, 0.10, 0.30, 0.50)
  out <- data.frame()
  for (pi_frac in pi_grid) {
    set.seed(SEED_BASE + 300L + round(pi_frac * 100))
    estimates <- numeric(R)
    intercepts <- numeric(R)
    # We pre-simulate a large population, then sample with target case fraction.
    n_pop <- max(n, round(n * pi_frac / p_Y) + n)  # ensure enough cases available
    for (r in 1:R) {
      H <- rnorm(n_pop, 0, sqrt(rho2))
      G <- rnorm(n_pop)
      Y_star <- H + gamma * G + rnorm(n_pop, 0, sqrt(1 - rho2 - gamma^2))
      Y      <- as.numeric(Y_star > qnorm(1 - p_Y))
      cases   <- which(Y == 1)
      ctrls   <- which(Y == 0)
      n_case_target <- round(n * pi_frac)
      n_ctrl_target <- n - n_case_target
      if (length(cases) < n_case_target || length(ctrls) < n_ctrl_target) {
        estimates[r]  <- NA
        intercepts[r] <- NA
        next
      }
      kept <- c(sample(cases, n_case_target), sample(ctrls, n_ctrl_target))
      fit  <- glm(Y[kept] ~ G[kept] + H[kept], family = binomial())
      estimates[r]  <- coef(fit)[2]
      intercepts[r] <- coef(fit)[1]
    }
    out <- rbind(out, data.frame(
      pi = pi_frac, pi_over_pY = pi_frac / p_Y,
      beta_G_mean = mean(estimates, na.rm = TRUE),
      SE          = sd(estimates,   na.rm = TRUE) / sqrt(R),
      alpha       = mean(intercepts, na.rm = TRUE),
      theory      = lambda_cond_gauss(p_Y, rho2) * gamma
    ))
  }
  cat("\n=== Supp Table S3: Prentice-Pyke invariance ===\n")
  .print_rounded(out, 3)
  invisible(out)
}

# ============================================================================= 
# Supp Table S4: bias-slope recovery across estimators
# ============================================================================= 

#' Simplified bias-slope recovery test: fit log-OR GWAS regressed on SNP-covariate
#' associations and recover the bias slope on the log-OR scale.
table_bias_slope_recovery <- function(R = 200, K = 200, n = 50000,
                                                                                  b_lin = 0.50, gamma_sd = 0.025,
                                                                                  beta_GH_sd = 0.05) {
  cells <- list(c(0.10, 0.30), c(0.20, 0.30), c(0.10, 0.40))
  out <- data.frame()
  for (cell in cells) {
    p_Y <- cell[1]; rho2 <- cell[2]
    set.seed(SEED_BASE + 400L + round(p_Y * 100) + round(rho2 * 100))
    lc <- lambda_cond_gauss(p_Y, rho2)
    b_logOR_true <- lc * b_lin

    # Three estimators (all reduce to ordinary-regression slope at first order)
    est_OLS  <- numeric(R)
    est_SH   <- numeric(R)        # Slope-Hunter (we approximate by trimmed mean)
    est_cML  <- numeric(R)        # MVMR-cML (we approximate by IV regression)
    for (r in 1:R) {
      gamma_k     <- rnorm(K, 0, gamma_sd)
      beta_GH_k   <- rnorm(K, 0, beta_GH_sd)
      se_C        <- 1 / sqrt(n * p_Y * (1 - p_Y))
      beta_C      <- lc * (gamma_k + b_lin * beta_GH_k) + rnorm(K, 0, se_C)
      beta_GH_obs <- beta_GH_k + rnorm(K, 0, 1 / sqrt(n))

      est_OLS[r] <- coef(lm(beta_C ~ beta_GH_obs))[2]

      # Slope-Hunter approximation: cluster on near-zero gamma_k SNPs (those with
      # |adj-gwas - lc * b_lin * beta_GH_obs| small) and re-fit
      nh_idx <- order(abs(beta_C - est_OLS[r] * beta_GH_obs))[1:max(20, K %/% 4)]
      est_SH[r] <- coef(lm(beta_C[nh_idx] ~ beta_GH_obs[nh_idx]))[2]

      # MVMR-cML approximation: weighted regression with constraint
      w <- 1 / se_C^2
      est_cML[r] <- sum(w * beta_C * beta_GH_obs) / sum(w * beta_GH_obs^2)
    }
    out <- rbind(out,
      data.frame(p_Y = p_Y, rho2 = rho2, lambda_cond = lc, true_b_logOR = b_logOR_true,
                              estimator = "Ordinary regression",
                              b_hat = mean(est_OLS), mcse = sd(est_OLS) / sqrt(R)),
      data.frame(p_Y = p_Y, rho2 = rho2, lambda_cond = lc, true_b_logOR = b_logOR_true,
                              estimator = "Slope-Hunter (approx)",
                              b_hat = mean(est_SH),  mcse = sd(est_SH)  / sqrt(R)),
      data.frame(p_Y = p_Y, rho2 = rho2, lambda_cond = lc, true_b_logOR = b_logOR_true,
                              estimator = "MVMR-cML (approx)",
                              b_hat = mean(est_cML), mcse = sd(est_cML) / sqrt(R))
    )
  }
  cat("\n=== Supp Table S4: bias-slope recovery across estimators ===\n")
  .print_rounded(out, 3)
  invisible(out)
}

# ============================================================================= 
# Supp Table S5: bias-slope robustness to sigma_beta
# ============================================================================= 

table_sigma_beta <- function(R = 800, K = 200, n = 50000, b_lin = 0.50,
                                                              p_Y = 0.10, rho2 = 0.30) {
  sigma_beta_grid <- c(0.04, 0.06, 0.08, 0.12, 0.20)
  lc <- lambda_cond_gauss(p_Y, rho2)
  b_logOR_true <- lc * b_lin
  out <- data.frame()
  for (sb in sigma_beta_grid) {
    set.seed(SEED_BASE + 500L + round(sb * 1000))
    est <- numeric(R)
    for (r in 1:R) {
      gamma_k     <- rnorm(K, 0, 0.025)
      beta_GH_k   <- rnorm(K, 0, sb)
      se_C        <- 1 / sqrt(n * p_Y * (1 - p_Y))
      beta_C      <- lc * (gamma_k + b_lin * beta_GH_k) + rnorm(K, 0, se_C)
      beta_GH_obs <- beta_GH_k + rnorm(K, 0, 1 / sqrt(n))
      est[r]      <- coef(lm(beta_C ~ beta_GH_obs))[2]
    }
    out <- rbind(out, data.frame(
      sigma_beta = sb,
      b_hat      = mean(est),
      mcse       = sd(est) / sqrt(R),
      rel_bias   = (mean(est) - b_logOR_true) / b_logOR_true
    ))
  }
  cat("\n=== Supp Table S5: bias-slope robustness to sigma_beta ===\n")
  .print_rounded(out, 4)
  invisible(out)
}

# ============================================================================= 
# Supp Table S6: multi-covariate adjustment
# ============================================================================= 

table_multi_covar <- function(R = 800, n = 20000, p_Y = 0.10, gamma = 0.05,
                                                          total_rho2 = 0.30) {
  q_grid <- c(1, 3, 5)
  out <- data.frame()
  for (q in q_grid) {
    set.seed(SEED_BASE + 700L + q)
    rho2_each <- total_rho2 / q  # split equally
    lc <- lambda_cond_gauss(p_Y, total_rho2)
    lm <- lambda_marg(p_Y)

    rel_marg <- rel_score <- numeric(R)
    cov_marg <- cov_score <- logical(R)
    for (r in 1:R) {
      H <- mvtnorm::rmvnorm(n, sigma = diag(rho2_each, q))
      G <- rnorm(n)
      M <- rowSums(H)                # delta = 1 per dim
      res_var <- max(1 - total_rho2 - gamma^2, 1e-6)
      Y_star  <- M + gamma * G + rnorm(n, 0, sqrt(res_var))
      Y       <- as.numeric(Y_star > qnorm(1 - p_Y))
      df      <- data.frame(Y = Y, G = G, H = H)
      fit     <- glm(Y ~ ., data = df, family = binomial())
      bC      <- coef(fit)["G"]
      seC     <- summary(fit)$coefficients["G", "Std. Error"]

      bm     <- bC / lm; sm <- seC / lm
      bs     <- bC / lc; ss <- seC / lc
      rel_marg[r]  <- bm / gamma - 1
      rel_score[r] <- bs / gamma - 1
      cov_marg[r]  <- abs(bm - gamma) < 1.96 * sm
      cov_score[r] <- abs(bs - gamma) < 1.96 * ss
    }
    out <- rbind(out, data.frame(
      q = q,
      score_bias  = mean(rel_score),
      marg_bias   = mean(rel_marg),
      score_cov   = mean(cov_score) * 100,
      marg_cov    = mean(cov_marg)  * 100
    ))
  }
  cat("\n=== Supp Table S6: multi-covariate adjustment ===\n")
  .print_rounded(out, 3)
  invisible(out)
}

# ============================================================================= 
# Supp Table S9: MVMR under non-IVW estimators
# ============================================================================= 

#' Wrapper around MVMRcML::MVmr_cML_DP from the published MVMR-cML package
#' (Lin et al. 2023, AJHG 110:592-605). MVMRcML is GitHub-only:
#'   remotes::install_github("ZhaotongL/MVMR-cML")
#'
#' API note: MVMRcML's MVmr_cML_DP requires:
#'   - b_exp:     K x p matrix of exposure effects
#'   - b_out:     K x 1 column matrix (NOT a vector!) of outcome effects
#'   - se_bx:     K x p matrix of exposure SEs
#'   - Sig_inv_l: list of K (p+1) x (p+1) inverse covariance matrices,
#'                built via MVMRcML::invcov_mvmr(se_bx, se_by, rho_mat).
#'                For two-sample MVMR with independent samples, rho_mat is
#'                the identity matrix of dimension (p+1).
#'   - K_vec:     range of candidate invalid-IV counts. Must be small relative
#'                to K (the README uses 0:20 for K=132). Setting K_vec too
#'                large (e.g., 0:50 for K=80) drives the optimizer into
#'                degenerate matrix algebra and produces cryptic Armadillo
#'                errors.
#'
#' Argument names for invcov_mvmr vary across package versions
#' (se_bx/se_by vs se_exp/se_out), so we introspect formals().
#'
#' Set verbose = TRUE for one-shot diagnostics.
.mvmr_cml_one <- function(X_mat, SE_mat, beta_Y, se_Y, n_X, verbose = FALSE) {
  if (!requireNamespace("MVMRcML", quietly = TRUE))
    return(rep(NA_real_, ncol(X_mat)))

  pp1 <- ncol(SE_mat) + 1L

  # --- Build SigmaInv via invcov_mvmr (introspect arg names by version) ---
  inv_args <- names(formals(MVMRcML::invcov_mvmr))
  call_args <- list()
  if      ("se_bx"  %in% inv_args) call_args$se_bx  <- SE_mat
  else if ("se_exp" %in% inv_args) call_args$se_exp <- SE_mat
  if      ("se_by"  %in% inv_args) call_args$se_by  <- as.numeric(se_Y)
  else if ("se_out" %in% inv_args) call_args$se_out <- as.numeric(se_Y)
  # rho_mat (per-SNP correlation between exposure and outcome stats):
  # for independent-sample two-sample MVMR, this is the identity.
  if ("rho_mat" %in% inv_args) call_args$rho_mat <- diag(pp1)
  if (length(call_args) < 2)
    return(rep(NA_real_, ncol(X_mat)))

  Sigma_inv <- tryCatch(
    do.call(MVMRcML::invcov_mvmr, call_args),
    error = function(e) {
      if (verbose) cat("    invcov_mvmr error: ", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(Sigma_inv)) return(rep(NA_real_, ncol(X_mat)))

  # Reshape to list of (p+1) x (p+1) matrices if needed
  if (is.list(Sigma_inv) && length(Sigma_inv) > 0 &&
            !is.matrix(Sigma_inv[[1]])) {
    Sigma_inv <- lapply(Sigma_inv, function(v) matrix(v, pp1, pp1))
  } else if (is.array(Sigma_inv) && length(dim(Sigma_inv)) == 3) {
    d <- dim(Sigma_inv)
    Sigma_inv <- if (d[3] == nrow(SE_mat))
      lapply(seq_len(d[3]), function(i) Sigma_inv[, , i])
    else
      lapply(seq_len(d[1]), function(i) Sigma_inv[i, , ])
  }
  # Strip class/dimnames per Rcpp coercion requirements
  Sigma_inv <- lapply(Sigma_inv, function(m) {
    m <- unclass(m); dimnames(m) <- NULL
    if (!is.matrix(m)) m <- matrix(m, pp1, pp1)
    m
  })

  # --- Call MVmr_cML_DP exactly as in the package README example ---
  # Critical: b_out must be a column matrix (NOT a numeric vector); K_vec
  # must be small relative to K (the README uses 0:20 for K=132).
  X_clean  <- matrix(as.numeric(X_mat),  nrow = nrow(X_mat),  ncol = ncol(X_mat))
  SE_clean <- matrix(as.numeric(SE_mat), nrow = nrow(SE_mat), ncol = ncol(SE_mat))

  fit <- tryCatch(
    MVMRcML::MVmr_cML_DP(
      b_exp     = X_clean,
      b_out     = as.matrix(as.numeric(beta_Y)),
      se_bx     = SE_clean,
      Sig_inv_l = Sigma_inv,
      n         = n_X,
      num_pert  = 100,
      K_vec     = 0:as.integer(min(floor(nrow(X_mat) / 4), 20))
    ),
    error = function(e) {
      if (verbose) cat("    MVmr_cML_DP error: ", conditionMessage(e), "\n")
      else        message("  [MVMRcML] error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) return(rep(NA_real_, ncol(X_mat)))

  # The DP variant returns BIC_DP_theta as the canonical estimate (Lin et al.
  # AJHG 2023 supplement). Some versions name it differently; try several.
  out <- NULL
  for (fld in c("BIC_DP_theta", "theta", "BIC_theta", "theta_DP")) {
    if (!is.null(fit[[fld]])) { out <- fit[[fld]]; break }
  }
  if (is.null(out)) return(rep(NA_real_, ncol(X_mat)))
  as.numeric(out)
}

#' Wrapper around MendelianRandomization::mr_mvmedian (weighted-median MVMR).
#' MendelianRandomization returns S4 objects; we use slotNames-based access
#' which works whether the field is an S3 list element or an S4 slot.
.mvmr_median_one <- function(X_mat, SE_mat, beta_Y, se_Y) {
  obj <- tryCatch(
    MendelianRandomization::mr_mvinput(
      bx     = X_mat,
      bxse   = SE_mat,
      by     = as.numeric(beta_Y),
      byse   = as.numeric(se_Y)
    ),
    error = function(e) NULL
  )
  if (is.null(obj)) return(rep(NA_real_, ncol(X_mat)))
  fit <- tryCatch(MendelianRandomization::mr_mvmedian(obj),
                                  error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, ncol(X_mat)))
  est <- if (isVirtualClass(class(fit)) || isS4(fit)) tryCatch(fit@Estimate, error = function(e) NULL)
              else fit$Estimate
  if (is.null(est)) return(rep(NA_real_, ncol(X_mat)))
  as.numeric(est)
}

#' Wrapper around MendelianRandomization::mr_mvegger (multivariable MR-Egger).
#' Returns slope estimates and the pleiotropy intercept.
.mvmr_egger_one <- function(X_mat, SE_mat, beta_Y, se_Y) {
  obj <- tryCatch(
    MendelianRandomization::mr_mvinput(
      bx     = X_mat,
      bxse   = SE_mat,
      by     = as.numeric(beta_Y),
      byse   = as.numeric(se_Y)
    ),
    error = function(e) NULL
  )
  na_out <- list(theta = rep(NA_real_, ncol(X_mat)), intercept = NA_real_)
  if (is.null(obj)) return(na_out)

  fit <- tryCatch(MendelianRandomization::mr_mvegger(obj, orientate = 1),
                                  error = function(e) NULL)
  if (is.null(fit)) return(na_out)

  # MR-Egger return object is S4; use slot access.
  est <- if (isS4(fit)) tryCatch(fit@Estimate, error = function(e) NULL)
              else fit$Estimate
  ix  <- if (isS4(fit)) tryCatch(fit@Intercept, error = function(e) NULL)
              else fit$Intercept

  if (is.null(est)) est <- rep(NA_real_, ncol(X_mat))
  if (is.null(ix))  ix  <- NA_real_
  list(theta = as.numeric(est),
              intercept = as.numeric(ix)[1])
}

# ============================================================================= 
# Supp Table S9: MVMR under non-IVW estimators
# ============================================================================= 
#
# Reproduces Table S9 by calling the published packages MVMRcML and 
# MendelianRandomization. Falls back to IVW-only output with installation
# instructions if those packages are missing.

table_mvmr_non_ivw <- function(R = 300, K = 80, n_X = 20000) {
  has_MR  <- requireNamespace("MendelianRandomization", quietly = TRUE)
  has_cML <- requireNamespace("MVMRcML",                 quietly = TRUE)

  if (!has_MR || !has_cML) {
    cat("\nNote: external MVMR packages not installed.\n")
    cat("      For full Supp Table S9 reproduction install:\n")
    cat("        install.packages('MendelianRandomization')\n")
    cat("        # MVMRcML available from https://github.com/ZhaotongL/MVMR-cML\n")
    cat("        # remotes::install_github('ZhaotongL/MVMR-cML')\n")
    cat("      Reporting IVW baseline only.\n")
  }

  set.seed(SEED_BASE + 1000L)
  theta <- c(0.4, -0.3)

  # Storage
  ivw_marg_t1   <- numeric(R)   # IVW under marginal scaling
  ivw_score_t1  <- numeric(R)   # IVW under score scaling
  cml_marg_t1   <- numeric(R)   # MVMRcML under marginal scaling
  cml_score_t1  <- numeric(R)   # MVMRcML under score scaling
  med_marg_t1   <- numeric(R)
  med_score_t1  <- numeric(R)
  egg_marg_t1   <- numeric(R)
  egg_score_t1  <- numeric(R)
  egg_marg_int  <- numeric(R)
  egg_score_int <- numeric(R)

  for (r in 1:R) {
    rep <- sim_mvmr_one(K = K, n_X = n_X,
                                            p1 = 0.10, rho2_1 = 0.30,
                                            p2 = 0.10, rho2_2 = 0.30,
                                            theta = theta)
    raw <- rep$raw

    # IVW (always available)
    ivw_marg_t1[r]   <- rep$marg$theta[1]
    ivw_score_t1[r]  <- rep$score$theta[1]

    if (has_cML) {
      cml_marg_t1[r]  <- .mvmr_cml_one(raw$X_marg,  raw$SE_marg,
                                                                          raw$beta_Y, raw$se_Y, raw$n_X)[1]
      cml_score_t1[r] <- .mvmr_cml_one(raw$X_score, raw$SE_score,
                                                                          raw$beta_Y, raw$se_Y, raw$n_X)[1]
    } else {
      cml_marg_t1[r] <- cml_score_t1[r] <- NA
    }

    if (has_MR) {
      med_marg_t1[r]  <- .mvmr_median_one(raw$X_marg,  raw$SE_marg,
                                                                                  raw$beta_Y, raw$se_Y)[1]
      med_score_t1[r] <- .mvmr_median_one(raw$X_score, raw$SE_score,
                                                                                  raw$beta_Y, raw$se_Y)[1]
      egg_m <- .mvmr_egger_one(raw$X_marg,  raw$SE_marg,
                                                          raw$beta_Y, raw$se_Y)
      egg_s <- .mvmr_egger_one(raw$X_score, raw$SE_score,
                                                          raw$beta_Y, raw$se_Y)
      egg_marg_t1[r]   <- egg_m$theta[1];  egg_marg_int[r]  <- egg_m$intercept
      egg_score_t1[r]  <- egg_s$theta[1];  egg_score_int[r] <- egg_s$intercept
    } else {
      med_marg_t1[r] <- med_score_t1[r] <- NA
      egg_marg_t1[r] <- egg_score_t1[r] <- NA
      egg_marg_int[r] <- egg_score_int[r] <- NA
    }

    if (r %% 50 == 0) cat(sprintf("  ... replicate %d / %d\n", r, R))
  }

  out <- data.frame(
    estimator = c("IVW", "MVMR-cML", "MVMR-Median", "MVMR-Egger"),
    marg_rel_bias_pct = c(
      (mean(ivw_marg_t1) / theta[1] - 1) * 100,
      (mean(cml_marg_t1, na.rm = TRUE) / theta[1] - 1) * 100,
      (mean(med_marg_t1, na.rm = TRUE) / theta[1] - 1) * 100,
      (mean(egg_marg_t1, na.rm = TRUE) / theta[1] - 1) * 100
    ),
    score_rel_bias_pct = c(
      (mean(ivw_score_t1) / theta[1] - 1) * 100,
      (mean(cml_score_t1, na.rm = TRUE) / theta[1] - 1) * 100,
      (mean(med_score_t1, na.rm = TRUE) / theta[1] - 1) * 100,
      (mean(egg_score_t1, na.rm = TRUE) / theta[1] - 1) * 100
    ),
    intercept_est = c(NA, NA, NA, mean(egg_score_int, na.rm = TRUE))
  )

  cat("\n=== Supp Table S9: MVMR under non-IVW estimators ===\n")
  if (!has_MR || !has_cML) {
    cat("    (rows for missing packages will show NaN)\n")
  }
  .print_rounded(out, 3)
  invisible(out)
}

# ============================================================================= 
# Master driver
# ============================================================================= 

run_quick <- function() {
  # Closed-form-only tables, < 1 minute
  table_kappa()
  table_ad_apoe()
  table_pl_gap()
  table_nongaussian()
}

#' Debug helper: verify that MVMRcML and MendelianRandomization are correctly
#' installed and the wrapper functions return finite values on a single replicate.
#' Run once before table_mvmr_non_ivw() to surface installation issues.
#' @export
diagnose_mvmr_packages <- function() {
  cat("=== MVMR external-package diagnostic ===\n")
  has_MR  <- requireNamespace("MendelianRandomization", quietly = TRUE)
  has_cML <- requireNamespace("MVMRcML",                 quietly = TRUE)
  cat("MendelianRandomization installed:", has_MR, "\n")
  cat("MVMRcML installed:               ", has_cML, "\n")
  if (!has_MR || !has_cML) {
    cat("\nInstall the missing package(s) before running table_mvmr_non_ivw().\n")
    return(invisible(NULL))
  }

  # Check exported function names
  cat("\nMVMRcML exports:\n")
  exports <- getNamespaceExports("MVMRcML")
  cat("  ", paste(exports, collapse = ", "), "\n")

  # Run a single replicate
  cat("\nRunning a single MVMR simulation replicate...\n")
  set.seed(SEED_BASE + 1000L)
  rep <- sim_mvmr_one(K = 80, n_X = 20000,
                                          p1 = 0.10, rho2_1 = 0.30,
                                          p2 = 0.10, rho2_2 = 0.30,
                                          theta = c(0.4, -0.3))
  raw <- rep$raw

  cat("\n[IVW under score scaling] theta1 =",
          round(rep$score$theta[1], 4), "\n")

  cat("\n[MVMR-cML score]\n")
  cml <- tryCatch(.mvmr_cml_one(raw$X_score, raw$SE_score,
                                                                      raw$beta_Y, raw$se_Y, raw$n_X,
                                                                      verbose = TRUE),
                                  error = function(e) {
                                    cat("  ERROR:", conditionMessage(e), "\n"); NA
                                  })
  cat("  theta_hat:", round(cml, 4), "\n")
  if (any(is.na(cml))) {
    cat("  --> MVMRcML returned NA. The verbose output above shows the\n")
    cat("      introspected function signatures and any error message.\n")
    cat("      Edit the .mvmr_cml_one wrapper to match if needed.\n")
  }

  cat("\n[MVMR-Median score]\n")
  med <- .mvmr_median_one(raw$X_score, raw$SE_score,
                                                  raw$beta_Y, raw$se_Y)
  cat("  theta_hat:", round(med, 4), "\n")

  cat("\n[MVMR-Egger score]\n")
  egg <- .mvmr_egger_one(raw$X_score, raw$SE_score,
                                                  raw$beta_Y, raw$se_Y)
  cat("  theta_hat:", round(egg$theta, 4), "\n")
  cat("  intercept:", round(egg$intercept, 4), "\n")

  cat("\nDiagnostic complete.\n")
  invisible(list(cml = cml, med = med, egg = egg))
}

run_all <- function() {
  cat("=================================================================\n")
  cat("Running all replication simulations. Expected total time: 30-60 min.\n")
  cat("=================================================================\n")

  # Closed-form
  table_kappa()
  table_ad_apoe()
  table_pl_gap()
  table_nongaussian()

  # Per-SNP simulations
  table_bias_coverage()
  table_rho_sens()
  table_rho_sens_n20k()
  table_cc()

  # Bias-slope simulations
  table_bias_slope_recovery()
  table_sigma_beta()
  table_multi_covar()

  # MR simulations
  table_t2d_cad()
  table_theta_sens()
  table_residual_decomp()
  table_mvmr()
  table_mvmr_strength()
  table_mvmr_non_ivw()

  # Downstream-analysis demonstrations
  table_heritability()
  table_pgs()

  cat("\n=================================================================\n")
  cat("All tables complete.\n")
  cat("=================================================================\n")
}

# Auto-print closed-form results when sourced (sanity check)
if (interactive()) {
  cat("Loaded cbcGWAS_replication.R\n")
  cat("Run run_quick() for closed-form tables (< 1 min)\n")
  cat("Run run_all()   for all tables (30-60 min)\n")
  cat("Run table_<name>() for specific tables; type ls() to see options.\n")
  cat("Run diagnose_mvmr_packages() to check MVMR-cML / MendelianRandomization setup.\n")
}
