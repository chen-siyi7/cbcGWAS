.print_rounded <- function(df, digits = 3) {
    if (is.data.frame(df)) {
        num_cols <- sapply(df, is.numeric)
        df[, num_cols] <- round(df[, num_cols], digits)
        print(df)
    }
    else {
        print(round(df, digits))
    }
}

sim_one_no_collider <- function(n, p_Y, rho2, gamma = 0.05) {
    H <- rnorm(n, 0, sqrt(rho2))
    G <- rnorm(n)
    res_var <- max(1 - rho2 - gamma^2, 1e-08)
    Y_star <- H + gamma * G + rnorm(n, 0, sqrt(res_var))
    t_Y <- qnorm(1 - p_Y)
    Y <- as.numeric(Y_star > t_Y)
    fit <- glm(Y ~ G + H, family = binomial())
    beta_C <- coef(fit)["G"]
    se_C <- summary(fit)$coefficients["G", "Std. Error"]
    list(beta_C = beta_C, se_C = se_C)
}

table_bias_coverage <- function(R = 500, n = 20000, gamma = 0.05) {
    cells <- list(c(0.1, 0.05), c(0.1, 0.1), c(0.2, 0.1), c(0.05, 0.3), c(0.1, 0.3), c(0.2, 0.3), 
        c(0.5, 0.3), c(0.5, 0.5))
    out <- data.frame(p_Y = numeric(), rho2 = numeric(), kappa = numeric(), marg_bias = numeric(), 
        marg_mcse = numeric(), score_bias = numeric(), score_mcse = numeric(), marg_cov = numeric(), 
        score_cov = numeric())
    for (k in seq_along(cells)) {
        set.seed(20260101L + 100L + k)
        p_Y <- cells[[k]][1]
        rho2 <- cells[[k]][2]
        kap <- kappa_score(p_Y, rho2)
        lm <- lambda_marg(p_Y)
        lc <- lm * kap
        rel_marg <- numeric(R)
        rel_score <- numeric(R)
        cov_marg <- logical(R)
        cov_score <- logical(R)
        for (r in 1:R) {
            rep <- sim_one_no_collider(n, p_Y, rho2, gamma)
            beta_marg <- rep$beta_C/lm
            beta_score <- rep$beta_C/lc
            se_marg <- rep$se_C/lm
            se_score <- rep$se_C/lc
            rel_marg[r] <- beta_marg/gamma - 1
            rel_score[r] <- beta_score/gamma - 1
            cov_marg[r] <- abs(beta_marg - gamma) < 1.96 * se_marg
            cov_score[r] <- abs(beta_score - gamma) < 1.96 * se_score
        }
        out[k, ] <- list(p_Y, rho2, kap, mean(rel_marg), sd(rel_marg)/sqrt(R), mean(rel_score), 
            sd(rel_score)/sqrt(R), mean(cov_marg) * 100, mean(cov_score) * 100)
    }
    cat("\n=== Population simulation: bias and coverage ===\n")
    .print_rounded(out, 3)
    invisible(out)
}

table_cc <- function(R = 200, n = 20000, p_Y = 0.1, rho2 = 0.3, gamma = 0.05) {
    pi_grid <- c(0.05, 0.1, 0.3, 0.5)
    out <- data.frame()
    for (pi_frac in pi_grid) {
        set.seed(20260101L + 300L + round(pi_frac * 100))
        estimates <- numeric(R)
        intercepts <- numeric(R)
        n_pop <- max(n, round(n * pi_frac/p_Y) + n)
        for (r in 1:R) {
            H <- rnorm(n_pop, 0, sqrt(rho2))
            G <- rnorm(n_pop)
            Y_star <- H + gamma * G + rnorm(n_pop, 0, sqrt(1 - rho2 - gamma^2))
            Y <- as.numeric(Y_star > qnorm(1 - p_Y))
            cases <- which(Y == 1)
            ctrls <- which(Y == 0)
            n_case_target <- round(n * pi_frac)
            n_ctrl_target <- n - n_case_target
            if (length(cases) < n_case_target || length(ctrls) < n_ctrl_target) {
                estimates[r] <- NA
                intercepts[r] <- NA
                next
            }
            kept <- c(sample(cases, n_case_target), sample(ctrls, n_ctrl_target))
            fit <- glm(Y[kept] ~ G[kept] + H[kept], family = binomial())
            estimates[r] <- coef(fit)[2]
            intercepts[r] <- coef(fit)[1]
        }
        out <- rbind(out, data.frame(pi = pi_frac, pi_over_pY = pi_frac/p_Y, beta_G_mean = mean(estimates, 
            na.rm = TRUE), SE = sd(estimates, na.rm = TRUE)/sqrt(R), alpha = mean(intercepts, na.rm = TRUE), 
            theory = lambda_working(p_Y, rho2, pi_frac)$lambda * gamma))
    }
    cat("\n=== Case-control simulation ===\n")
    .print_rounded(out, 3)
    invisible(out)
}

sim_t2d_cad_one <- function(K = 120, n_T2D = 80000, n_BMI = 80000, n_CAD = 1e+05, p_T2D = 0.1, p_CAD = 0.06, 
    rho2 = 0.3, theta = 0.4, b_lin = 0.5, gamma_sd = 0.025, beta_GH_sd = 0.02, known_b = FALSE) {
    gamma_k <- rnorm(K, 0, gamma_sd)
    beta_GH_k <- rnorm(K, 0, beta_GH_sd)
    lc_T2D <- lambda_cond_gauss(p_T2D, rho2)
    lm_T2D <- lambda_marg(p_T2D)
    lc_CAD <- lambda_marg(p_CAD)
    beta_C_T2D <- lc_T2D * (gamma_k + b_lin * beta_GH_k)
    se_C_T2D <- 1/sqrt(n_T2D * p_T2D * (1 - p_T2D))
    beta_C_T2D <- beta_C_T2D + rnorm(K, 0, se_C_T2D)
    beta_GH_obs <- beta_GH_k + rnorm(K, 0, 1/sqrt(n_BMI))
    se_GH <- rep(1/sqrt(n_BMI), K)
    beta_GY_CAD <- lc_CAD * theta * gamma_k + rnorm(K, 0, 1/sqrt(n_CAD * p_CAD * (1 - p_CAD)))
    se_GY_CAD <- rep(1/sqrt(n_CAD * p_CAD * (1 - p_CAD)), K)
    beta_Y_liab <- beta_GY_CAD/lc_CAD
    se_Y_liab <- se_GY_CAD/lc_CAD
    b_logOR_true <- lc_T2D * b_lin
    if (known_b) {
        b_hat <- b_logOR_true
    }
    else {
        b_hat <- coef(lm(beta_C_T2D ~ beta_GH_obs))[2]
    }
    A <- beta_C_T2D - b_hat * beta_GH_obs
    beta_X_marg <- A/lm_T2D
    se_X_marg <- se_C_T2D/lm_T2D
    beta_X_score <- A/lc_T2D
    se_X_score <- se_C_T2D/lc_T2D
    beta_X_uncorr <- A
    se_X_uncorr <- se_C_T2D
    marg <- ivw_mr(beta_X_marg, se_X_marg, beta_Y_liab, se_Y_liab)
    score <- ivw_mr(beta_X_score, se_X_score, beta_Y_liab, se_Y_liab)
    uncorr <- ivw_mr(beta_X_uncorr, se_X_uncorr, beta_Y_liab, se_Y_liab)
    list(marg = marg, score = score, uncorr = uncorr)
}

run_t2d_cad <- function(R = 100, rho2_grid = c(0.2, 0.3, 0.4), theta = 0.4, K = 120, n_T2D = 80000, 
    known_b = FALSE) {
    out <- data.frame()
    for (rho2 in rho2_grid) {
        set.seed(20260101L + 600L + round(rho2 * 100))
        kap <- kappa_score(0.1, rho2)
        em <- es <- eu <- numeric(R)
        sm <- ss <- su <- numeric(R)
        for (r in 1:R) {
            rep <- sim_t2d_cad_one(K = K, n_T2D = n_T2D, p_T2D = 0.1, p_CAD = 0.06, rho2 = rho2, 
                theta = theta, known_b = known_b)
            em[r] <- rep$marg$theta
            es[r] <- rep$score$theta
            eu[r] <- rep$uncorr$theta
            sm[r] <- rep$marg$se
            ss[r] <- rep$score$se
            su[r] <- rep$uncorr$se
        }
        cov_m <- mean(abs(em - theta) < 1.96 * sm) * 100
        cov_s <- mean(abs(es - theta) < 1.96 * ss) * 100
        cov_u <- mean(abs(eu - theta) < 1.96 * su) * 100
        out <- rbind(out, data.frame(rho2 = rho2, kappa = kap, method = "Marginal-scaled", theta_hat = mean(em), 
            SE = mean(sm), rel_bias = mean(em)/theta - 1, coverage = cov_m), data.frame(rho2 = rho2, 
            kappa = kap, method = "Conditional-scaled", theta_hat = mean(es), SE = mean(ss), rel_bias = mean(es)/theta - 
                1, coverage = cov_s), data.frame(rho2 = rho2, kappa = kap, method = "Uncorrected log-OR", 
            theta_hat = mean(eu), SE = mean(su), rel_bias = mean(eu)/theta - 1, coverage = cov_u))
    }
    out
}

sim_mvmr_one <- function(K = 30, n_X = 8000, p1 = 0.1, rho2_1 = 0.3, p2 = 0.1, rho2_2 = 0.3, theta = c(0.4, 
    -0.3), gamma_sd = 0.06) {
    gamma_1 <- rnorm(K, 0, gamma_sd)
    gamma_2 <- rnorm(K, 0, gamma_sd)
    lc1 <- lambda_cond_gauss(p1, rho2_1)
    lc2 <- lambda_cond_gauss(p2, rho2_2)
    lm1 <- lambda_marg(p1)
    lm2 <- lambda_marg(p2)
    se_X1 <- 1/sqrt(n_X * p1 * (1 - p1))
    se_X2 <- 1/sqrt(n_X * p2 * (1 - p2))
    A1 <- lc1 * gamma_1 + rnorm(K, 0, se_X1)
    A2 <- lc2 * gamma_2 + rnorm(K, 0, se_X2)
    beta_Y <- theta[1] * gamma_1 + theta[2] * gamma_2 + rnorm(K, 0, 0.005)
    se_Y <- rep(0.005, K)
    X_naive <- cbind(A1, A2)
    X_marg <- cbind(A1/lm1, A2/lm2)
    X_score <- cbind(A1/lc1, A2/lc2)
    SE_naive <- cbind(rep(se_X1, K), rep(se_X2, K))
    SE_marg <- cbind(rep(se_X1/lm1, K), rep(se_X2/lm2, K))
    SE_score <- cbind(rep(se_X1/lc1, K), rep(se_X2/lc2, K))
    list(naive = mvmr_ivw(X_naive, beta_Y, se_Y), marg = mvmr_ivw(X_marg, beta_Y, se_Y), score = mvmr_ivw(X_score, 
        beta_Y, se_Y), raw = list(X_naive = X_naive, SE_naive = SE_naive, X_marg = X_marg, SE_marg = SE_marg, 
        X_score = X_score, SE_score = SE_score, beta_Y = beta_Y, se_Y = se_Y, n_X = n_X))
}

table_mvmr <- function(R = 500, K = 30, n_X = 8000, theta = c(0.4, -0.3)) {
    scenarios <- list(list(p1 = 0.1, r1 = 0.3, p2 = 0.1, r2 = 0.3, label = "S1: equal kappa"), list(p1 = 0.1, 
        r1 = 0.3, p2 = 0.05, r2 = 0.1, label = "S2: unequal kappa"))
    out <- data.frame()
    for (s in scenarios) {
        set.seed(20260101L + 800L + which(sapply(scenarios, function(x) x$label) == s$label))
        k1 <- kappa_score(s$p1, s$r1)
        k2 <- kappa_score(s$p2, s$r2)
        naive <- marg <- score <- matrix(NA_real_, nrow = R, ncol = 2)
        for (r in 1:R) {
            rep <- sim_mvmr_one(K = K, n_X = n_X, p1 = s$p1, rho2_1 = s$r1, p2 = s$p2, rho2_2 = s$r2, 
                theta = theta)
            naive[r, ] <- rep$naive$theta
            marg[r, ] <- rep$marg$theta
            score[r, ] <- rep$score$theta
        }
        out <- rbind(out, data.frame(scenario = s$label, kappa1 = k1, kappa2 = k2, method = "Naive (log-OR)", 
            theta1 = mean(naive[, 1]), theta1_mcse = sd(naive[, 1])/sqrt(R), theta2 = mean(naive[, 
                2]), theta2_mcse = sd(naive[, 2])/sqrt(R)), data.frame(scenario = s$label, kappa1 = k1, 
            kappa2 = k2, method = "Marginal", theta1 = mean(marg[, 1]), theta1_mcse = sd(marg[, 
                1])/sqrt(R), theta2 = mean(marg[, 2]), theta2_mcse = sd(marg[, 2])/sqrt(R)), data.frame(scenario = s$label, 
            kappa1 = k1, kappa2 = k2, method = "Conditional-scaled", theta1 = mean(score[, 1]), theta1_mcse = sd(score[, 
                1])/sqrt(R), theta2 = mean(score[, 2]), theta2_mcse = sd(score[, 2])/sqrt(R)))
    }
    cat("\n=== MVMR with two binary exposures ===\n")
    .print_rounded(out, 3)
    invisible(out)
}

table_bias_slope_recovery <- function(R = 200, K = 200, n = 50000, b_lin = 0.5, gamma_sd = 0.025, 
    beta_GH_sd = 0.05) {
    cells <- list(c(0.1, 0.3), c(0.2, 0.3), c(0.1, 0.4))
    out <- data.frame()
    for (cell in cells) {
        p_Y <- cell[1]
        rho2 <- cell[2]
        set.seed(20260101L + 400L + round(p_Y * 100) + round(rho2 * 100))
        lc <- lambda_cond_gauss(p_Y, rho2)
        b_logOR_true <- lc * b_lin
        est_OLS <- numeric(R)
        est_SH <- numeric(R)
        est_cML <- numeric(R)
        for (r in 1:R) {
            gamma_k <- rnorm(K, 0, gamma_sd)
            beta_GH_k <- rnorm(K, 0, beta_GH_sd)
            se_C <- 1/sqrt(n * p_Y * (1 - p_Y))
            beta_C <- lc * (gamma_k + b_lin * beta_GH_k) + rnorm(K, 0, se_C)
            beta_GH_obs <- beta_GH_k + rnorm(K, 0, 1/sqrt(n))
            est_OLS[r] <- coef(lm(beta_C ~ beta_GH_obs))[2]
            nh_idx <- order(abs(beta_C - est_OLS[r] * beta_GH_obs))[1:max(20, K%/%4)]
            est_SH[r] <- coef(lm(beta_C[nh_idx] ~ beta_GH_obs[nh_idx]))[2]
            w <- 1/se_C^2
            est_cML[r] <- sum(w * beta_C * beta_GH_obs)/sum(w * beta_GH_obs^2)
        }
        out <- rbind(out, data.frame(p_Y = p_Y, rho2 = rho2, lambda_cond = lc, true_b_logOR = b_logOR_true, 
            estimator = "Ordinary regression", b_hat = mean(est_OLS), mcse = sd(est_OLS)/sqrt(R)), 
            data.frame(p_Y = p_Y, rho2 = rho2, lambda_cond = lc, true_b_logOR = b_logOR_true, estimator = "Residual-trimmed regression", 
                b_hat = mean(est_SH), mcse = sd(est_SH)/sqrt(R)), data.frame(p_Y = p_Y, rho2 = rho2, 
                lambda_cond = lc, true_b_logOR = b_logOR_true, estimator = "Precision-weighted regression", 
                b_hat = mean(est_cML), mcse = sd(est_cML)/sqrt(R)))
    }
    cat("\n=== Bias-slope regression comparison ===\n")
    .print_rounded(out, 3)
    invisible(out)
}

table_multi_covar <- function(R = 300, n = 20000, p_Y = 0.1, gamma = 0.05, total_rho2 = 0.3) {
    q_grid <- c(1, 3, 5)
    out <- data.frame()
    for (q in q_grid) {
        set.seed(20260101L + 700L + q)
        rho2_each <- total_rho2/q
        lc <- lambda_cond_gauss(p_Y, total_rho2)
        lm <- lambda_marg(p_Y)
        rel_marg <- rel_score <- numeric(R)
        cov_marg <- cov_score <- logical(R)
        for (r in 1:R) {
            H <- mvtnorm::rmvnorm(n, sigma = diag(rho2_each, q))
            G <- rnorm(n)
            M <- rowSums(H)
            res_var <- max(1 - total_rho2 - gamma^2, 1e-06)
            Y_star <- M + gamma * G + rnorm(n, 0, sqrt(res_var))
            Y <- as.numeric(Y_star > qnorm(1 - p_Y))
            df <- data.frame(Y = Y, G = G, H = H)
            fit <- glm(Y ~ ., data = df, family = binomial())
            bC <- coef(fit)["G"]
            seC <- summary(fit)$coefficients["G", "Std. Error"]
            bm <- bC/lm
            sm <- seC/lm
            bs <- bC/lc
            ss <- seC/lc
            rel_marg[r] <- bm/gamma - 1
            rel_score[r] <- bs/gamma - 1
            cov_marg[r] <- abs(bm - gamma) < 1.96 * sm
            cov_score[r] <- abs(bs - gamma) < 1.96 * ss
        }
        out <- rbind(out, data.frame(q = q, score_bias = mean(rel_score), marg_bias = mean(rel_marg), 
            score_cov = mean(cov_score) * 100, marg_cov = mean(cov_marg) * 100))
    }
    cat("\n=== Multiple adjusted covariates ===\n")
    .print_rounded(out, 3)
    invisible(out)
}

# Corrected explicit-collider model from revision_checks.R (S1 Text S4.2).
# Residual variance includes the structural covariance terms so liability variance is one.
#' Simulate the corrected unit-variance collider model
#'
#' Evaluate marginal, conditional-prevalence, and working-logistic conversion
#' after subtracting a known first-order collider contribution.
#' @param R Integer number of replicates per design cell, at least two.
#' @param n Integer sample size per replicate, at least 100.
#' @param seed Integer base seed; each design cell uses \code{seed + cell_index}.
#' @details Six cells vary disease prevalence and the variant-covariate effect.
#'   The residual variance includes all structural covariance terms so that
#'   liability has unit variance. Calibration and the subtracted bias slope
#'   are fixed at their generating values. Small runs only check execution.
#' @return A data frame with design inputs, conversion method, mean effect,
#'   relative bias, Monte Carlo SE, and percent interval coverage.
#' @examples
#' x <- simulate_collider(R = 2, n = 1000)
#' head(x)
#' @export
simulate_collider <- function(R = 300L, n = 30000L, seed = 2026091400L) {
  restore_rng <- .preserve_rng()
  on.exit(restore_rng(), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  .check_count(seed, "seed", 0L)
  if (seed > .Machine$integer.max - 6L) stop("seed is too large.", call. = FALSE)
  .check_count(R, "R", 2L)
  .check_count(n, "n", 100L)
  collider <- expand.grid(p = c(0.1, 0.2), a = c(0.02, 0.06, 0.1))
  gamma <- 0.05
  alpha <- 0.5
  c_u <- 0.4
  vz <- 0.45
  delta <- 0.4
  v0 <- c_u^2 + vz
  b_lin <- -alpha * c_u/v0
  rho0 <- (delta + alpha * c_u/v0)^2 * v0
  out <- list()
  for (k in seq_len(nrow(collider))) {
      p <- collider$p[k]
      a <- collider$a[k]
      sigma2 <- 1 - (gamma + delta * a)^2 - (delta * c_u + alpha)^2 - delta^2 * vz
      stopifnot(sigma2 > 0)
      lam <- lambda_working(p, rho0)$lambda
      approx <- lambda_approx(p, rho0)
      marg <- lambda_marginal(p)
      set.seed(seed + k)
      estimates <- ses <- numeric(R)
      for (j in seq_len(R)) {
          G <- rnorm(n)
          U <- rnorm(n)
          H <- a * G + c_u * U + rnorm(n, sd = sqrt(vz))
          L <- gamma * G + delta * H + alpha * U + rnorm(n, sd = sqrt(sigma2))
          Y <- as.numeric(L > qnorm(1 - p))
          fit <- glm(Y ~ G + H, family = binomial())
          estimates[j] <- coef(fit)[2] - lam * b_lin * a
          ses[j] <- sqrt(vcov(fit)[2, 2])
      }
      for (method in c("Marginal", "Conditional approximation", "Working logistic")) {
          scale <- switch(method, Marginal = marg, `Conditional approximation` = approx, `Working logistic` = lam)
          est <- estimates/scale
          out[[length(out) + 1]] <- data.frame(p = p, beta_GH = a, rho2_null = rho0, method = method, 
              n = n, R = R, gamma = gamma, mean = mean(est), relative_bias = mean(est)/gamma - 1, 
              mcse_bias = sd(est/gamma)/sqrt(R), coverage = 100 * mean(abs(est - gamma) <= qnorm(0.975) * 
                  ses/scale))
      }
      cat("Collider cell", k, "completed\n")
  }
  do.call(rbind, out)
}
