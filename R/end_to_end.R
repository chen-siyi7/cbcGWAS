#' Simulate estimated bias correction and scale conversion
#'
#' Run the individual-level population and case-control experiments combining
#' CWLS and Slope-Hunter bias estimates with estimated conversion factors.
#' @param R Integer replicates per sampling design, at least one.
#' @param cores Integer worker count. Use one on Windows.
#' @details The full design has 300 variants, 50,000 GWAS individuals, independent
#'   covariate and outcome samples of 100,000 each, and a calibration cohort of
#'   20,000. The two designs use population sampling or a sample case fraction
#'   of one half. The population prevalence is one tenth.
#'
#'   This experiment needs the optional \pkg{indexevent}, \pkg{SlopeHunter}, and
#'   \pkg{mr.raps} packages. Full runs are computationally expensive; multiple
#'   workers also require substantial memory. Replicate seeds are fixed and use
#'   the L'Ecuyer-CMRG generator. The original output key \code{Dudbridge}
#'   denotes CWLS. Marginal SEs omit covariance between the bias-subtracted
#'   association and conversion factor, and cross-variant covariance.
#' @return A data frame of per-replicate estimates and metrics for each sampling
#'   design, bias estimator, and conversion factor. No files are written.
#' @examples
#' \dontrun{
#' x <- simulate_end_to_end(R = 200, cores = 1)
#' summarize_end_to_end(x)
#' }
#' @export
simulate_end_to_end <- function(R = 200L, cores = 1L) {
  .check_count(R, "R")
  .check_count(cores, "cores")
  restore_rng <- .preserve_rng()
  on.exit(restore_rng(), add = TRUE)
  needed <- c('statmod', 'indexevent', 'SlopeHunter', 'mr.raps')
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Packages required: ', paste(missing, collapse = ', '))
  if (.Platform$OS.type == 'windows' && cores > 1L) stop('Use cores = 1 on Windows.')
  # Working-logistic factor for a population (pi = p) or case-control sample,
  # solved from the null logistic score equations (as in revision_checks.R).
  gh_rule <- function(n = 200) {
    a <- statmod::gauss.quad(n, 'hermite')
    list(z = sqrt(2) * a$nodes, w = a$weights / sqrt(pi))
  }
  lambda_working <- function(p, rho2, pi = p, n = 200) {
    if (rho2 <= 0) return(lambda_marg(p))
    g <- gh_rule(n)
    X <- cbind(1, g$z)
    risk <- pnorm((sqrt(rho2) * g$z - qnorm(1 - p)) / sqrt(1 - rho2))
    dotrisk <- dnorm((sqrt(rho2) * g$z - qnorm(1 - p)) / sqrt(1 - rho2)) / sqrt(1 - rho2)
    r <- pi * (1 - p) / (p * (1 - pi))
    sw <- 1 + (r - 1) * risk
    srisk <- r * risk / sw
    w <- g$w * sw / sum(g$w * sw)
    par <- c(qlogis(pi), 1)
    for (i in 1:100) {
      q <- plogis(drop(X %*% par))
      grad <- drop(crossprod(X, w * (srisk - q)))
      if (max(abs(grad)) < 1e-13) break
      par <- par + solve(crossprod(X, (w * q * (1 - q)) * X), grad)
    }
    q <- plogis(drop(X %*% par))
    sum(g$w * dotrisk * (r - (r - 1) * q)) / sum(g$w * sw * q * (1 - q))
  }
  dlambda <- function(f, rho2, h = 1e-4) (f(rho2 + h) - f(rho2 - h)) / (2 * h)
  
  # Design (fixed across replicates; SNP effects are redrawn in each replicate).
  DESIGN <- list(K_H = 100, K_Y = 50, K_both = 50, K_null = 100,
                 sd_beta = 0.0365, sd_gamma = 0.03,
                 a_U = 0.5, a_UH = 0.4, v_zeta = 0.45, delta = 0.4,
                 p = 0.10, n_gwas = 50000, n_H = 100000, n_W = 100000, n_cal = 20000,
                 theta = 0.4)
  K <- with(DESIGN, K_H + K_Y + K_both + K_null)
  
  draw_architecture <- function() with(DESIGN, {
    cls <- rep(c('H_only', 'Y_only', 'both', 'null'), c(K_H, K_Y, K_both, K_null))
    beta <- ifelse(cls %in% c('H_only', 'both'), rnorm(K, 0, sd_beta), 0)
    gamma <- ifelse(cls %in% c('Y_only', 'both'), rnorm(K, 0, sd_gamma), 0)
    maf <- runif(K, 0.05, 0.5)
    # Unit liability variance, including all structural covariance terms.
    v_eta <- 1 - sum((gamma + delta * beta)^2) - (delta * a_UH + a_U)^2 - delta^2 * v_zeta
    stopifnot(v_eta > 0)
    # R_Y is the non-focal part of Y* - delta*H; R_H is the non-focal part of H.
    covRYRH <- a_U * a_UH + (sum(gamma * beta) - gamma * beta)
    varRH <- (sum(beta^2) - beta^2) + a_UH^2 + v_zeta
    b_lin <- -covRYRH / varRH
    delta_eff <- delta - b_lin
    rho2 <- delta_eff^2 * varRH
    list(cls = cls, beta = beta, gamma = gamma, maf = maf, v_eta = v_eta,
         b_lin = b_lin, rho2 = rho2)
  })
  
  simulate_people <- function(n, A) with(DESIGN, {
    g <- matrix(rbinom(n * K, 2, rep(A$maf, each = n)), n, K)
    G <- sweep(sweep(g, 2, 2 * A$maf), 2, sqrt(2 * A$maf * (1 - A$maf)), '/')
    U <- rnorm(n)
    H <- drop(G %*% A$beta) + a_UH * U + rnorm(n, 0, sqrt(v_zeta))
    L <- drop(G %*% A$gamma) + delta * H + a_U * U + rnorm(n, 0, sqrt(A$v_eta))
    list(G = G, H = H, Y = as.numeric(L > qnorm(1 - p)), Gdirect = drop(G %*% A$gamma))
  })
  
  case_control_sample <- function(n, pi, A) {
    need <- c(case = round(n * pi), control = n - round(n * pi))
    keep <- list(); got <- c(case = 0, control = 0)
    while (any(got < need)) {
      s <- simulate_people(100000, A)
      for (st in c('case', 'control')) {
        idx <- which(s$Y == (st == 'case'))
        idx <- head(idx, need[st] - got[st])
        if (length(idx)) {
          keep[[length(keep) + 1]] <- list(G = s$G[idx, , drop = FALSE], H = s$H[idx], Y = s$Y[idx])
          got[st] <- got[st] + length(idx)
        }
      }
    }
    list(G = do.call(rbind, lapply(keep, `[[`, 'G')), H = unlist(lapply(keep, `[[`, 'H')),
         Y = unlist(lapply(keep, `[[`, 'Y')))
  }
  
  # Adjusted logistic GWAS of Y on (1, G_k, H), fitted by Newton-Raphson for all
  # SNPs simultaneously; returns ML estimates and model-based SEs of the SNP term.
  logistic_gwas <- function(Y, G, H) {
    n <- length(Y); K <- ncol(G)
    f0 <- glm.fit(cbind(1, H), Y, family = binomial())
    a <- rep(f0$coefficients[1], K); b <- rep(0, K); d <- rep(f0$coefficients[2], K)
    for (it in 1:25) {
      eta <- G * rep(b, each = n) + outer(H, d) + rep(a, each = n)
      P <- plogis(eta); rm(eta)
      Rr <- Y - P; W <- P * (1 - P); rm(P)
      ga <- colSums(Rr); gb <- colSums(G * Rr); gd <- colSums(H * Rr); rm(Rr)
      GW <- G * W
      Iaa <- colSums(W); Iab <- colSums(GW); Iad <- colSums(H * W)
      Ibb <- colSums(G * GW); Ibd <- colSums(H * GW); Idd <- colSums(H^2 * W)
      rm(GW, W)
      # 3x3 inverse by cofactors, vectorized over SNPs.
      C11 <- Ibb * Idd - Ibd^2; C12 <- -(Iab * Idd - Ibd * Iad); C13 <- Iab * Ibd - Ibb * Iad
      C22 <- Iaa * Idd - Iad^2; C23 <- -(Iaa * Ibd - Iab * Iad); C33 <- Iaa * Ibb - Iab^2
      det <- Iaa * C11 + Iab * C12 + Iad * C13
      da <- (C11 * ga + C12 * gb + C13 * gd) / det
      db <- (C12 * ga + C22 * gb + C23 * gd) / det
      dd <- (C13 * ga + C23 * gb + C33 * gd) / det
      a <- a + da; b <- b + db; d <- d + dd
      if (max(abs(c(da, db, dd))) < 1e-10) break
    }
    stopifnot(it < 25)
    list(beta = b, se = sqrt(C22 / det))
  }
  
  linear_gwas <- function(Z, G) {
    n <- length(Z)
    Gc <- sweep(G, 2, colMeans(G)); Zc <- Z - mean(Z)
    sxx <- colSums(Gc^2); b <- drop(crossprod(Gc, Zc)) / sxx
    s2 <- (sum(Zc^2) - b^2 * sxx) / (n - 2)
    list(beta = b, se = sqrt(s2 / sxx))
  }
  
  # Calibration inputs.
  # Conditional delta SE: empirical var(H) is treated as fixed.
  probit_rho2 <- function(Y, H) {
    fit <- glm(Y ~ H, family = binomial(link = 'probit'))
    u <- coef(fit)[2]; s2 <- var(H); V <- u^2 * s2
    grad <- 2 * u * s2 / (1 + V)^2
    c(rho2 = unname(V / (1 + V)), se = unname(abs(grad) * sqrt(vcov(fit)[2, 2])))
  }
  fitted_variance_factor <- function(Y, H, B = 40) {
    est <- function(y, h) {
      q <- glm.fit(cbind(1, h), y, family = binomial())$fitted.values
      dnorm(qnorm(1 - mean(y))) / mean(q * (1 - q))
    }
    boot <- replicate(B, { i <- sample.int(length(Y), replace = TRUE); est(Y[i], H[i]) })
    c(lambda = est(Y, H), se = sd(boot))
  }
  
  # Bias-slope estimation by the published implementations.
  # Legacy key 'Dudbridge' denotes the CWLS method; figures label it CWLS.
  bias_dudbridge <- function(bx, sx, by, sy) {
    f <- indexevent::indexevent(bx, sx, by, sy, method = 'CWLS')
    list(b = f$b, b_se = f$b.se, A = f$ybeta.adj, A_se = f$yse.adj)
  }
  bias_slopehunter <- function(bx, sx, by, sy, seed) {
    dat <- data.frame(SNP = paste0('snp', seq_along(bx)), BETA.incidence = bx, SE.incidence = sx,
                      BETA.prognosis = by, SE.prognosis = sy)
    out <- NULL
    invisible(capture.output(suppressWarnings(suppressMessages(
      out <- SlopeHunter::hunt(dat, xp_thresh = 0.001, Bootstrapping = TRUE, M = 100,
                               seed = seed, Plot = FALSE, show_adjustments = TRUE)))))
    est <- out$est[match(dat$SNP, out$est$SNP), ]
    list(b = out$b, b_se = out$bse, A = est$ybeta_adj, A_se = est$yse_adj)
  }
  bias_oracle <- function(b_true, bx, sx, by, sy)
    list(b = b_true, b_se = 0, A = by - b_true * bx, A_se = sqrt(sy^2 + b_true^2 * sx^2))
  
  mr_ivw <- function(bx, bw, sw) {
    w <- 1 / sw^2
    c(est = sum(w * bx * bw) / sum(w * bx^2), se = sqrt(1 / sum(w * bx^2)))
  }
  mr_raps <- function(bx, sx, bw, sw) {
    # Profile-score MR-RAPS (Zhao et al. 2020), without over-dispersion.
    f <- mr.raps::mr.raps.mle(bx, bw, sx, sw, over.dispersion = FALSE, loss.function = 'l2',
                              se.method = 'sandwich', suppress.warning = TRUE)
    c(est = f$beta.hat, se = f$beta.se)
  }
  
  one_replicate <- function(r, scenario) {
    set.seed(20260928L + 1000L * scenario$id + r)
    A <- draw_architecture()
    p <- DESIGN$p
    gw <- if (scenario$sampling == 'population') simulate_people(DESIGN$n_gwas, A) else
      case_control_sample(DESIGN$n_gwas, scenario$pi, A)
    Ysum <- logistic_gwas(gw$Y, gw$G, gw$H)
    fv <- if (scenario$sampling == 'population') fitted_variance_factor(gw$Y, gw$H) else NULL
    rm(gw); gc(FALSE)
    hs <- simulate_people(DESIGN$n_H, A); Hsum <- linear_gwas(hs$H, hs$G); rm(hs)
    ws <- simulate_people(DESIGN$n_W, A)
    # Proportional-association target; W is generated from the direct genetic component.
    W <- DESIGN$theta * ws$Gdirect + rnorm(DESIGN$n_W, 0, sqrt(1 - DESIGN$theta^2 * sum(A$gamma^2)))
    Wsum <- linear_gwas(W, ws$G); rm(ws, W)
    cal <- simulate_people(DESIGN$n_cal, A); pr <- probit_rho2(cal$Y, cal$H); rm(cal); gc(FALSE)
  
    pi_s <- if (scenario$sampling == 'population') p else scenario$pi
    rho2_true <- mean(A$rho2)
    lambda_true <- lambda_working(p, rho2_true, pi_s)
    b_true <- lambda_true * mean(A$b_lin)
  
    f_cond <- function(x) lambda_cond_gauss(p, x)
    f_pi <- function(x) lambda_working(p, x, pi_s)
    factors <- list(Marginal = c(lambda_marg(p), 0),
                    Conditional = c(f_cond(pr['rho2']), abs(dlambda(f_cond, pr['rho2'])) * pr['se']))
    if (scenario$sampling == 'population') {
      factors$`Fitted variance` <- unname(fv)
    } else {
      factors$`Sampled factor` <- c(f_pi(pr['rho2']), abs(dlambda(f_pi, pr['rho2'])) * pr['se'])
    }
    factors$`True factor` <- c(lambda_true, 0)
  
    bias <- list(`Known slope` = bias_oracle(b_true, Hsum$beta, Hsum$se, Ysum$beta, Ysum$se),
                 Dudbridge = bias_dudbridge(Hsum$beta, Hsum$se, Ysum$beta, Ysum$se),
                 `Slope-Hunter` = bias_slopehunter(Hsum$beta, Hsum$se, Ysum$beta, Ysum$se, seed = r))
  
    h2_true <- sum(A$gamma^2)
    sig <- Ysum$beta * 0
    out <- list()
    for (bn in names(bias)) {
      bs <- bias[[bn]]
      selected <- which(2 * pnorm(-abs(bs$A / bs$A_se)) < 5e-8)
      for (fn in names(factors)) {
        lam <- factors[[fn]][1]; lam_se <- factors[[fn]][2]
        est <- bs$A / lam
        # Marginal approximation: Cov(A_i, lambda) and cross-SNP covariance omitted.
        se <- sqrt(bs$A_se^2 / lam^2 + bs$A^2 * lam_se^2 / lam^4)
        cover <- abs(est - A$gamma) <= qnorm(0.975) * se
        nz <- A$gamma != 0
        ivw <- if (length(selected) >= 3) mr_ivw(est[selected], Wsum$beta[selected], Wsum$se[selected]) else c(NA, NA)
        raps <- tryCatch(mr_raps(est, se, Wsum$beta, Wsum$se), error = function(e) c(NA, NA))
        out[[length(out) + 1]] <- data.frame(
          scenario = scenario$name, replicate = r, bias_method = bn, factor = fn,
          rho2_true = rho2_true, rho2_hat = unname(pr['rho2']), lambda_true = lambda_true,
          lambda_used = lam, b_true = b_true, b_hat = bs$b, b_se = bs$b_se,
          scale_slope = sum(A$gamma[nz] * est[nz]) / sum(A$gamma[nz]^2),
          coverage_direct = mean(cover[nz]), coverage_H_only = mean(cover[A$cls == 'H_only']),
          coverage_null = mean(cover[A$cls == 'null']),
          h2_true = h2_true, h2_hat = sum(est^2 - se^2),
          n_instruments = length(selected),
          ivw = unname(ivw[1]), ivw_se = unname(ivw[2]),
          raps = unname(raps[1]), raps_se = unname(raps[2]))
      }
    }
    do.call(rbind, out)
  }
  
  scenarios <- list(list(id = 1L, name = 'Population sample', sampling = 'population'),
                    list(id = 2L, name = 'Case-control sample', sampling = 'case-control', pi = .5))
  RNGkind("L'Ecuyer-CMRG")
  results <- list()
  for (sc in scenarios) {
    reps <- parallel::mclapply(seq_len(R), one_replicate, scenario = sc,
                               mc.cores = cores, mc.preschedule = FALSE)
    failed <- vapply(reps, inherits, logical(1), 'try-error')
    if (any(failed)) stop('Replicate failures in ', sc$name, ': ', paste(which(failed), collapse = ', '),
                          '\n', as.character(reps[[which(failed)[1]]]))
    results[[sc$name]] <- do.call(rbind, reps)
    message(sc$name, ': ', R, ' replicates completed.')
  }
  do.call(rbind, results)
}

# Summarize replicate data returned by simulate_end_to_end; matches the S18 export.
#' Summarize the individual-level validation simulations
#'
#' Aggregate the per-replicate output of \code{simulate_end_to_end()}.
#' @param x Data frame returned by \code{simulate_end_to_end()}.
#' @param theta Specified MR association-ratio target used in the simulation.
#' @return A data frame grouped by sampling scenario, bias method, and factor,
#'   with scale bias, coverage, MR metrics, and Monte Carlo SEs. SEs require
#'   at least two nonmissing replicate values in each group.
#' @seealso \code{\link{simulate_end_to_end}}
#' @export
summarize_end_to_end <- function(x, theta = .4) {
  mcse <- function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))
  keys <- unique(x[, c('scenario', 'bias_method', 'factor')])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    d <- merge(keys[i, ], x)
    cover <- function(est, se) 100 * mean(abs(est - theta) <= qnorm(0.975) * se, na.rm = TRUE)
    data.frame(
      keys[i, ], replicates = nrow(d),
      rho2_true = mean(d$rho2_true), rho2_hat = mean(d$rho2_hat),
      lambda_true = mean(d$lambda_true), lambda_used = mean(d$lambda_used),
      bias_slope_ratio = mean(d$b_hat / d$b_true), bias_slope_ratio_mcse = mcse(d$b_hat / d$b_true),
      effect_relative_bias_pct = 100 * (mean(d$scale_slope) - 1),
      effect_relative_bias_mcse = 100 * mcse(d$scale_slope),
      coverage_direct_pct = 100 * mean(d$coverage_direct),
      coverage_H_only_pct = 100 * mean(d$coverage_H_only),
      coverage_null_pct = 100 * mean(d$coverage_null),
      h2_relative_bias_pct = 100 * mean(d$h2_hat / d$h2_true - 1),
      h2_relative_bias_mcse = 100 * mcse(d$h2_hat / d$h2_true),
      ivw_instruments = mean(d$n_instruments),
      ivw_relative_bias_pct = 100 * (mean(d$ivw, na.rm = TRUE) / theta - 1),
      ivw_relative_bias_mcse = 100 * mcse(d$ivw) / theta,
      ivw_coverage_pct = cover(d$ivw, d$ivw_se),
      raps_relative_bias_pct = 100 * (mean(d$raps, na.rm = TRUE) / theta - 1),
      raps_relative_bias_mcse = 100 * mcse(d$raps) / theta,
      raps_coverage_pct = cover(d$raps, d$raps_se),
      check.names = FALSE)
  })
  out <- do.call(rbind, rows)
  ord <- order(match(out$scenario, c('Population sample', 'Case-control sample')),
               match(out$bias_method, c('Known slope', 'Dudbridge', 'Slope-Hunter')),
               match(out$factor, c('Marginal', 'Conditional', 'Fitted variance', 'Sampled factor', 'True factor')))
  out <- out[ord, ]
  rownames(out) <- NULL
  out
}
