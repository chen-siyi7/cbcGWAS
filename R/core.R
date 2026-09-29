gh_rule <- function(n = 200) {
    .check_count(n, "n", 2L)
    a <- statmod::gauss.quad(n, "hermite")
    list(z = sqrt(2) * a$nodes, w = a$weights/sqrt(pi))
}

#' Liability-to-logit conversion factors
#'
#' Compute the marginal factor, a conditional-prevalence approximation, or the
#' working-logistic factor under a Gaussian liability model.
#' @param p,p_Y Population disease prevalence. The marginal functions also accept a vector.
#' @param rho2 Scalar covariate-explained liability variance, in \eqn{[0,1)}.
#' @param n,n_quad Integer Gauss-Hermite quadrature order, at least two.
#' @param pi Sample case fraction. Set to \code{p} for population sampling.
#' @details
#' The direct and collider contributions share a conversion factor to first order
#' in small genetic effects. Divide a bias-subtracted adjusted log-odds coefficient
#' by the factor to express it in liability units. Bias identification is a separate
#' step and depends on the chosen estimator's assumptions.
#'
#' \code{lambda_approx()} and \code{lambda_cond_gauss()} compute the
#' conditional-prevalence approximation; they do not solve the logistic score
#' equations. \code{lambda_working()} solves the null working-logistic projection
#' and accounts for outcome-dependent sampling through \code{pi}. The Gaussian
#' population risk model and its prevalence remain inputs in case-control samples.
#'
#' \code{lambda_marg()} is a compatibility name with argument \code{p_Y}.
#' \code{kappa_cond()} is the conditional-to-marginal factor ratio;
#' \code{kappa_score()} retains its original code name.
#' @return The marginal and approximate functions return positive numeric factors;
#'   the kappa functions return a numeric ratio. \code{lambda_working()} returns a
#'   list with \code{lambda} and \code{score_error}, plus \code{intercept} and
#'   \code{slope} when \code{rho2 > 0}.
#' @examples
#' lambda_marginal(0.10)
#' lambda_approx(0.10, 0.30)
#' lambda_working(0.10, 0.30)$lambda
#' lambda_working(0.10, 0.30, pi = 0.50)$lambda
#' kappa_cond(0.10, 0.30)
#' @name conversion_factors
#' @export
lambda_marginal <- function(p) {
    .check_probability(p, "p", scalar = FALSE)
    dnorm(qnorm(1 - p))/(p * (1 - p))
}

#' @rdname conversion_factors
#' @export
lambda_approx <- function(p, rho2, n = 200) {
    .check_probability(p, "p")
    .check_rho2(rho2)
    g <- gh_rule(n)
    risk <- pnorm((sqrt(rho2) * g$z - qnorm(1 - p))/sqrt(1 - rho2))
    dnorm(qnorm(1 - p))/sum(g$w * risk * (1 - risk))
}

#' @rdname conversion_factors
#' @export
lambda_working <- function(p, rho2, pi = p, n = 200) {
    .check_probability(p, "p")
    .check_rho2(rho2)
    .check_probability(pi, "pi")
    .check_count(n, "n", 2L)
    if (rho2 == 0) 
        return(list(lambda = lambda_marginal(p), score_error = 0))
    g <- gh_rule(n)
    X <- cbind(1, g$z)
    risk <- pnorm((sqrt(rho2) * g$z - qnorm(1 - p))/sqrt(1 - rho2))
    dotrisk <- dnorm((sqrt(rho2) * g$z - qnorm(1 - p))/sqrt(1 - rho2))/sqrt(1 - rho2)
    odds_sampling <- pi * (1 - p)/(p * (1 - pi))
    sw <- 1 + (odds_sampling - 1) * risk
    sample_risk <- odds_sampling * risk/sw
    w <- g$w * sw/sum(g$w * sw)
    par <- c(qlogis(pi), 1)
    softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
    objective <- function(a) sum(w * (softplus(drop(X %*% a)) - sample_risk * drop(X %*% a)))
    for (i in seq_len(100)) {
        q <- plogis(drop(X %*% par))
        gradient <- drop(crossprod(X, w * (sample_risk - q)))
        information <- crossprod(X, (w * q * (1 - q)) * X)
        if (max(abs(gradient)) < 1e-13) 
            break
        step <- solve(information, gradient)
        rate <- 1
        while (objective(par + rate * step) > objective(par) + 1e-14 && rate > 1e-08) rate <- rate/2
        par <- par + rate * step
    }
    q <- plogis(drop(X %*% par))
    num <- sum(g$w * dotrisk * (odds_sampling - (odds_sampling - 1) * q))
    den <- sum(g$w * sw * q * (1 - q))
    err <- max(abs(crossprod(X, w * (sample_risk - q))))
    stopifnot(err < 1e-09)
    list(lambda = num/den, score_error = err, intercept = par[1], slope = par[2])
}

# Original function names
# These definitions are retained for compatibility with the original analyses.
# p_Y is population prevalence and n_quad is the quadrature order.
# lambda_cond_gauss equals lambda_approx; it is an approximation to lambda_working.
# kappa_score is the original name for the ratio called kappa_cond in the revision.

#' @rdname conversion_factors
#' @export
lambda_marg <- function(p_Y) {
    .check_probability(p_Y, "p_Y", scalar = FALSE)
    t_Y <- qnorm(1 - p_Y)
    dnorm(t_Y)/(p_Y * (1 - p_Y))
}

#' @rdname conversion_factors
#' @export
lambda_cond_gauss <- function(p_Y, rho2, n_quad = 200) {
    .check_probability(p_Y, "p_Y")
    .check_rho2(rho2)
    .check_count(n_quad, "n_quad", 2L)
    if (rho2 <= 0) 
        return(lambda_marg(p_Y))
    t_Y <- qnorm(1 - p_Y)
    tau <- sqrt(1 - rho2)
    rho <- sqrt(rho2)
    gh <- statmod::gauss.quad(n_quad, kind = "hermite")
    M <- sqrt(2) * rho * gh$nodes
    p_H <- pnorm((M - t_Y)/tau)
    D <- sum(gh$weights * p_H * (1 - p_H))/sqrt(pi)
    dnorm(t_Y)/D
}

#' @rdname conversion_factors
#' @export
kappa_score <- function(p_Y, rho2, n_quad = 200) {
    lambda_cond_gauss(p_Y, rho2, n_quad)/lambda_marg(p_Y)
}

# Association estimates
# Y is a binary outcome, G a numeric genotype matrix (individuals by variants),
# and H one numeric adjustment covariate. logistic_gwas returns beta and se.
# linear_gwas fits separate unadjusted linear regressions of Z on columns of G.
# The returned coefficients use the units of the supplied genotype matrix.

#' Fit variant-wise association models
#'
#' Fit separate logistic regressions adjusted for one covariate, or separate
#' unadjusted linear regressions, for columns of a genotype matrix.
#' @param Y Numeric binary outcome vector containing both zero and one.
#' @param G Numeric matrix with individuals in rows and variants in columns.
#' @param H Numeric vector of one adjustment covariate.
#' @param Z Numeric continuous outcome vector for the linear regressions.
#' @details Inputs must be complete and finite. Coefficients use the units of the
#'   supplied genotype matrix. The logistic implementation uses Newton-Raphson
#'   updates and requires finite maximum-likelihood estimates and convergence.
#' @return A list with vectors \code{beta} and \code{se}, one element per variant.
#' @examples
#' set.seed(24)
#' G <- matrix(rnorm(1200), ncol = 3)
#' H <- rnorm(nrow(G))
#' Y <- rbinom(nrow(G), 1, plogis(-1 + 0.2 * G[, 1] + 0.3 * H))
#' logistic_gwas(Y, G, H)
#' linear_gwas(H, G)
#' @name association_models
#' @export
logistic_gwas <- function(Y, G, H) {
    .check_individual_data(Y, H)
    .check_genotype_matrix(G, length(Y))
    n <- length(Y)
    K <- ncol(G)
    f0 <- glm.fit(cbind(1, H), Y, family = binomial())
    a <- rep(unname(f0$coefficients[1]), K)
    b <- rep(0, K)
    d <- rep(unname(f0$coefficients[2]), K)
    for (it in 1:25) {
        eta <- G * rep(b, each = n) + outer(H, d) + rep(a, each = n)
        P <- plogis(eta)
        rm(eta)
        Rr <- Y - P
        W <- P * (1 - P)
        rm(P)
        ga <- colSums(Rr)
        gb <- colSums(G * Rr)
        gd <- colSums(H * Rr)
        rm(Rr)
        GW <- G * W
        Iaa <- colSums(W)
        Iab <- colSums(GW)
        Iad <- colSums(H * W)
        Ibb <- colSums(G * GW)
        Ibd <- colSums(H * GW)
        Idd <- colSums(H^2 * W)
        rm(GW, W)
        C11 <- Ibb * Idd - Ibd^2
        C12 <- -(Iab * Idd - Ibd * Iad)
        C13 <- Iab * Ibd - Ibb * Iad
        C22 <- Iaa * Idd - Iad^2
        C23 <- -(Iaa * Ibd - Iab * Iad)
        C33 <- Iaa * Ibb - Iab^2
        det <- Iaa * C11 + Iab * C12 + Iad * C13
        da <- (C11 * ga + C12 * gb + C13 * gd)/det
        db <- (C12 * ga + C22 * gb + C23 * gd)/det
        dd <- (C13 * ga + C23 * gb + C33 * gd)/det
        a <- a + da
        b <- b + db
        d <- d + dd
        if (max(abs(c(da, db, dd))) < 1e-10) 
            break
    }
    stopifnot(it < 25)
    list(beta = b, se = sqrt(C22/det))
}

#' @rdname association_models
#' @export
linear_gwas <- function(Z, G) {
    if (!is.numeric(Z) || any(!is.finite(Z)) || length(Z) < 3L)
        stop("Z must be a finite numeric vector with at least three observations.", call. = FALSE)
    .check_genotype_matrix(G, length(Z))
    n <- length(Z)
    Gc <- sweep(G, 2, colMeans(G))
    Zc <- Z - mean(Z)
    sxx <- colSums(Gc^2)
    b <- drop(crossprod(Gc, Zc))/sxx
    s2 <- (sum(Zc^2) - b^2 * sxx)/(n - 2)
    list(beta = b, se = sqrt(s2/sxx))
}

# Population calibration
# probit_rho2 estimates explained liability variance for one numeric H; its
# delta-method SE treats empirical var(H) as fixed.
# fitted_variance_factor uses population data and the fitted null logistic model.
# It returns lambda and a bootstrap SE; B is the number of bootstrap samples.
# Call set.seed() before fitted_variance_factor for reproducible bootstrap results.
# These calibration fits assume a representative population sample.
# dlambda uses a central difference; both rho2 +/- h must be in the valid domain.

#' Estimate calibration quantities from population data
#'
#' Estimate covariate-explained liability variance from a probit model or estimate
#' the working-logistic factor from the fitted null-model Bernoulli variance.
#' @param Y Numeric binary outcome vector containing both zero and one.
#' @param H Numeric vector of one adjustment covariate.
#' @param B Integer number of bootstrap samples, at least two.
#' @details These functions assume a representative population sample. The probit
#'   delta-method SE treats the empirical variance of \code{H} as fixed.
#'   The fitted-variance factor uses the observed population disease fraction.
#'   It should not be applied directly to an ascertained case-control sample.
#'   Set the random seed before calling \code{fitted_variance_factor()} for
#'   reproducible bootstrap results.
#' @return A named numeric vector: \code{rho2} and \code{se} for
#'   \code{probit_rho2()}, or \code{lambda} and \code{se} for
#'   \code{fitted_variance_factor()}.
#' @examples
#' set.seed(25)
#' H <- rnorm(1000)
#' Y <- rbinom(1000, 1, pnorm(-1 + 0.4 * H))
#' probit_rho2(Y, H)
#' fitted_variance_factor(Y, H, B = 10)
#' @name population_calibration
#' @export
probit_rho2 <- function(Y, H) {
    .check_individual_data(Y, H)
    fit <- glm(Y ~ H, family = binomial(link = "probit"))
    u <- coef(fit)[2]
    s2 <- var(H)
    V <- u^2 * s2
    grad <- 2 * u * s2/(1 + V)^2
    c(rho2 = unname(V/(1 + V)), se = unname(abs(grad) * sqrt(vcov(fit)[2, 2])))
}

#' @rdname population_calibration
#' @export
fitted_variance_factor <- function(Y, H, B = 40) {
    .check_individual_data(Y, H)
    .check_count(B, "B", 2L)
    est <- function(y, h) {
        q <- glm.fit(cbind(1, h), y, family = binomial())$fitted.values
        dnorm(qnorm(1 - mean(y)))/mean(q * (1 - q))
    }
    boot <- replicate(B, {
        i <- sample.int(length(Y), replace = TRUE)
        est(Y[i], H[i])
    })
    c(lambda = est(Y, H), se = sd(boot))
}

dlambda <- function(f, rho2, h = 1e-04) (f(rho2 + h) - f(rho2 - h))/(2 * h)

# Published bias-estimation wrappers
# bx/sx: variant-covariate coefficients and SEs; by/sy: adjusted log-odds coefficients and SEs.
# These wrappers use one adjusted covariate and the identifying assumptions of each method.
# bias_dudbridge is the retained function name for indexevent CWLS (indexevent 0.2.0 used).
# bias_slopehunter uses SlopeHunter (1.1.0 used); seed controls its bootstrap.
# Both return b, b_se, A (bias-subtracted associations), and A_se.
# The external packages are required only when the corresponding wrapper is called.

#' Estimate collider-bias slopes with published methods
#'
#' Call the corrected weighted least-squares (CWLS) method in \pkg{indexevent}
#' or the clustering method in \pkg{SlopeHunter} for one adjusted covariate.
#' @param bx,sx Numeric vectors of variant-covariate associations and their SEs.
#' @param by,sy Numeric vectors of adjusted outcome log-odds associations and their SEs.
#' @param seed Integer random seed for the Slope-Hunter bootstrap.
#' @details All vectors must refer to the same aligned variants. These wrappers
#'   require the identifying assumptions of the chosen method and adjusted
#'   log-odds units appropriate to the target GWAS. \code{bias_cwls()} and
#'   \code{bias_dudbridge()} call the same CWLS implementation; the latter name
#'   is retained for compatibility. Slope-Hunter uses 100 bootstrap samples.
#'   The optional packages must be installed before calling their wrappers.
#' @return A list with \code{b} (bias slope), \code{b_se}, \code{A}
#'   (bias-subtracted coefficients), and \code{A_se}.
#' @seealso \code{\link{conversion_factors}}, \code{\link{simulate_end_to_end}}
#' @examples
#' \dontrun{
#' # Supply aligned GWAS summary vectors, with the optional package installed.
#' corrected <- bias_cwls(bx, sx, by, sy)
#' liability_effect <- corrected$A / lambda_working(0.10, 0.30)$lambda
#' }
#' @name bias_estimators
#' @export
bias_dudbridge <- function(bx, sx, by, sy) {
    .require_packages("indexevent")
    f <- indexevent::indexevent(bx, sx, by, sy, method = "CWLS")
    list(b = f$b, b_se = f$b.se, A = f$ybeta.adj, A_se = f$yse.adj)
}

#' @rdname bias_estimators
#' @export
bias_slopehunter <- function(bx, sx, by, sy, seed) {
    .require_packages("SlopeHunter")
    dat <- data.frame(SNP = paste0("snp", seq_along(bx)), BETA.incidence = bx, SE.incidence = sx, 
        BETA.prognosis = by, SE.prognosis = sy)
    out <- NULL
    invisible(capture.output(suppressWarnings(suppressMessages(out <- SlopeHunter::hunt(dat, xp_thresh = 0.001, 
        Bootstrapping = TRUE, M = 100, seed = seed, Plot = FALSE, show_adjustments = TRUE)))))
    est <- out$est[match(dat$SNP, out$est$SNP), ]
    list(b = out$b, b_se = out$bse, A = est$ybeta_adj, A_se = est$yse_adj)
}

# Univariable MR
# bx/sx: exposure coefficients and SEs; bw/sw: outcome coefficients and SEs.
# mr_ivw uses inverse outcome-variance weights and returns est and se.
# mr_raps requires mr.raps (0.4.3 used) and uses the profile-score estimator
# without overdispersion, with L2 loss and sandwich SEs.

#' Summary-statistic Mendelian randomization helpers
#'
#' Fit inverse-variance weighted (IVW) univariable or multivariable models, or
#' call the profile-score MR-RAPS estimator.
#' @param bx,sx Exposure association coefficients and their SEs.
#' @param bw,sw Outcome association coefficients and their SEs.
#' @param beta_X,se_X,beta_Y,se_Y Compatibility arguments for \code{ivw_mr()}.
#'   \code{se_X} is retained in the interface but is not used by IVW.
#'   \code{se_Y} also supplies outcome SEs for \code{mvmr_ivw()}.
#' @param B_X Numeric matrix of exposure associations, variants by exposures.
#' @param B_Y Numeric vector of outcome associations.
#' @details IVW uses inverse outcome-variance weights and treats exposure
#'   associations as fixed. Its SE is not adjusted for residual overdispersion.
#'   The multivariable exposure matrix must have full column rank. Coefficients
#'   and SEs must be harmonized and on the intended scales before use.
#'   MR-RAPS requires the optional \pkg{mr.raps} package and uses L2 loss,
#'   no overdispersion, and sandwich SEs.
#' @return \code{mr_ivw()} and \code{mr_raps()} return a named vector with
#'   \code{est} and \code{se}. \code{ivw_mr()} returns a list with \code{theta}
#'   and \code{se}. \code{mvmr_ivw()} returns \code{theta} and covariance
#'   matrix \code{cov}.
#' @examples
#' bx <- c(0.05, 0.10, 0.15, 0.20)
#' by <- c(0.02, 0.05, 0.06, 0.07)
#' mr_ivw(bx, by, rep(0.01, 4))
#' @name mr_helpers
#' @export
mr_ivw <- function(bx, bw, sw) {
    w <- 1/sw^2
    c(est = sum(w * bx * bw)/sum(w * bx^2), se = sqrt(1/sum(w * bx^2)))
}

#' @rdname mr_helpers
#' @export
mr_raps <- function(bx, sx, bw, sw) {
    .require_packages("mr.raps")
    f <- mr.raps::mr.raps.mle(bx, bw, sx, sw, over.dispersion = FALSE, loss.function = "l2", se.method = "sandwich", 
        suppress.warning = TRUE)
    c(est = f$beta.hat, se = f$beta.se)
}

# Original IVW helpers
# ivw_mr returns theta and se; se_X is retained in its interface but is not used.
# mvmr_ivw takes a matrix B_X (variants by exposures), vector B_Y, and outcome SEs se_Y.
# It returns theta and its covariance matrix cov, treating exposure estimates as fixed.

#' @rdname mr_helpers
#' @export
ivw_mr <- function(beta_X, se_X, beta_Y, se_Y) {
    w <- 1/se_Y^2
    num <- sum(w * beta_X * beta_Y)
    den <- sum(w * beta_X^2)
    theta <- num/den
    se <- sqrt(1/den)
    list(theta = theta, se = se)
}

#' @rdname mr_helpers
#' @export
mvmr_ivw <- function(B_X, B_Y, se_Y) {
    W <- diag(1/se_Y^2)
    inv <- solve(t(B_X) %*% W %*% B_X)
    list(theta = as.numeric(inv %*% t(B_X) %*% W %*% B_Y), cov = inv)
}


#' @rdname conversion_factors
#' @export
kappa_cond <- function(p_Y, rho2, n_quad = 200) {
    kappa_score(p_Y, rho2, n_quad)
}

#' @rdname bias_estimators
#' @export
bias_cwls <- function(bx, sx, by, sy) {
    bias_dudbridge(bx, sx, by, sy)
}
