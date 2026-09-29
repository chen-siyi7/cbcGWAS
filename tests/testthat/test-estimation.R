test_that("variant-wise associations agree with independent glm and lm fits", {
    set.seed(421)
    G <- matrix(rnorm(1800), ncol = 3)
    H <- rnorm(nrow(G))
    Y <- rbinom(nrow(G), 1, plogis(-1 + .2 * G[, 1] + .4 * H))
    Z <- .2 * G[, 1] + H
    logistic <- logistic_gwas(Y, G, H)
    linear <- linear_gwas(Z, G)
    for (j in seq_len(ncol(G))) {
        logfit <- glm(Y ~ G[, j] + H, family = binomial())
        linfit <- lm(Z ~ G[, j])
        expect_equal(logistic$beta[j], unname(coef(logfit)[2]), tolerance = 1e-7)
        expect_equal(logistic$se[j], unname(sqrt(vcov(logfit)[2, 2])), tolerance = 1e-6)
        expect_equal(linear$beta[j], unname(coef(linfit)[2]), tolerance = 1e-10)
        expect_equal(linear$se[j], unname(sqrt(vcov(linfit)[2, 2])), tolerance = 1e-10)
    }
    expect_error(logistic_gwas(rep(0, 600), G, H), "binary")
    expect_error(linear_gwas(Z, G[-1, ]), "matching rows")
    expect_error(logistic_gwas(Y, G, rep(1, 600)), "nonconstant")
})

test_that("population calibration agrees with its defining regression", {
    set.seed(501)
    H <- rnorm(1000)
    Y <- rbinom(1000, 1, pnorm(-1 + .5 * H))
    fit <- glm(Y ~ H, family = binomial(link = "probit"))
    explained <- unname(coef(fit)[2]^2 * var(H))
    expect_equal(unname(probit_rho2(Y, H)["rho2"]), explained / (1 + explained))
    null <- glm(Y ~ H, family = binomial())$fitted.values
    expected <- dnorm(qnorm(1 - mean(Y))) / mean(null * (1 - null))
    set.seed(19)
    actual <- fitted_variance_factor(Y, H, B = 5)
    expect_equal(unname(actual["lambda"]), expected)
    expect_true(is.finite(actual["se"]) && actual["se"] > 0)
    expect_error(fitted_variance_factor(Y, H, B = 1), "integer")
})

test_that("MR estimators recover exact linear relations and obey scale transformations", {
    bx <- c(.05, .1, .15, .2)
    by <- .4 * bx
    se <- rep(.01, 4)
    fit <- mr_ivw(bx, by, se)
    expect_equal(unname(fit["est"]), .4)
    expect_equal(mr_ivw(bx / 2, by, se), 2 * fit)
    legacy <- ivw_mr(bx, se, by, se)
    expect_equal(legacy$theta, unname(fit["est"]))
    expect_equal(legacy$se, unname(fit["se"]))
    B <- cbind(bx, c(.2, .05, -.1, .1))
    theta <- c(.4, -.3)
    expect_equal(mvmr_ivw(B, drop(B %*% theta), se)$theta, theta)
})
