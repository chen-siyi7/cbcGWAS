test_that("conditional approximation agrees with independent adaptive integration", {
    for (p in c(0.01, 0.1, 0.5)) for (rho2 in c(0.05, 0.3, 0.5)) {
        risk_variance <- function(z) {
            risk <- pnorm((sqrt(rho2) * z - qnorm(1 - p)) / sqrt(1 - rho2))
            risk * (1 - risk) * dnorm(z)
        }
        expected <- dnorm(qnorm(1 - p)) /
            integrate(risk_variance, -Inf, Inf, rel.tol = 1e-10)$value
        expect_equal(lambda_approx(p, rho2), expected, tolerance = 1e-8)
        expect_equal(lambda_cond_gauss(p, rho2), expected, tolerance = 1e-8)
    }
})

test_that("zero covariate variance gives the marginal factor under either sampling scheme", {
    for (p in c(0.01, 0.1, 0.5)) {
        expected <- dnorm(qnorm(1 - p)) / (p * (1 - p))
        expect_equal(lambda_marginal(p), expected)
        expect_equal(lambda_approx(p, 0), expected)
        expect_equal(kappa_cond(p, 0), 1)
        for (pi in c(0.05, 0.5, 0.8))
            expect_equal(lambda_working(p, 0, pi)$lambda, expected)
    }
})

test_that("working factors retain the reviewed population and case-control values", {
    # Reference values from the reviewed S04 calculation, before packaging.
    expected <- c(2.23664657486938, 2.25726187446102, 2.31313966427218, 2.35862766784216)
    for (i in seq_along(expected)) {
        fit <- lambda_working(0.1, 0.3, pi = c(.05, .10, .30, .50)[i])
        expect_equal(fit$lambda, expected[i], tolerance = 1e-10)
        expect_lt(fit$score_error, 1e-9)
    }
    expect_equal(lambda_marg(c(.05, .1)), lambda_marginal(c(.05, .1)))
    expect_equal(kappa_score(.1, .3), kappa_cond(.1, .3))
})

test_that("invalid factor inputs fail clearly", {
    for (p in list(0, 1, NA_real_, Inf, "0.1"))
        expect_error(lambda_marginal(p), "probabilities")
    for (rho2 in list(-.1, 1, NA_real_, c(.2, .3)))
        expect_error(lambda_working(.1, rho2), "rho2")
    expect_error(lambda_approx(c(.1, .2), .3), "probabilities")
    expect_error(lambda_working(.1, .3, pi = 1), "probabilities")
    expect_error(lambda_approx(.1, .3, n = 1.5), "integer")
})
