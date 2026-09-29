#' Liability-scale conversion for covariate-adjusted binary GWAS
#'
#' Conversion factors, calibration, bias-estimation interfaces, and simulations
#' under a Gaussian liability model. Start with [lambda_working()] or
#' [run_simulations()].
#' @importFrom stats binomial coef dnorm glm glm.fit lm plogis pnorm qlogis qnorm rbinom rnorm runif sd var vcov
#' @importFrom utils capture.output head sessionInfo write.csv
#' @keywords internal
"_PACKAGE"

# These symbols are evaluated within the explicit, local DESIGN list in
# simulate_end_to_end(), through with(DESIGN, ...).
utils::globalVariables(c("a_U", "a_UH", "delta", "K_both", "K_H", "K_null", "K_Y",
                         "sd_beta", "sd_gamma", "v_zeta"))
