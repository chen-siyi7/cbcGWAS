#' Run the conversion and Mendelian randomization simulations
#'
#' Run one experiment, the complete set, or a small installation check.
#' Results are returned as data frames and can optionally be saved as CSV files.
#' @param mode Experiment to run. The default, \code{"smoke"}, uses three
#'   replicates and smaller individual-level samples. See Details.
#' @param replicates Optional integer number of replicates, at least two.
#'   \code{NULL} uses the experiment's original settings, or three in smoke mode.
#' @param output_dir Optional directory for CSV files, run settings, and session
#'   information. \code{NULL} writes no files. Existing files of the same names
#'   in this directory are replaced.
#' @param cores Positive integer worker count for the end-to-end experiment.
#'   Use one on Windows.
#' @details
#' The modes are:
#' \describe{
#'   \item{\code{factors}}{Population and case-control conversion factors.}
#'   \item{\code{population}}{A no-collider experiment; 500 replicates of 20,000 individuals.}
#'   \item{\code{casecontrol}}{Four sample case fractions; 200 replicates of 20,000 individuals.}
#'   \item{\code{collider}}{The corrected unit-variance collider model; 300 replicates of 30,000 individuals.}
#'   \item{\code{mr}}{Univariable and multivariable summary-statistic experiments;
#'     100 and 500 replicates, respectively. The disease labels describe simulated traits.}
#'   \item{\code{multicovariate}}{One, three, and five covariates; 300 replicates of 20,000 individuals.}
#'   \item{\code{bias}}{Residual-trimmed and precision-weighted bias-slope regressions;
#'     200 replicates. These are illustrative regressions, not the published CWLS
#'     or Slope-Hunter estimators.}
#'   \item{\code{end-to-end}}{Estimated CWLS and Slope-Hunter correction and MR;
#'     200 full-size replicates per sampling design. See \code{simulate_end_to_end()}.}
#'   \item{\code{all}}{All experiments, including end-to-end validation.}
#'   \item{\code{smoke}}{All modes except end-to-end, using three replicates and
#'     2,000 individuals in the individual-level experiments. This checks execution.}
#' }
#' The optional \pkg{indexevent}, \pkg{SlopeHunter}, and \pkg{mr.raps} packages
#' are needed for \code{"end-to-end"} and \code{"all"}. No packages or data
#' are downloaded by this function. Full runs can take hours on one core.
#'
#' Simulations use fixed seeds and restore the caller's RNG state. Column names
#' containing \code{score} retain their original meaning: the conditional-prevalence
#' approximation, not the working-logistic solution. Monte Carlo results depend
#' on replicate count and dependency versions. Numeric output is unrounded.
#' @return A named list of data frames. End-to-end mode returns both replicate
#'   data and a summary; smoke mode returns nine tables.
#' @examples
#' factors <- run_simulations("factors")
#' head(factors$population_factors)
#' \dontrun{
#' small_run <- run_simulations("smoke", output_dir = "smoke_results")
#' validation <- run_simulations("end-to-end", output_dir = "results", cores = 1)
#' }
#' @seealso \code{\link{simulate_collider}}, \code{\link{simulate_end_to_end}}
#' @export
run_simulations <- function(mode = c("smoke", "population", "casecontrol", "collider", "mr",
                                     "multicovariate", "bias", "factors", "end-to-end", "all"),
                            replicates = NULL, output_dir = NULL, cores = 1L) {
    mode <- match.arg(mode)
    if (!is.null(replicates)) .check_count(replicates, "replicates", 2L)
    .check_count(cores, "cores")
    if (!is.null(output_dir) && (!is.character(output_dir) || length(output_dir) != 1L ||
        is.na(output_dir) || !nzchar(output_dir)))
        stop("output_dir must be NULL or one nonempty directory path.", call. = FALSE)
    if (mode %in% c("all", "end-to-end")) {
        .require_packages(c("indexevent", "SlopeHunter", "mr.raps"))
        if (.Platform$OS.type == "windows" && cores > 1L)
            stop("Use cores = 1 on Windows.", call. = FALSE)
    }
    if (!is.null(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        if (!dir.exists(output_dir)) stop("Cannot create output_dir.", call. = FALSE)
    }
    restore_rng <- .preserve_rng()
    on.exit(restore_rng(), add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    nr <- function(default) if (!is.null(replicates)) as.integer(replicates)
                           else if (mode == "smoke") 3L else default
    nn <- function(default) if (mode == "smoke") 2000L else default
    run <- function(which) mode %in% c("smoke", "all", which)
    results <- list()
    if (run("factors")) {
        grid <- expand.grid(p = c(.50, .20, .10, .05, .01),
                            rho2 = c(.05, .10, .20, .30, .40, .50))
        grid$lambda_approx <- mapply(lambda_approx, grid$p, grid$rho2)
        grid$kappa <- grid$lambda_approx / lambda_marginal(grid$p)
        grid$lambda_working <- mapply(function(p, r) lambda_working(p, r)$lambda,
                                      grid$p, grid$rho2)
        grid$gap_percent <- 100 * (grid$lambda_approx / grid$lambda_working - 1)
        results$population_factors <- grid
        cc <- expand.grid(p = .10, rho2 = .30, pi = c(.05, .10, .30, .50))
        cc$lambda_working <- mapply(function(p, r, pi) lambda_working(p, r, pi)$lambda,
                                   cc$p, cc$rho2, cc$pi)
        results$casecontrol_factors <- cc
    }
    if (run("population")) results$no_collider_simulation <-
        .quiet_call(table_bias_coverage, R = nr(500L), n = nn(20000L))
    if (run("casecontrol")) results$casecontrol_simulation <-
        .quiet_call(table_cc, R = nr(200L), n = nn(20000L))
    if (run("collider")) results$explicit_collider_simulation <-
        simulate_collider(R = nr(300L), n = nn(30000L))
    if (run("mr")) {
        results$t2d_cad_simulation <- .quiet_call(run_t2d_cad, R = nr(100L))
        results$mvmr_simulation <- .quiet_call(table_mvmr, R = nr(500L))
    }
    if (run("multicovariate")) results$multiple_covariate_simulation <-
        .quiet_call(table_multi_covar, R = nr(300L), n = nn(20000L))
    if (run("bias")) results$bias_slope_regressions <-
        .quiet_call(table_bias_slope_recovery, R = nr(200L))
    if (mode %in% c("all", "end-to-end")) {
        results$end_to_end_replicates <- simulate_end_to_end(R = nr(200L), cores = cores)
        results$end_to_end_summary <- summarize_end_to_end(results$end_to_end_replicates)
    }
    if (!is.null(output_dir)) {
        for (name in names(results))
            write.csv(results[[name]], file.path(output_dir, paste0(name, ".csv")), row.names = FALSE)
        writeLines(capture.output(sessionInfo()), file.path(output_dir, paste0("session_", mode, ".txt")))
        writeLines(c(paste("mode:", mode),
                     paste("replicates:", if (is.null(replicates)) "default" else replicates),
                     paste("cores:", cores)), file.path(output_dir, paste0("run_", mode, ".txt")))
    }
    results
}
