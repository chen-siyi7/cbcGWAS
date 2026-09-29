test_that("smoke simulations are reproducible and preserve the caller's RNG state", {
    set.seed(741)
    seed <- .Random.seed
    kind <- RNGkind()
    x <- run_simulations("smoke", replicates = 2)
    expect_identical(.Random.seed, seed)
    expect_identical(RNGkind(), kind)
    expect_length(x, 9)
    expect_equal(vapply(x, nrow, integer(1)),
                 c(population_factors = 30L, casecontrol_factors = 4L,
                   no_collider_simulation = 8L, casecontrol_simulation = 4L,
                   explicit_collider_simulation = 18L, t2d_cad_simulation = 9L,
                   mvmr_simulation = 6L, multiple_covariate_simulation = 3L,
                   bias_slope_regressions = 9L))
    for (table in x) for (column in table[vapply(table, is.numeric, logical(1))])
        expect_true(all(is.finite(column)))
    # A different incoming RNG kind must not change the experiment.
    on.exit(do.call(RNGkind, as.list(kind)), add = TRUE)
    RNGkind("L'Ecuyer-CMRG")
    y <- run_simulations("smoke", replicates = 2)
    expect_identical(x, y)
})

test_that("factor tables are written only to the requested directory", {
    out <- tempfile("cbcGWAS-test-")
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    x <- run_simulations("factors", output_dir = out)
    expect_setequal(list.files(out), c("population_factors.csv", "casecontrol_factors.csv",
                                      "session_factors.txt", "run_factors.txt"))
    expect_equal(as.matrix(read.csv(file.path(out, "population_factors.csv"))),
                 as.matrix(x$population_factors),
                 tolerance = 1e-12)
    expect_error(run_simulations("unknown"), "arg")
    expect_error(run_simulations(replicates = 1), "integer")
    expect_error(run_simulations(cores = 0), "integer")
    expect_error(run_simulations(output_dir = NA_character_), "output_dir")
    expect_error(simulate_collider(R = 1), "integer")
    expect_error(simulate_end_to_end(R = 0), "integer")
})
