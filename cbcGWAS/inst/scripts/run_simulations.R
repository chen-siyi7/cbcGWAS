# Command-line interface to the installed cbcGWAS package.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || args[1] %in% c("--help", "-h", "help")) {
    cat("Usage: Rscript run_simulations.R MODE [REPLICATES] [OUTPUT_DIR] [CORES]\n",
        "Modes: smoke, population, casecontrol, collider, mr, multicovariate,\n",
        "       bias, factors, end-to-end, all\n",
        "Use 'default' for the original replicate count. Start with 'smoke'.\n",
        "Example: Rscript run_simulations.R smoke default smoke_results\n", sep = "")
    quit(status = 0)
}
if (length(args) > 4L) stop("Expected at most four arguments. Use --help.")
replicates <- if (length(args) >= 2L && args[2] != "default") as.numeric(args[2]) else NULL
output_dir <- if (length(args) >= 3L) args[3] else if (args[1] == "smoke") "smoke_results" else "results"
cores <- if (length(args) >= 4L) as.numeric(args[4]) else 1L
result <- cbcGWAS::run_simulations(args[1], replicates, output_dir, cores)
cat("Wrote", length(result), "tables to", normalizePath(output_dir), "\n")
