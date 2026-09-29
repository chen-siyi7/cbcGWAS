.check_count <- function(x, name, minimum = 1L) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x < minimum || x > .Machine$integer.max || x != floor(x))
        stop(name, " must be an integer >= ", minimum, ".", call. = FALSE)
    invisible(x)
}

.check_probability <- function(x, name, scalar = TRUE) {
    if (!is.numeric(x) || !length(x) || (scalar && length(x) != 1L) ||
        any(!is.finite(x)) || any(x <= 0 | x >= 1))
        stop(name, " must contain finite probabilities strictly between zero and one.", call. = FALSE)
    invisible(x)
}

.check_rho2 <- function(x) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x >= 1)
        stop("rho2 must be a finite scalar in [0, 1).", call. = FALSE)
    invisible(x)
}

.check_individual_data <- function(Y, H) {
    if (!is.numeric(Y) || any(!is.finite(Y)) || !all(Y %in% c(0, 1)) ||
        length(unique(Y)) != 2L || length(Y) < 4L)
        stop("Y must be a complete numeric binary vector containing both zero and one.", call. = FALSE)
    if (!is.numeric(H) || length(H) != length(Y) || any(!is.finite(H)) || var(H) <= 0)
        stop("H must be a finite, nonconstant numeric vector of the same length as Y.", call. = FALSE)
    invisible(NULL)
}

.check_genotype_matrix <- function(G, n) {
    if (!is.matrix(G) || !is.numeric(G) || nrow(G) != n || ncol(G) < 1L ||
        any(!is.finite(G)) || any(apply(G, 2, var) <= 0))
        stop("G must be a finite numeric matrix with matching rows and nonconstant columns.", call. = FALSE)
    invisible(NULL)
}

.require_packages <- function(packages) {
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop("Install the optional package(s) first: ",
                             paste(missing, collapse = ", "), ".", call. = FALSE)
    invisible(NULL)
}

.quiet_call <- function(fun, ...) {
    invisible(capture.output(result <- fun(...)))
    result
}

.preserve_rng <- function() {
    kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    function() {
        do.call(RNGkind, as.list(kind))
        if (had_seed) assign(".Random.seed", seed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
        invisible(NULL)
    }
}
