# =============================================================================
#  realdata_helpers.R
#  Shared utilities for the real-data pipeline:
#    - prepare_sumstats(): normalize column names from common GWAS formats
#    - sanity_check_sumstats(): inspect a file before committing to munge/MR
#    - preflight(): verify packages, reference panels, and config paths
#  =============================================================================

# data.table is needed for the I/O functions but not for inspecting the
# OVERRIDES_* dictionaries. We attempt to load it but allow source() to
# succeed (for inspection) if it's missing. The functions that need it will
# error informatively at call time.
if (requireNamespace("data.table", quietly = TRUE)) {
  suppressPackageStartupMessages(library(data.table))
} else {
  warning("data.table not installed. Functions that read files will error ",
          "until you install it: install.packages('data.table').")
}

# -----------------------------------------------------------------------------
# Column-name dictionaries for common public GWAS releases.
# Each entry maps the standard column we want (left) to a vector of possible
# names found in the wild (right).
#
# When you point prepare_sumstats() at a new file, it tries each candidate
# name in turn. If none match, it errors and prints what it found, so you
# can extend the dictionary.
# -----------------------------------------------------------------------------

COLUMN_DICTIONARY <- list(
  SNP   = c("SNP", "snp", "rsid", "variant_id", "RSID", "rs_id", "MarkerName",
            "ID", "id", "rs_number"),
  A1    = c("A1", "a1", "EFFECT_ALLELE", "effect_allele", "Allele1",
            "ALT", "alt", "Tested_Allele"),
  A2    = c("A2", "a2", "OTHER_ALLELE", "other_allele", "Allele2",
            "REF", "ref", "Other_Allele"),
  FRQ   = c("FRQ", "frq", "EAF", "eaf", "MAF", "maf",
            "effect_allele_frequency", "Freq_Tested_Allele",
            "FreqAllele1HapMapCEU", "Freq.Allele1.HapMapCEU"),
  BETA  = c("BETA", "beta", "Effect", "effect", "logOR", "log_OR",
            "b", "Z"),
  OR    = c("OR", "or", "odds_ratio", "OddsRatio"),
  SE    = c("SE", "se", "StdErr", "stderr", "standard_error", "SD"),
  P     = c("P", "p", "PVAL", "pval", "p_value", "P_VALUE", "Pvalue", "P-value"),
  N     = c("N", "n", "samplesize", "SampleSize", "sample_size",
            "N_total", "TotalN", "Neff", "N_eff", "neff"),
  CHR   = c("CHR", "chr", "Chromosome", "chromosome", "#CHROM", "chrom"),
  BP    = c("BP", "bp", "Position", "POSITION", "POS", "pos",
            "base_pair_location")
)

# -----------------------------------------------------------------------------
# Sanity-check a raw summary-statistics file:
# - Read first 1000 rows
# - Print column names with mapping suggestions
# - Show numeric ranges for typical sanity flags (P in [0,1], MAF/FRQ in
#   [0,1], BETA centered near 0, etc.)
# -----------------------------------------------------------------------------

sanity_check_sumstats <- function(filepath, sep = "\t", n_preview = 1000) {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  cat(sprintf("=== Sanity-checking %s ===\n", filepath))
  cat(sprintf("File size: %.1f MB\n",
              file.info(filepath)$size / (1024 * 1024)))

  # Read a preview
  preview <- tryCatch({
    fread(filepath, sep = sep, nrows = n_preview, header = TRUE)
  }, error = function(e) {
    cat("ERROR reading file:", conditionMessage(e), "\n")
    cat("Try a different separator (',' for csv, ' ' for whitespace).\n")
    return(NULL)
  })
  if (is.null(preview)) return(invisible(NULL))

  cat(sprintf("Columns (%d): %s\n", ncol(preview),
              paste(colnames(preview), collapse = ", ")))
  cat(sprintf("Preview rows: %d\n", nrow(preview)))

  # Try to detect each standard column
  cat("\n--- Column mapping (auto-detected) ---\n")
  detected <- list()
  for (std_name in names(COLUMN_DICTIONARY)) {
    candidates <- COLUMN_DICTIONARY[[std_name]]
    found <- intersect(candidates, colnames(preview))
    if (length(found) > 0) {
      detected[[std_name]] <- found[1]
      cat(sprintf("  %-6s -> %s\n", std_name, found[1]))
    } else {
      cat(sprintf("  %-6s -> NOT FOUND (tried: %s)\n",
                  std_name, paste(candidates[1:3], collapse = ", ")))
    }
  }

  # Sanity-check numeric ranges
  cat("\n--- Numeric range sanity checks ---\n")
  if (!is.null(detected$P)) {
    p_vals <- as.numeric(preview[[detected$P]])
    p_vals <- p_vals[is.finite(p_vals)]
    if (length(p_vals) > 0) {
      cat(sprintf("  P range:    [%.2e, %.2e]  (expect [0, 1])\n",
                  min(p_vals), max(p_vals)))
      if (max(p_vals) > 1.01) {
        cat("    WARNING: P > 1 detected. Column may not be a p-value.\n")
      }
    }
  }
  if (!is.null(detected$FRQ)) {
    frq_vals <- as.numeric(preview[[detected$FRQ]])
    frq_vals <- frq_vals[is.finite(frq_vals)]
    if (length(frq_vals) > 0) {
      cat(sprintf("  FRQ range:  [%.3f, %.3f]  (expect [0, 1])\n",
                  min(frq_vals), max(frq_vals)))
    }
  }
  if (!is.null(detected$BETA)) {
    b_vals <- as.numeric(preview[[detected$BETA]])
    b_vals <- b_vals[is.finite(b_vals)]
    if (length(b_vals) > 0) {
      cat(sprintf("  BETA range: [%.3f, %.3f]  median=%.4f\n",
                  min(b_vals), max(b_vals), median(b_vals)))
      if (median(b_vals) > 0.5) {
        cat("    WARNING: median |BETA| seems large. Could be an OR column.\n")
      }
    }
  }
  if (!is.null(detected$OR) && is.null(detected$BETA)) {
    or_vals <- as.numeric(preview[[detected$OR]])
    or_vals <- or_vals[is.finite(or_vals)]
    if (length(or_vals) > 0) {
      cat(sprintf("  OR range:   [%.3f, %.3f]  median=%.3f\n",
                  min(or_vals), max(or_vals), median(or_vals)))
      cat("    Note: BETA absent but OR present. Pipeline will use log(OR).\n")
    }
  }
  if (!is.null(detected$N)) {
    n_vals <- as.numeric(preview[[detected$N]])
    n_vals <- n_vals[is.finite(n_vals)]
    if (length(n_vals) > 0) {
      cat(sprintf("  N range:    [%d, %d]\n",
                  as.integer(min(n_vals)), as.integer(max(n_vals))))
    }
  } else {
    cat("  N: NOT FOUND -- you must pass N as an argument to munge()\n")
  }

  invisible(detected)
}

# -----------------------------------------------------------------------------
# Prepare a normalized summary-statistics file from a raw GWAS download.
# Reads the raw file, renames columns to standard names (SNP, A1, A2, FRQ,
# BETA, SE, P, N, optional CHR, BP), converts OR->BETA if necessary, applies
# basic QC filters, and writes a tab-separated gz-compressed output.
# -----------------------------------------------------------------------------

prepare_sumstats <- function(input_file, output_file,
                             column_overrides = list(),
                             N_constant       = NULL,
                             apply_qc         = TRUE,
                             maf_threshold    = 0.01,
                             info_threshold   = 0.9,
                             info_col         = NULL,
                             sep              = "\t") {
  cat(sprintf("=== prepare_sumstats: %s -> %s ===\n", input_file, output_file))
  if (!file.exists(input_file)) stop("Input file not found: ", input_file)

  # Read full file
  d <- fread(input_file, sep = sep, header = TRUE,
             showProgress = TRUE)
  cat(sprintf("Read %d rows, %d cols\n", nrow(d), ncol(d)))

  # Auto-detect columns, with overrides
  detected <- list()
  for (std_name in names(COLUMN_DICTIONARY)) {
    if (!is.null(column_overrides[[std_name]])) {
      detected[[std_name]] <- column_overrides[[std_name]]
    } else {
      candidates <- COLUMN_DICTIONARY[[std_name]]
      found <- intersect(candidates, colnames(d))
      if (length(found) > 0) detected[[std_name]] <- found[1]
    }
  }

  # Required columns
  required <- c("SNP", "A1", "A2", "P")
  missing_required <- setdiff(required, names(detected))
  if (length(missing_required) > 0) {
    stop("Missing required columns: ",
         paste(missing_required, collapse = ", "),
         "\nPass them via column_overrides = list(SNP = 'your_col_name', ...)")
  }

  # Need at least BETA or OR
  if (is.null(detected$BETA) && is.null(detected$OR)) {
    stop("Neither BETA nor OR found. Pass via column_overrides.")
  }

  # Build output
  out <- data.table(
    SNP  = as.character(d[[detected$SNP]]),
    A1   = toupper(as.character(d[[detected$A1]])),
    A2   = toupper(as.character(d[[detected$A2]])),
    P    = as.numeric(d[[detected$P]])
  )

  if (!is.null(detected$BETA)) {
    out$BETA <- as.numeric(d[[detected$BETA]])
  } else {
    cat("  Converting OR -> BETA via log(OR)\n")
    out$BETA <- log(as.numeric(d[[detected$OR]]))
  }

  if (!is.null(detected$SE)) {
    out$SE <- as.numeric(d[[detected$SE]])
  } else {
    stop("SE column required. Pass via column_overrides$SE.")
  }

  if (!is.null(detected$FRQ)) {
    out$FRQ <- as.numeric(d[[detected$FRQ]])
  } else {
    cat("  WARNING: no FRQ/EAF column found. Setting FRQ = NA.\n")
    cat("  This may cause issues with TwoSampleMR clump_data().\n")
    out$FRQ <- NA_real_
  }

  if (!is.null(detected$N)) {
    out$N <- as.numeric(d[[detected$N]])
  } else if (!is.null(N_constant)) {
    cat(sprintf("  Setting constant N = %.0f from N_constant argument\n",
                N_constant))
    out$N <- N_constant
  } else {
    stop("No N column found and N_constant not provided.")
  }

  if (!is.null(detected$CHR)) out$CHR <- as.character(d[[detected$CHR]])
  if (!is.null(detected$BP))  out$BP  <- as.numeric(d[[detected$BP]])

  cat(sprintf("Built standardized table: %d rows\n", nrow(out)))

  # QC filters
  if (apply_qc) {
    n_before <- nrow(out)
    out <- out[is.finite(BETA) & is.finite(SE) & is.finite(P) & SE > 0]
    out <- out[A1 %in% c("A", "C", "G", "T") & A2 %in% c("A", "C", "G", "T")]
    out <- out[!duplicated(SNP)]
    if (!all(is.na(out$FRQ))) {
      out <- out[is.na(FRQ) | (FRQ >= maf_threshold & FRQ <= 1 - maf_threshold)]
    }
    if (!is.null(info_col) && info_col %in% colnames(d)) {
      info_vals <- as.numeric(d[[info_col]])[match(out$SNP,
                                                   d[[detected$SNP]])]
      out <- out[info_vals >= info_threshold | is.na(info_vals)]
    }
    cat(sprintf("After QC: %d rows (%.1f%% kept)\n",
                nrow(out), 100 * nrow(out) / n_before))
  }

  # Write
  fwrite(out, output_file, sep = "\t", compress = "gzip")
  cat(sprintf("Wrote %s\n", output_file))
  invisible(out)
}

# -----------------------------------------------------------------------------
# Preflight check: verify packages, reference panels, and configured paths
# before launching the heavy pipelines.
# -----------------------------------------------------------------------------

preflight <- function(ld_path, hm3_snplist,
                      sumstats_files = list(),
                      check_packages = TRUE) {
  cat("=== Preflight check ===\n\n")
  all_ok <- TRUE

  if (check_packages) {
    cat("Required R packages:\n")
    for (pkg in c("data.table", "GenomicSEM", "TwoSampleMR")) {
      ok <- requireNamespace(pkg, quietly = TRUE)
      cat(sprintf("  %-15s %s\n", pkg, if (ok) "OK" else "MISSING"))
      if (!ok) all_ok <- FALSE
    }
    cat("\n")
  }

  cat("LDSC reference files:\n")
  if (dir.exists(ld_path)) {
    n_files <- length(list.files(ld_path, pattern = "\\.gz$"))
    cat(sprintf("  %-40s OK (%d .gz files)\n", ld_path, n_files))
    if (n_files < 22) {
      cat("    WARNING: expected 22+ chromosome-level files; got", n_files, "\n")
      all_ok <- FALSE
    }
  } else {
    cat(sprintf("  %-40s MISSING\n", ld_path))
    all_ok <- FALSE
  }
  if (file.exists(hm3_snplist)) {
    n_lines <- length(readLines(hm3_snplist, n = 5))
    cat(sprintf("  %-40s OK\n", hm3_snplist))
  } else {
    cat(sprintf("  %-40s MISSING\n", hm3_snplist))
    all_ok <- FALSE
  }
  cat("\n")

  if (length(sumstats_files) > 0) {
    cat("Summary-statistics files:\n")
    for (label in names(sumstats_files)) {
      f <- sumstats_files[[label]]
      if (file.exists(f)) {
        size_mb <- file.info(f)$size / (1024 * 1024)
        cat(sprintf("  %-12s %-40s OK (%.0f MB)\n", label, f, size_mb))
      } else {
        cat(sprintf("  %-12s %-40s MISSING\n", label, f))
        all_ok <- FALSE
      }
    }
    cat("\n")
  }

  if (all_ok) {
    cat("=== ALL CHECKS PASSED -- ready to run pipelines ===\n")
  } else {
    cat("=== CHECKS FAILED -- fix the issues above before running ===\n")
  }
  invisible(all_ok)
}

# -----------------------------------------------------------------------------
# Quick reference: column overrides for common public GWAS releases
# Copy/paste into your pipeline scripts to avoid auto-detection failures
# -----------------------------------------------------------------------------

OVERRIDES_WIGHTMAN_2021 <- list(
  SNP   = "SNP",
  A1    = "EA",
  A2    = "NEA",
  FRQ   = "EAF",
  BETA  = "BETA",
  SE    = "SE",
  P     = "P",
  N     = "N",
  CHR   = "CHR",
  BP    = "BP"
)

OVERRIDES_BELLENGUEZ_2022 <- list(
  SNP   = "variant_id",
  A1    = "effect_allele",
  A2    = "other_allele",
  FRQ   = "effect_allele_frequency",
  BETA  = "beta",
  SE    = "standard_error",
  P     = "p_value",
  N     = "n_total",
  CHR   = "chromosome",
  BP    = "base_pair_location"
)

OVERRIDES_ARAGAM_2022_CAD <- list(
  SNP   = "MarkerName",
  A1    = "Allele1",
  A2    = "Allele2",
  FRQ   = "Freq1",
  BETA  = "Effect",
  SE    = "StdErr",
  P     = "P-value"
  # N: pass as N_constant = 1165690 (cases + controls)
)

OVERRIDES_YENGO_2018_BMI <- list(
  SNP   = "SNP",
  A1    = "Tested_Allele",
  A2    = "Other_Allele",
  FRQ   = "Freq_Tested_Allele",
  BETA  = "BETA",
  SE    = "SE",
  P     = "P",
  N     = "N",
  CHR   = "CHR",
  BP    = "POS"
)

OVERRIDES_MAHAJAN_2018_T2D_BMIADJ <- list(
  SNP   = "SNP",
  A1    = "EA",
  A2    = "NEA",
  FRQ   = "EAF",
  BETA  = "Beta",
  SE    = "SE",
  P     = "Pvalue"
  # Mahajan files vary; check column names with sanity_check_sumstats() first
)

cat("Loaded realdata_helpers.R\n")
cat("Run sanity_check_sumstats('/path/to/file.txt.gz') to inspect a raw file\n")
cat("Run preflight(ld_path = '...', hm3_snplist = '...') before the pipelines\n")
