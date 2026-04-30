# =============================================================================
#  realdata_mr.R
#  Real-data IVW Mendelian randomization with the conditional-prevalence
#  correction. Uses LOCAL plink clumping (no OpenGWAS API dependency).
#
#  Analysis: BMI -> T2D (BMI as exposure, BMI-adjusted T2D as outcome)
#    Yengo et al. 2018 BMI (continuous exposure)
#    Mahajan et al. 2018 BMI-adjusted T2D (binary outcome)
#
#  IVW is run on the log-OR outcome scale, then the constant-rescaling
#  identity is applied:
#    theta_marg  = theta_logOR / lambda_marg
#    theta_score = theta_logOR / lambda_cond  =  theta_marg / kappa_score
# =============================================================================
#
#  PREREQUISITES (one-time)
#  ------------------------
#  R packages:
#    install.packages(c("data.table", "remotes"))
#    remotes::install_github("MRCIEU/TwoSampleMR")
#    remotes::install_github("MRCIEU/genetics.binaRies")  # bundled plink
#    remotes::install_github("MRCIEU/ieugwasr")           # local clumping
#
#  LD reference (one-time, ~1 GB unpacked):
#    cd ld_reference
#    curl -L -o 1kg.v3.tgz 'http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz'
#    tar -xzf 1kg.v3.tgz
#    # produces EUR.bed/.bim/.fam (and AFR, EAS, AMR, SAS variants)
#
#  USAGE
#  -----
#  Step 1: edit CONFIG below to point at your local files
#  Step 2: run preflight first:
#            Rscript R/realdata/bmi_t2d_mr.R 2>&1 | tee realdata_mr_preflight.log
#  Step 3: change CONFIG$phase to "RUN" and execute:
#            Rscript R/realdata/bmi_t2d_mr.R 2>&1 | tee realdata_mr.log
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source("R/realdata/helpers.R")
source("R/sims/cbcGWAS_replication.R")

# -----------------------------------------------------------------------------
# Column-name overrides verified by preflight diagnostics
# -----------------------------------------------------------------------------

OVERRIDES_YENGO_2018_BMI <- list(
  SNP   = "SNP",
  A1    = "Tested_Allele",
  A2    = "Other_Allele",
  FRQ   = "Freq_Tested_Allele_in_HRS",
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
  P     = "Pvalue",
  N     = "Neff",
  CHR   = "Chr",
  BP    = "Pos"
)

# -----------------------------------------------------------------------------
# CONFIG
# -----------------------------------------------------------------------------
CONFIG <- list(
  phase = "PREFLIGHT",

  # ----- Input GWAS files -----
  bmi_sumstats_raw     = "data/yengo2018_bmi.txt.gz",
  bmi_overrides        = OVERRIDES_YENGO_2018_BMI,
  bmi_sumstats_norm    = "data/bmi_normalized.txt.gz",
  bmi_pthresh          = 5e-8,

  t2d_sumstats_raw     = "data/mahajan2018_t2d_bmi_adjusted.txt.gz",
  t2d_overrides        = OVERRIDES_MAHAJAN_2018_T2D_BMIADJ,
  t2d_N_constant       = NULL,
  t2d_sumstats_norm    = "data/t2d_bmiadj_normalized.txt.gz",
  t2d_pop_prev         = 0.10,

  # ----- LD reference (local plink clumping) -----
  # Path to the 1000G EUR PLINK fileset (without the .bed/.bim/.fam extension)
  ld_bfile             = "/Users/bai-zhu/Downloads/realdata/ldsc_data/1kg/EUR",

  # Path to the PLINK 1.9 binary (overrides genetics.binaRies search)
  plink_bin            = "/Users/bai-zhu/Downloads/realdata/tools/plink",

  # Clumping parameters (TwoSampleMR defaults)
  clump_r2             = 0.001,
  clump_kb             = 10000,

  # ----- rho^2 plug-in grid for T2D-on-BMI -----
  t2d_rho2_grid        = c(0.20, 0.30, 0.40, 0.49, 0.64),

  # ----- Output -----
  out_results          = "realdata_mr_results.rds"
)

# =============================================================================
#  PREFLIGHT phase
# =============================================================================

if (CONFIG$phase == "PREFLIGHT") {
  cat("\n>>> PREFLIGHT phase. Edit CONFIG$phase = 'RUN' to execute.\n\n")

  cat("Checking required files:\n")
  files_to_check <- list(
    BMI = CONFIG$bmi_sumstats_raw,
    T2D = CONFIG$t2d_sumstats_raw
  )
  for (label in names(files_to_check)) {
    f <- files_to_check[[label]]
    if (file.exists(f)) {
      cat(sprintf("  %-5s %-50s OK\n", label, f))
    } else {
      cat(sprintf("  %-5s %-50s MISSING\n", label, f))
    }
  }

  cat("\nLD reference panel:\n")
  for (ext in c(".bed", ".bim", ".fam")) {
    p <- paste0(CONFIG$ld_bfile, ext)
    if (file.exists(p)) {
      cat(sprintf("  %-50s OK (%.0f MB)\n", p,
                  file.info(p)$size / (1024 * 1024)))
    } else {
      cat(sprintf("  %-50s MISSING\n", p))
    }
  }

  cat("\nRequired R packages:\n")
  for (pkg in c("data.table", "TwoSampleMR", "ieugwasr",
                "genetics.binaRies")) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    cat(sprintf("  %-20s %s\n", pkg, if (ok) "OK" else "MISSING"))
  }

  # Test plink binary availability
  if (requireNamespace("genetics.binaRies", quietly = TRUE)) {
    plink_bin <- tryCatch(genetics.binaRies::get_plink_binary(),
                          error = function(e) NA_character_)
    if (!is.na(plink_bin) && file.exists(plink_bin)) {
      cat(sprintf("\nPLINK binary: %s OK\n", plink_bin))
    } else {
      cat("\nPLINK binary: NOT FOUND via genetics.binaRies\n")
    }
  }

  for (label in c("BMI", "T2D")) {
    f <- files_to_check[[label]]
    if (file.exists(f)) {
      cat(sprintf("\n>>> Inspecting %s:\n", label))
      sanity_check_sumstats(f)
    }
  }

  cat("\n>>> If everything above looks correct, set CONFIG$phase = 'RUN'.\n")
  quit(save = "no", status = 0)
}

# =============================================================================
#  RUN phase
# =============================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
})

# Use the configured PLINK path (avoids dependency on genetics.binaRies)
PLINK_BIN <- CONFIG$plink_bin
if (!file.exists(PLINK_BIN)) {
  stop("PLINK binary not found at CONFIG$plink_bin = ", PLINK_BIN)
}
cat(sprintf("Using PLINK at: %s\n", PLINK_BIN))
cat(sprintf("Using LD reference: %s\n", CONFIG$ld_bfile))

# -----------------------------------------------------------------------------
# Step 1: Normalize columns
# -----------------------------------------------------------------------------
cat("\n=== Step 1: normalizing columns ===\n")

normalize_if_needed <- function(raw, norm, overrides, N_constant = NULL) {
  if (file.exists(norm)) {
    cat(sprintf("  %s already exists; skipping\n", norm))
    return(invisible())
  }
  prepare_sumstats(input_file       = raw,
                   output_file      = norm,
                   column_overrides = overrides,
                   N_constant       = N_constant,
                   apply_qc         = TRUE)
}

normalize_if_needed(CONFIG$bmi_sumstats_raw, CONFIG$bmi_sumstats_norm,
                    CONFIG$bmi_overrides)
normalize_if_needed(CONFIG$t2d_sumstats_raw, CONFIG$t2d_sumstats_norm,
                    CONFIG$t2d_overrides, CONFIG$t2d_N_constant)

# -----------------------------------------------------------------------------
# Helper: load harmonized exposure-outcome instrument set (LOCAL clumping)
# -----------------------------------------------------------------------------

load_instruments_local <- function(exposure_file, outcome_file, pthresh,
                                   exposure_label, outcome_label,
                                   ld_bfile, plink_bin,
                                   clump_r2 = 0.001, clump_kb = 10000) {
  cat(sprintf("\n  loading %s -> %s instruments (p < %.0e):\n",
              exposure_label, outcome_label, pthresh))

  exposure_dat <- read_exposure_data(
    filename            = exposure_file,
    sep                 = "\t",
    snp_col             = "SNP",
    beta_col            = "BETA",
    se_col              = "SE",
    effect_allele_col   = "A1",
    other_allele_col    = "A2",
    eaf_col             = "FRQ",
    pval_col            = "P",
    samplesize_col      = "N",
    chr_col             = "CHR",
    pos_col             = "BP"
  )
  exposure_dat <- exposure_dat[exposure_dat$pval.exposure < pthresh, ]
  cat(sprintf("    %d SNPs at p < %.0e\n", nrow(exposure_dat), pthresh))

  # Local plink clumping via ieugwasr
  cat("    running local PLINK clumping (this may take a few minutes)...\n")
  clump_input <- data.frame(
    rsid = exposure_dat$SNP,
    pval = exposure_dat$pval.exposure,
    id   = "exposure"
  )
  clumped <- ld_clump(
    dat       = clump_input,
    clump_r2  = clump_r2,
    clump_kb  = clump_kb,
    plink_bin = plink_bin,
    bfile     = ld_bfile
  )
  exposure_dat <- exposure_dat[exposure_dat$SNP %in% clumped$rsid, ]
  cat(sprintf("    %d independent instruments after clumping\n",
              nrow(exposure_dat)))

  # ---------------------------------------------------------------------------
  # SNP-ID format reconciliation: BMI uses rsIDs but Mahajan 2018 BMI-adjusted
  # T2D uses chr:pos. We join the two files on chr:pos, then attach the rsID
  # from the exposure side so that harmonise_data() works on a unified key.
  # ---------------------------------------------------------------------------
  cat("    looking up exposure SNPs in outcome file by chr:pos...\n")
  outcome_raw <- data.table::fread(outcome_file)
  exp_dt <- data.table::data.table(
    SNP = exposure_dat$SNP,
    CHR = exposure_dat$chr.exposure,
    BP  = exposure_dat$pos.exposure
  )
  setkey(exp_dt, CHR, BP)
  data.table::setDT(outcome_raw)
  setkey(outcome_raw, CHR, BP)

  joined <- outcome_raw[exp_dt, on = c("CHR", "BP"), nomatch = NULL]
  joined[, SNP := i.SNP]
  joined[, i.SNP := NULL]
  cat(sprintf("    %d instruments found in outcome (chr:pos join)\n",
              nrow(joined)))

  # Write the joined outcome to a temp file in TwoSampleMR-readable format
  tmp_outcome <- tempfile(fileext = ".tsv")
  data.table::fwrite(joined, tmp_outcome, sep = "\t")

  outcome_dat <- read_outcome_data(
    snps                = exposure_dat$SNP,
    filename            = tmp_outcome,
    sep                 = "\t",
    snp_col             = "SNP",
    beta_col            = "BETA",
    se_col              = "SE",
    effect_allele_col   = "A1",
    other_allele_col    = "A2",
    eaf_col             = "FRQ",
    pval_col            = "P",
    samplesize_col      = "N"
  )
  unlink(tmp_outcome)
  cat(sprintf("    %d instruments parsed by read_outcome_data\n",
              nrow(outcome_dat)))

  harmonised <- harmonise_data(exposure_dat, outcome_dat)
  k_keep <- sum(harmonised$mr_keep)
  cat(sprintf("    %d instruments after harmonisation (mr_keep == TRUE)\n",
              k_keep))
  if (k_keep < 5) {
    warning("Fewer than 5 valid instruments; results may be unreliable")
  }
  harmonised
}

# -----------------------------------------------------------------------------
# Helper: IVW with marginal-scaled and score-corrected liability conversion
# -----------------------------------------------------------------------------
# theta_marg  = theta_logOR / lambda_marg
# theta_score = theta_logOR / lambda_cond  =  theta_marg / kappa_score

ivw_with_corrections <- function(harmonised, p_Y, rho2, label) {
  d <- harmonised[harmonised$mr_keep == TRUE, ]
  k <- nrow(d)
  beta_X <- d$beta.exposure
  beta_Y <- d$beta.outcome
  se_Y   <- d$se.outcome

  w           <- 1 / se_Y^2
  theta_logOR <- sum(w * beta_X * beta_Y) / sum(w * beta_X^2)
  se_logOR    <- sqrt(1 / sum(w * beta_X^2))

  lm <- lambda_marg(p_Y)
  lc <- lambda_cond_gauss(p_Y, rho2)
  kappa <- lc / lm

  theta_marg  <- theta_logOR / lm
  se_marg     <- se_logOR    / lm
  theta_score <- theta_logOR / lc
  se_score    <- se_logOR    / lc

  data.frame(
    label                  = label,
    rho2                   = rho2,
    n_instruments          = k,
    theta_logOR            = theta_logOR,
    se_logOR               = se_logOR,
    lambda_marg            = lm,
    lambda_cond            = lc,
    kappa_score            = kappa,
    theta_marg             = theta_marg,
    se_marg                = se_marg,
    theta_score            = theta_score,
    se_score               = se_score,
    pct_diff_score_vs_marg = 100 * (theta_score - theta_marg) / theta_marg
  )
}

# =============================================================================
#  Analysis: BMI -> T2D
# =============================================================================
cat("\n========================================================\n")
cat("Analysis: BMI -> T2D (BMI-adjusted T2D summary statistics)\n")
cat("========================================================\n")

harm <- load_instruments_local(
  exposure_file  = CONFIG$bmi_sumstats_norm,
  outcome_file   = CONFIG$t2d_sumstats_norm,
  pthresh        = CONFIG$bmi_pthresh,
  exposure_label = "BMI",
  outcome_label  = "T2D-bmi-adj",
  ld_bfile       = CONFIG$ld_bfile,
  plink_bin      = PLINK_BIN,
  clump_r2       = CONFIG$clump_r2,
  clump_kb       = CONFIG$clump_kb
)

cat("\n--- Results ---\n")
res <- do.call(rbind, lapply(CONFIG$t2d_rho2_grid, function(rho2)
  ivw_with_corrections(harm, CONFIG$t2d_pop_prev, rho2,
                       label = sprintf("BMI->T2D, rho2=%.2f", rho2))
))
print(res, row.names = FALSE, digits = 4)

# Print key numbers in a manuscript-friendly format
cat("\n=== Key numbers for the manuscript ===\n")
cat(sprintf("  Number of independent instruments: K = %d\n",
            res$n_instruments[1]))
cat(sprintf("  IVW estimate on log-OR scale:      %.4f (SE %.4f)\n",
            res$theta_logOR[1], res$se_logOR[1]))
cat(sprintf("  Marginal-scaled liability MR:      theta_marg  = %.4f (SE %.4f)\n",
            res$theta_marg[1], res$se_marg[1]))
cat(sprintf("\n  At rho^2 = 0.30 (primary plug-in):\n"))
i <- which(abs(res$rho2 - 0.30) < 1e-6)
if (length(i) == 1) {
  cat(sprintf("    kappa_score      = %.4f\n", res$kappa_score[i]))
  cat(sprintf("    theta_score      = %.4f (SE %.4f)\n",
              res$theta_score[i], res$se_score[i]))
  cat(sprintf("    pct diff vs marg = %.2f%%\n",
              res$pct_diff_score_vs_marg[i]))
}
cat(sprintf("\n  Range across rho^2 plug-ins (%s):\n",
            paste(CONFIG$t2d_rho2_grid, collapse = ", ")))
cat(sprintf("    theta_score range: %.4f to %.4f\n",
            min(res$theta_score), max(res$theta_score)))
cat(sprintf("    pct diff range:    %.2f%% to %.2f%%\n",
            min(res$pct_diff_score_vs_marg),
            max(res$pct_diff_score_vs_marg)))

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
saveRDS(list(
  config     = CONFIG,
  res        = res,
  n_instr    = sum(harm$mr_keep),
  timestamp  = Sys.time()
), file = CONFIG$out_results)

cat(sprintf("\n=== saved %s ===\n", CONFIG$out_results))
cat("Send this file (or the printed key numbers above) back to fill in §4.3.\n")
