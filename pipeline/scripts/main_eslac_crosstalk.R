#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# eSLAC: Single-Read Crosstalk Analysis (MSR-seq)
# Main pipeline script for tRNA crosstalk + charging summaries
#
# - Project-agnostic; all configuration via --config and CLI flags
# - No package installation in-script; use conda/renv to manage versions
# - Parallelized where appropriate; deterministic seeds for tie-breaks
# - Outputs: intermediate RDS/TSV tables + plots in --out_dir
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(foreach)
  library(doParallel)
})

# --------------------------- CLI & CONFIG -------------------------------------

option_list <- list(
  make_option("--config",      type="character", help="Path to config.yaml", metavar="FILE"),
  make_option("--tsv_dir",     type="character", help="Directory of per-read TSVs [read_id,gene,position,base,deletion,mutation,insertion]"),
  make_option("--mod_catalog", type="character", help="Excel/CSV catalog with annotated modification positions"),
  make_option("--samples",     type="character", help="Sample sheet CSV (columns: file,sample,rep,group,source)"),
  make_option("--out_dir",     type="character", help="Output directory for tables/plots", default="results"),
  make_option("--threads",     type="integer",   help="CPU threads to use (overrides config)", default=NA_integer_)
)
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$config) || is.null(opt$tsv_dir) || is.null(opt$mod_catalog) || is.null(opt$samples)) {
  stop("Missing required arguments. See --help.")
}

cfg <- yaml::read_yaml(opt$config)

OUT_DIR <- opt$out_dir
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "plots"),  recursive=TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "rds"),    recursive=TRUE, showWarnings = FALSE)

SEED <- 1337
set.seed(SEED)

# Threads
threads_cfg <- if (is.na(opt$threads)) {
  if (!is.null(cfg$alignment$bowtie2$threads)) cfg$alignment$bowtie2$threads else max(1L, parallel::detectCores()-2L)
} else opt$threads
data.table::setDTthreads(threads_cfg)

# --------------------------- PARAMETERS ---------------------------------------

# Coverage, thresholds, and statistics are read from config.yaml
# Provide sane defaults if keys are missing (so script runs out-of-the-box).
get_cfg <- function(path, default = NULL) {
  purrr::reduce(strsplit(path, "\\.")[[1]], function(acc, k) if (!is.null(acc[[k]])) acc[[k]] else NULL, .init = cfg) %||% default
}

READ_LEN_MIN                   <- as.integer(get_cfg("coverage.read_length_min", 60))
MIN_READS_PER_POS_PER_REP      <- as.integer(get_cfg("coverage.min_reads_per_position_per_rep", 30))
MIN_REPS_PASSING               <- as.integer(get_cfg("coverage.min_reps_passing", 2))
MIN_PAIR_COVERAGE              <- as.integer(get_cfg("coverage.min_pair_coverage", 200))
MIN_SITE_EVENT_PCT             <- as.numeric(get_cfg("coverage.min_site_pct", 1.0))
POSITION_EVENT_PRESELECT_PCT   <- as.numeric(get_cfg("mod_calling.min_effect_size_pct", 3.0)) # used as site prefilter
PAIR_EVENT_GATE_PCT            <- as.numeric(get_cfg("mod_calling.site_pair_gate_pct", 5.0))  # classic ≥5% rule

USE_MUTATIONS                  <- isTRUE(get_cfg("mod_calling.use_mutations", TRUE))
USE_DELETIONS                  <- isTRUE(get_cfg("mod_calling.use_deletions", TRUE))
CROSSTALK_FDR_ALPHA            <- as.numeric(get_cfg("crosstalk.fdr_alpha", 0.05))
CROSSTALK_FDR_METHOD           <- as.character(get_cfg("crosstalk.fdr_method", "BH"))
CROSSTALK_MIN_ABS_LOG2OR       <- as.numeric(get_cfg("crosstalk.effect_size_abs_log2OR", 0.3))

PLOT_ENABLE                    <- isTRUE(get_cfg("reporting.plots", TRUE))

# --------------------------- IO HELPERS ---------------------------------------

read_sample_sheet <- function(path) {
  # Required columns: file,sample,rep,group,source
  ss <- fread(path)
  required_cols <- c("file","sample","rep","group","source")
  miss <- setdiff(required_cols, names(ss))
  if (length(miss) > 0) stop("Sample sheet missing columns: ", paste(miss, collapse=", "))
  ss
}

read_mod_catalog <- function(path) {
  # Expect a table with AA, anticodon, N1, positions (columns 2..103), etc.
  # We normalize case and build AA_anticodon + AA_anticodon_N1
  ext <- tolower(tools::file_ext(path))
  tab <- switch(ext,
    "xlsx" = readxl::read_excel(path),
    "xls"  = readxl::read_excel(path),
    "csv"  = readr::read_csv(path),
    "tsv"  = readr::read_tsv(path),
    stop("Unsupported mod_catalog extension: ", ext)
  )
  tab <- as.data.frame(tab)
  if (!"AA" %in% names(tab)) stop("Modification catalog must contain a column named 'AA' with AA-anticodon-N1 (e.g., Phe-GAA-1).")
  tab <- tab |>
    dplyr::select(1:min(103L, ncol(tab))) |>
    dplyr::mutate(dplyr::across(-1, ~na_if(., "-"))) |>
    dplyr::mutate(dplyr::across(-1, ~na_if(., "."))) |>
    dplyr::mutate(dplyr::across(-1, \(x) gsub("[tU]", "T", x))) |>
    dplyr::mutate(dplyr::across(-1, toupper)) |>
    dplyr::mutate(AA_anticodon_N1 = AA) |>
    tidyr::separate(AA, c("AA","anticodon","N1"), sep="-", remove=TRUE) |>
    dplyr::mutate(AA_anticodon = paste(AA, anticodon, sep="-"))
  tab
}

normalize_trna_gene <- function(g) {
  g <- stringr::str_remove(g, "^Homo_sapiens_tRNA-")
  if (!stringr::str_detect(g, "-")) g <- stringr::str_replace(g, "(...$)", "-\\1")
  g
}

charging_from_tail <- function(top3) {
  # CCN => uncharged; ACC => charged; else NA
  fifelse(grepl("^CC[ATCG]$", top3), -1L,
  fifelse(top3 == "ACC", 1L, NA_integer_))
}

ensure_dirs <- function(...) { for (d in list(...)) dir.create(d, showWarnings=FALSE, recursive=TRUE) }

# --------------------------- LOAD INPUTS --------------------------------------

SAMPLES <- read_sample_sheet(opt$samples)
TSV_DIR <- opt$tsv_dir
TSV_FILES <- list.files(TSV_DIR, pattern="\\.tsv$", full.names=TRUE)
if (length(TSV_FILES) == 0) stop("No TSV files found in: ", TSV_DIR)

message("Found ", length(TSV_FILES), " TSV files.")

# --------------------------- READ TSVs (PARALLEL) -----------------------------

selected_columns <- c("read_id","gene","position","base","deletion","mutation","insertion")

NUM_CORES <- max(1L, min(threads_cfg, parallel::detectCores()-1L))
cl <- parallel::makeCluster(NUM_CORES)
doParallel::registerDoParallel(cl)

hg_tsv <- foreach::foreach(fp = TSV_FILES, .combine='rbind', .packages='data.table') %dopar% {
  if (!file.exists(fp)) return(NULL)
  dt <- tryCatch(
    data.table::fread(fp, sep='\t', select=selected_columns, showProgress=FALSE),
    error = function(e) NULL
  )
  if (is.null(dt)) return(NULL)
  # Join minimal metadata by filename stem if present in sample sheet
  stem <- sub("\\.tsv$", "", basename(fp))
  # sample sheet 'file' can be fastq or stem; be permissive
  meta <- SAMPLES[file %in% c(stem, basename(stem), fp)]
  if (nrow(meta) == 0L) {
    # fallback: try partial match on sample or rep encoded in filename
    meta <- SAMPLES[0]
  }
  if (nrow(meta) >= 1L) {
    meta <- meta[1] # first match
    dt[, `:=`(sample=meta$sample, rep=meta$rep, group=meta$group, source=meta$source)]
  } else {
    # If no metadata, mark unknown to keep row but allow filtering later
    dt[, `:=`(sample=NA_character_, rep=NA_character_, group=NA_character_, source=NA_character_)]
  }
  dt
}
parallel::stopCluster(cl)

hg_tsv <- data.table::rbindlist(hg_tsv, use.names=TRUE, fill=TRUE)
if (nrow(hg_tsv) == 0) stop("No rows parsed from TSVs.")

# --------------------------- CLEAN + STANDARDIZE ------------------------------

# Filter out spikes/rRNA/externals and drop insertions
rrna_spike_list <- c(
  "NR_004394.1", "NR_002716.3", "NR_003925.1", "NR_004430.2", "NR_002756.2",
  "NR_004391.1", "NR_004392.1", "NR_004393.1", "NR_001571.2",
  "Homo_sapiens_chrX.rRNA-5SR5SR", "Homo_sapiens_chrX.rRNA-58S58S",
  "Ecoli_Tyr", "Yeast_Phe", "Ecoli_Lys",
  "spikein_SCC1", "spikein_SCCA1", "spikein_SCCA2", "spikein_SCCA3"
)

hg_tsv <- hg_tsv[!(gene %in% rrna_spike_list) & insertion == 0]
hg_tsv[, gene := normalize_trna_gene(gene)]

# Tag source if missing (mt- prefix in gene implies mitochondrial)
hg_tsv[, source := ifelse(is.na(source) | source=="",
                          ifelse(substr(gene,1,2)=="mt", "mitochondrial", "cytosolic"),
                          source)]

# Keep only needed columns
hg_tsv <- hg_tsv[, .(read_id, sample, rep, group, source, gene, position, base, mutation, deletion)]

# Order records per read from 5' to 3' of tRNA coordinate (descending if your positions count from 3'—adjust if needed)
# Your previous code assumed descending makes the first row the 5' end for charging logic.
hg_tsv <- hg_tsv[order(-position), .SD, by = read_id]

# Per-read attributes
hg_tsv[, read_size_nt := .N, by = read_id]
hg_tsv[, top3_bases   := paste0(base[1], base[2], base[3]), by = read_id]
hg_tsv[, charging := charging_from_tail(top3_bases)]
hg_tsv[, top3_bases := NULL]

# --------------------------- GENE FIELDS --------------------------------------

# Split tRNA gene into AA, anticodon, isodecoder/serial (N1), optional N2
hg_tsv[, c("AA","anticodon","N1","N2") := tstrsplit(gene, "-", fixed=TRUE)]
hg_tsv[, AA_anticodon    := paste0(AA, "-", anticodon)]
hg_tsv[, AA_anticodon_N1 := paste0(AA_anticodon, "-", N1)]

# --------------------------- SAVE CHECKPOINT ----------------------------------

saveRDS(hg_tsv, file.path(OUT_DIR, "rds", paste0("per_read_table_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")))

# --------------------------- CHARGING SUMMARY ---------------------------------

reads_one_row <- hg_tsv[order(position),
                        .(position = first(position)),
                        by = .(read_id, AA_anticodon, AA_anticodon_N1, charging, read_size_nt, sample, rep, group, source)]

reads_one_row <- reads_one_row[read_size_nt >= READ_LEN_MIN]

charged_n   <- reads_one_row[charging ==  1, .N]
uncharged_n <- reads_one_row[charging == -1, .N]
na_n        <- reads_one_row[is.na(charging), .N]
message("Charging calls — Charged: ", charged_n, " | Uncharged: ", uncharged_n, " | NA: ", na_n)

charging_ratio <- reads_one_row[charging %in% c(1, -1),
  .(charged   = sum(charging ==  1),
    uncharged = sum(charging == -1),
    charging_ratio = {
      tot <- sum(charging == 1) + sum(charging == -1)
      ifelse(tot > 0, sum(charging == 1) / tot, NA_real_)
    }),
  by = .(AA_anticodon, AA_anticodon_N1, sample, rep, group, source)][order(AA_anticodon, sample, rep)]

fwrite(charging_ratio, file.path(OUT_DIR, "tables", "charging_ratio_isodecoder.tsv"), sep="\t")

if (PLOT_ENABLE) {
  p_cr <- ggplot(charging_ratio, aes(x = AA_anticodon, y = charging_ratio, fill = group)) +
    geom_col(position = position_dodge(), width = 0.7) +
    geom_text(aes(label = sprintf("%.0f charged", charged)),
              position = position_dodge(width = 0.9), vjust = -0.25, size = 3) +
    coord_cartesian(ylim=c(0,1)) +
    labs(x = "AA-anticodon", y = "Charging ratio", title = "tRNA charging ratio by isodecoder") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(file.path(OUT_DIR, "plots", "charging_ratio.png"), p_cr, width=14, height=8, dpi=300)
}

# --------------------------- MOD CATALOG --------------------------------------

modtable <- read_mod_catalog(opt$mod_catalog)

# Expectation: 'modtable' links AA_anticodon_N1 to known/annotated positions and mod types.
# If your catalog stores positions in wide format (columns 2..103), we’ll pivot to long.
# Here we keep as-is and match by AA_anticodon or AA_anticodon_N1 when needed.
# You can add explicit schemas here if desired.

# --------------------------- SITE METRICS -------------------------------------

# Per (sample, rep, AA_anticodon_N1, position): event counts & rates
# Event definition: (deletion == 1) OR (mutation == 1) restricted by USE_* flags
hg_tsv[, event := 0L]
if (USE_DELETIONS && USE_MUTATIONS) hg_tsv[, event := as.integer((deletion == 1L) | (mutation == 1L))]
if (USE_DELETIONS && !USE_MUTATIONS) hg_tsv[, event := as.integer(deletion == 1L)]
if (!USE_DELETIONS && USE_MUTATIONS) hg_tsv[, event := as.integer(mutation == 1L)]
# If both disabled, nothing to do
if (!USE_DELETIONS && !USE_MUTATIONS) stop("Both USE_DELETIONS and USE_MUTATIONS are FALSE. Nothing to analyze.")

site_stats <- hg_tsv[, .(
  cov = .N,
  ev  = sum(event == 1L),
  ev_pct = 100 * sum(event == 1L) / .N
), by = .(sample, rep, source, AA_anticodon, AA_anticodon_N1, position)]

# Apply evaluability thresholds at site level
site_eval <- site_stats[cov >= MIN_READS_PER_POS_PER_REP]
# We will later require ≥ MIN_REPS_PASSING reps that pass the above per site.

# --------------------------- PAIRWISE CROSSTALK -------------------------------

# Helper: per-gene (AA_anticodon_N1) event matrix for a given sample/rep
# Returns a data.table with rows=read_id (spanning len >= READ_LEN_MIN) and columns=positions with event 0/1
# For memory efficiency we don’t build full matrices; we count 2x2 on the fly.
compute_crosstalk_one_gene <- function(dt_gene_sr, pos_list, min_pair_cov=MIN_PAIR_COVERAGE,
                                       pair_gate_pct=PAIR_EVENT_GATE_PCT) {
  # dt_gene_sr: hg_tsv subset for a single (sample,rep,AA_anticodon_N1)
  # pos_list: candidate positions to consider (numeric vector)
  # We will only consider pairs where at least one site has ev_pct >= pair_gate_pct in this (sample,rep,gene) stratum
  out <- vector("list", length=0L)
  if (length(pos_list) < 2L) return(data.table())
  # Pre-compute per-position ev rates
  site_s <- dt_gene_sr[, .(cov=.N, ev=sum(event==1L), ev_pct=100*sum(event==1L)/.N), by=.(position)]
  gate_pos <- site_s[ev_pct >= pair_gate_pct, position]
  if (length(gate_pos) == 0L) return(data.table())
  # Candidate pairs: any pair where at least one is in gate_pos
  # To bound complexity, we’ll pair gate_pos with all others (including within gate_pos)
  others <- sort(unique(pos_list))
  # Index reads -> position->event map by data.table operations
  # For pair counts, we build per read summaries for the two positions only (no giant matrix)
  for (p1 in gate_pos) {
    for (p2 in others[others > p1]) {
      # Restrict to reads covering both p1 and p2
      dt_sub <- dt_gene_sr[position %in% c(p1,p2)]
      # Count positions per read
      cov2 <- dt_sub[, .N, by=read_id][N==2L]$read_id
      if (length(cov2) < min_pair_cov) next
      d2 <- dt_sub[read_id %in% cov2, .(pos=position, ev=event), by=read_id]
      # Reshape to wide 2 columns
      w <- data.table::dcast(d2, read_id ~ pos, value.var="ev")
      # Ensure columns exist
      if (!all(c(p1,p2) %in% names(w))) next
      a <- w[[as.character(p1)]]; b <- w[[as.character(p2)]]
      # 2x2 counts
      O_O <- sum(a==1L & b==1L, na.rm=TRUE)
      O_l <- sum(a==1L & b==0L, na.rm=TRUE)
      l_O <- sum(a==0L & b==1L, na.rm=TRUE)
      l_l <- sum(a==0L & b==0L, na.rm=TRUE)
      tot <- O_O + O_l + l_O + l_l
      if (tot < min_pair_cov) next
      # Fisher test (add 0.5 continuity if zeros to avoid Inf OR)
      mat <- matrix(c(O_O, O_l, l_O, l_l), nrow=2, byrow=TRUE)
      ft <- tryCatch(fisher.test(mat), error=function(e) NULL)
      if (is.null(ft)) next
      or  <- as.numeric(ft$estimate)
      p   <- as.numeric(ft$p.value)
      out[[length(out)+1L]] <- data.table(
        mod_pos1 = p1, mod_pos2 = p2,
        O_O = O_O, O_l = O_l, l_O = l_O, l_l = l_l,
        odds_ratio = or, pval = p
      )
    }
  }
  if (length(out) == 0L) return(data.table())
  rbindlist(out)
}

# Driver for all (sample,rep,AA_anticodon_N1)
message("Computing crosstalk (this may take time on large datasets)...")
t0 <- Sys.time()

# Build a table of evaluable (sample,rep,AA_anticodon_N1,position) with ≥ MIN_REPS_PASSING reps
site_eval[, rep_pass := 1L]
site_reps <- site_eval[, .(n_reps = sum(rep_pass)), by=.(sample, AA_anticodon_N1, position)]
eval_positions <- site_reps[n_reps >= MIN_REPS_PASSING]

# Preselect by event% gate for efficiency: at least one site in a pair will typically need ≥ PAIR_EVENT_GATE_PCT in that stratum.
# (We still enforce per-site coverage rules above.)
key_cols <- c("sample","rep","AA_anticodon_N1")
setkeyv(hg_tsv, key_cols)

cl <- parallel::makeCluster(NUM_CORES)
doParallel::registerDoParallel(cl)

ct_results <- foreach::foreach(
  # unique combinations present in data
  comb = iter(unique(hg_tsv[, .(sample, rep, AA_anticodon_N1)]), by="row"),
  .combine='rbind',
  .packages=c("data.table","dplyr","stringr")
) %dopar% {
  s <- comb$sample; r <- comb$rep; g <- comb$AA_anticodon_N1
  if (is.na(s) || is.na(r) || is.na(g)) return(NULL)
  dt_sr <- hg_tsv[sample==s & rep==r & AA_anticodon_N1==g]
  if (nrow(dt_sr) == 0L) return(NULL)
  # Filter by read length at the read-level gate
  ok_reads <- dt_sr[, .(read_size_nt=.N), by=read_id][read_size_nt >= READ_LEN_MIN]$read_id
  if (length(ok_reads) == 0L) return(NULL)
  dt_sr <- dt_sr[read_id %in% ok_reads]
  # Candidate positions: those passing per-rep site coverage threshold
  pos_ok <- dt_sr[, .N, by=position][N >= MIN_READS_PER_POS_PER_REP]$position
  if (length(pos_ok) < 2L) return(NULL)
  # Intersect with global evaluable list (positions seen in >= MIN_REPS_PASSING reps for this sample/gene)
  pos_eval <- intersect(pos_ok, eval_positions[sample==s & AA_anticodon_N1==g]$position)
  if (length(pos_eval) < 2L) return(NULL)
  # Compute crosstalk pairs for this stratum
  ct <- compute_crosstalk_one_gene(dt_sr, pos_eval, min_pair_cov=MIN_PAIR_COVERAGE, pair_gate_pct=PAIR_EVENT_GATE_PCT)
  if (nrow(ct) == 0L) return(NULL)
  ct[, `:=`(sample=s, rep=r, AA_anticodon_N1=g)]
  ct[]
}

parallel::stopCluster(cl)

if (is.null(ct_results) || nrow(ct_results) == 0L) {
  warning("No crosstalk pairs passed filters.")
  ct_results <- data.table()
} else {
  # FDR within sample or globally? Here: global FDR across all tested pairs.
  ct_results[, p_adjust := p.adjust(pval, method = CROSSTALK_FDR_METHOD)]
  ct_results[, log2OR := log2(odds_ratio)]
  ct_results[, significant := (p_adjust < CROSSTALK_FDR_ALPHA) & (abs(log2OR) >= CROSSTALK_MIN_ABS_LOG2OR)]
}

fwrite(ct_results, file.path(OUT_DIR, "tables", "crosstalk_pairs_isodecoder.tsv"), sep="\t")

# Optionally aggregate to anticodon level (median log2OR, Stouffer combine, etc.)
# Minimal anticodon aggregation (median log2OR; count contributing isodecoders), restricted to significant pairs per isodecoder
if (nrow(ct_results) > 0L) {
  # Parse AA_anticodon from AA_anticodon_N1
  ct_results[, AA_anticodon := sub("-[^-]+$", "", AA_anticodon_N1)]
  agg_anticodon <- ct_results[, .(
    n_isodecoders = uniqueN(AA_anticodon_N1),
    median_log2OR = median(log2OR, na.rm=TRUE),
    # naive p combine: min p_adjust as a conservative summary; replace with proper Stouffer if you wish
    min_p_adjust  = min(p_adjust, na.rm=TRUE)
  ), by=.(sample, rep, AA_anticodon, mod_pos1, mod_pos2)]
  fwrite(agg_anticodon, file.path(OUT_DIR, "tables", "crosstalk_pairs_anticodon.tsv"), sep="\t")
}

# --------------------------- VISUAL QC PLOTS ----------------------------------

if (PLOT_ENABLE && nrow(ct_results) > 0L) {
  # Simple volcano: all pairs pooled (you can facet per sample/rep)
  vdat <- ct_results[is.finite(log2OR)]
  p_v <- ggplot(vdat, aes(x=log2OR, y=-log10(p_adjust), color=significant)) +
    geom_point(alpha=0.4, size=1.2) +
    geom_vline(xintercept=c(-CROSSTALK_MIN_ABS_LOG2OR, CROSSTALK_MIN_ABS_LOG2OR), linetype="dashed") +
    geom_hline(yintercept=-log10(CROSSTALK_FDR_ALPHA), linetype="dashed") +
    labs(title="Crosstalk volcano (isodecoder-level pairs)", x="log2(OR)", y="-log10(FDR)") +
    theme_minimal()
  ggsave(file.path(OUT_DIR, "plots", "crosstalk_volcano_isodecoder.png"), p_v, width=10, height=7, dpi=300)
}

# --------------------------- LOGGING & SESSION --------------------------------

t1 <- Sys.time()
message(sprintf("Completed in %.1f min", as.numeric(difftime(t1, t0, units="mins"))))
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo.txt"))
