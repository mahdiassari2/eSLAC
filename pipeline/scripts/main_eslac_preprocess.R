#!/usr/bin/env Rscript
# ==============================================================================
# eSLAC Preprocessing & Analysis (MSR-seq Input / Global Summaries)
# Author: Mahdi Assari
# Repo:   https://github.com/mahdiassari2/eSLAC
# License: MIT
# ------------------------------------------------------------------------------
# What this script does
#   1) Loads per-read basewise TSVs and gene-level count files in parallel.
#   2) Standardizes tRNA gene identifiers (AA-anticodon-N1[-N2]).
#   3) Builds abundance summaries at isodecoder and anticodon levels.
#   4) Computes charging summaries from per-read tails (CCA vs CCN) when available.
#   5) Produces a volcano plot for a user-specified contrast (group1 vs group2).
#
# All parameters come from --config (YAML) and CLI flags. No in-script installs.
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(foreach)
  library(doParallel)
  library(iterators)
  library(scales)
  library(writexl)
})

# ------------------------------ CLI -------------------------------------------

option_list <- list(
  make_option("--config",     type="character", help="Path to config.yaml", metavar="FILE"),
  make_option("--tsv_dir",    type="character", help="Directory with basewise TSVs"),
  make_option("--count_dir",  type="character", help="Directory with gene-level count files"),
  make_option("--samples",    type="character", help="Sample sheet CSV: file,sample,rep,group,source"),
  make_option("--out_dir",    type="character", help="Output directory (tables/plots/rds)", default="results"),
  make_option("--volcano_group1", type="character", help="Group name (e.g., input)", default=NA_character_),
  make_option("--volcano_group2", type="character", help="Group name (e.g., poly)",  default=NA_character_),
  make_option("--threads",    type="integer",   help="Threads override", default=NA_integer_)
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$config) || is.null(opt$tsv_dir) || is.null(opt$count_dir) || is.null(opt$samples)) {
  stop("Missing required arguments. See --help.")
}

cfg <- yaml::read_yaml(opt$config)

OUT_DIR <- opt$out_dir
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR,"tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR,"plots"),  recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR,"rds"),    recursive=TRUE, showWarnings=FALSE)

# ------------------------------ Helpers ---------------------------------------

`%nin%` <- Negate(`%in%`)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n")
`%||%` <- function(a,b) if (is.null(a)) b else a

get_cfg <- function(path, default=NULL) {
  purrr::reduce(strsplit(path, "\\.")[[1]], function(acc,k) if(!is.null(acc[[k]])) acc[[k]] else NULL, .init=cfg) %||% default
}

READ_LEN_MIN          <- as.integer(get_cfg("coverage.read_length_min", 60))
DROP_INSERTIONS       <- isTRUE(get_cfg("preprocess.drop_insertions", TRUE))
REMOVE_RRNA_SPIKES    <- isTRUE(get_cfg("preprocess.remove_rrna_spikes", TRUE))
PLOT_ENABLE           <- isTRUE(get_cfg("reporting.plots", TRUE))

threads_cfg <- if (is.na(opt$threads)) {
  if (!is.null(cfg$alignment$bowtie2$threads)) cfg$alignment$bowtie2$threads else max(1L, parallel::detectCores()-2L)
} else opt$threads
data.table::setDTthreads(threads_cfg)

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

read_sample_sheet <- function(path) {
  ss <- data.table::fread(path)
  req <- c("file","sample","rep","group","source")
  miss <- setdiff(req, names(ss))
  if (length(miss)) stop("Sample sheet missing columns: ", paste(miss, collapse=", "))
  ss
}

# ------------------------------ Load Sample Sheet -----------------------------

SAMPLES <- read_sample_sheet(opt$samples)

# ------------------------------ Enumerate Files --------------------------------

TSV_FILES   <- list.files(opt$tsv_dir,   pattern="\\.tsv$", full.names=TRUE)
COUNT_FILES <- list.files(opt$count_dir, pattern="\\.(tsv|txt|csv)$", ignore.case=TRUE, full.names=TRUE)

if (!length(TSV_FILES))   stop("No TSV files found in: ", opt$tsv_dir)
if (!length(COUNT_FILES)) stop("No count files found in: ", opt$count_dir)

msg("Basewise TSVs: ", length(TSV_FILES))
msg("Count files:   ", length(COUNT_FILES))

# ------------------------------ Read TSVs (parallel) -------------------------

selected_cols_tsv <- c("read_id","gene","position","base","pileup","deletion","mutation","insertion")

NUM_CORES <- max(1L, min(threads_cfg, parallel::detectCores()-1L))
cl <- parallel::makeCluster(NUM_CORES)
doParallel::registerDoParallel(cl)

read_tsv <- function(fp) {
  dt <- tryCatch(data.table::fread(fp, sep="\t", showProgress=FALSE), error=function(e) NULL)
  if (is.null(dt)) return(NULL)
  keep <- intersect(selected_cols_tsv, names(dt))
  if (!length(keep)) return(NULL)
  dt <- dt[, ..keep]
  stem <- sub("\\.tsv$","",basename(fp))
  meta <- SAMPLES[file %in% c(stem, basename(stem), fp)]
  if (nrow(meta)) meta <- meta[1] else meta <- SAMPLES[0]
  if (nrow(meta)) {
    dt[, `:=`(sample=meta$sample, rep=meta$rep, group=meta$group, source=meta$source)]
  } else {
    dt[, `:=`(sample=NA_character_, rep=NA_character_, group=NA_character_, source=NA_character_)]
  }
  dt
}

tsv_parts <- foreach::foreach(fp = TSV_FILES, .packages="data.table") %dopar% read_tsv(fp)
parallel::stopCluster(cl)

TSV <- data.table::rbindlist(tsv_parts, use.names=TRUE, fill=TRUE)
rm(tsv_parts); gc()
if (!nrow(TSV)) stop("Failed to parse TSVs.")

# ------------------------------ Read Counts -----------------------------------

# Expect at least two columns: name/gene + count (handle common headers)
read_counts <- function(fp){
  dt <- tryCatch(data.table::fread(fp, sep="\t", showProgress=FALSE), error=function(e) NULL)
  if (is.null(dt)) dt <- tryCatch(readr::read_csv(fp, show_col_types = FALSE) |> as.data.table(), error=function(e) NULL)
  if (is.null(dt)) return(NULL)

  nms <- tolower(names(dt))
  gene_col  <- which(nms %in% c("name","gene","target","id"))[1]
  count_col <- which(nms %in% c("count","reads","n","nreads"))[1]
  if (is.na(gene_col) || is.na(count_col)) return(NULL)

  dt <- dt[, .(gene = as.character(dt[[gene_col]]),
               count = as.numeric(dt[[count_col]]))]
  stem <- sub("\\.(tsv|txt|csv)$","",basename(fp), ignore.case=TRUE)
  meta <- SAMPLES[file %in% c(stem, basename(stem), fp)]
  if (nrow(meta)) meta <- meta[1] else meta <- SAMPLES[0]
  if (nrow(meta)) {
    dt[, `:=`(sample=meta$sample, rep=meta$rep, group=meta$group, source=meta$source)]
  } else {
    dt[, `:=`(sample=NA_character_, rep=NA_character_, group=NA_character_, source=NA_character_)]
  }
  dt
}

count_list <- lapply(COUNT_FILES, read_counts)
COUNTS <- data.table::rbindlist(count_list, use.names=TRUE, fill=TRUE)
rm(count_list); gc()
if (!nrow(COUNTS)) stop("Failed to parse count files.")

# ------------------------------ Clean / Normalize -----------------------------

# Remove spikes / rRNA (optional)
if (REMOVE_RRNA_SPIKES) {
  rrna_spike <- c(
    "NR_004394.1","NR_002716.3","NR_003925.1","NR_004430.2","NR_002756.2",
    "NR_004391.1","NR_004392.1","NR_004393.1","NR_001571.2",
    "Homo_sapiens_chrX.rRNA-5SR5SR","Homo_sapiens_chrX.rRNA-58S58S",
    "Ecoli_Tyr","Yeast_Phe","Ecoli_Lys",
    "spikein_SCC1","spikein_SCCA1","spikein_SCCA2","spikein_SCCA3"
  )
  TSV    <- TSV[gene %nin% rrna_spike]
  COUNTS <- COUNTS[gene %nin% rrna_spike]
}

# Drop insertions if requested (we use pileup, deletion, mutation for site rates)
if (DROP_INSERTIONS && "insertion" %in% names(TSV)) {
  TSV <- TSV[insertion == 0 | is.na(insertion)]
}

# Normalize gene ID
TSV[,    gene := normalize_trna_gene(gene)]
COUNTS[, gene := normalize_trna_gene(gene)]

# Source fallback
TSV[, source := ifelse(is.na(source) | source=="",
                       ifelse(substr(gene,1,2)=="mt","mitochondrial","cytosolic"),
                       source)]
COUNTS[, source := ifelse(is.na(source) | source=="",
                          ifelse(substr(gene,1,2)=="mt","mitochondrial","cytosolic"),
                          source)]

# Split gene into fields
TSV[,    c("AA","anticodon","N1","N2") := tstrsplit(gene, "-", fixed=TRUE)]
TSV[,    AA_anticodon    := paste0(AA,"-",anticodon)]
TSV[,    AA_anticodon_N1 := paste0(AA_anticodon,"-",N1)]

COUNTS[, c("AA","anticodon","N1","N2") := tstrsplit(gene, "-", fixed=TRUE)]
COUNTS[, AA_anticodon    := paste0(AA,"-",anticodon)]
COUNTS[, AA_anticodon_N1 := paste0(AA_anticodon,"-",N1)]

# ------------------------------ Per-read features -----------------------------

# Order records per read from 5' to 3' (descending index in your convention)
TSV <- TSV[order(-position), .SD, by=read_id]

TSV[, read_size_nt := .N, by=read_id]
TSV[, tail3 := paste0(base[1], base[2], base[3]), by=read_id]
TSV[, charging := charging_from_tail(tail3)]
TSV[, tail3 := NULL]

saveRDS(TSV,    file.path(OUT_DIR,"rds", paste0("preprocess_perread_",format(Sys.time(),"%Y%m%d_%H%M%S"),".rds")))
saveRDS(COUNTS, file.path(OUT_DIR,"rds", paste0("preprocess_counts_" ,format(Sys.time(),"%Y%m%d_%H%M%S"),".rds")))

# ------------------------------ Abundance Summaries ---------------------------

# Isodecoder counts per (sample,rep)
iso_rep <- COUNTS[, .(reads = sum(count, na.rm=TRUE)), by=.(AA_anticodon_N1, AA_anticodon, sample, rep, group, source)]
fwrite(iso_rep, file.path(OUT_DIR,"tables","abundance_isodecoder_rep.tsv"), sep="\t")

# Isodecoder totals per (sample)
iso_smp <- iso_rep[, .(reads = sum(reads, na.rm=TRUE)), by=.(AA_anticodon_N1, AA_anticodon, sample, group, source)]
fwrite(iso_smp, file.path(OUT_DIR,"tables","abundance_isodecoder_sample.tsv"), sep="\t")

# Anticodon totals per (sample)
anti_smp <- iso_smp[, .(reads = sum(reads, na.rm=TRUE)), by=.(AA_anticodon, sample, group, source)]
fwrite(anti_smp, file.path(OUT_DIR,"tables","abundance_anticodon_sample.tsv"), sep="\t")

# RPM within sample (anticodon level)
anti_smp[, total_reads := sum(reads, na.rm=TRUE), by=.(sample)]
anti_smp[, rpm := 1e6 * reads / pmax(1, total_reads)]
fwrite(anti_smp[, .(AA_anticodon, sample, group, source, reads, rpm)],
       file.path(OUT_DIR,"tables","abundance_anticodon_sample_rpm.tsv"), sep="\t")

# ------------------------------ Charging Summaries (per-read) -----------------

# One row per read_id for charging analysis; keep adequate read length
reads_1 <- TSV[read_size_nt >= READ_LEN_MIN][order(position),
           .(position = first(position)),
           by=.(read_id, AA_anticodon, AA_anticodon_N1, charging, sample, rep, group, source)]

charging_iso <- reads_1[charging %in% c(1,-1),
  .(charged   = sum(charging == 1),
    uncharged = sum(charging == -1),
    charging_ratio = {
      tot <- sum(charging == 1) + sum(charging == -1)
      ifelse(tot > 0, sum(charging == 1) / tot, NA_real_)
    }),
  by=.(AA_anticodon_N1, AA_anticodon, sample, rep, group, source)]

fwrite(charging_iso, file.path(OUT_DIR,"tables","charging_ratio_isodecoder.tsv"), sep="\t")

if (PLOT_ENABLE) {
  p_cr <- ggplot(charging_iso, aes(x=AA_anticodon, y=charging_ratio, fill=group)) +
    geom_col(position = position_dodge(), width = 0.7) +
    coord_cartesian(ylim=c(0,1)) +
    labs(x="AA-anticodon", y="Charging ratio", title="Charging by isodecoder") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
  ggsave(file.path(OUT_DIR,"plots","charging_isodecoder.png"), p_cr, width=14, height=8, dpi=300)
}

# ------------------------------ Volcano Plot (group1 vs group2) --------------

g1 <- opt$volcano_group1
g2 <- opt$volcano_group2

if (!is.na(g1) && !is.na(g2)) {
  msg("Volcano contrast: ", g1, " vs ", g2, " (anticodon-level RPM)")

  # Use anticodon-level RPM table; need per-sample reps for a t-test-like comparison.
  # If you have true biological replicates, ensure the sample sheet encodes them distinctly.
  vul <- anti_smp[, .(AA_anticodon, sample, group, rpm = 1e6*reads/pmax(1, total_reads)), by=.(AA_anticodon, sample, group, total_reads)]
  vul <- vul[group %in% c(g1, g2)]

  # Filter to anticodons present in both groups with ≥2 samples per group (basic guard)
  keep <- vul[, .N, by=.(AA_anticodon, group)][, .N, by=.(AA_anticodon)][N == 2]$AA_anticodon
  vul  <- vul[AA_anticodon %in% keep]

  # Compute per-anticodon stats (Welch test on log2 RPM+1)
  stats <- vul[, {
      x <- log2(rpm[group==g1] + 1)
      y <- log2(rpm[group==g2] + 1)
      if (length(x) >= 2 && length(y) >= 2) {
        tt <- tryCatch(t.test(x, y), error=function(e) NULL)
        pv <- if (!is.null(tt)) tt$p.value else NA_real_
        lfc <- mean(x) - mean(y)
        .(log2FC = lfc, p = pv)
      } else .(log2FC = NA_real_, p = NA_real_)
    }, by=AA_anticodon]

  stats[, padj := p.adjust(p, method="BH")]
  fwrite(stats, file.path(OUT_DIR,"tables", sprintf("volcano_stats_%s_vs_%s.tsv", g1, g2)), sep="\t")

  if (PLOT_ENABLE) {
    p_v <- ggplot(stats, aes(x=log2FC, y=-log10(padj))) +
      geom_point(alpha=0.6) +
      geom_vline(xintercept=c(-1,1), linetype="dashed") +
      geom_hline(yintercept=-log10(0.05), linetype="dashed") +
      labs(title=sprintf("Volcano: %s vs %s (anticodon RPM)", g1, g2),
           x="log2FC", y="-log10(FDR)") +
      theme_bw()
    ggsave(file.path(OUT_DIR,"plots", sprintf("volcano_%s_vs_%s.png", g1, g2)),
           p_v, width=7, height=5, dpi=300)
  }
} else {
  msg("Volcano contrast not requested (use --volcano_group1 and --volcano_group2).")
}

# ------------------------------ Outputs & Session -----------------------------

capture.output(sessionInfo(), file=file.path(OUT_DIR,"sessionInfo_preprocess.txt"))

msg("Done. Tables at: ", file.path(OUT_DIR,"tables"))
msg("Plots  at: ", file.path(OUT_DIR,"plots"))
