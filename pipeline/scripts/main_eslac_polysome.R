#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# eSLAC: Polysome tRNA Profiling (MSR-seq)
# - Abundance shifts (polysome vs input) at isodecoder and anticodon levels
# - Charging ratios per condition and Δcharging (poly - input)
# - Ribosomal-site proxy: E-site (C75) vs A/P-site (A76) comparisons of event rates
# Configuration via --config + CLI flags. No in-script installs (use conda/renv).
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(yaml)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

# --------------------------- CLI & CONFIG -------------------------------------

option_list <- list(
  make_option("--config",   type="character", help="Path to config.yaml", metavar="FILE"),
  make_option("--tsv_dir",  type="character", help="Directory of per-read TSVs"),
  make_option("--samples",  type="character", help="Sample sheet CSV: file,sample,rep,group,source"),
  make_option("--out_dir",  type="character", help="Output directory", default="results"),
  make_option("--threads",  type="integer",   help="Threads override", default=NA_integer_)
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$config) || is.null(opt$tsv_dir) || is.null(opt$samples)) {
  stop("Missing required arguments. See --help.")
}

cfg <- yaml::read_yaml(opt$config)

OUT_DIR <- opt$out_dir
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "plots"),  recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "rds"),    recursive=TRUE, showWarnings=FALSE)

SEED <- 20251028
set.seed(SEED)

# Threading
threads_cfg <- if (is.na(opt$threads)) {
  if (!is.null(cfg$alignment$bowtie2$threads)) cfg$alignment$bowtie2$threads else max(1L, parallel::detectCores()-2L)
} else opt$threads
data.table::setDTthreads(threads_cfg)

# --------------------------- PARAMS & HELPERS ---------------------------------

`%nin%` <- Negate(`%in%`)
msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "-", ..., "\n")

get_cfg <- function(path, default=NULL) {
  purrr::reduce(strsplit(path, "\\.")[[1]], function(acc,k) if(!is.null(acc[[k]])) acc[[k]] else NULL, .init=cfg) %||% default
}

READ_LEN_MIN              <- as.integer(get_cfg("coverage.read_length_min", 60))
MIN_READS_PER_CELL        <- as.integer(get_cfg("coverage.min_reads_per_position_per_rep", 30)) # for site rates
USE_MUTATIONS             <- isTRUE(get_cfg("mod_calling.use_mutations", TRUE))
USE_DELETIONS             <- isTRUE(get_cfg("mod_calling.use_deletions", TRUE))
PLOT_ENABLE               <- isTRUE(get_cfg("reporting.plots", TRUE))

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

# --------------------------- LOAD INPUTS --------------------------------------

SAMPLES <- read_sample_sheet(opt$samples)
TSV_DIR <- opt$tsv_dir
TSV_FILES <- list.files(TSV_DIR, pattern="\\.tsv$", full.names=TRUE)
if (!length(TSV_FILES)) stop("No TSVs in: ", TSV_DIR)

msg("TSVs found: ", length(TSV_FILES))

# --------------------------- READ + STANDARDIZE -------------------------------

selected_cols <- c("read_id","gene","position","base","deletion","mutation","insertion")

dt_list <- lapply(TSV_FILES, function(fp){
  dt <- tryCatch(data.table::fread(fp, sep="\t", showProgress=FALSE), error=function(e) NULL)
  if (is.null(dt)) return(NULL)
  dt <- dt[, intersect(selected_cols, names(dt)), with=FALSE]
  stem <- sub("\\.tsv$","",basename(fp))
  meta <- SAMPLES[file %in% c(stem, basename(stem), fp)]
  if (nrow(meta)) meta <- meta[1] else meta <- SAMPLES[0]
  if (nrow(meta)) {
    dt[, `:=`(sample=meta$sample, rep=meta$rep, group=meta$group, source=meta$source)]
  } else {
    dt[, `:=`(sample=NA_character_, rep=NA_character_, group=NA_character_, source=NA_character_)]
  }
  dt
})

DT <- data.table::rbindlist(dt_list, use.names=TRUE, fill=TRUE)
rm(dt_list); gc()

if (!nrow(DT)) stop("No data parsed.")

# Remove spikes/rRNA/external; drop insertions
rrna_spike <- c(
  "NR_004394.1","NR_002716.3","NR_003925.1","NR_004430.2","NR_002756.2",
  "NR_004391.1","NR_004392.1","NR_004393.1","NR_001571.2",
  "Homo_sapiens_chrX.rRNA-5SR5SR","Homo_sapiens_chrX.rRNA-58S58S",
  "Ecoli_Tyr","Yeast_Phe","Ecoli_Lys",
  "spikein_SCC1","spikein_SCCA1","spikein_SCCA2","spikein_SCCA3"
)
DT <- DT[gene %nin% rrna_spike & insertion == 0]

# Normalize gene and source
DT[, gene := normalize_trna_gene(gene)]
DT[, source := ifelse(is.na(source) | source=="",
                      ifelse(substr(gene,1,2)=="mt","mitochondrial","cytosolic"),
                      source)]

# Per-read ordering (descending if that makes the first row 5' end for your numbering)
DT <- DT[order(-position), .SD, by=read_id]

# Per-read features
DT[, read_size_nt := .N, by=read_id]
DT[, top3_bases   := paste0(base[1], base[2], base[3]), by=read_id]
DT[, charging     := charging_from_tail(top3_bases)]
DT[, top3_bases := NULL]

# Gene fields
DT[, c("AA","anticodon","N1","N2") := tstrsplit(gene, "-", fixed=TRUE)]
DT[, AA_anticodon    := paste0(AA,"-",anticodon)]
DT[, AA_anticodon_N1 := paste0(AA_anticodon,"-",N1)]

# Save checkpoint
saveRDS(DT, file.path(OUT_DIR,"rds", paste0("polysome_per_read_",format(Sys.time(),"%Y%m%d_%H%M%S"),".rds")))

# --------------------------- ABUNDANCE (POLY vs INPUT) -----------------------

# We expect 'group' to label conditions like "poly" and "input"
# Guard: drop NA samples and small reads
DTq <- DT[!is.na(sample) & read_size_nt >= READ_LEN_MIN]

# Count reads per stratum
counts_iso <- DTq[, .N, by=.(AA_anticodon_N1, AA_anticodon, sample, rep, group, source)]
data.table::setnames(counts_iso, "N", "reads")

# Summaries at isodecoder
iso_sum <- counts_iso[, .(reads = sum(reads)), by=.(AA_anticodon_N1, AA_anticodon, group, source)]
# Wide to compute log2FC polysome/input
iso_w <- data.table::dcast(iso_sum, AA_anticodon_N1 + AA_anticodon + source ~ group, value.var="reads", fill=0)
grp_levels <- setdiff(names(iso_w), c("AA_anticodon_N1","AA_anticodon","source"))
if (!all(c("poly","input") %in% grp_levels)) {
  msg("WARNING: 'poly'/'input' labels not both present in sample sheet 'group'. Found: ", paste(grp_levels, collapse=", "))
}
iso_w[, log2FC_poly_vs_input := log2( (poly + 1) / (input + 1) ) ]
fwrite(iso_w, file.path(OUT_DIR,"tables","abundance_isodecoder_poly_vs_input.tsv"), sep="\t")

# Summaries at anticodon
anti_sum <- counts_iso[, .(reads = sum(reads)), by=.(AA_anticodon, group, source)]
anti_w <- data.table::dcast(anti_sum, AA_anticodon + source ~ group, value.var="reads", fill=0)
anti_w[, log2FC_poly_vs_input := log2( (poly + 1) / (input + 1) ) ]
fwrite(anti_w, file.path(OUT_DIR,"tables","abundance_anticodon_poly_vs_input.tsv"), sep="\t")

# Plot abundance (anticodon)
if (PLOT_ENABLE) {
  p_ab <- ggplot(anti_w, aes(x=reorder(AA_anticodon, log2FC_poly_vs_input),
                             y=log2FC_poly_vs_input, fill=source)) +
    geom_col() +
    coord_flip() +
    labs(title="Polysome enrichment (log2FC poly/input) by anticodon",
         x="AA-anticodon", y="log2FC (poly/input)") +
    theme_minimal()
  ggsave(file.path(OUT_DIR,"plots","abundance_anticodon_log2FC.png"), p_ab, width=10, height=12, dpi=300)
}

# --------------------------- CHARGING RATIOS ----------------------------------

# One row per read
reads_1 <- DTq[order(position),
               .(position = first(position)),
               by=.(read_id, AA_anticodon, AA_anticodon_N1, charging, sample, rep, group, source)]

# Charging per isodecoder x group
charge_iso <- reads_1[charging %in% c(1,-1),
  .(charged=sum(charging==1), uncharged=sum(charging==-1)),
  by=.(AA_anticodon_N1, AA_anticodon, group, source)]

charge_iso[, charging_ratio := charged / pmax(1, (charged + uncharged))]
fwrite(charge_iso, file.path(OUT_DIR,"tables","charging_ratio_isodecoder_by_group.tsv"), sep="\t")

# Δcharging (poly - input) per isodecoder
charge_w <- data.table::dcast(charge_iso, AA_anticodon_N1 + AA_anticodon + source ~ group,
                              value.var="charging_ratio", fill=NA_real_)
if (all(c("poly","input") %in% names(charge_w))) {
  charge_w[, delta_charging_poly_minus_input := poly - input]
} else {
  charge_w[, delta_charging_poly_minus_input := NA_real_]
}
fwrite(charge_w, file.path(OUT_DIR,"tables","charging_ratio_isodecoder_delta.tsv"), sep="\t")

if (PLOT_ENABLE && "delta_charging_poly_minus_input" %in% names(charge_w)) {
  p_dc <- ggplot(charge_w, aes(x=reorder(AA_anticodon, delta_charging_poly_minus_input),
                               y=delta_charging_poly_minus_input, fill=source)) +
    geom_col() + coord_flip() +
    labs(title="Δ charging (poly - input) by isodecoder (aggregated to anticodon on x)",
         x="AA-anticodon", y="Δ charging ratio") +
    theme_minimal()
  ggsave(file.path(OUT_DIR,"plots","delta_charging_isodecoder.png"), p_dc, width=10, height=12, dpi=300)
}

# --------------------------- RIBOSOMAL SITE PROXY -----------------------------

# Define per-read terminal base (positioned relative to 3' tail)
# We already derived charging from top3; now classify 3′ end as A76 (A/P sites) vs C75 (E site proxy)
DTq[, tail3 := NA_character_]
DTq[, tail3 := paste0(base[1], base[2], base[3]), by=read_id]
DTq[, ribo_site := fifelse(tail3 == "ACC", "AP_site",
                     fifelse(grepl("^CC[ATCG]$", tail3), "E_site", NA_character_))]
DTq[, tail3 := NULL]

# Event definition for rates
DTq[, event := 0L]
if (USE_DELETIONS && USE_MUTATIONS) DTq[, event := as.integer((deletion==1L) | (mutation==1L))]
if (USE_DELETIONS && !USE_MUTATIONS) DTq[, event := as.integer(deletion==1L)]
if (!USE_DELETIONS && USE_MUTATIONS) DTq[, event := as.integer(mutation==1L)]
if (!USE_DELETIONS && !USE_MUTATIONS) stop("Both USE_DELETIONS and USE_MUTATIONS are FALSE; no event channel to evaluate.")

# Site-level event rates per ribo_site
site_rates <- DTq[!is.na(ribo_site),
  .(cov=.N, ev=sum(event==1L), ev_pct = 100*sum(event==1L)/.N),
  by=.(AA_anticodon_N1, AA_anticodon, source, sample, rep, ribo_site, position)]

# Filter for minimal coverage for robust comparisons
site_rates <- site_rates[cov >= MIN_READS_PER_CELL]

# E vs AP comparison: fold-change of ev_pct (add small epsilon to avoid /0)
eps <- 1e-6
cmp <- data.table::dcast(site_rates,
                         AA_anticodon_N1 + AA_anticodon + source + sample + rep + position ~ ribo_site,
                         value.var = "ev_pct", fill=NA_real_)
if (!all(c("E_site","AP_site") %in% names(cmp))) {
  msg("WARNING: Could not find both E_site and AP_site categories in data.")
  cmp[, log2FC_E_vs_AP := NA_real_]
} else {
  cmp[, log2FC_E_vs_AP := log2( (E_site + eps) / (AP_site + eps) ) ]
}
fwrite(cmp, file.path(OUT_DIR,"tables","eventrate_E_vs_AP_by_position.tsv"), sep="\t")

# Optional plot: distribution of log2FC across positions
if (PLOT_ENABLE && "log2FC_E_vs_AP" %in% names(cmp)) {
  p_ev <- ggplot(cmp[is.finite(log2FC_E_vs_AP)],
                 aes(x=log2FC_E_vs_AP, fill=source)) +
    geom_histogram(bins=80, alpha=0.8, position="identity") +
    geom_vline(xintercept=0, linetype="dashed") +
    labs(title="Event-rate log2FC (E-site vs A/P-site)", x="log2FC(E/AP)", y="counts") +
    theme_minimal()
  ggsave(file.path(OUT_DIR,"plots","eventrate_E_vs_AP_hist.png"), p_ev, width=10, height=6, dpi=300)
}

# --------------------------- POLYSOME vs INPUT: SITE RATES --------------------

# Compare event rates per position between groups (poly vs input), per isodecoder
site_rates_grp <- DTq[, .(cov=.N, ev=sum(event==1L), ev_pct=100*sum(event==1L)/.N),
                      by=.(AA_anticodon_N1, AA_anticodon, source, group, position)]
site_rates_grp <- site_rates_grp[cov >= MIN_READS_PER_CELL]

site_cmp <- data.table::dcast(site_rates_grp,
                              AA_anticodon_N1 + AA_anticodon + source + position ~ group,
                              value.var="ev_pct", fill=NA_real_)
if (all(c("poly","input") %in% names(site_cmp))) {
  site_cmp[, log2FC_poly_vs_input := log2( (poly + eps) / (input + eps) ) ]
} else {
  site_cmp[, log2FC_poly_vs_input := NA_real_]
}
fwrite(site_cmp, file.path(OUT_DIR,"tables","eventrate_poly_vs_input_by_position.tsv"), sep="\t")

# --------------------------- WRAP UP ------------------------------------------

msg("Writing sessionInfo and finishing…")
capture.output(sessionInfo(), file = file.path(OUT_DIR, "sessionInfo_polysome.txt"))
msg("Done.")
