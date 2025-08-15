#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# -------------------- Helpers --------------------

collapse_field <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (length(x) == 0) return(NA_character_)
  if (inherits(x, "CompressedList") || is.list(x)) {
    return(sapply(x, function(e) paste(as.character(e), collapse = ",")))
  }
  as.character(x)
}

parse_ann_character <- function(ann_char) {
  if (is.na(ann_char) || ann_char == "") return(NULL)
  recs <- unlist(strsplit(ann_char, ",", fixed = TRUE), use.names = FALSE)
  if (length(recs) == 0) return(NULL)
  rows <- lapply(recs, function(rec) {
    f <- strsplit(rec, "\\|", perl = TRUE)[[1]]
    if (length(f) >= 18) f <- f[1:18] else length(f) <- 18
    f
  })
  do.call(rbind, rows)
}

split_ratio_cols <- function(x) {
  m <- regexec("^(\\d+)/(\\d+)$", as.character(x))
  parts <- regmatches(as.character(x), m)
  tibble(
    num   = suppressWarnings(as.integer(ifelse(lengths(parts)>0, sapply(parts, `[`, 2), NA))),
    total = suppressWarnings(as.integer(ifelse(lengths(parts)>0, sapply(parts, `[`, 3), NA)))
  )
}

first_numeric <- function(s) {
  s <- as.character(s)
  s <- sub(",", ".", s, fixed = TRUE)
  m <- regexpr("[0-9]+(\\.[0-9]+)?", s, perl = TRUE)
  res <- ifelse(m > 0, substr(s, m, m + attr(m, "match.length") - 1), NA_character_)
  suppressWarnings(as.numeric(res))
}

protect_slash <- function(x) gsub("/", "\u2215", as.character(x), fixed = TRUE)

ann_colnames_18 <- c(
  "Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID",
  "Feature_Type","Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p",
  "cDNA.pos.length","CDS.pos.length","AA.pos.length","Distance","Errors",
  "Warnings","Info"
)

impact_order <- c("HIGH","MODERATE","LOW","MODIFIER")

# Xuất .xlsx – không gọi library(openxlsx) (tránh lỗi); dùng openxlsx::*
write_xlsx_text <- function(df, path, text_cols) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    message("ℹ Gợi ý: cài 'openxlsx' để xuất .xlsx (install.packages(\"openxlsx\"))")
    return(invisible(FALSE))
  }
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Sheet1")
  openxlsx::writeData(wb, "Sheet1", df)

  textStyle <- openxlsx::createStyle(numFmt = "@")
  cols_to_text <- intersect(text_cols, names(df))
  if (length(cols_to_text) > 0) {
    idx <- match(cols_to_text, names(df))
    for (j in idx) {
      openxlsx::addStyle(wb, "Sheet1", textStyle,
                         rows = 1:(nrow(df) + 1), cols = j, gridExpand = TRUE)
    }
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(TRUE)
}

# -------------------- Core --------------------

# Xử lý 1 VCF -> xuất bảng vào out_dir
process_one_vcf <- function(sample_id, vcf_file, genome = "hg38", out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  message("▶ Đọc: ", vcf_file)
  vcf <- readVcf(vcf_file, genome = genome)
  nvar <- nrow(vcf)
  if (is.null(nvar) || nvar == 0) stop("VCF không có biến thể: ", vcf_file)

  rr <- rowRanges(vcf); fx <- fixed(vcf)

  CHROM  <- as.character(seqnames(rr))
  POS    <- start(rr)
  ID     <- names(rr)
  REF    <- as.character(fx$REF)
  ALT    <- vapply(as.list(fx$ALT), function(x) paste(as.character(x), collapse = ","), character(1))
  QUAL   <- as.numeric(fx$QUAL)
  FILTER <- as.character(fx$FILTER)

  variant_df <- tibble(VAR_IDX = seq_len(nvar), CHROM, POS, ID, REF, ALT, QUAL, FILTER)

  inf  <- info(vcf)
  grab <- function(key) if (key %in% names(inf)) collapse_field(inf[[key]]) else rep(NA_character_, nvar)
  info_df <- tibble(
    VAR_IDX = seq_len(nvar),
    AF      = grab("AF"),
    AC      = grab("AC"),
    AN      = grab("AN"),
    CLNSIG  = grab("CLNSIG"),
    CLNDN   = grab("CLNDN"),
    LOF     = grab("LOF"),
    NMD     = grab("NMD")
  ) %>% mutate(AF_num = first_numeric(AF))

  ann_list <- if ("ANN" %in% names(inf)) inf$ANN else NULL
  if (is.null(ann_list)) {
    primary_only <- variant_df %>% left_join(info_df, by = "VAR_IDX")
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_primary.csv")))
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.csv")))
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_primary_excel.csv")))
    message("✅ Xuất (no-ANN): ", out_dir)
    return(invisible(TRUE))
  }

  ann_strings <- collapse_field(ann_list)
  rows <- vector("list", nvar)
  for (i in seq_len(nvar)) {
    mat <- parse_ann_character(ann_strings[i])
    if (is.null(mat)) next
    df_i <- as.data.frame(mat, stringsAsFactors = FALSE)
    colnames(df_i) <- ann_colnames_18
    df_i$VAR_IDX <- i
    rows[[i]] <- df_i
  }
  ann_expanded <- bind_rows(rows)

  if (is.null(ann_expanded) || nrow(ann_expanded) == 0) {
    primary_only <- variant_df %>% left_join(info_df, by = "VAR_IDX")
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_primary.csv")))
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.csv")))
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_primary_excel.csv")))
    message("✅ Xuất (empty-ANN): ", out_dir)
    return(invisible(TRUE))
  }

  expanded <- ann_expanded %>%
    left_join(variant_df, by = "VAR_IDX") %>%
    left_join(info_df,    by = "VAR_IDX") %>%
    relocate(VAR_IDX, CHROM, POS, ID, REF, ALT, QUAL, FILTER)

  # Thêm Rank_num/Rank_total nếu chưa có
  if ("Rank" %in% names(expanded) && !("Rank_num" %in% names(expanded))) {
    rank_split <- split_ratio_cols(expanded$Rank)
    names(rank_split) <- c("Rank_num", "Rank_total")
    expanded <- bind_cols(expanded, rank_split)
  }

  # Chọn annotation đại diện → primary
  expanded$Annotation_Impact <- toupper(expanded$Annotation_Impact)
  expanded$imp_rank <- match(expanded$Annotation_Impact, impact_order)
  expanded$imp_rank[is.na(expanded$imp_rank)] <- length(impact_order) + 1

  primary <- expanded %>%
    group_by(VAR_IDX) %>%
    arrange(imp_rank, desc(!is.na(`HGVS.p`)), desc(!is.na(`HGVS.c`))) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(-imp_rank)

  # KHÔNG bind lại Rank_num/Rank_total ở primary (đã có từ expanded)

  # Xuất
  write_csv(expanded, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.csv")))
  write_csv(primary,  file.path(out_dir, paste0(sample_id, "_variants_primary.csv")))

  ratio_cols <- c("Rank","cDNA.pos.length","CDS.pos.length","AA.pos.length")
  expanded_excel <- expanded; primary_excel <- primary
  for (cc in ratio_cols) {
    if (cc %in% names(expanded_excel)) expanded_excel[[cc]] <- protect_slash(expanded_excel[[cc]])
    if (cc %in% names(primary_excel))  primary_excel[[cc]]  <- protect_slash(primary_excel[[cc]])
  }
  write_csv(expanded_excel, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded_excel.csv")))
  write_csv(primary_excel,  file.path(out_dir, paste0(sample_id, "_variants_primary_excel.csv")))

  text_cols <- unique(c(ratio_cols,"Gene_Name","Feature_ID","HGVS.c","HGVS.p","ID","REF","ALT","FILTER","CLNSIG","CLNDN"))
  write_xlsx_text(expanded, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.xlsx")), text_cols)
  write_xlsx_text(primary,  file.path(out_dir, paste0(sample_id, "_variants_primary.xlsx")),        text_cols)

  message("✅ Xuất xong: ", out_dir)
  invisible(TRUE)
}

# Tìm VCF (gatk/dv) cho 1 sample theo chế độ which
discover_sample_vcfs <- function(project_dir, sample_id, which = c("both","gatk","dv")) {
  which <- match.arg(which)
  base <- file.path(project_dir, "results", sample_id, "snpeff")
  cand <- list(
    gatk_gz = file.path(base, paste0(sample_id, "_final_annotated.vcf.gz")),
    gatk    = file.path(base, paste0(sample_id, "_final_annotated.vcf")),
    dv_gz   = file.path(base, paste0(sample_id, "_dv_final_annotated.vcf.gz")),
    dv      = file.path(base, paste0(sample_id, "_dv_final_annotated.vcf"))
  )
  have <- cand[file.exists(unlist(cand))]
  out <- list()
  if (which %in% c("both","gatk")) {
    if (!is.null(have$gatk_gz)) out <- c(out, list(list(tag="gatk", path=have$gatk_gz)))
    else if (!is.null(have$gatk)) out <- c(out, list(list(tag="gatk", path=have$gatk)))
  }
  if (which %in% c("both","dv")) {
    if (!is.null(have$dv_gz)) out <- c(out, list(list(tag="dv", path=have$dv_gz)))
    else if (!is.null(have$dv)) out <- c(out, list(list(tag="dv", path=have$dv)))
  }
  out
}

# -------------------- Main --------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Cách dùng:\n",
      "  PROJECT_DIR=/path/to/project Rscript brca_annot_to_table.R [--which gatk|dv|both] SAMPLE1 [SAMPLE2 ...]\n",
      "  (mặc định --which=both)\n", sep = "")
  quit(status = 1)
}

which_mode <- "both"
if (grepl("^--which=", args[1])) {
  which_mode <- sub("^--which=", "", args[1])
  args <- args[-1]
}

project_dir <- Sys.getenv("PROJECT_DIR", ".")  # set từ shell

for (sid in args) {
  vcfs <- discover_sample_vcfs(project_dir, sid, which = which_mode)
  if (length(vcfs) == 0) {
    message("❌ Không tìm thấy VCF cho ", sid, " (", which_mode, ")")
    next
  }
  for (v in vcfs) {
    out_dir_name <- if (v$tag == "dv") "tables_dv" else "tables_gatk"
    out_dir <- file.path(project_dir, "results", sid, out_dir_name)
    try({
      process_one_vcf(sid, v$path, genome = "hg38", out_dir = out_dir)
    }, silent = FALSE)
  }
}
