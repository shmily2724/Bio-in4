#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# ---------- Helpers ----------
collapse_field <- function(x) {
  # Ép mọi kiểu list/CompressedList/AtomicList về chuỗi "a,b,c"
  if (is.null(x)) return(NA_character_)
  if (length(x) == 0) return(NA_character_)
  if (inherits(x, "CompressedList") || is.list(x)) {
    return(sapply(x, function(e) paste(as.character(e), collapse = ",")))
  }
  as.character(x)
}

parse_ann_character <- function(ann_char) {
  # ann_char có thể là NA hoặc "" hoặc "A|missense|...,..."
  if (is.na(ann_char) || ann_char == "") return(NULL)
  recs <- unlist(strsplit(ann_char, ",", fixed = TRUE), use.names = FALSE)
  if (length(recs) == 0) return(NULL)

  rows <- lapply(recs, function(rec) {
    fields <- strsplit(rec, "\\|", perl = TRUE)[[1]]
    # Chuẩn SnpEff tối đa 18 trường → pad/truncate về 18
    if (length(fields) >= 18) {
      fields <- fields[1:18]
    } else {
      length(fields) <- 18
    }
    fields
  })
  do.call(rbind, rows)
}

ann_colnames_18 <- c(
  "Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID",
  "Feature_Type","Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p",
  "cDNA.pos.length","CDS.pos.length","AA.pos.length","Distance","Errors",
  "Warnings","Info"
)

impact_order <- c("HIGH","MODERATE","LOW","MODIFIER")

# ---------- Core per-sample ----------
process_one_sample <- function(sample_id,
                               project_dir = ".",
                               genome = "hg38",
                               out_dir_name = "tables") {

  # B1. Tìm file VCF (ưu tiên .vcf.gz, fallback .vcf)
  vcf_gz <- file.path(project_dir, "results", sample_id, "snpeff",
                      paste0(sample_id, "_final_annotated.vcf.gz"))
  vcf_plain <- file.path(project_dir, "results", sample_id, "snpeff",
                         paste0(sample_id, "_final_annotated.vcf"))

  if (file.exists(vcf_gz)) {
    vcf_file <- vcf_gz
  } else if (file.exists(vcf_plain)) {
    vcf_file <- vcf_plain
  } else {
    stop("Không tìm thấy file VCF: ", vcf_gz, " hoặc ", vcf_plain)
  }

  message("▶ Đọc: ", vcf_file)
  vcf <- readVcf(vcf_file, genome = genome)

  nvar <- nrow(vcf)
  if (is.null(nvar) || nvar == 0) {
    stop("VCF không có biến thể nào: ", vcf_file)
  }

  # B2. Fixed (REF/ALT/QUAL/FILTER) + vị trí
  rr <- rowRanges(vcf)
  fx <- fixed(vcf)

  CHROM  <- as.character(seqnames(rr))
  POS    <- start(rr)
  ID     <- names(rr)  # có thể NA nếu ID trống
  REF    <- collapse_field(fx$REF)
  ALT    <- sapply(as.list(fx$ALT), function(x) paste(as.character(x), collapse = ","))
  QUAL   <- as.numeric(fx$QUAL)
  FILTER <- collapse_field(fx$FILTER)

  variant_df <- tibble(
    VAR_IDX = seq_len(nvar),
    CHROM, POS, ID, REF, ALT, QUAL, FILTER
  )

  # B3. INFO fields quan trọng (có thì lấy, không có thì NA)
  inf <- info(vcf)
  grab <- function(key)
    if (key %in% names(inf)) collapse_field(inf[[key]]) else rep(NA_character_, nvar)

  info_df <- tibble(
    VAR_IDX = seq_len(nvar),
    AF      = grab("AF"),
    AC      = grab("AC"),
    AN      = grab("AN"),
    CLNSIG  = grab("CLNSIG"),
    CLNDN   = grab("CLNDN"),
    LOF     = grab("LOF"),
    NMD     = grab("NMD")
  )

  # B4. Parse ANN → bảng expanded (mỗi transcript/annotation 1 dòng)
  ann_list <- if ("ANN" %in% names(inf)) inf$ANN else NULL

  ann_expanded <- NULL
  if (!is.null(ann_list)) {
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
  }

  # Nếu không có ANN, vẫn tạo bảng primary = variant + info
  out_dir <- file.path(project_dir, "results", sample_id, out_dir_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(ann_expanded) || nrow(ann_expanded) == 0) {
    warning("VCF không chứa trường ANN hoặc ANN rỗng — chỉ tạo bảng biến thể.")
    primary <- variant_df %>% left_join(info_df, by = "VAR_IDX")
    write_csv(primary, file.path(out_dir, paste0(sample_id, "_variants_primary.csv")))
    write_csv(primary, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.csv")))
    message("✅ Xuất (no-ANN) xong: ", out_dir)
    return(invisible(TRUE))
  }

  # B5. Expanded table = variant + info + từng ANN
  expanded <- ann_expanded %>%
    left_join(variant_df, by = "VAR_IDX") %>%
    left_join(info_df,    by = "VAR_IDX") %>%
    relocate(VAR_IDX, CHROM, POS, ID, REF, ALT, QUAL, FILTER)

  # B6. Chọn annotation “đại diện” cho mỗi biến thể → primary
  expanded$Annotation_Impact <- toupper(expanded$Annotation_Impact)
  expanded$imp_rank <- match(expanded$Annotation_Impact, impact_order)
  expanded$imp_rank[is.na(expanded$imp_rank)] <- length(impact_order) + 1

  primary <- expanded %>%
    group_by(VAR_IDX) %>%
    arrange(imp_rank,
            desc(!is.na(`HGVS.p`)), desc(!is.na(`HGVS.c`))) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(-imp_rank)

  # B7. Ghi file
  write_csv(expanded, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.csv")))
  write_csv(primary,  file.path(out_dir, paste0(sample_id, "_variants_primary.csv")))

  message("✅ Xuất xong: ", out_dir)
  invisible(TRUE)
}

# ---------- Main: nhận nhiều SAMPLE_ID ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Cách dùng: Rscript brca_annot_to_table.R <SAMPLE_ID_1> [SAMPLE_ID_2 ...]")
}

for (sid in args) {
  try({
    process_one_sample(sid, project_dir = ".", genome = "hg38")
  }, silent = FALSE)
}
