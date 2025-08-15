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

# Ép mọi kiểu list/CompressedList về vector ký tự, ghép bằng ","
collapse_field <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (length(x) == 0) return(NA_character_)
  if (inherits(x, "CompressedList") || is.list(x)) {
    return(sapply(x, function(e) paste(as.character(e), collapse = ",")))
  }
  as.character(x)
}

# Parse một chuỗi ANN (có thể chứa nhiều entry ngăn bởi ",")
# Mỗi entry tách theo "|" và pad/truncate về đúng 18 trường SnpEff
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

# Tách "n/m" -> 2 cột số
split_ratio_cols <- function(x) {
  m <- regexec("^(\\d+)/(\\d+)$", as.character(x))
  parts <- regmatches(as.character(x), m)
  tibble(
    num   = suppressWarnings(as.integer(ifelse(lengths(parts)>0, sapply(parts, `[`, 2), NA))),
    total = suppressWarnings(as.integer(ifelse(lengths(parts)>0, sapply(parts, `[`, 3), NA)))
  )
}

# Lấy số thực đầu tiên trong chuỗi (dùng cho AF nếu có nhiều giá trị)
first_numeric <- function(s) {
  s <- as.character(s)
  s <- sub(",", ".", s, fixed = TRUE) # trong trường hợp dùng dấu phẩy thập phân
  m <- regexpr("[0-9]+(\\.[0-9]+)?", s, perl = TRUE)
  res <- ifelse(m > 0, substr(s, m, m + attr(m, "match.length") - 1), NA_character_)
  suppressWarnings(as.numeric(res))
}

# Đổi "/" thành "∕" (U+2215) để Excel/Sheets không ép thành ngày
protect_slash <- function(x) gsub("/", "\u2215", as.character(x), fixed = TRUE)

ann_colnames_18 <- c(
  "Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID",
  "Feature_Type","Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p",
  "cDNA.pos.length","CDS.pos.length","AA.pos.length","Distance","Errors",
  "Warnings","Info"
)

impact_order <- c("HIGH","MODERATE","LOW","MODIFIER")

# Xuất .xlsx với một số cột ép Text (nếu có openxlsx)
write_xlsx_text <- function(df, path, text_cols) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    message("ℹ Gợi ý: cài 'openxlsx' để xuất .xlsx (install.packages(\"openxlsx\"))")
    return(invisible(FALSE))
  }
  openxlsx::library(openxlsx)
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Sheet1")
  openxlsx::writeData(wb, "Sheet1", df)
  textStyle <- openxlsx::createStyle(numFmt = "@")
  for (nm in intersect(text_cols, names(df))) {
    j <- match(nm, names(df))
    openxlsx::addStyle(wb, "Sheet1", textStyle, rows = 1:(nrow(df)+1), cols = j, gridExpand = TRUE)
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(TRUE)
}

# -------------------- Core per-sample --------------------

process_one_sample <- function(sample_id,
                               project_dir = ".",
                               genome = "hg38",
                               out_dir_name = "tables") {

  # Ưu tiên .vcf.gz, fallback .vcf
  vcf_gz    <- file.path(project_dir, "results", sample_id, "snpeff",
                         paste0(sample_id, "_final_annotated.vcf.gz"))
  vcf_plain <- file.path(project_dir, "results", sample_id, "snpeff",
                         paste0(sample_id, "_final_annotated.vcf"))

  if (file.exists(vcf_gz)) {
    vcf_file <- vcf_gz
  } else if (file.exists(vcf_plain)) {
    vcf_file <- vcf_plain
  } else {
    stop("Không tìm thấy VCF: ", vcf_gz, " hoặc ", vcf_plain)
  }

  message("▶ Đọc: ", vcf_file)
  vcf <- readVcf(vcf_file, genome = genome)
  nvar <- nrow(vcf)
  if (is.null(nvar) || nvar == 0) stop("VCF không có biến thể: ", vcf_file)

  # ---- Vị trí & Fixed ----
  rr <- rowRanges(vcf)
  fx <- fixed(vcf)

  CHROM  <- as.character(seqnames(rr))
  POS    <- start(rr)
  ID     <- names(rr)
  REF    <- as.character(fx$REF)
  ALT    <- vapply(as.list(fx$ALT), function(x) paste(as.character(x), collapse = ","), character(1))
  QUAL   <- as.numeric(fx$QUAL)
  FILTER <- as.character(fx$FILTER)

  variant_df <- tibble(
    VAR_IDX = seq_len(nvar),
    CHROM, POS, ID, REF, ALT, QUAL, FILTER
  )

  # ---- INFO (mềm: có thì lấy, không có thì NA) ----
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
  ) %>%
    mutate(
      AF_num = first_numeric(AF)  # tiện lọc số
    )

  # ---- Parse ANN → expanded (mỗi transcript/allele 1 dòng) ----
  ann_list <- if ("ANN" %in% names(inf)) inf$ANN else NULL
  out_dir  <- file.path(project_dir, "results", sample_id, out_dir_name)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(ann_list)) {
    warning("Không có trường ANN — tạo bảng primary từ variant + info.")
    primary_only <- variant_df %>% left_join(info_df, by = "VAR_IDX")
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_primary.csv")))
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.csv")))
    # Excel-safe CSV
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
    warning("ANN rỗng — tạo bảng primary từ variant + info.")
    primary_only <- variant_df %>% left_join(info_df, by = "VAR_IDX")
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_primary.csv")))
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.csv")))
    write_csv(primary_only, file.path(out_dir, paste0(sample_id, "_variants_primary_excel.csv")))
    message("✅ Xuất (empty-ANN): ", out_dir)
    return(invisible(TRUE))
  }

  # ---- Expanded = variant + info + ANN ----
  expanded <- ann_expanded %>%
    left_join(variant_df, by = "VAR_IDX") %>%
    left_join(info_df,    by = "VAR_IDX") %>%
    relocate(VAR_IDX, CHROM, POS, ID, REF, ALT, QUAL, FILTER)

  # Thêm Rank_num/Rank_total (hữu ích cho lọc)
  if ("Rank" %in% names(expanded)) {
    rank_split <- split_ratio_cols(expanded$Rank)
    names(rank_split) <- c("Rank_num", "Rank_total")
    expanded <- bind_cols(expanded, rank_split)
  }

  # ---- Chọn annotation "đại diện" → primary ----
  expanded$Annotation_Impact <- toupper(expanded$Annotation_Impact)
  expanded$imp_rank <- match(expanded$Annotation_Impact, impact_order)
  expanded$imp_rank[is.na(expanded$imp_rank)] <- length(impact_order) + 1

  primary <- expanded %>%
    group_by(VAR_IDX) %>%
    arrange(imp_rank, desc(!is.na(`HGVS.p`)), desc(!is.na(`HGVS.c`))) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(-imp_rank)

  if ("Rank" %in% names(primary)) {
    rank_split_p <- split_ratio_cols(primary$Rank)
    names(rank_split_p) <- c("Rank_num", "Rank_total")
    primary <- bind_cols(primary, rank_split_p)
  }

  # -------------------- Xuất file --------------------

  # CSV "gốc"
  write_csv(expanded, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.csv")))
  write_csv(primary,  file.path(out_dir, paste0(sample_id, "_variants_primary.csv")))

  # CSV "Excel-safe": đổi "/" thành "∕" ở 4 cột dễ bị hiểu là ngày
  ratio_cols <- c("Rank", "cDNA.pos.length", "CDS.pos.length", "AA.pos.length")
  expanded_excel <- expanded
  primary_excel  <- primary
  for (cc in ratio_cols) {
    if (cc %in% names(expanded_excel)) expanded_excel[[cc]] <- protect_slash(expanded_excel[[cc]])
    if (cc %in% names(primary_excel))  primary_excel[[cc]]  <- protect_slash(primary_excel[[cc]])
  }
  write_csv(expanded_excel, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded_excel.csv")))
  write_csv(primary_excel,  file.path(out_dir, paste0(sample_id, "_variants_primary_excel.csv")))

  # XLSX (nếu có openxlsx): ép Text cho các cột hay bị "méo"
  text_cols <- unique(c(
    ratio_cols,
    "Gene_Name","Feature_ID","HGVS.c","HGVS.p","ID","REF","ALT","FILTER","CLNSIG","CLNDN"
  ))
  write_xlsx_text(expanded, file.path(out_dir, paste0(sample_id, "_variants_ANN_expanded.xlsx")), text_cols)
  write_xlsx_text(primary,  file.path(out_dir, paste0(sample_id, "_variants_primary.xlsx")),        text_cols)

  message("✅ Xuất xong: ", out_dir)
  invisible(TRUE)
}

# -------------------- Main --------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Cách dùng: Rscript brca_annot_to_table.R <SAMPLE1> [SAMPLE2 ...]")
}

project_dir <- Sys.getenv("PROJECT_DIR", ".")  # có thể set PROJECT_DIR từ shell
for (sid in args) {
  try({
    process_one_sample(sid, project_dir = project_dir, genome = "hg38")
  }, silent = FALSE)
}
