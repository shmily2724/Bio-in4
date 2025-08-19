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

# ======================= Helpers =======================

impact_order <- c("HIGH","MODERATE","LOW","MODIFIER")

ann_colnames_18 <- c(
  "Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID",
  "Feature_Type","Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p",
  "cDNA.pos.length","CDS.pos.length","AA.pos.length","Distance","Errors",
  "Warnings","Info"
)

ensure_len <- function(x, n, fill = NA_character_) {
  if (is.null(x)) return(rep(fill, n))
  lx <- length(x)
  if (lx == n) return(x)
  if (lx == 0) return(rep(fill, n))
  if (n %% lx == 0) return(rep(x, length.out = n))
  rep(fill, n)
}

collapse_field <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (length(x) == 0) return(NA_character_)
  if (inherits(x, "CompressedList") || is.list(x)) {
    return(vapply(x, function(e) paste(as.character(e), collapse = ","), character(1)))
  }
  as.character(x)
}

collapse_field_n <- function(x, n) {
  if (is.null(x) || length(x) == 0) return(rep(NA_character_, n))
  if (inherits(x, "CompressedList") || is.list(x)) {
    v <- vapply(x, function(e) paste(as.character(e), collapse=","), character(1))
    return(ensure_len(as.character(v), n))
  }
  ensure_len(as.character(x), n)
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

# Đổi "/" thường thành "∕" (U+2215) để Excel không hiểu là ngày
protect_slash <- function(x) gsub("/", "\u2215", as.character(x), fixed = TRUE)

# ---- Safe extractors ----
safe_num_vec <- function(x, n) {
  if (is.null(x)) return(rep(NA_real_, n))
  if (inherits(x, "CompressedList") || is.list(x)) {
    return(vapply(x, function(e) suppressWarnings(as.numeric(e)[1]), numeric(1)))
  }
  out <- suppressWarnings(as.numeric(x))
  if (length(out) != n) out <- rep_len(out, n)
  out
}

safe_chr_vec <- function(x, n) {
  if (is.null(x)) return(rep(NA_character_, n))
  if (inherits(x, "CompressedList") || is.list(x)) {
    return(vapply(x, function(e) paste(as.character(e), collapse=","), character(1)))
  }
  out <- as.character(x)
  if (length(out) != n) out <- rep_len(out, n)
  out
}

safe_filter_vec <- function(f, n) {
  if (inherits(f, "CompressedList") || is.list(f)) {
    return(vapply(f, function(e) paste(as.character(e), collapse=";"), character(1)))
  }
  out <- as.character(f)
  if (length(out) != n) out <- rep_len(out, n)
  out
}

safe_alt_parts <- function(ALT_field, n) {
  ALT_LIST <- tryCatch(lapply(as.list(ALT_field), function(x) as.character(x)),
                       error = function(e) replicate(n, character(0), simplify = FALSE))
  ALT_STR  <- tryCatch(vapply(ALT_LIST, function(x) paste(x, collapse=","), character(1)),
                       error = function(e) rep(NA_character_, n))
  list(ALT_LIST = ALT_LIST, ALT_STR = ALT_STR)
}

pick_af_for_alt <- function(af_str, alt_selected, alt_list) {
  if (is.na(af_str) || is.na(alt_selected) || is.null(alt_list)) return(NA_real_)
  parts <- unlist(strsplit(af_str, ",", fixed=TRUE), use.names = FALSE)
  if (length(parts) == 0) return(NA_real_)
  idx <- match(alt_selected, alt_list)
  if (is.na(idx) || idx < 1 || idx > length(parts)) return(NA_real_)
  suppressWarnings(as.numeric(parts[idx]))
}

vaf_from_AD <- function(count_vec, alt_idx) {
  if (is.null(count_vec) || length(count_vec) < (alt_idx + 1)) return(NA_real_)
  counts <- suppressWarnings(as.numeric(count_vec))
  if (any(is.na(counts))) return(NA_real_)
  tot <- sum(counts)
  if (!is.finite(tot) || tot <= 0) return(NA_real_)
  alt <- counts[alt_idx + 1]  # pos 1 là REF
  alt / tot
}

maf_from_af_vec <- function(af) {
  x <- suppressWarnings(as.numeric(af))
  ifelse(is.na(x), NA_real_, pmin(x, 1 - x))
}

grab_any <- function(info_list, keys, nvar) {
  for (k in keys) {
    if (k %in% names(info_list)) {
      val <- collapse_field(info_list[[k]])
      return(ensure_len(val, nvar))
    }
  }
  rep(NA_character_, nvar)
}

# Xuất .xlsx với định dạng Text cho các cột dễ bị Excel auto-format
write_xlsx_text <- function(df, path, text_cols, ratio_cols = character()) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    message("ℹ Chưa cài 'openxlsx' → chỉ xuất CSV. Cài: install.packages('openxlsx')")
    return(invisible(FALSE))
  }
  df2 <- df
  for (cc in intersect(ratio_cols, names(df2))) df2[[cc]] <- protect_slash(df2[[cc]])
  for (cc in intersect(text_cols, names(df2)))  df2[[cc]] <- as.character(df2[[cc]])
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "variants")
  openxlsx::writeData(wb, "variants", df2, keepNA = TRUE)
  textStyle <- openxlsx::createStyle(numFmt = "@")
  idx <- match(intersect(text_cols, names(df2)), names(df2))
  if (length(idx) > 0) for (j in idx) openxlsx::addStyle(wb, "variants", textStyle,
                                                          rows = 1:(nrow(df2)+1), cols = j, gridExpand = TRUE)
  openxlsx::addFilter(wb, "variants", rows = 1, cols = 1:ncol(df2))
  openxlsx::freezePane(wb, "variants", firstActiveRow = 2, firstActiveCol = 1)
  suppressWarnings(openxlsx::setColWidths(wb, "variants", cols = 1:ncol(df2), widths = "auto"))
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  invisible(TRUE)
}

# Ẩn cảnh báo "duplicate keys in header"
readVcf_quiet <- function(path, genome) {
  withCallingHandlers(
    readVcf(path, genome = genome),
    warning = function(w) {
      if (grepl("duplicate keys in header", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

# ======================= Core =======================

extract_and_write <- function(sample_id, vcf_file, tag = c("gatk","dv"),
                              genome = "hg38", out_dir) {
  tag <- match.arg(tag)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  message("▶ Đọc VCF (", tag, "): ", vcf_file)
  vcf <- readVcf_quiet(vcf_file, genome = genome)
  nvar <- nrow(vcf)
  if (is.null(nvar) || nvar == 0) stop("VCF không có biến thể: ", vcf_file)

  rr <- rowRanges(vcf); fx <- fixed(vcf)
  var_idx <- seq_len(nvar)

  # ALT (mọi kiểu) → ALT_LIST + ALT_STR
  alt_parts <- safe_alt_parts(fx$ALT, nvar)
  ALT_LIST <- alt_parts$ALT_LIST
  ALT_STR  <- alt_parts$ALT_STR

  # ID/REF/QUAL/FILTER
  ID <- safe_chr_vec(names(rr), nvar); ID[ID %in% c(".", "")] <- NA_character_
  REF    <- safe_chr_vec(fx$REF,    nvar)
  QUAL   <- safe_num_vec(fx$QUAL,   nvar)
  FILTER <- safe_filter_vec(fx$FILTER, nvar)

  variant_df <- tibble(
    VAR_IDX = var_idx,
    CHROM   = as.character(seqnames(rr)),
    POS     = as.integer(start(rr)),
    ID      = ID,
    REF     = REF,
    ALT     = ALT_STR,
    QUAL    = QUAL,
    FILTER  = FILTER
  )

  # -------- INFO (bỏ AF tổng quát, chỉ chọn theo allele) --------
  inf <- info(vcf)
  AC        <- collapse_field_n(if ("AC" %in% names(inf)) inf$AC else NULL, nvar)
  AN        <- collapse_field_n(if ("AN" %in% names(inf)) inf$AN else NULL, nvar)
  GNOMAD_AF <- grab_any(inf, c("GNOMAD_AF","gnomAD_AF","gnomad_AF","gnomADg_AF","gnomadg_AF"), nvar)
  KG_AF     <- grab_any(inf, c("KG_AF","1000G_AF","1KG_AF","1000g2015aug_all"), nvar)
  CLNSIG    <- grab_any(inf, c("CLNSIG","CLIN_SIG","CLNSIGCONF"), nvar)
  CLNDN     <- grab_any(inf, c("CLNDN","CLNDISDB","CLNDISDBINCL"), nvar)
  LOF       <- collapse_field_n(if ("LOF" %in% names(inf)) inf$LOF else NULL, nvar)
  NMD       <- collapse_field_n(if ("NMD" %in% names(inf)) inf$NMD else NULL, nvar)

  info_df <- tibble(VAR_IDX = var_idx, AC, AN, GNOMAD_AF, KG_AF, CLNSIG, CLNDN, LOF, NMD)

  # -------- ANN --------
  ann_list <- if ("ANN" %in% names(inf)) inf$ANN else NULL

  # Các cột cần Text trong Excel
  excel_text_cols <- c(
    "CHROM","ID","REF","ALT","FILTER",
    "Gene_Name","Annotation","Annotation_Impact","Transcript_BioType",
    "Rank","HGVS.c","HGVS.p","CLNSIG","CLNDN","AF_source"
  )
  excel_ratio_cols <- c("Rank","cDNA.pos.length","CDS.pos.length","AA.pos.length")

  # Hàm xuất (duy nhất)
  export_csv_xlsx <- function(df_core, note_tag) {
    suffix <- if (tag == "dv") "_dv_variants_full" else "_variants_full"
    out_csv  <- file.path(out_dir, paste0(sample_id, suffix, ".csv"))
    out_xlsx <- file.path(out_dir, paste0(sample_id, suffix, ".xlsx"))
    write_csv(df_core, out_csv)
    write_xlsx_text(df_core, out_xlsx, text_cols = excel_text_cols, ratio_cols = excel_ratio_cols)
    message("✅ Xuất (", note_tag, "): ", out_csv, " & .xlsx")
  }

  # --- Không có ANN ---
  if (is.null(ann_list)) {
    out_min <- variant_df %>%
      left_join(info_df, by = "VAR_IDX") %>%
      transmute(
        CHROM, POS, REF, ALT, ID, QUAL, FILTER,
        Gene_Name = NA_character_, Annotation = NA_character_, Annotation_Impact = NA_character_,
        Transcript_BioType = NA_character_,
        Rank = NA_character_, `HGVS.c` = NA_character_, `HGVS.p` = NA_character_,
        GNOMAD_AF, KG_AF,
        GNOMAD_AF_sel = NA_real_, KG_AF_sel = NA_real_,
        AF_source = "none",
        VAF = NA_real_, MAF = NA_real_,
        CLNSIG, CLNDN, LOF, NMD
      )
    export_csv_xlsx(out_min, "no-ANN")
    return(invisible(TRUE))
  }

  # --- Có ANN: bung & chọn annotation đại diện ---
  ann_strings <- ensure_len(collapse_field(ann_list), nvar)
  rows <- vector("list", nvar)
  for (i in seq_len(nvar)) {
    mat <- parse_ann_character(ann_strings[i])
    if (!is.null(mat)) {
      df_i <- as.data.frame(mat, stringsAsFactors = FALSE)
      colnames(df_i) <- ann_colnames_18
      df_i$VAR_IDX <- i
      rows[[i]] <- df_i
    }
  }
  ann_expanded <- bind_rows(rows)

  if (is.null(ann_expanded) || nrow(ann_expanded) == 0) {
    out_min <- variant_df %>%
      left_join(info_df, by = "VAR_IDX") %>%
      transmute(
        CHROM, POS, REF, ALT, ID, QUAL, FILTER,
        Gene_Name = NA_character_, Annotation = NA_character_, Annotation_Impact = NA_character_,
        Transcript_BioType = NA_character_,
        Rank = NA_character_, `HGVS.c` = NA_character_, `HGVS.p` = NA_character_,
        GNOMAD_AF, KG_AF,
        GNOMAD_AF_sel = NA_real_, KG_AF_sel = NA_real_,
        AF_source = "none",
        VAF = NA_real_, MAF = NA_real_,
        CLNSIG, CLNDN, LOF, NMD
      )
    export_csv_xlsx(out_min, "empty-ANN")
    return(invisible(TRUE))
  }

  expanded <- ann_expanded %>%
    left_join(variant_df, by = "VAR_IDX") %>%
    left_join(info_df,    by = "VAR_IDX") %>%
    relocate(VAR_IDX, CHROM, POS, ID, REF, ALT, QUAL, FILTER)

  if ("Rank" %in% names(expanded) && !("Rank_num" %in% names(expanded))) {
    rank_split <- split_ratio_cols(expanded$Rank)
    names(rank_split) <- c("Rank_num", "Rank_total")
    expanded <- bind_cols(expanded, rank_split)
  }

  expanded$Annotation_Impact <- toupper(expanded$Annotation_Impact)
  expanded$imp_rank <- match(expanded$Annotation_Impact, impact_order)
  expanded$imp_rank[is.na(expanded$imp_rank)] <- length(impact_order) + 1

  primary <- expanded %>%
    group_by(VAR_IDX) %>%
    arrange(imp_rank, desc(!is.na(`HGVS.p`)), desc(!is.na(`HGVS.c`))) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(-imp_rank)

  # ALT list cho từng biến thể
  primary$ALT_list <- ALT_LIST[primary$VAR_IDX]

  # AF theo đúng allele
  primary$GNOMAD_AF_sel <- purrr::pmap_dbl(
    list(primary$GNOMAD_AF, primary$Allele, primary$ALT_list),
    function(gn, al, alts) pick_af_for_alt(gn, al, alts)
  )
  primary$KG_AF_sel <- purrr::pmap_dbl(
    list(primary$KG_AF, primary$Allele, primary$ALT_list),
    function(kg, al, alts) pick_af_for_alt(kg, al, alts)
  )

  # VAF từ FORMAT/AD
  vaf_vec <- rep(NA_real_, nvar)
  ad <- geno(vcf)$AD
  if (!is.null(ad)) {
    if (is.array(ad) && length(dim(ad)) >= 3) {
      for (i in seq_len(nvar)) {
        counts <- suppressWarnings(as.numeric(ad[i, 1, ]))  # sample 1
        alts <- primary$ALT_list[[i]]
        alt_idx <- match(primary$Allele[i], alts)
        if (!is.na(alt_idx) && length(counts) >= (alt_idx + 1)) vaf_vec[i] <- vaf_from_AD(counts, alt_idx)
      }
    } else if (inherits(ad, "CompressedIntegerList") || inherits(ad, "IntegerList") || is.list(ad)) {
      ad_list <- as.list(ad)
      for (i in seq_len(min(length(ad_list), nvar))) {
        counts <- suppressWarnings(as.numeric(unlist(ad_list[[i]], recursive = TRUE, use.names = FALSE)))
        alts <- primary$ALT_list[[i]]
        alt_idx <- match(primary$Allele[i], alts)
        if (!is.na(alt_idx) && length(counts) >= (alt_idx + 1)) vaf_vec[i] <- vaf_from_AD(counts, alt_idx)
      }
    } else if (is.matrix(ad)) {
      for (i in seq_len(nvar)) {
        counts <- suppressWarnings(as.numeric(ad[i, ]))
        alts <- primary$ALT_list[[i]]
        alt_idx <- match(primary$Allele[i], alts)
        if (!is.na(alt_idx) && length(counts) >= (alt_idx + 1)) vaf_vec[i] <- vaf_from_AD(counts, alt_idx)
      }
    }
  }
  primary$VAF <- vaf_vec

  # MAF + AF_source
  primary <- primary %>%
    mutate(
      MAF = ifelse(
        !is.na(GNOMAD_AF_sel),
        maf_from_af_vec(GNOMAD_AF_sel),
        ifelse(!is.na(KG_AF_sel), maf_from_af_vec(KG_AF_sel), NA_real_)
      ),
      AF_source = ifelse(!is.na(GNOMAD_AF_sel), "gnomad",
                         ifelse(!is.na(KG_AF_sel), "1000g", "none"))
    )

  # Output (CSV + XLSX)
  out_min <- primary %>%
    transmute(
      CHROM, POS, REF, ALT, ID, QUAL, FILTER,
      Gene_Name, Annotation, Annotation_Impact, Transcript_BioType,
      Rank, `HGVS.c`, `HGVS.p`,
      GNOMAD_AF, KG_AF,
      GNOMAD_AF_sel, KG_AF_sel,
      AF_source, VAF, MAF,
      CLNSIG, CLNDN, LOF, NMD
    )

  export_csv_xlsx(out_min, note_tag = tag)
  invisible(TRUE)
}

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

# ======================= Main =======================

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
project_dir <- Sys.getenv("PROJECT_DIR", ".")

for (sid in args) {
  vcfs <- discover_sample_vcfs(project_dir, sid, which = which_mode)
  if (length(vcfs) == 0) {
    message("❌ Không tìm thấy VCF cho ", sid, " (", which_mode, ")")
    next
  }
  for (v in vcfs) {
    out_dir <- file.path(project_dir, "results", sid,
                         if (v$tag == "dv") "tables_dv" else "tables_gatk")
    try({
      extract_and_write(sid, v$path, tag = v$tag, genome = "hg38", out_dir = out_dir)
    }, silent = FALSE)
  }
}
