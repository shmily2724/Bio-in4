#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(fs)
  library(openxlsx)
})

# ======================= CONFIG =======================
BASE <- "/mnt/sarek_storage/sarek_results"

SAREK_DIRS <- list(
  gatk    = file.path(BASE, "table_gatk"),
  dv      = file.path(BASE, "table_dv"),
  mutect2 = file.path(BASE, "table_mt2")
)
MANUAL_DIRS <- list(
  gatk    = file.path(BASE, "results_tables_migrated/manual_table_gatk"),
  dv      = file.path(BASE, "results_tables_migrated/manual_table_dv"),
  mutect2 = file.path(BASE, "results_tables_migrated/manual_table_mt2")
)

OUT_BASE <- file.path(BASE, "compare")

# Thứ tự cột output cuối cùng
REQ_COLS <- c(
  "CHROM","POS","REF","ALT",
  "ID",
  "QUAL_sarek","QUAL_manual",
  "FILTER","Gene_Name","Gene_ID",
  "Annotation","Annotation_Impact","Transcript_BioType",
  "Rank_sarek","Rank_manual",
  "HGVS.c","HGVS.p",
  "VAF_sarek","VAF_manual",
  "LOF"
)

VAF_CANDIDATES <- c("VAF","Variant_AF","ALT_AF","AF","VAF(%)","TumorVAF","Allele_Fraction","AlleleFraction")
# ======================================================

# -------- Helpers --------
norm_col <- function(x) x |>
  str_trim() |>
  str_replace_all("\\s+", "") |>
  str_replace_all("[^0-9A-Za-z_\\.]+", "") |>
  tolower()
no_udot <- function(x) gsub("[_\\.]", "", x)

charify <- function(v) {
  v <- as.character(v); v[is.na(v)] <- ""; v
}
normalize_chrom <- function(v) { str_replace_all(charify(v), "(?i)^chr", "") }
to_int_safe <- function(v) suppressWarnings(as.integer(as.numeric(v)))

# -------- Column detection --------
detect_columns <- function(df) {
  nms <- names(df)
  nm_norm <- sapply(nms, norm_col)
  nm_noundot <- sapply(nm_norm, no_udot)
  map_exact <- setNames(nms, nm_norm)
  map_loose <- setNames(nms, nm_noundot)
  find_name <- function(targets, loosen_starts=NULL){
    for (t in targets){
      t1 <- norm_col(t); t2 <- no_udot(t1)
      if (t1 %in% names(map_exact)) return(map_exact[[t1]])
      if (t2 %in% names(map_loose)) return(map_loose[[t2]])
    }
    if (!is.null(loosen_starts)){
      p1 <- norm_col(loosen_starts); hits1 <- names(map_exact)[startsWith(names(map_exact), p1)]
      if (length(hits1)>0) return(map_exact[[hits1[1]]])
      p2 <- no_udot(p1); hits2 <- names(map_loose)[startsWith(names(map_loose), p2)]
      if (length(hits2)>0) return(map_loose[[hits2[1]]])
    }
    NULL
  }
  # phát hiện cột chính
  res <- list(
    chrom = find_name(c("chrom","chr","chromosome"), "chrom"),
    pos   = find_name(c("pos","position","start"),   "pos"),
    ref   = find_name(c("ref","refallele")),
    alt   = find_name(c("alt","altallele")),
    id    = find_name(c("id","variantid","rsid","dbsnpid","ID")),
    qual  = find_name(c("qual","quality","QUAL")),
    filter= find_name(c("filter","filters","pass","FILTER")),
    gene_name = find_name(c("gene_name","genename","gene","symbol")),
    gene_id   = find_name(c("gene_id","geneid","ensembl_gene_id","ensemblgeneid","ensg","entrezid")),
    annotation= find_name(c("annotation","ann","effect","consequence")),
    annotation_impact = find_name(c("annotation_impact","annotationimpact","impact")),
    transcript_biotype = find_name(c("transcript_biotype","transcriptbiotype","biotype")),
    rank = find_name(c("rank","Rank")),
    "hgvs.c" = find_name(c("hgvs.c","hgvsc","hgvs_c","cdna","hgvscdna")),
    "hgvs.p" = find_name(c("hgvs.p","hgvsp","hgvs_p","protein","hgvspprotein")),
    lof = find_name(c("lof","loss_of_function","lofpred","LOF"))
  )
  # VAF
  vf <- find_name(VAF_CANDIDATES)
  if (is.null(vf)) {
    nm_all <- c(names(map_exact), names(map_loose))
    if (any(grepl("vaf|allelefraction|af$", nm_all)))
      vf <- c(map_exact, map_loose)[[ nm_all[grepl("vaf|allelefraction|af$", nm_all)][1] ]]
  }
  res$vaf <- vf
  res
}

# -------- Load & standardize one CSV --------
load_and_prepare <- function(path) {
  if (!file_exists(path)) stop("File not found: ", path)
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  if (ncol(df) == 0) stop("Empty/invalid CSV: ", path)

  map <- detect_columns(df)
  pick <- function(key, default="") {
    col <- map[[key]]
    if (!is.null(col) && col %in% names(df)) df[[col]] else default
  }

  out <- tibble::tibble(
    CHROM = normalize_chrom(pick("chrom")),
    POS   = to_int_safe(pick("pos")),
    REF   = charify(pick("ref")),
    ALT   = charify(pick("alt")),
    ID    = charify(pick("id")),
    QUAL  = charify(pick("qual")),              # giữ dạng text để an toàn kiểu
    FILTER= charify(pick("filter")),
    Gene_Name = charify(pick("gene_name")),
    Gene_ID   = charify(pick("gene_id")),
    Annotation = charify(pick("annotation")),
    Annotation_Impact = charify(pick("annotation_impact")),
    Transcript_BioType = charify(pick("transcript_biotype")),
    Rank = charify(pick("rank")),
    `HGVS.c` = charify(pick("hgvs.c")),
    `HGVS.p` = charify(pick("hgvs.p")),
    LOF = charify(pick("lof"))
  ) %>% distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)

  list(data = out, mapping = map)
}

# -------- Attach VAF theo khóa --------
attach_vaf <- function(base_df, full_csv, vaf_col, new_name) {
  if (is.null(vaf_col)) { base_df[[new_name]] <- NA_real_; return(base_df) }
  full <- suppressMessages(readr::read_csv(full_csv, show_col_types = FALSE))
  m <- detect_columns(full)

  keys <- tibble::tibble(
    CHROM = normalize_chrom(if (!is.null(m$chrom)) full[[m$chrom]] else ""),
    POS   = to_int_safe(if (!is.null(m$pos))   full[[m$pos]]   else NA),
    REF   = charify(if (!is.null(m$ref)) full[[m$ref]] else ""),
    ALT   = charify(if (!is.null(m$alt)) full[[m$alt]] else "")
  )
  vaf_raw <- if (vaf_col %in% names(full)) full[[vaf_col]] else NA
  vaf_num <- suppressWarnings(as.numeric(gsub(",", ".", as.character(vaf_raw))))
  keys[[new_name]] <- vaf_num

  base_df %>% left_join(keys, by = c("CHROM","POS","REF","ALT"))
}

# -------- Write CSV + XLSX --------
write_dual <- function(df, csv_path, xlsx_path) {
  dir_create(dirname(csv_path), recurse = TRUE)
  readr::write_csv(df, csv_path, na = "")
  wb <- createWorkbook(); addWorksheet(wb, "Sheet1"); writeData(wb, "Sheet1", df)
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
}

# -------- Resolve filename by pipeline --------
resolve_variants_path <- function(base_dir, sample, pipeline) {
  cands <- switch(
    pipeline,
    "gatk"   = c(paste0(sample, "_variants_full.csv")),
    "dv"     = c(paste0(sample, "_dv_variants_full.csv"), paste0(sample, "_variants_full.csv")),
    "mutect2"= c(paste0(sample, "_mt2_variants_full.csv"), paste0(sample, "_variants_full.csv")),
    c(paste0(sample, "_variants_full.csv"))
  )
  for (fn in cands) {
    p <- file.path(base_dir, sample, fn)
    if (file_exists(p)) return(p)
  }
  file.path(base_dir, sample, cands[[1]])
}

# -------- Compare one sample --------
compare_one_sample <- function(pipeline, sample, sarek_csv, manual_csv, out_dir) {
  s <- load_and_prepare(sarek_csv); ds <- s$data
  m <- load_and_prepare(manual_csv); dm <- m$data

  # VAF
  vaf_s <- detect_columns(suppressMessages(readr::read_csv(sarek_csv, show_col_types = FALSE)))[["vaf"]]
  vaf_m <- detect_columns(suppressMessages(readr::read_csv(manual_csv, show_col_types = FALSE)))[["vaf"]]
  ds <- attach_vaf(ds, sarek_csv, vaf_s, "VAF_sarek")
  dm <- attach_vaf(dm, manual_csv, vaf_m, "VAF_manual")

  keys <- c("CHROM","POS","REF","ALT")

  # BOTH = biến thể chung
  joined <- ds %>% inner_join(dm, by = keys, suffix = c("_s","_m"))

  co_char <- function(a, b) dplyr::coalesce(charify(a), charify(b))

  # Lưu ý: cột có dấu chấm sau join sẽ thành "HGVS.c_s" / "HGVS.c_m"
  both <- tibble::tibble(
    CHROM = joined$CHROM, POS=joined$POS, REF=joined$REF, ALT=joined$ALT,

    ID = co_char(joined$ID_m, joined$ID_s),
    QUAL_sarek  = charify(joined$QUAL_s),
    QUAL_manual = charify(joined$QUAL_m),
    FILTER = co_char(joined$FILTER_m, joined$FILTER_s),
    Gene_Name = co_char(joined$Gene_Name_m, joined$Gene_Name_s),
    Gene_ID   = co_char(joined$Gene_ID_m,   joined$Gene_ID_s),

    Annotation        = co_char(joined$Annotation_m,        joined$Annotation_s),
    Annotation_Impact = co_char(joined$Annotation_Impact_m, joined$Annotation_Impact_s),
    Transcript_BioType= co_char(joined$Transcript_BioType_m, joined$Transcript_BioType_s),

    Rank_sarek  = charify(joined$Rank_s),
    Rank_manual = charify(joined$Rank_m),

    `HGVS.c` = co_char(joined[["HGVS.c_m"]], joined[["HGVS.c_s"]]),
    `HGVS.p` = co_char(joined[["HGVS.p_m"]], joined[["HGVS.p_s"]]),

    VAF_sarek  = joined$VAF_sarek,
    VAF_manual = joined$VAF_manual,

    LOF = co_char(joined$LOF_m, joined$LOF_s)
  ) %>% select(all_of(REQ_COLS))

  # ONLY SAREK
  only_sarek_idx <- anti_join(ds %>% select(all_of(keys)), dm %>% select(all_of(keys)), by = keys)
  only_sarek <- ds %>% semi_join(only_sarek_idx, by = keys)
  only_sarek_out <- tibble::tibble(
    CHROM=only_sarek$CHROM, POS=only_sarek$POS, REF=only_sarek$REF, ALT=only_sarek$ALT,
    ID=only_sarek$ID,
    QUAL_sarek=only_sarek$QUAL, QUAL_manual="",
    FILTER=only_sarek$FILTER,
    Gene_Name=only_sarek$Gene_Name, Gene_ID=only_sarek$Gene_ID,
    Annotation=only_sarek$Annotation, Annotation_Impact=only_sarek$Annotation_Impact,
    Transcript_BioType=only_sarek$Transcript_BioType,
    Rank_sarek=only_sarek$Rank, Rank_manual="",
    `HGVS.c`=only_sarek$`HGVS.c`, `HGVS.p`=only_sarek$`HGVS.p`,
    VAF_sarek=only_sarek$VAF_sarek, VAF_manual=NA_real_,
    LOF=only_sarek$LOF
  ) %>% select(all_of(REQ_COLS))

  # ONLY MANUAL
  only_manual_idx <- anti_join(dm %>% select(all_of(keys)), ds %>% select(all_of(keys)), by = keys)
  only_manual <- dm %>% semi_join(only_manual_idx, by = keys)
  only_manual_out <- tibble::tibble(
    CHROM=only_manual$CHROM, POS=only_manual$POS, REF=only_manual$REF, ALT=only_manual$ALT,
    ID=only_manual$ID,
    QUAL_sarek="", QUAL_manual=only_manual$QUAL,
    FILTER=only_manual$FILTER,
    Gene_Name=only_manual$Gene_Name, Gene_ID=only_manual$Gene_ID,
    Annotation=only_manual$Annotation, Annotation_Impact=only_manual$Annotation_Impact,
    Transcript_BioType=only_manual$Transcript_BioType,
    Rank_sarek="", Rank_manual=only_manual$Rank,
    `HGVS.c`=only_manual$`HGVS.c`, `HGVS.p`=only_manual$`HGVS.p`,
    VAF_sarek=NA_real_, VAF_manual=only_manual$VAF_manual,
    LOF=only_manual$LOF
  ) %>% select(all_of(REQ_COLS))

  # Ghi bảng
  dir_create(out_dir, recurse = TRUE)
  write_dual(both,             file.path(out_dir, paste0(sample, "_both.csv")),   file.path(out_dir, paste0(sample, "_both.xlsx")))
  write_dual(only_manual_out,  file.path(out_dir, paste0(sample, "_manual.csv")), file.path(out_dir, paste0(sample, "_manual.xlsx")))
  write_dual(only_sarek_out,   file.path(out_dir, paste0(sample, "_sarek.csv")),  file.path(out_dir, paste0(sample, "_sarek.xlsx")))

  # Summary cho sample
  sample_sum <- tibble::tibble(
    sample = sample, pipeline = pipeline,
    n_sarek = nrow(ds), n_manual = nrow(dm),
    n_common = nrow(both),
    n_only_sarek = nrow(only_sarek_out),
    n_only_manual = nrow(only_manual_out)
  )
  write_dual(sample_sum,
             file.path(out_dir, paste0("summary_sample_", pipeline, ".csv")),
             file.path(out_dir, paste0("summary_sample_", pipeline, ".xlsx")))

  c(n_sarek=nrow(ds), n_manual=nrow(dm), n_common=nrow(both),
    n_only_sarek=nrow(only_sarek_out), n_only_manual=nrow(only_manual_out))
}

# -------- Discover samples --------
find_samples <- function(dirpath){
  if (!dir_exists(dirpath)) return(character())
  entries <- dir_ls(dirpath, type = "directory", recurse = FALSE)
  sort(path_file(entries))
}

# ------------------ MAIN ------------------
main <- function(){
  cat("=== Compare Sarek vs Manual per pipeline (R) ===\n")
  for (pipeline in c("gatk","dv","mutect2")){
    cat("\n--- Pipeline:", pipeline, "---\n")
    sdir <- SAREK_DIRS[[pipeline]]
    mdir <- MANUAL_DIRS[[pipeline]]
    out_pipeline <- file.path(OUT_BASE, pipeline)

    samples_s <- find_samples(sdir)
    samples_m <- find_samples(mdir)
    common_samples <- intersect(samples_s, samples_m)
    if (length(common_samples) == 0){ cat("  ! No common samples between", sdir, "and", mdir, "\n"); next }

    rows <- list()
    for (smp in common_samples){
      sarek_csv  <- resolve_variants_path(sdir, smp, pipeline)
      manual_csv <- resolve_variants_path(mdir, smp, pipeline)
      out_dir <- file.path(out_pipeline, smp)

      missing <- c(if (!file_exists(sarek_csv)) sarek_csv, if (!file_exists(manual_csv)) manual_csv)
      if (length(missing) > 0){ cat("  ! Missing files for", smp, ":", paste(missing, collapse = ", "), "-> skipped\n"); next }

      res <- try(compare_one_sample(pipeline, smp, sarek_csv, manual_csv, out_dir), silent = TRUE)
      if (inherits(res, "try-error")){
        cat("  ! Error on sample", smp, ":", as.character(res), "\n")
      } else {
        rows[[length(rows)+1]] <- c(sample = smp, res)
        cat(sprintf("  ✓ %s: S=%s, M=%s, ∩=%s, S\\M=%s, M\\S=%s\n",
                    smp, res["n_sarek"], res["n_manual"], res["n_common"],
                    res["n_only_sarek"], res["n_only_manual"]))
      }
    }

    # summary gộp theo pipeline
    if (length(rows) > 0){
      df <- bind_rows(lapply(rows, as.data.frame.list)) %>%
        mutate(across(c(n_sarek,n_manual,n_common,n_only_sarek,n_only_manual), ~as.integer(.))) %>%
        select(sample, n_sarek, n_manual, n_common, n_only_sarek, n_only_manual)
      sum_dir <- file.path(out_pipeline, "sample")
      dir_create(sum_dir, recurse = TRUE)
      write_dual(df,
                 file.path(sum_dir, paste0("summary_", pipeline, ".csv")),
                 file.path(sum_dir, paste0("summary_", pipeline, ".xlsx")))
      cat("  → Wrote summary:", file.path(sum_dir, paste0("summary_", pipeline, ".{csv,xlsx}")), "\n")
    } else {
      cat("  ! No rows to summarize for this pipeline.\n")
    }
  }
  cat("\nDone.\n")
}

invisible(main())
