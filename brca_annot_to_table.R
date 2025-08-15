#!/usr/bin/env Rscript

# ========== BƯỚC 1: Đọc đối số ==========
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Vui lòng truyền vào SAMPLE_ID (ví dụ: ERR2228152)")
}
SAMPLE_ID <- args[1]

# ========== BƯỚC 2: Load thư viện ==========
suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(dplyr)
  library(readr)
})

# ========== BƯỚC 3: Xác định đường dẫn VCF (ưu tiên .vcf.gz) ==========
vcf_file <- file.path(
  "results", SAMPLE_ID, "snpeff",
  paste0(SAMPLE_ID, "_final_annotated.vcf.gz")
)

if (!file.exists(vcf_file)) {
  stop("❌ Không tìm thấy file .vcf.gz: ", vcf_file)
}

# ========== BƯỚC 4: Đọc file VCF ==========
vcf <- readVcf(vcf_file, "hg38")

# ========== BƯỚC 5: Lấy thông tin vị trí biến thể ==========
variant_range <- rowRanges(vcf)
variant_df <- data.frame(
  CHROM = as.character(seqnames(variant_range)),
  POS   = start(variant_range),
  REF   = as.character(variant_range@elementMetadata$REF),
  ALT   = sapply(variant_range@elementMetadata$ALT, function(x) paste(as.character(x), collapse = ","))
)

# ========== BƯỚC 6: Lấy ANN annotation ==========
ann_list <- info(vcf)$ANN
ann_df_list <- lapply(ann_list, function(x) {
  strsplit(x, ",")[[1]] %>%
    lapply(function(record) strsplit(record, "\\|")[[1]]) %>%
    do.call(what = rbind) %>%
    as.data.frame(stringsAsFactors = FALSE)
})

# Ghép tất cả annotation lại
ann_df <- bind_rows(ann_df_list)

# Gán tên cột (16 cột phổ biến nhất)
colnames(ann_df) <- c(
  "Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID",
  "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c",
  "HGVS.p", "cDNA_pos_len", "CDS_pos_len", "AA_pos_len", "Distance", "Info"
)

# ========== BƯỚC 7: Lấy thông tin từ ClinVar ==========
clinvar_info <- info(vcf)
CLNSIG <- as.character(unlist(clinvar_info$CLNSIG))
CLNDN  <- as.character(unlist(clinvar_info$CLNDN))
max_len <- max(length(CLNSIG), length(CLNDN), nrow(ann_df))

clinvar_df <- data.frame(
  CLNSIG = rep(CLNSIG, length.out = max_len),
  CLNDN  = rep(CLNDN, length.out = max_len),
  stringsAsFactors = FALSE
)

# ========== BƯỚC 8: Gộp tất cả vào bảng kết quả ==========
final_df <- cbind(variant_df, ann_df[1:nrow(variant_df), ], clinvar_df[1:nrow(variant_df), ])

# ========== BƯỚC 9: Ghi ra file ==========
output_file <- file.path("results", SAMPLE_ID, "snpeff", "vcf_annotated_output.csv")
write_csv(final_df, output_file)

cat("✅ Xuất file CSV thành công:", output_file, "\n")
