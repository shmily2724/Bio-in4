#!/bin/bash

# Tự động dừng script khi có lỗi
set -e

# CẢI TIẾN: Kiểm tra tham số đầu vào
if [ -z "$1" ]; then
    echo "Lỗi: Vui lòng cung cấp tên mẫu (SAMPLE_NAME)."
    echo "Cách dùng: $0 <sample_name>"
    exit 1
fi

SAMPLE_NAME=$1

# --- TẠO CÁC BIẾN ĐƯỜNG DẪN ---
PROJECT_DIR="/media/shmily/writable/BRCA_project"
REF="/media/shmily/writable/BRCA_project/reference/Homo_sapiens_assembly38.fasta"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
BWA_DIR="${SAMPLE_DIR}/Bwa_alignments"
RECAL_DIR="${SAMPLE_DIR}/recal"

# Các file VCF chứa các biến thể đã biết (known sites)
KNOWN_SNP="/media/shmily/writable/BRCA_project/reference/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_INDEL="/media/shmily/writable/BRCA_project/reference/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS_1000G_INDEL="/media/shmily/writable/BRCA_project/reference/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
# Tạo thư mục output nếu chưa tồn tại
mkdir -p "${RECAL_DIR}"

# File input
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam"
BAM_BAI="${BWA_DIR}/${SAMPLE_NAME}_dedup.bai"
# File output
RECAL_DATA_TABLE="${RECAL_DIR}/${SAMPLE_NAME}_recal_data.table"
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"
RECAL_BAI="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bai"
#  Kiểm tra sự tồn tại của các file input quan trọng
for f in "${BAM_DEDUP}" "${BAM_BAI}" "${REF}" "${KNOWN_SNP}" "${KNOWN_INDEL}" "${MILLS_1000G_INDEL}"; do
    if [ ! -f "$f" ]; then
        echo "Lỗi: Không tìm thấy file input hoặc index cần thiết: $f"
        exit 1
    fi
done

# Kích hoạt môi trường Conda
eval "$(conda shell.bash hook)"
conda activate GATK

# --- BƯỚC 1: TẠO BẢNG RECALIBRATION VỚI BaseRecalibrator ---
echo "Bắt đầu BaseRecalibrator cho mẫu: ${SAMPLE_NAME}"


gatk BaseRecalibrator \
    -I "${BAM_DEDUP}" \
    -R "${REF}" \
    --known-sites "${KNOWN_SNP}" \
    --known-sites "${KNOWN_INDEL}" \
    --known-sites "${MILLS_1000G_INDEL}" \
    -O "${RECAL_DATA_TABLE}"

echo "BaseRecalibrator hoàn tất."

# --- BƯỚC 2: ÁP DỤNG BQSR VỚI ApplyBQSR ---
echo "Bắt đầu ApplyBQSR cho mẫu: ${SAMPLE_NAME}"

gatk ApplyBQSR \
    -R "${REF}" \
    -I "${BAM_DEDUP}" \
    --bqsr-recal-file "${RECAL_DATA_TABLE}" \
    -O "${RECAL_BAM}"

echo "ApplyBQSR hoàn tất."

# --- BƯỚC 3: ĐÁNH INDEX CHO FILE BAM CUỐI CÙNG ---
echo "Bắt đầu đánh index cho file BAM đã recalibrate..."

samtools index "${RECAL_BAM}"

echo "Pipeline BQSR hoàn tất cho mẫu: ${SAMPLE_NAME}!"
