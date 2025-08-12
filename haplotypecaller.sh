#!/bin/bash

# Tự động dừng script khi có lỗi
set -e

# Kiểm tra tham số đầu vào
if [ -z "$1" ]; then
    echo "Lỗi: Vui lòng cung cấp tên mẫu (SAMPLE_NAME)."
    echo "Cách dùng: $0 <sample_name>"
    exit 1
fi

SAMPLE_NAME=$1

# --- TẠO CÁC BIẾN ĐƯỜNG DẪN ---
PROJECT_DIR="/media/shmily/writable/BRCA_project"
REF="/media/shmily/writable/BRCA_project/reference/Homo_sapiens_assembly38.fasta"
TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed" #link tải ở đây https://support.illumina.com/downloads/enrichment-bed-files-hg38.html
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
RECAL_DIR="${SAMPLE_DIR}/recal"
HAPLO_DIR="${SAMPLE_DIR}/haplotypecaller"

# Tạo các thư mục output
mkdir -p "${HAPLO_DIR}"

# File input
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"

# File output
HAPLO_GVCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.g.vcf.gz"
HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"

# Kích hoạt môi trường Conda (chỉ một lần)
eval "$(conda shell.bash hook)"
conda activate GATK

# --- BƯỚC 1: GỌI BIẾN THỂ VỚI HaplotypeCaller ---
echo "Bắt đầu HaplotypeCaller cho mẫu: ${SAMPLE_NAME}"
# Nếu chỉ phân tích 1 mẫu, có thể bỏ -ERC GVCF để ra thẳng VCF.
gatk HaplotypeCaller \
    -R "${REF}" \
    -I "${RECAL_BAM}" \
    -O "${HAPLO_GVCF}" \
    -L "${TARGET_BED}" \
    -ERC GVCF

echo "HaplotypeCaller hoàn tất."

# --- BƯỚC 2: CHUYỂN GVCF SANG VCF VỚI GenotypeGVCFs ---
echo "Bắt đầu GenotypeGVCFs cho mẫu: ${SAMPLE_NAME}"

gatk GenotypeGVCFs \
    -R "${REF}" \
    -V "${HAPLO_GVCF}" \
    -O "${HAPLO_VCF}" \
    -L "${TARGET_BED}"

echo "GenotypeGVCFs hoàn tất."

# --- BƯỚC 3: ĐÁNH INDEX CHO FILE VCF ---
echo "Đánh index cho file VCF..."
tabix -f -p vcf "${HAPLO_VCF}" # -f để overwrite tránh dừng script tại đây khi chạy lại
echo "Index hoàn tất cho mẫu: ${SAMPLE_NAME}!"

