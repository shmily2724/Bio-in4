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
ANN_DIR="${SAMPLE_DIR}/snpeff"
SNPEFF_STATS="${ANN_DIR}/${SAMPLE_NAME}" 

# Tạo các thư mục output
mkdir -p "${HAPLO_DIR}"
mkdir -p "${ANN_DIR}"

# File input
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"

# File output
HAPLO_GVCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.g.vcf.gz"
HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"
ANN_VCF="${ANN_DIR}/${SAMPLE_NAME}_gatk_annotated.vcf"

# Các biến cho SnpEff
SNPEFF_JAR="/home/shmily/miniconda/envs/GATK/share/snpeff-5.2-1/snpEff.jar"
SNPEFF_CONFIG="/home/shmily/miniconda/envs/GATK/share/snpeff-5.2-1/snpEff.config"
SNPEFF_DB="GRCh38.86"

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
echo "Index hoàn tất."

# --- BƯỚC 4: CHÚ GIẢI VỚI SnpEff ---
echo "Bắt đầu chú giải với SnpEff..."

# Kiểm tra xem database của SnpEff đã tồn tại chưa
SNPEFF_DATA_DIR=$(dirname "$SNPEFF_JAR")/data
if [ ! -d "${SNPEFF_DATA_DIR}/${SNPEFF_DB}" ]; then
    echo "Thông báo: Database ${SNPEFF_DB} của SnpEff chưa tồn tại. Bắt đầu tải về..."
    java -jar "$SNPEFF_JAR" download "$SNPEFF_DB" -c "$SNPEFF_CONFIG"
fi

# SỬA LỖI: Xóa dấu \ ở cuối lệnh
java -Xmx6g -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \
    -stats "${SNPEFF_STATS}" \
    "${HAPLO_VCF}" > "${ANN_VCF}"

echo "Pipeline hoàn tất cho mẫu: ${SAMPLE_NAME}!"

