#!/bin/bash

# Tự động dừng script khi có lỗi
set -e

# Kiểm tra xem tham số đầu vào đã được cung cấp chưa
if [ -z "$1" ]; then
    echo "Lỗi: Vui lòng cung cấp tên mẫu (SAMPLE_NAME)."
    echo "Cách dùng: $0 <sample_name>"
    exit 1
fi

SAMPLE_NAME=$1

# --- TẠO CÁC BIẾN ĐƯỜNG DẪN ---
PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
RECAL_DIR="${SAMPLE_DIR}/recal"
COVERAGE_DIR="${SAMPLE_DIR}/coverage"
TARGET_BED="/media/shmily/writable/BRCA_project/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Tạo thư mục output nếu chưa tồn tại
mkdir -p "${COVERAGE_DIR}"

# --- Định nghĩa file Input và Output ---
# File input là file BAM đã qua BQSR
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"

# Tên tiền tố cho các file output của mosdepth
COVERAGE_PREFIX="${COVERAGE_DIR}/${SAMPLE_NAME}"

# CẢI TIẾN: Kiểm tra xem file BAM đầu vào có tồn tại không
if [ ! -f "${RECAL_BAM}" ]; then
    echo "Lỗi: Không tìm thấy file input ${RECAL_BAM}."
    echo "Hãy chắc chắn rằng bạn đã chạy script BQSR thành công trước."
    exit 1
fi

# --- KÍCH HOẠT MÔI TRƯỜNG CONDA ---
# Giả sử mosdepth được cài trong môi trường GATK
eval "$(conda shell.bash hook)"
conda activate BRCA

# --- CHẠY LỆNH MOSDEPTH ---
echo "Bắt đầu tính toán độ phủ với Mosdepth cho mẫu: ${SAMPLE_NAME}"

mosdepth --threads 8 \
    -n \
    --by ${TARGET_BED} \
    "${COVERAGE_PREFIX}" \
    "${RECAL_BAM}"

echo "Mosdepth hoàn tất cho mẫu ${SAMPLE_NAME}!"
echo "File tóm tắt kết quả: ${COVERAGE_PREFIX}.mosdepth.summary.txt"
