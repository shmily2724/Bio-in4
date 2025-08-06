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

# Thư mục output cho báo cáo MultiQC
MULTIQC_DIR="${SAMPLE_DIR}/multiqc_report"

# Tạo thư mục output nếu chưa tồn tại
mkdir -p "${MULTIQC_DIR}"

# --- CHẠY LỆNH MultiQC ---
# Không cần kích hoạt Conda vì multiqc đã được cài toàn hệ thống
echo "Bắt đầu chạy MultiQC trên toàn bộ kết quả của mẫu: ${SAMPLE_NAME}"

multiqc "${SAMPLE_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --title "Báo cáo QC cho mẫu ${SAMPLE_NAME}" \
    --force

echo "MultiQC hoàn tất!"
echo "Báo cáo đã được lưu tại: ${MULTIQC_DIR}/multiqc_report.html"
