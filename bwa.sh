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
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"
BWA_DIR="${SAMPLE_DIR}/Bwa_alignments"

# Tạo thư mục output nếu chưa tồn tại
mkdir -p "$BWA_DIR"

# File input
TRIMMED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"

# File output
BAM="${BWA_DIR}/${SAMPLE_NAME}_aligned.bam"
SORTED_BAM="${BWA_DIR}/${SAMPLE_NAME}_sorted.bam"
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam"
DEDUP_METRICS="${BWA_DIR}/${SAMPLE_NAME}_dedup_metrics.txt"

# Đường dẫn đến Picard Jar
PICARD_JAR="/home/shmily/miniconda/envs/BRCA/share/picard-2.20.4-0/picard.jar"

# Định nghĩa thông tin Read Group
RG_ID="${SAMPLE_NAME}_RG"
PLATFORM="Illumina"
LIBRARY_ID="Lib1"
PLATFORM_UNIT="${SAMPLE_NAME}_${PLATFORM}_${LIBRARY_ID}"

# Kích hoạt môi trường Conda (chỉ một lần)
eval "$(conda shell.bash hook)"
conda activate BRCA

# --- BƯỚC 1: ALIGNMENT VỚI BWA-MEM ---
echo "Bắt đầu BWA-MEM cho mẫu: ${SAMPLE_NAME}"

# CẢI TIẾN: Kiểm tra sự tồn tại của BWA index trước khi chạy
if [ ! -f "${REF}.64.bwt" ]; then
    echo "Lỗi: Không tìm thấy BWA index cho file tham chiếu. Hãy chạy 'bwa index ${REF}' trước."
    exit 1
fi

# KHỐI LỆNH ĐÚNG (đã dùng khoảng trắng thông thường)
bwa mem -Y \
    -K 100000000 \
    -t 16 \
    -R "@RG\tID:${RG_ID}\tSM:${SAMPLE_NAME}\tPL:${PLATFORM}\tLB:${LIBRARY_ID}\tPU:${PLATFORM_UNIT}" \
    "$REF" \
    "${TRIMMED_R1}" \
    "${TRIMMED_R2}" | \
samtools view -Shb -o "${BAM}" -

echo "BWA-MEM hoàn tất."

# --- BƯỚC 2: SORT BAM VỚI PICARD ---
echo "Bắt đầu SortSam cho mẫu: ${SAMPLE_NAME}"

# KHỐI LỆNH ĐÚNG (đã dùng khoảng trắng thông thường)
java -Djava.awt.headless=true -jar "${PICARD_JAR}" SortSam \
    I="${BAM}" \
    O="${SORTED_BAM}" \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=2000000

echo "SortSam hoàn tất."

# --- BƯỚC 3: MARK DUPLICATES VỚI PICARD ---
echo "Bắt đầu MarkDuplicates cho mẫu: ${SAMPLE_NAME}"

# KHỐI LỆNH ĐÚNG (đã dùng khoảng trắng thông thường)
java -Djava.awt.headless=true -jar "${PICARD_JAR}" MarkDuplicates \
    I="${SORTED_BAM}" \
    O="${BAM_DEDUP}" \
    M="${DEDUP_METRICS}" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=2000000

echo "Pipeline hoàn tất cho mẫu: ${SAMPLE_NAME}!"
    
