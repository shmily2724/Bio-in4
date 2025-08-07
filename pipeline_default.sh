#!/bin/bash

# ===================================================================
# ==         PIPELINE PHÂN TÍCH DỮ LIỆU BRCA TOÀN DIỆN            ==
# ===================================================================
#
# Script này thực hiện toàn bộ quy trình từ dữ liệu FASTQ thô
# đến file VCF đã chú giải và báo cáo chất lượng tổng hợp.
#
# CÁCH DÙNG:
#   bash run_brca_pipeline.sh <sample_name>
#
# CÁC BƯỚC CHÍNH:
#   1. QC (FastQC) và Trimming (Trimmomatic)
#   2. Gióng hàng (BWA), Sắp xếp & Đánh dấu trùng lặp (Picard)
#   3. Hiệu chỉnh điểm chất lượng Base (GATK BQSR)
#   4. Gọi biến thể (GATK HaplotypeCaller) & Chú giải (SnpEff)
#   5. Tính toán độ phủ (Mosdepth)
#   6. Tổng hợp báo cáo (MultiQC)
#   7. Dọn dẹp các file trung gian (Tùy chọn)
#
# ===================================================================


# --- BƯỚC 0: KIỂM TRA ĐẦU VÀO VÀ THIẾT LẬP BAN ĐẦU ---

# Tự động dừng script khi có bất kỳ lệnh nào thất bại
set -e

# Kiểm tra tham số đầu vào (tên mẫu)
if [ -z "$1" ]; then
    echo "LỖI: Vui lòng cung cấp tên mẫu (SAMPLE_NAME)."
    echo "Cách dùng: $0 <sample_name>"
    exit 1
fi

echo "---=== BẮT ĐẦU PIPELINE CHO MẪU: $1 ===---"

# --- KHAI BÁO BIẾN TOÀN CỤC ---

# Cài đặt chung
SAMPLE_NAME=$1
THREADS=8 # Số luồng CPU sử dụng
CLEANUP=false # Đặt là 'false' nếu bạn muốn giữ lại tất cả các file trung gian

# Đường dẫn chính và file tham chiếu
PROJECT_DIR="/media/shmily/writable/BRCA_project"
REF="/media/shmily/writable/BRCA_project/reference/Homo_sapiens_assembly38.fasta"
ADAPTER_FILE="/home/shmily/miniconda/envs/BRCA/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"
TARGET_BED="/media/shmily/writable/BRCA_project/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Đường dẫn đến các công cụ
PICARD_JAR="/home/shmily/miniconda/envs/BRCA/share/picard-2.20.4-0/picard.jar"
SNPEFF_JAR="/home/shmily/miniconda/envs/GATK/share/snpeff-5.2-1/snpEff.jar"
SNPEFF_CONFIG="/home/shmily/miniconda/envs/GATK/share/snpeff-5.2-1/snpEff.config"
SNPEFF_DB="GRCh38.86"

# Các file VCF chứa các biến thể đã biết (known sites) cho BQSR
KNOWN_SNP="/media/shmily/writable/BRCA_project/reference/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_INDEL="/media/shmily/writable/BRCA_project/reference/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS_1000G_INDEL="/media/shmily/writable/BRCA_project/reference/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# --- KHAI BÁO CÁC THƯ MỤC VÀ FILE OUTPUT ---

# Thư mục chính
RAW_DIR="${PROJECT_DIR}/raw_data"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"

# Thư mục cho các bước
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"
FASTQC_RAW_DIR="${SAMPLE_DIR}/fastqc_raw"
FASTQC_TRIM_DIR="${SAMPLE_DIR}/fastqc_trimmed"
BWA_DIR="${SAMPLE_DIR}/Bwa_alignments"
RECAL_DIR="${SAMPLE_DIR}/recal"
HAPLO_DIR="${SAMPLE_DIR}/haplotypecaller"
ANN_DIR="${SAMPLE_DIR}/snpeff"
COVERAGE_DIR="${SAMPLE_DIR}/coverage"
MULTIQC_DIR="${SAMPLE_DIR}/multiqc_report"

# Tạo tất cả các thư mục output nếu chưa tồn tại
mkdir -p "${TRIM_DIR}" "${FASTQC_RAW_DIR}" "${FASTQC_TRIM_DIR}" "${BWA_DIR}" \
         "${RECAL_DIR}" "${HAPLO_DIR}" "${ANN_DIR}" "${COVERAGE_DIR}" "${MULTIQC_DIR}"

# File Input
READ1="${RAW_DIR}/${SAMPLE_NAME}_1.fastq.gz"
READ2="${RAW_DIR}/${SAMPLE_NAME}_2.fastq.gz"

# Các file trung gian và cuối cùng
# B1: Trimming
TRIMMED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"
UNPAIRED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_unpaired.fastq.gz"
UNPAIRED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_unpaired.fastq.gz"
# B2: Alignment
BAM="${BWA_DIR}/${SAMPLE_NAME}_aligned.bam"
SORTED_BAM="${BWA_DIR}/${SAMPLE_NAME}_sorted.bam"
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam"
BAM_DEDUP_BAI="${BWA_DIR}/${SAMPLE_NAME}_dedup.bai"
DEDUP_METRICS="${BWA_DIR}/${SAMPLE_NAME}_dedup_metrics.txt"
STATS_FILE="${BWA_DIR}/${SAMPLE_NAME}_samtools_stats.txt"
FLAGSTAT_FILE="${BWA_DIR}/${SAMPLE_NAME}_samtools_flagstat.txt"
# B3: BQSR
RECAL_DATA_TABLE="${RECAL_DIR}/${SAMPLE_NAME}_recal_data.table"
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"
# B4: Variant Calling
HAPLO_GVCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.g.vcf.gz"
HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"
ANN_VCF="${ANN_DIR}/${SAMPLE_NAME}_gatk_annotated.vcf"
SNPEFF_STATS="${ANN_DIR}/${SAMPLE_NAME}"
# B5: Coverage
COVERAGE_PREFIX="${COVERAGE_DIR}/${SAMPLE_NAME}"
# B6: MultiQC
MULTIQC_FILENAME="${SAMPLE_NAME}_report.html"

# --- KÍCH HOẠT CONDA ---
# Kích hoạt hook của shell một lần duy nhất
eval "$(conda shell.bash hook)"


# ===================================================================
# --- BƯỚC 1: QC (FASTQC) & TRIMMING (TRIMMOMATIC) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 1: Bắt đầu QC và Trimming ===---"
conda activate BRCA

echo "[1.1] Chạy FastQC trên dữ liệu thô..."
fastqc --threads ${THREADS} -o "${FASTQC_RAW_DIR}" "${READ1}" "${READ2}"
echo "FastQC trên dữ liệu thô hoàn tất."

echo "[1.2] Chạy Trimmomatic để loại bỏ adapter và read chất lượng thấp..."
trimmomatic PE -threads ${THREADS} -phred33 \
    "${READ1}" "${READ2}" \
    "${TRIMMED_R1}" "${UNPAIRED_R1}" \
    "${TRIMMED_R2}" "${UNPAIRED_R2}" \
    ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
echo "Trimmomatic hoàn tất."

echo "[1.3] Chạy FastQC trên dữ liệu đã trim..."
fastqc --threads ${THREADS} -o "${FASTQC_TRIM_DIR}" "${TRIMMED_R1}" "${TRIMMED_R2}"
echo "FastQC trên dữ liệu đã trim hoàn tất."


# ===================================================================
# --- BƯỚC 2: GIÓNG HÀNG (BWA), SẮP XẾP & ĐÁNH DẤU TRÙNG LẶP ---
# ===================================================================
echo ""
echo "---=== BƯỚC 2: Bắt đầu Gióng hàng, Sắp xếp và Đánh dấu trùng lặp ===---"
# Môi trường BRCA vẫn đang được kích hoạt

# Định nghĩa thông tin Read Group
RG_ID="${SAMPLE_NAME}_RG"
PLATFORM="Illumina"
LIBRARY_ID="Lib1"
PLATFORM_UNIT="${SAMPLE_NAME}_${PLATFORM}_${LIBRARY_ID}"

echo "[2.1] Gióng hàng với BWA-MEM..."
bwa mem -Y -K 100000000 -t ${THREADS} \
    -R "@RG\tID:${RG_ID}\tSM:${SAMPLE_NAME}\tPL:${PLATFORM}\tLB:${LIBRARY_ID}\tPU:${PLATFORM_UNIT}" \
    "$REF" \
    "${TRIMMED_R1}" \
    "${TRIMMED_R2}" | \
samtools view -Shb -o "${BAM}" -
echo "BWA-MEM hoàn tất."

echo "[2.2] Sắp xếp file BAM theo tọa độ với Picard SortSam..."
java -Djava.awt.headless=true -jar "${PICARD_JAR}" SortSam \
    I="${BAM}" O="${SORTED_BAM}" \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
echo "SortSam hoàn tất."

echo "[2.3] Đánh dấu các read trùng lặp với Picard MarkDuplicates..."
java -Djava.awt.headless=true -jar "${PICARD_JAR}" MarkDuplicates \
    I="${SORTED_BAM}" O="${BAM_DEDUP}" M="${DEDUP_METRICS}" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
echo "MarkDuplicates hoàn tất."

echo "[2.4] Tạo file thống kê cho MultiQC..."
samtools stats "${BAM_DEDUP}" > "${STATS_FILE}"
samtools flagstat "${BAM_DEDUP}" > "${FLAGSTAT_FILE}"
echo "Tạo file thống kê hoàn tất."


# ===================================================================
# --- BƯỚC 3: HIỆU CHỈNH ĐIỂM CHẤT LƯỢNG BASE (GATK BQSR) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 3: Bắt đầu Hiệu chỉnh điểm chất lượng Base (BQSR) ===---"
conda activate GATK

echo "[3.1] Tạo bảng hiệu chỉnh với BaseRecalibrator..."
gatk BaseRecalibrator \
    -I "${BAM_DEDUP}" \
    -R "${REF}" \
    --known-sites "${KNOWN_SNP}" \
    --known-sites "${KNOWN_INDEL}" \
    --known-sites "${MILLS_1000G_INDEL}" \
    -O "${RECAL_DATA_TABLE}"
echo "BaseRecalibrator hoàn tất."

echo "[3.2] Áp dụng BQSR để tạo file BAM mới..."
gatk ApplyBQSR \
    -R "${REF}" \
    -I "${BAM_DEDUP}" \
    --bqsr-recal-file "${RECAL_DATA_TABLE}" \
    -O "${RECAL_BAM}"
echo "ApplyBQSR hoàn tất."

echo "[3.3] Đánh index cho file BAM cuối cùng..."
samtools index "${RECAL_BAM}"
echo "Đánh index hoàn tất."


# ===================================================================
# --- BƯỚC 4: GỌI BIẾN THỂ (GATK) & CHÚ GIẢI (SNPEFF) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 4: Bắt đầu Gọi và Chú giải biến thể ===---"
# Môi trường GATK vẫn đang được kích hoạt

echo "[4.1] Gọi biến thể với HaplotypeCaller..."
gatk HaplotypeCaller \
    -R "${REF}" \
    -I "${RECAL_BAM}" \
    -O "${HAPLO_GVCF}" \
    -L "${TARGET_BED}" \
    -ERC GVCF
echo "HaplotypeCaller hoàn tất."

echo "[4.2] Chuyển gVCF sang VCF với GenotypeGVCFs..."
gatk GenotypeGVCFs \
    -R "${REF}" \
    -V "${HAPLO_GVCF}" \
    -O "${HAPLO_VCF}" \
    -L "${TARGET_BED}"
echo "GenotypeGVCFs hoàn tất."

echo "[4.3] Đánh index cho file VCF..."
tabix -f -p vcf "${HAPLO_VCF}"
echo "Index VCF hoàn tất."

echo "[4.4] Chú giải biến thể với SnpEff..."
# Kiểm tra và tải database của SnpEff nếu cần
SNPEFF_DATA_DIR=$(dirname "$SNPEFF_JAR")/data
if [ ! -d "${SNPEFF_DATA_DIR}/${SNPEFF_DB}" ]; then
    echo "Thông báo: Database ${SNPEFF_DB} của SnpEff chưa tồn tại. Bắt đầu tải về..."
    java -jar "$SNPEFF_JAR" download "$SNPEFF_DB" -c "$SNPEFF_CONFIG"
fi

# Chạy SnpEff
java -Xmx6g -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \
    -stats "${SNPEFF_STATS}.html" \
    "${HAPLO_VCF}" > "${ANN_VCF}"
echo "SnpEff hoàn tất."


# ===================================================================
# --- BƯỚC 5: TÍNH TOÁN ĐỘ PHỦ (MOSDEPTH) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 5: Bắt đầu Tính toán độ phủ ===---"
conda activate BRCA

echo "[5.1] Chạy Mosdepth trên vùng mục tiêu..."
mosdepth --threads ${THREADS} -n --by ${TARGET_BED} "${COVERAGE_PREFIX}" "${RECAL_BAM}"
echo "Mosdepth hoàn tất. File tóm tắt: ${COVERAGE_PREFIX}.mosdepth.summary.txt"


# ===================================================================
# --- BƯỚC 6: TỔNG HỢP BÁO CÁO (MULTIQC) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 6: Bắt đầu Tổng hợp báo cáo với MultiQC ===---"
conda activate MQC

echo "[6.1] Chạy MultiQC trên toàn bộ kết quả của mẫu..."
multiqc "${SAMPLE_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --title "Báo cáo QC cho mẫu ${SAMPLE_NAME}" \
    --filename "${MULTIQC_FILENAME}" \
    --force
echo "MultiQC hoàn tất."

# ===================================================================
# --- BƯỚC 7: DỌN DẸP FILE TRUNG GIAN ---
# ===================================================================

if [ "$CLEANUP" = true ]; then
    echo ""
    echo "---=== BƯỚC 7: Bắt đầu dọn dẹp các file trung gian ===---"
    
    rm -f \
        "${TRIMMED_R1}" \
        "${UNPAIRED_R1}" \
        "${TRIMMED_R2}" \
        "${UNPAIRED_R2}" \
        "${BAM}" \
        "${SORTED_BAM}" \
        "${BAM_DEDUP}" \
        "${BAM_DEDUP_BAI}" \
        "${RECAL_DATA_TABLE}" \
        "${HAPLO_GVCF}" \
        "${HAPLO_GVCF}.tbi" \
        "${HAPLO_VCF}" \
        "${HAPLO_VCF}.tbi"
        
    echo "Dọn dẹp hoàn tất."
else
    echo ""
    echo "---=== BƯỚC 7: Bỏ qua bước dọn dẹp file trung gian (CLEANUP=false) ===---"
fi


# --- HOÀN TẤT ---
echo ""
echo "==================================================================="
echo "==      TOÀN BỘ PIPELINE ĐÃ HOÀN TẤT THÀNH CÔNG CHO MẪU: ${SAMPLE_NAME}     =="
echo "== Kết quả cuối cùng và các báo cáo đã được lưu trong thư mục:"
echo "== ${SAMPLE_DIR}"
echo "== Báo cáo MultiQC tổng hợp: ${MULTIQC_DIR}/${MULTIQC_FILENAME}"
echo "==================================================================="

