#!/bin/bash
# ===================================================================
# ==         PIPELINE PHÂN TÍCH DỮ LIỆU BRCA TOÀN DIỆN            ==
# ===================================================================

set -euo pipefail

# ========== CẤU HÌNH CHUNG ==========
SAMPLE_NAME="${1:?Vui lòng truyền SAMPLE_NAME}"
THREADS="${THREADS:-8}"            # override bằng biến môi trường nếu muốn
CLEANUP="${CLEANUP:-true}"         # true|false – xóa file trung gian
TIMESTAMP() { date '+%Y-%m-%d %H:%M:%S'; }

# Conda envs (đổi theo máy bạn)
ENV_BRCA="BRCA"
ENV_GATK="GATK"
ENV_MQC="MQC"

# Logging (tuỳ chọn)
PROJECT_DIR="/media/shmily/writable/BRCA_project"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
LOG="${SAMPLE_DIR}/${SAMPLE_NAME}_pipeline.log"
mkdir -p "${SAMPLE_DIR}"
# exec > >(tee -i "$LOG") 2>&1    # bật nếu muốn ghi log toàn bộ

# ========== THAM CHIẾU & VÙNG MỤC TIÊU ==========
REF="${PROJECT_DIR}/reference/Homo_sapiens_assembly38.fasta"
TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Known sites cho GATK BQSR
KNOWN_SNP="${PROJECT_DIR}/reference/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_INDEL="${PROJECT_DIR}/reference/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS_1000G_INDEL="${PROJECT_DIR}/reference/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# ========== TOOL PATHS ==========
# Binaries (đảm bảo đã có trong PATH hoặc chỉnh tuyệt đối)
FASTQC_BIN="fastqc"
TRIMMOMATIC_BIN="trimmomatic"
BWA_BIN="bwa"
SAMTOOLS_BIN="samtools"
GATK_BIN="gatk"
MOSDEPTH_BIN="mosdepth"
MULTIQC_BIN="multiqc"

# JAR / cấu hình
PICARD_JAR="/home/shmily/miniconda/envs/BRCA/share/picard-2.20.4-0/picard.jar"

SNPEFF_HOME="${PROJECT_DIR}/snpEff"
SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_HOME}/SnpSift.jar"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
SNPEFF_DB="GRCh38.86"
SNPEFF_DATA_DIR="${SNPEFF_HOME}/data"

# Java opts (tuỳ chọn tinh chỉnh GC)
JAVA_OPTS_SNPEFF="-Xmx6g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_SNPSIFT="-Xmx4g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_PICARD="-Xmx4g -Djava.awt.headless=true -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"

# ========== EXTERNAL RESOURCES (ANNOTATION) ==========
GNOMAD_VCF="${PROJECT_DIR}/reference/resources/gnomad.v4.1.panel.merged.vcf.gz"
CLINVAR_VCF="${PROJECT_DIR}/reference/resources/clinvar_20250810.vcf.gz"
THOUSANDG_VCF="${PROJECT_DIR}/reference/resources/1000g.panel.merged.vcf.gz"

# ========== DỮ LIỆU ĐẦU VÀO ==========
RAW_DIR="${PROJECT_DIR}/raw_data"
READ1="${RAW_DIR}/${SAMPLE_NAME}_1.fastq.gz"
READ2="${RAW_DIR}/${SAMPLE_NAME}_2.fastq.gz"
ADAPTER_FILE="/home/shmily/miniconda/envs/BRCA/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa"

# ========== THƯ MỤC KẾT QUẢ THEO BƯỚC ==========
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"
FASTQC_RAW_DIR="${SAMPLE_DIR}/fastqc_raw"
FASTQC_TRIM_DIR="${SAMPLE_DIR}/fastqc_trimmed"
BWA_DIR="${SAMPLE_DIR}/Bwa_alignments"
RECAL_DIR="${SAMPLE_DIR}/recal"
HAPLO_DIR="${SAMPLE_DIR}/haplotypecaller"
ANN_DIR="${SAMPLE_DIR}/snpeff"
COVERAGE_DIR="${SAMPLE_DIR}/coverage"
MULTIQC_DIR="${SAMPLE_DIR}/multiqc_report"
mkdir -p "${TRIM_DIR}" "${FASTQC_RAW_DIR}" "${FASTQC_TRIM_DIR}" "${BWA_DIR}" \
         "${RECAL_DIR}" "${HAPLO_DIR}" "${ANN_DIR}" "${COVERAGE_DIR}" "${MULTIQC_DIR}"

# ========== FILE TRUNG GIAN & ĐẦU RA ==========
# Trimming
TRIMMED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"
UNPAIRED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_unpaired.fastq.gz"
UNPAIRED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_unpaired.fastq.gz"

# Alignment & BAM
BAM="${BWA_DIR}/${SAMPLE_NAME}_aligned.bam"
SORTED_BAM="${BWA_DIR}/${SAMPLE_NAME}_sorted.bam"
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam"
BAM_DEDUP_BAI="${BWA_DIR}/${SAMPLE_NAME}_dedup.bai"
DEDUP_METRICS="${BWA_DIR}/${SAMPLE_NAME}_dedup_metrics.txt"
STATS_FILE="${BWA_DIR}/${SAMPLE_NAME}_samtools_stats.txt"
FLAGSTAT_FILE="${BWA_DIR}/${SAMPLE_NAME}_samtools_flagstat.txt"

# BQSR
RECAL_DATA_TABLE="${RECAL_DIR}/${SAMPLE_NAME}_recal_data.table"
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"

# GATK calling
HAPLO_GVCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.g.vcf.gz"
HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"

# Annotation outputs
SNPEFF_STATS="${ANN_DIR}/${SAMPLE_NAME}"
ANN_VCF_FINAL="${ANN_DIR}/${SAMPLE_NAME}_final_annotated.vcf"
ANN_VCF_SNPEFF="${ANN_DIR}/${SAMPLE_NAME}_snpeff.vcf"
ANN_VCF_GNOMAD="${ANN_DIR}/${SAMPLE_NAME}_gnomad.vcf"
ANN_VCF_CLINVAR="${ANN_DIR}/${SAMPLE_NAME}_clinvar.vcf"

# Map tên file tạm dùng ở bước 4
TMP1="${ANN_VCF_SNPEFF}"
TMP2="${ANN_VCF_GNOMAD}"
TMP3="${ANN_VCF_CLINVAR}"
ANN_VCF="${ANN_VCF_FINAL}"

# Coverage & Reporting
COVERAGE_PREFIX="${COVERAGE_DIR}/${SAMPLE_NAME}"
MULTIQC_FILENAME="${SAMPLE_NAME}_report.html"

# ========== READ GROUP (BWA-MEM) ==========
RG_ID="${SAMPLE_NAME}_RG"
PLATFORM="Illumina"
LIBRARY_ID="Lib1"
PLATFORM_UNIT="${SAMPLE_NAME}_${PLATFORM}_${LIBRARY_ID}"
RG_STR="@RG\tID:${RG_ID}\tSM:${SAMPLE_NAME}\tPL:${PLATFORM}\tLB:${LIBRARY_ID}\tPU:${PLATFORM_UNIT}"

# --- KÍCH HOẠT CONDA ---
eval "$(conda shell.bash hook)"

# --- Tiền kiểm tra input quan trọng ---
for f in "${READ1}" "${READ2}" "${ADAPTER_FILE}" "${REF}" "${TARGET_BED}"; do
  [ -f "$f" ] || { echo "❌ Missing: $f"; exit 1; }
done

# ===================================================================
# --- BƯỚC 1: QC (FASTQC) & TRIMMING (TRIMMOMATIC) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 1: Bắt đầu QC và Trimming ===---"
conda activate "${ENV_BRCA}"

echo "[1.1] Chạy FastQC trên dữ liệu thô..."
${FASTQC_BIN} --threads ${THREADS} -o "${FASTQC_RAW_DIR}" "${READ1}" "${READ2}"
echo "FastQC trên dữ liệu thô hoàn tất."

echo "[1.2] Chạy Trimmomatic để loại bỏ adapter và read chất lượng thấp..."
${TRIMMOMATIC_BIN} PE -threads ${THREADS} -phred33 \
    "${READ1}" "${READ2}" \
    "${TRIMMED_R1}" "${UNPAIRED_R1}" \
    "${TRIMMED_R2}" "${UNPAIRED_R2}" \
    ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
echo "Trimmomatic hoàn tất."

echo "[1.3] Chạy FastQC trên dữ liệu đã trim..."
${FASTQC_BIN} --threads ${THREADS} -o "${FASTQC_TRIM_DIR}" "${TRIMMED_R1}" "${TRIMMED_R2}"
echo "FastQC trên dữ liệu đã trim hoàn tất."

# ===================================================================
# --- BƯỚC 2: GIÓNG HÀNG (BWA), SẮP XẾP & ĐÁNH DẤU TRÙNG LẶP ---
# ===================================================================
echo ""
echo "---=== BƯỚC 2: Bắt đầu Gióng hàng, Sắp xếp và Đánh dấu trùng lặp ===---"

echo "[2.1] Gióng hàng với BWA-MEM..."
${BWA_BIN} mem -Y -K 100000000 -t ${THREADS} \
    -R "${RG_STR}" \
    "$REF" \
    "${TRIMMED_R1}" \
    "${TRIMMED_R2}" | \
${SAMTOOLS_BIN} view -b -o "${BAM}" -
echo "BWA-MEM hoàn tất."

echo "[2.2] Sắp xếp file BAM theo tọa độ với Picard SortSam..."
java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" SortSam \
    I="${BAM}" O="${SORTED_BAM}" \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
echo "SortSam hoàn tất."

echo "[2.3] Đánh dấu các read trùng lặp với Picard MarkDuplicates..."
java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" MarkDuplicates \
    I="${SORTED_BAM}" O="${BAM_DEDUP}" M="${DEDUP_METRICS}" \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
echo "MarkDuplicates hoàn tất."

echo "[2.4] Tạo file thống kê cho MultiQC..."
${SAMTOOLS_BIN} stats "${BAM_DEDUP}" > "${STATS_FILE}"
${SAMTOOLS_BIN} flagstat "${BAM_DEDUP}" > "${FLAGSTAT_FILE}"
echo "Tạo file thống kê hoàn tất."

# ===================================================================
# --- BƯỚC 3: HIỆU CHỈNH ĐIỂM CHẤT LƯỢNG BASE (GATK BQSR) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 3: Bắt đầu Hiệu chỉnh điểm chất lượng Base (BQSR) ===---"
conda activate "${ENV_GATK}"

echo "[3.1] Tạo bảng hiệu chỉnh với BaseRecalibrator..."
${GATK_BIN} BaseRecalibrator \
    -I "${BAM_DEDUP}" \
    -R "${REF}" \
    --known-sites "${KNOWN_SNP}" \
    --known-sites "${KNOWN_INDEL}" \
    --known-sites "${MILLS_1000G_INDEL}" \
    -O "${RECAL_DATA_TABLE}"
echo "BaseRecalibrator hoàn tất."

echo "[3.2] Áp dụng BQSR để tạo file BAM mới..."
${GATK_BIN} ApplyBQSR \
    -R "${REF}" \
    -I "${BAM_DEDUP}" \
    --bqsr-recal-file "${RECAL_DATA_TABLE}" \
    -O "${RECAL_BAM}"
echo "ApplyBQSR hoàn tất."

echo "[3.3] Đánh index cho file BAM cuối cùng..."
${SAMTOOLS_BIN} index "${RECAL_BAM}"
echo "Đánh index hoàn tất."

# ===================================================================
# --- BƯỚC 4: GỌI BIẾN THỂ (GATK) & CHÚ GIẢI (SNPEFF + SNPSIFT) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 4: Bắt đầu Gọi và Chú giải biến thể ===---"

echo "[4.1] Gọi biến thể với HaplotypeCaller..."
${GATK_BIN} HaplotypeCaller \
    -R "${REF}" \
    -I "${RECAL_BAM}" \
    -O "${HAPLO_GVCF}" \
    -L "${TARGET_BED}" \
    -ERC GVCF
echo "HaplotypeCaller hoàn tất."

echo "[4.1b] Index gVCF..."
tabix -f -p vcf "${HAPLO_GVCF}"

echo "[4.2] Chuyển gVCF sang VCF với GenotypeGVCFs..."
${GATK_BIN} GenotypeGVCFs \
    -R "${REF}" \
    -V "${HAPLO_GVCF}" \
    -O "${HAPLO_VCF}" \
    -L "${TARGET_BED}"
echo "GenotypeGVCFs hoàn tất."

echo "[4.3] Đánh index cho file VCF..."
tabix -f -p vcf "${HAPLO_VCF}"
echo "Index VCF hoàn tất."

echo "[4.4] Chú giải biến thể với SnpEff + SnpSift..."

# Đảm bảo DB SnpEff sẵn sàng
if [ ! -d "${SNPEFF_DATA_DIR}/${SNPEFF_DB}" ]; then
    echo "⚠️ $(TIMESTAMP) Tải database SnpEff ${SNPEFF_DB}..."
    java ${JAVA_OPTS_SNPEFF} -jar "$SNPEFF_JAR" download "$SNPEFF_DB" -c "$SNPEFF_CONFIG"
fi

# [1] SnpEff
echo "🔬 [1] Annotating với SnpEff..."
java ${JAVA_OPTS_SNPEFF} -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \
    -stats "${SNPEFF_STATS}" \
    "${HAPLO_VCF}" > "${TMP1}"

# [2] gnomAD (thêm AF) → ĐỔI TÊN AF → GNOMAD_AF (nếu có bcftools)
echo "📊 [2] Annotating với gnomAD..."
[ -f "${GNOMAD_VCF}" ] || { echo "❌ Không tìm thấy gnomAD VCF: ${GNOMAD_VCF}"; exit 1; }
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -info AF "${GNOMAD_VCF}" "${TMP1}" > "${TMP2}"
#  (KHÔNG dùng -id: sẽ ghép theo CHROM:POS:REF:ALT, an toàn & đầy đủ hơn)

if command -v bcftools >/dev/null 2>&1; then
  echo "📝 Đổi INFO/AF → GNOMAD_AF..."
  echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="Allele frequency from gnomAD">' > "${ANN_DIR}/gn_hdr.hdr"
  bcftools annotate -h "${ANN_DIR}/gn_hdr.hdr" \
    -c INFO/GNOMAD_AF:=INFO/AF -x INFO/AF \
    -O v -o "${TMP2}.renamed" "${TMP2}"
  mv -f "${TMP2}.renamed" "${TMP2}"
  rm -f "${ANN_DIR}/gn_hdr.hdr"
else
  echo "ℹ️ Không có bcftools → GIỮ nguyên INFO/AF sau gnomAD."
fi

# [3] ClinVar
echo "🧬 [3] Annotating với ClinVar..."
[ -f "${CLINVAR_VCF}" ] || { echo "❌ Không tìm thấy ClinVar VCF: ${CLINVAR_VCF}"; exit 1; }
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate \
    -info CLNSIG,CLNDN "${CLINVAR_VCF}" "${TMP2}" > "${TMP3}"

# [4] 1000 Genomes (thêm AF), rồi ĐỔI TÊN AF -> KG_AF (nếu có bcftools)
echo "🌍 [4] Annotating với 1000 Genomes..."
[ -f "${THOUSANDG_VCF}" ] || { echo "❌ Không tìm thấy 1000 Genomes VCF: ${THOUSANDG_VCF}"; exit 1; }
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate \
    -info AF "${THOUSANDG_VCF}" "${TMP3}" > "${ANN_VCF}"

if command -v bcftools >/dev/null 2>&1; then
  echo "📝 Đổi INFO/AF → KG_AF..."
  echo '##INFO=<ID=KG_AF,Number=A,Type=Float,Description="Allele frequency from 1000 Genomes Project">' > "${ANN_DIR}/kg_hdr.hdr"
  bcftools annotate -h "${ANN_DIR}/kg_hdr.hdr" \
    -c INFO/KG_AF:=INFO/AF -x INFO/AF \
    -O v -o "${ANN_VCF}.tmp" "${ANN_VCF}"
  mv -f "${ANN_VCF}.tmp" "${ANN_VCF}"
  rm -f "${ANN_DIR}/kg_hdr.hdr"
else
  echo "ℹ️ Không có bcftools → GIỮ nguyên INFO/AF sau 1000G."
fi

echo "✅ Hoàn tất tạo annotated VCF (đã tách GNOMAD_AF & KG_AF nếu có bcftools): ${ANN_VCF}"

# (Tuỳ chọn) Nén + index VCF cuối để truy vấn nhanh
if command -v bgzip >/dev/null 2>&1; then
  bgzip -f "${ANN_VCF}"
  tabix -f -p vcf "${ANN_VCF}.gz"
  echo "📦 Đã nén & index: ${ANN_VCF}.gz"
fi

# ===================================================================
# --- BƯỚC 5: TÍNH TOÁN ĐỘ PHỦ (MOSDEPTH) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 5: Bắt đầu Tính toán độ phủ ===---"
conda activate "${ENV_BRCA}"

echo "[5.1] Chạy Mosdepth trên vùng mục tiêu..."
${MOSDEPTH_BIN} --threads ${THREADS} -n --by ${TARGET_BED} "${COVERAGE_PREFIX}" "${RECAL_BAM}"
echo "Mosdepth hoàn tất. File tóm tắt: ${COVERAGE_PREFIX}.mosdepth.summary.txt"

# ===================================================================
# --- BƯỚC 6: TỔNG HỢP BÁO CÁO (MULTIQC) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 6: Bắt đầu Tổng hợp báo cáo với MultiQC ===---"
conda activate "${ENV_MQC}"

echo "[6.1] Chạy MultiQC trên toàn bộ kết quả của mẫu..."
${MULTIQC_BIN} "${SAMPLE_DIR}" \
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
        "${ANN_DIR}/gn_hdr.hdr" "${ANN_DIR}/kg_hdr.hdr" \
        "${TMP2}.renamed" "${ANN_VCF}.tmp" \
        "${TMP1}" "${TMP2}" "${TMP3}"
    echo "Dọn dẹp hoàn tất."
else
    echo ""
    echo "---=== BƯỚC 7: Bỏ qua bước dọn dẹp file trung gian (CLEANUP=false) ===---"
fi

# --- HOÀN TẤT ---
echo ""
echo "==================================================================="
echo "==  PIPELINE ĐÃ HOÀN TẤT THÀNH CÔNG CHO MẪU: ${SAMPLE_NAME}"
echo "==  Kết quả & báo cáo tại: ${SAMPLE_DIR}"
echo "==  Báo cáo MultiQC: ${MULTIQC_DIR}/${MULTIQC_FILENAME}"
echo "==================================================================="
