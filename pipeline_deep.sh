#!/bin/bash
# ===================================================================
# ==         PIPELINE PHÂN TÍCH DỮ LIỆU BRCA (DeepVariant)        ==
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
DEEP_DIR="${SAMPLE_DIR}/Deepvariants"
ANN_DIR="${SAMPLE_DIR}/snpeff"
COVERAGE_DIR="${SAMPLE_DIR}/coverage"
MULTIQC_DIR="${SAMPLE_DIR}/multiqc_report"
mkdir -p "${TRIM_DIR}" "${FASTQC_RAW_DIR}" "${FASTQC_TRIM_DIR}" "${BWA_DIR}" \
         "${RECAL_DIR}" "${DEEP_DIR}" "${ANN_DIR}" "${COVERAGE_DIR}" "${MULTIQC_DIR}"

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

# DeepVariant outputs (thay cho GATK calling)
DV_VCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.vcf.gz"
DV_GVCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.g.vcf.gz"
DV_DOCKER_IMAGE="${DV_DOCKER_IMAGE:-google/deepvariant:1.9.0}"

# Annotation outputs (đặt nhãn _dv_ để không trùng với GATK)
SNPEFF_STATS="${ANN_DIR}/${SAMPLE_NAME}_dv"
ANN_VCF_FINAL="${ANN_DIR}/${SAMPLE_NAME}_dv_final_annotated.vcf"
ANN_VCF_SNPEFF="${ANN_DIR}/${SAMPLE_NAME}_dv_snpeff.vcf"
ANN_VCF_GNOMAD="${ANN_DIR}/${SAMPLE_NAME}_dv_gnomad.vcf"
ANN_VCF_CLINVAR="${ANN_DIR}/${SAMPLE_NAME}_dv_clinvar.vcf"

# Map tên file tạm dùng ở bước annotate
TMP1="${ANN_VCF_SNPEFF}"
TMP2="${ANN_VCF_GNOMAD}"
TMP3="${ANN_VCF_CLINVAR}"
ANN_VCF="${ANN_VCF_FINAL}"

# Coverage & Reporting
COVERAGE_PREFIX="${COVERAGE_DIR}/${SAMPLE_NAME}"
MULTIQC_FILENAME="${SAMPLE_NAME}_dv_report.html"

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
# --- BƯỚC 4: GỌI BIẾN THỂ (DV) & CHÚ GIẢI (SnpEff + SnpSift) ---
# ===================================================================
echo ""
echo "---=== BƯỚC 4: Bắt đầu Gọi và Chú giải biến thể (DeepVariant) ===---"

echo "[4.1] DeepVariant (Docker) → VCF/GVCF..."
# (giữ nguyên khối run_deepvariant của bạn)
tabix -f -p vcf "${DV_VCF}"  || true
tabix -f -p vcf "${DV_GVCF}" || true

command -v bcftools >/dev/null 2>&1 || { echo "❌ Cần bcftools trong PATH"; exit 1; }
command -v bgzip    >/dev/null 2>&1 || { echo "❌ Cần bgzip trong PATH"; exit 1; }
command -v tabix    >/dev/null 2>&1 || { echo "❌ Cần tabix trong PATH"; exit 1; }

supports_atomize() { bcftools norm -h 2>&1 | grep -q -- '--atomize'; }

echo "[4.2] Chuẩn hoá + tách đa-allele..."
NORM_VCF_GZ="${ANN_DIR}/${SAMPLE_NAME}_dv.norm.vcf.gz"
bcftools norm -m -both -f "${REF}" -O z -o "${NORM_VCF_GZ}" "${DV_VCF}"
tabix -f "${NORM_VCF_GZ}"

echo "[4.2b] Atomize primitives..."
ATOM_VCF_GZ="${ANN_DIR}/${SAMPLE_NAME}_dv.atom.vcf.gz"
if supports_atomize; then
  bcftools norm --atomize -f "${REF}" -O z -o "${ATOM_VCF_GZ}" "${NORM_VCF_GZ}"
  tabix -f "${ATOM_VCF_GZ}"
elif command -v vt >/dev/null 2>&1; then
  bcftools view -Ov -o "${ANN_DIR}/tmp.dv.norm.vcf" "${NORM_VCF_GZ}"
  vt decompose -s "${ANN_DIR}/tmp.dv.norm.vcf" -o "${ANN_DIR}/tmp.dv.atom.vcf"
  bgzip -f "${ANN_DIR}/tmp.dv.atom.vcf"; tabix -f -p vcf "${ANN_DIR}/tmp.dv.atom.vcf.gz"
  mv -f "${ANN_DIR}/tmp.dv.atom.vcf.gz" "${ATOM_VCF_GZ}"
  mv -f "${ANN_DIR}/tmp.dv.atom.vcf.gz.tbi" "${ATOM_VCF_GZ}.tbi"
  rm -f "${ANN_DIR}/tmp.dv.norm.vcf"
else
  echo "⚠️ Không có --atomize hoặc vt → bỏ qua atomize."
  cp -f "${NORM_VCF_GZ}" "${ATOM_VCF_GZ}"; tabix -f "${ATOM_VCF_GZ}" || true
fi

echo "[4.3] Đồng bộ tiền tố chr với gnomAD..."
ADDCHR_MAP="$(mktemp)"; cat > "$ADDCHR_MAP" <<'EOF'
1 chr1
2 chr2
3 chr3
4 chr4
5 chr5
6 chr6
7 chr7
8 chr8
9 chr9
10 chr10
11 chr11
12 chr12
13 chr13
14 chr14
15 chr15
16 chr16
17 chr17
18 chr18
19 chr19
20 chr20
21 chr21
22 chr22
X chrX
Y chrY
MT chrM
EOF
RMCHR_MAP="$(mktemp)"; awk '{print $2"\t"$1}' "$ADDCHR_MAP" > "$RMCHR_MAP"

HARM_VCF_GZ="${ANN_DIR}/${SAMPLE_NAME}_dv.harm.vcf.gz"
if bcftools view -h "${ATOM_VCF_GZ}" | grep -m1 '^##contig' | grep -q 'ID=chr'; then
  if ! bcftools view -h "${GNOMAD_VCF}" | grep -m1 '^##contig' | grep -q 'ID=chr'; then
    bcftools annotate --rename-chrs "$RMCHR_MAP" -O z -o "${HARM_VCF_GZ}" "${ATOM_VCF_GZ}"
  else
    cp -f "${ATOM_VCF_GZ}" "${HARM_VCF_GZ}" && tabix -f "${HARM_VCF_GZ}" || true
  fi
else
  if bcftools view -h "${GNOMAD_VCF}" | grep -m1 '^##contig' | grep -q 'ID=chr'; then
    bcftools annotate --rename-chrs "$ADDCHR_MAP" -O z -o "${HARM_VCF_GZ}" "${ATOM_VCF_GZ}"
  else
    cp -f "${ATOM_VCF_GZ}" "${HARM_VCF_GZ}" && tabix -f "${HARM_VCF_GZ}" || true
  fi
fi
tabix -f "${HARM_VCF_GZ}"
rm -f "$ADDCHR_MAP" "$RMCHR_MAP"

echo "🔬 [4.4] SnpEff..."
java ${JAVA_OPTS_SNPEFF} -jar "$SNPEFF_JAR" ann -c "$SNPEFF_CONFIG" -v "$SNPEFF_DB" \
  -stats "${SNPEFF_STATS}" \
  "${HARM_VCF_GZ}" > "${TMP1}"

echo "📊 [4.5] gnomAD AF → GNOMAD_AF..."
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -info AF "${GNOMAD_VCF}" "${TMP1}" > "${TMP2}"
GN_HDR="${ANN_DIR}/gn_hdr.hdr"
echo '##INFO=<ID=GNOMAD_AF,Number=A,Type=Float,Description="Allele frequency from gnomAD">' > "$GN_HDR"
bcftools annotate -h "$GN_HDR" -c INFO/GNOMAD_AF:=INFO/AF -x INFO/AF -O v -o "${TMP2}.renamed" "${TMP2}"
mv -f "${TMP2}.renamed" "${TMP2}"
rm -f "$GN_HDR"

echo "🧬 [4.6] ClinVar (CLNSIG, CLNDN)..."
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -info CLNSIG,CLNDN "${CLINVAR_VCF}" "${TMP2}" > "${TMP3}"

echo "🌍 [4.7] 1000G AF → KG_AF..."
java ${JAVA_OPTS_SNPSIFT} -jar "$SNPSIFT_JAR" annotate -info AF "${THOUSANDG_VCF}" "${TMP3}" > "${ANN_VCF}"
KG_HDR="${ANN_DIR}/kg_hdr.hdr"
echo '##INFO=<ID=KG_AF,Number=A,Type=Float,Description="Allele frequency from 1000 Genomes Project">' > "$KG_HDR"
bcftools annotate -h "$KG_HDR" -c INFO/KG_AF:=INFO/AF -x INFO/AF -O v -o "${ANN_VCF}.tmp" "${ANN_VCF}"
mv -f "${ANN_VCF}.tmp" "${ANN_VCF}"
rm -f "$KG_HDR"

echo "📦 [4.8] Nén + index..."
bgzip -f "${ANN_VCF}"
tabix -f -p vcf "${ANN_VCF}.gz"
echo "✅ Annotated VCF (DV): ${ANN_VCF}.gz"
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
        "${DV_GVCF}" \
        "${DV_GVCF}.tbi" \
        "${ANN_DIR}/tmp.dv.norm.vcf" \
        "${ANN_DIR}/tmp.dv.atom.vcf" \
        "${ANN_DIR}/tmp.dv.atom.vcf.gz" \
        "${ANN_DIR}/tmp.dv.atom.vcf.gz.tbi" \
        "${ANN_DIR}/${SAMPLE_NAME}_dv.norm.vcf.gz" \
        "${ANN_DIR}/${SAMPLE_NAME}_dv.norm.vcf.gz.tbi" \
        "${ANN_DIR}/${SAMPLE_NAME}_dv.atom.vcf.gz" \
        "${ANN_DIR}/${SAMPLE_NAME}_dv.atom.vcf.gz.tbi" \
        "${ANN_DIR}/${SAMPLE_NAME}_dv.harm.vcf.gz" \
        "${ANN_DIR}/${SAMPLE_NAME}_dv.harm.vcf.gz.tbi" \
        "${TMP2}.renamed" \
        "${ANN_VCF}.tmp" \
        "${ANN_DIR}/gn_hdr.hdr" "${ANN_DIR}/kg_hdr.hdr" \
        "${TMP1}" "${TMP2}" "${TMP3}" \
        || true
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
echo "==  DeepVariant VCF: ${DV_VCF}"
echo "==  Annotated VCF  : ${ANN_VCF} (và ${ANN_VCF}.gz nếu đã nén)"
echo "==  Báo cáo MultiQC: ${MULTIQC_DIR}/${MULTIQC_FILENAME}"
echo "==================================================================="
