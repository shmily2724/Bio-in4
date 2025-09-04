#!/usr/bin/env bash
# ==============================================================
# SOMATIC (TUMOR-ONLY) — Mutect2 + Filter + SnpEff (khớp layout BRCA của bạn)
# ==============================================================

set -euo pipefail
SAMPLE_NAME="${1:?Vui lòng truyền SAMPLE_NAME}"
THREADS="${THREADS:-8}"

# ==== Kế thừa đúng biến/đường dẫn như script của bạn ====
PROJECT_DIR="${PROJECT_DIR:-/media/shmily/writable/BRCA_project}"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"

REF="${PROJECT_DIR}/reference/Homo_sapiens_assembly38.fasta"
TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Resources cần cho Mutect2 tumor-only
GERMLINE_AF="${PROJECT_DIR}/reference/resources/af-only-gnomad.hg38.vcf.gz"         # cần có
COMMON_SNP="${PROJECT_DIR}/reference/resources/small_exac_common_3.hg38.vcf.gz"    # cho contamination
PON_VCF="${PON_VCF:-}"   # để trống nếu chưa có PoN cùng hệ

# Conda envs theo bạn
ENV_GATK="${ENV_GATK:-GATK}"
ENV_ANN="${ENV_ANN:-BCF}"

# Tools (đồng bộ với script cũ)
GATK_BIN="${GATK_BIN:-gatk}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
TABIX_BIN="${TABIX_BIN:-tabix}"
BGZIP_BIN="${BGZIP_BIN:-bgzip}"

# snpEff/SnpSift (đồng bộ với script cũ)
SNPEFF_HOME="${PROJECT_DIR}/snpEff"
SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_HOME}/SnpSift.jar"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
SNPEFF_DB="${SNPEFF_DB:-GRCh38.86}"
SNPEFF_DATA_DIR="${SNPEFF_HOME}/data"
JAVA_OPTS_SNPEFF="-Xmx6g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_SNPSIFT="-Xmx4g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"

# ==== Input BAM đã BQSR từ pipeline germline hiện có ====
RECAL_DIR="${SAMPLE_DIR}/recal"
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"

# ==== Thư mục/đầu ra somatic (tách riêng) ====
SOM_DIR="${SAMPLE_DIR}/somatic"
ANN_DIR="${SOM_DIR}/snpeff"
mkdir -p "${SOM_DIR}" "${ANN_DIR}"

ts(){ date '+%Y-%m-%d %H:%M:%S'; }
ok(){ [[ -s "$1" ]]; }

# ==== Kiểm tra input tối thiểu ====
for f in "$RECAL_BAM" "$REF" "$TARGET_BED" "$GERMLINE_AF" "$COMMON_SNP"; do
  [[ -s "$f" ]] || { echo "❌ Thiếu: $f"; exit 1; }
done

# ==== Đảm bảo DB snpEff ====
if [[ ! -d "${SNPEFF_DATA_DIR}/${SNPEFF_DB}" ]]; then
  echo "[i] $(ts) Download snpEff DB ${SNPEFF_DB}"
  java -Xmx2g -jar "$SNPEFF_JAR" download "$SNPEFF_DB" -c "$SNPEFF_CONFIG"
fi

# -------------------- 1) Mutect2 --------------------
RAW_VCF="${SOM_DIR}/${SAMPLE_NAME}.mutect2.unfiltered.vcf.gz"
F1R2_TAR="${SOM_DIR}/${SAMPLE_NAME}.f1r2.tar.gz"

if ! ok "$RAW_VCF"; then
  conda activate "$ENV_GATK"
  CMD="${GATK_BIN} --java-options '-Xmx8g' Mutect2 \
      -R '${REF}' -I '${RECAL_BAM}' -tumor '${SAMPLE_NAME}' \
      --germline-resource '${GERMLINE_AF}' \
      --f1r2-tar-gz '${F1R2_TAR}' \
      -L '${TARGET_BED}' \
      -O '${RAW_VCF}'"
  [[ -n "$PON_VCF" ]] && CMD+=" --panel-of-normals '${PON_VCF}'"
  echo "[*] $(ts) Mutect2"; eval "$CMD"
  conda deactivate
else
  echo "[✓] Mutect2: đã có → bỏ qua"
fi

# -------------------- 2) LearnReadOrientationModel --------------------
ROM="${SOM_DIR}/${SAMPLE_NAME}.read-orientation-model.tar.gz"
if ! ok "$ROM"; then
  conda activate "$ENV_GATK"
  echo "[*] $(ts) LearnReadOrientationModel"
  ${GATK_BIN} LearnReadOrientationModel -I "${F1R2_TAR}" -O "${ROM}"
  conda deactivate
else
  echo "[✓] LearnReadOrientationModel: đã có → bỏ qua"
fi

# -------------------- 3) GetPileupSummaries --------------------
PUS="${SOM_DIR}/${SAMPLE_NAME}.pileups.table"
if ! ok "$PUS"; then
  conda activate "$ENV_GATK"
  echo "[*] $(ts) GetPileupSummaries"
  ${GATK_BIN} --java-options '-Xmx4g' GetPileupSummaries \
    -I "${RECAL_BAM}" \
    -V "${COMMON_SNP}" \
    -L "${COMMON_SNP}" \
    -L "${TARGET_BED}" \
    -O "${PUS}"
  conda deactivate
else
  echo "[✓] GetPileupSummaries: đã có → bỏ qua"
fi

# -------------------- 4) CalculateContamination --------------------
CONTAM="${SOM_DIR}/${SAMPLE_NAME}.contamination.table"
SEGS="${SOM_DIR}/${SAMPLE_NAME}.segments.table"
if ! ok "$CONTAM"; then
  conda activate "$ENV_GATK"
  echo "[*] $(ts) CalculateContamination"
  ${GATK_BIN} CalculateContamination -I "${PUS}" -O "${CONTAM}" --tumor-segmentation "${SEGS}"
  conda deactivate
else
  echo "[✓] CalculateContamination: đã có → bỏ qua"
fi

# -------------------- 5) FilterMutectCalls --------------------
FILT_VCF="${SOM_DIR}/${SAMPLE_NAME}.mutect2.filtered.vcf.gz"
if ! ok "$FILT_VCF"; then
  conda activate "$ENV_GATK"
  echo "[*] $(ts) FilterMutectCalls"
  ${GATK_BIN} --java-options '-Xmx6g' FilterMutectCalls \
    -R "${REF}" \
    -V "${RAW_VCF}" \
    --contamination-table "${CONTAM}" \
    --tumor-segmentation "${SEGS}" \
    --ob-priors "${ROM}" \
    -L "${TARGET_BED}" \
    -O "${FILT_VCF}"
  conda deactivate
else
  echo "[✓] FilterMutectCalls: đã có → bỏ qua"
fi

# -------------------- 6) Lấy PASS --------------------
PASS_VCF="${SOM_DIR}/${SAMPLE_NAME}.somatic.PASS.vcf.gz"
if ! ok "$PASS_VCF"; then
  conda activate "$ENV_ANN"
  echo "[*] $(ts) Extract PASS"
  ${BCFTOOLS_BIN} view -f PASS -Oz -o "${PASS_VCF}" "${FILT_VCF}"
  ${TABIX_BIN} -f -p vcf "${PASS_VCF}"
  conda deactivate
else
  echo "[✓] PASS VCF: đã có → bỏ qua"
fi

# -------------------- 7) Annotate (snpEff + SnpSift) --------------------
ANN_VCF="${ANN_DIR}/${SAMPLE_NAME}_somatic_annotated.vcf.gz"
if ! ok "$ANN_VCF"; then
  conda activate "$ENV_ANN"
  echo "[*] $(ts) snpEff"
  java ${JAVA_OPTS_SNPEFF} -jar "${SNPEFF_JAR}" ann -c "${SNPEFF_CONFIG}" -v "${SNPEFF_DB}" \
    "${PASS_VCF}" | bgzip -c > "${ANN_VCF}"
  ${TABIX_BIN} -f -p vcf "${ANN_VCF}"
  conda deactivate
else
  echo "[✓] snpEff: đã có → bỏ qua"
fi

echo -e "\n✅ $(ts) HOÀN TẤT somatic (tumor-only) cho ${SAMPLE_NAME}"
echo "  PASS VCF : ${PASS_VCF}"
echo "  Annotated: ${ANN_VCF}"
