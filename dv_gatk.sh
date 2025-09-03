#!/usr/bin/env bash
# ===================================================================
# BRCA pipeline (1 SAMPLE) chạy song song DeepVariant & GATK
# Bước: 1) FastQC raw → Trim → FastQC trimmed
#       2) BWA MEM (sort + markdup)
#       3) DeepVariant (dùng BAM dedup)
#       4) BQSR (cho nhánh GATK)
#       5) GATK HaplotypeCaller + GenotypeGVCFs
#       6) HARD FILTER (GATK) + ANNOTATION (SnpEff/SnpSift) cho DV & GATK
#       7) MultiQC
#       8) Cleanup
# ===================================================================
set -euo pipefail

# ======================= THAM SỐ & CẤU HÌNH =======================
SAMPLE_NAME="${1:?Vui lòng truyền SAMPLE_NAME}"
THREADS="${THREADS:-8}"
CLEANUP="${CLEANUP:-true}"
TIMESTAMP() { date '+%Y-%m-%d %H:%M:%S'; }

# Conda envs
ENV_BRCA="${ENV_BRCA:-BRCA}"     # fastqc, trimmomatic, bwa, samtools
ENV_GATK="${ENV_GATK:-GATK}"     # gatk
ENV_MQC="${ENV_MQC:-MQC}"        # multiqc
ENV_ANN="${ENV_ANN:-BCF}"        # bcftools, tabix/bgzip, java (snpEff/SnpSift)

# Project + IO
PROJECT_DIR="${PROJECT_DIR:-/media/shmily/writable/BRCA_project}"
SAMPLE_DIR="${PROJECT_DIR}/results/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_DIR}"

# Tham chiếu & targets
REF="${PROJECT_DIR}/reference/Homo_sapiens_assembly38.fasta"
TARGET_BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Known sites cho BQSR
KNOWN_SNP="${PROJECT_DIR}/reference/known_sites/Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_INDEL="${PROJECT_DIR}/reference/known_sites/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS_1000G_INDEL="${PROJECT_DIR}/reference/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# Tools
FASTQC_BIN="${FASTQC_BIN:-fastqc}"
TRIMMOMATIC_BIN="${TRIMMOMATIC_BIN:-trimmomatic}"
BWA_BIN="${BWA_BIN:-bwa}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"
GATK_BIN="${GATK_BIN:-gatk}"
MULTIQC_BIN="${MULTIQC_BIN:-multiqc}"
TABIX_BIN="${TABIX_BIN:-tabix}"
BGZIP_BIN="${BGZIP_BIN:-bgzip}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"

# DeepVariant (Docker)
DV_DOCKER_IMAGE="${DV_DOCKER_IMAGE:-google/deepvariant:1.9.0}"

# Picard + snpEff/SnpSift
PICARD_JAR="${PICARD_JAR:-/home/shmily/miniconda/envs/BRCA/share/picard-2.20.4-0/picard.jar}"
SNPEFF_HOME="${PROJECT_DIR}/snpEff"
SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
SNPSIFT_JAR="${SNPEFF_HOME}/SnpSift.jar"
SNPEFF_CONFIG="${SNPEFF_HOME}/snpEff.config"
SNPEFF_DB="${SNPEFF_DB:-GRCh38.86}"
SNPEFF_DATA_DIR="${SNPEFF_HOME}/data"

JAVA_OPTS_SNPEFF="-Xmx6g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_SNPSIFT="-Xmx4g -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"
JAVA_OPTS_PICARD="-Xmx4g -Djava.awt.headless=true -XX:+UseParallelGC -XX:ParallelGCThreads=${THREADS}"

# External annotation
GNOMAD_VCF="${PROJECT_DIR}/reference/resources/gnomad.v4.1.panel.merged.vcf.gz"
CLINVAR_VCF="${PROJECT_DIR}/reference/resources/clinvar_20250810.vcf.gz"
THOUSANDG_VCF="${PROJECT_DIR}/reference/resources/1000g.panel.merged.vcf.gz"

# Input FASTQ
RAW_DIR="${PROJECT_DIR}/raw_data"
READ1="${RAW_DIR}/${SAMPLE_NAME}_1.fastq.gz"
READ2="${RAW_DIR}/${SAMPLE_NAME}_2.fastq.gz"
ADAPTER_FILE="${ADAPTER_FILE:-/home/shmily/miniconda/envs/BRCA/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa}"

# Thư mục bước
TRIM_DIR="${SAMPLE_DIR}/trimmed_data"
FASTQC_RAW_DIR="${SAMPLE_DIR}/fastqc_raw"
FASTQC_TRIM_DIR="${SAMPLE_DIR}/fastqc_trimmed"
BWA_DIR="${SAMPLE_DIR}/Bwa_alignments"
RECAL_DIR="${SAMPLE_DIR}/recal"
HAPLO_DIR="${SAMPLE_DIR}/haplotypecaller"
DEEP_DIR="${SAMPLE_DIR}/Deepvariants"
ANN_DIR="${SAMPLE_DIR}/snpeff"
MULTIQC_DIR="${SAMPLE_DIR}/multiqc_report"
mkdir -p "${TRIM_DIR}" "${FASTQC_RAW_DIR}" "${FASTQC_TRIM_DIR}" \
         "${BWA_DIR}" "${RECAL_DIR}" "${HAPLO_DIR}" "${DEEP_DIR}" \
         "${ANN_DIR}" "${MULTIQC_DIR}"

# File trung gian & đầu ra
TRIMMED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_paired.fastq.gz"
TRIMMED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_paired.fastq.gz"
UNPAIRED_R1="${TRIM_DIR}/${SAMPLE_NAME}_1_unpaired.fastq.gz"
UNPAIRED_R2="${TRIM_DIR}/${SAMPLE_NAME}_2_unpaired.fastq.gz"

BAM="${BWA_DIR}/${SAMPLE_NAME}_aligned.bam"
SORTED_BAM="${BWA_DIR}/${SAMPLE_NAME}_sorted.bam"
BAM_DEDUP="${BWA_DIR}/${SAMPLE_NAME}_dedup.bam"
BAM_DEDUP_BAI="${BWA_DIR}/${SAMPLE_NAME}_dedup.bai"
DEDUP_METRICS="${BWA_DIR}/${SAMPLE_NAME}_dedup_metrics.txt"

RECAL_DATA_TABLE="${RECAL_DIR}/${SAMPLE_NAME}_recal_data.table"
RECAL_BAM="${RECAL_DIR}/${SAMPLE_NAME}_recalibrated.bam"

# DeepVariant outputs
DV_VCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.vcf.gz"
DV_GVCF="${DEEP_DIR}/${SAMPLE_NAME}_deepvariant.g.vcf.gz"

# GATK outputs
HAPLO_GVCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.g.vcf.gz"
HAPLO_VCF="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.vcf.gz"
HAPLO_VCF_FILTERED="${HAPLO_DIR}/${SAMPLE_NAME}_gatk.filtered.vcf.gz"

# Annotate outputs
ANN_DV_SNPEFF="${ANN_DIR}/${SAMPLE_NAME}_dv_snpeff.vcf"
ANN_DV_GNOMAD="${ANN_DIR}/${SAMPLE_NAME}_dv_gnomad.vcf"
ANN_DV_GNOMAD_R="${ANN_DIR}/${SAMPLE_NAME}_dv_gnomad_renamed.vcf"
ANN_DV_CLINVAR="${ANN_DIR}/${SAMPLE_NAME}_dv_clinvar.vcf"
ANN_DV_1KG="${ANN_DIR}/${SAMPLE_NAME}_dv_1kg.vcf"
ANN_DV_FINAL="${ANN_DIR}/${SAMPLE_NAME}_dv_final_annotated.vcf"

ANN_GATK_SNPEFF="${ANN_DIR}/${SAMPLE_NAME}_snpeff.vcf"
ANN_GATK_GNOMAD="${ANN_DIR}/${SAMPLE_NAME}_gnomad.vcf"
ANN_GATK_GNOMAD_R="${ANN_DIR}/${SAMPLE_NAME}_gnomad_renamed.vcf"
ANN_GATK_CLINVAR="${ANN_DIR}/${SAMPLE_NAME}_clinvar.vcf"
ANN_GATK_1KG="${ANN_DIR}/${SAMPLE_NAME}_1kg.vcf"
ANN_GATK_FINAL="${ANN_DIR}/${SAMPLE_NAME}_final_annotated.vcf"

# Tệp tạm cho chuẩn hoá
NORM_DV="${ANN_DIR}/${SAMPLE_NAME}_dv.norm.vcf.gz"
ATOM_DV="${ANN_DIR}/${SAMPLE_NAME}_dv.atom.vcf.gz"
HARM_DV="${ANN_DIR}/${SAMPLE_NAME}_dv.harm.vcf.gz"

NORM_GATK="${ANN_DIR}/${SAMPLE_NAME}.norm.vcf.gz"
ATOM_GATK="${ANN_DIR}/${SAMPLE_NAME}.atom.vcf.gz"
HARM_GATK="${ANN_DIR}/${SAMPLE_NAME}.harm.vcf.gz"

# MultiQC
MULTIQC_FILENAME="${SAMPLE_NAME}_report.html"

# ======================= HÀM PHỤ TRỢ ==============================
need_cmd()  { command -v "$1" >/dev/null 2>&1 || { echo "❌ Thiếu tool: $1"; exit 1; }; }
need_file() { [[ -f "$1" ]] || { echo "❌ Thiếu file: $1"; exit 1; }; }

supports_bcftools_atomize() {
  local exe; exe="$(command -v "${BCFTOOLS_BIN}" || true)"; [[ -n "$exe" ]] || return 1
  local help; help="$("$exe" norm -h 2>&1 || true)"
  case "$help" in *"--atomize"*) return 0 ;; *) return 1 ;; esac
}

has_chr() {
  "${BCFTOOLS_BIN}" view -h "$1" | grep -m1 '^##contig' | grep -q 'ID=chr' && return 0 || return 1
}

MKMAP_ADDCHR="$(mktemp)"; cat > "$MKMAP_ADDCHR" <<'EOF'
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
MKMAP_RMCHR="$(mktemp)"; awk '{print $2"\t"$1}' "$MKMAP_ADDCHR" > "$MKMAP_RMCHR"

harmonize_to_resource() {
  local sample="$1" resource="$2" out="$3"
  if has_chr "$sample" && ! has_chr "$resource"; then
    echo "[$(TIMESTAMP)] 🔁 Bỏ 'chr' để khớp $(basename "$resource")"
    "${BCFTOOLS_BIN}" annotate --rename-chrs "$MKMAP_RMCHR" -O z -o "$out" "$sample"
    "${TABIX_BIN}" -f -p vcf "$out"
  elif ! has_chr "$sample" && has_chr "$resource"; then
    echo "[$(TIMESTAMP)] 🔁 Thêm 'chr' để khớp $(basename "$resource")"
    "${BCFTOOLS_BIN}" annotate --rename-chrs "$MKMAP_ADDCHR" -O z -o "$out" "$sample"
    "${TABIX_BIN}" -f -p vcf "$out"
  else
    cp -f "$sample" "$out"
    "${TABIX_BIN}" -f -p vcf "$out" || true
  fi
}

rename_AF_to_new_tag() {
  local in="$1" newtag="$2" out="$3"
  local hdr tmp tsv_gz
  hdr="$(mktemp)"
  echo "##INFO=<ID=${newtag},Number=A,Type=Float,Description=\"Allele frequency from ${newtag}\">" > "$hdr"
  tmp="$(mktemp)"
  "${BCFTOOLS_BIN}" query -f'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$in" > "$tmp" || true
  tsv_gz="${tmp}.gz"; "${BGZIP_BIN}" -f -c "$tmp" > "$tsv_gz" && "${TABIX_BIN}" -f -s 1 -b 2 -e 2 "$tsv_gz"
  "${BCFTOOLS_BIN}" annotate -a "$tsv_gz" -c CHROM,POS,REF,ALT,INFO/"$newtag" -h "$hdr" -x INFO/AF -O v -o "$out" "$in"
  rm -f "$hdr" "$tmp" "$tsv_gz" "${tsv_gz}.tbi"
}

ensure_ref_index() {
  if [[ ! -f "${REF}.fai" ]]; then
    echo "[$(TIMESTAMP)] 🧩 Tạo index FASTA tham chiếu (.fai)..."
    "${SAMTOOLS_BIN}" faidx "$REF"
  fi
}

# ========================= KIỂM TRA INPUT ==========================
eval "$(conda shell.bash hook)"
for f in "${READ1}" "${READ2}" "${ADAPTER_FILE}" "${REF}" "${TARGET_BED}"; do need_file "$f"; done
ensure_ref_index
[[ -d "${SNPEFF_DATA_DIR}/${SNPEFF_DB}" ]] || { echo "[$(TIMESTAMP)] ⬇️  Tải SnpEff DB ${SNPEFF_DB}..."; java -jar "${SNPEFF_JAR}" download "${SNPEFF_DB}" -c "${SNPEFF_CONFIG}"; }

# ===================================================================
# 1) QC & TRIMMING
# ===================================================================
echo -e "\n---=== BƯỚC 1: QC & Trimming ===---"
conda activate "${ENV_BRCA}"
"${FASTQC_BIN}" --threads "${THREADS}" -o "${FASTQC_RAW_DIR}" "${READ1}" "${READ2}"
"${TRIMMOMATIC_BIN}" PE -threads "${THREADS}" -phred33 \
  "${READ1}" "${READ2}" \
  "${TRIMMED_R1}" "${UNPAIRED_R1}" \
  "${TRIMMED_R2}" "${UNPAIRED_R2}" \
  ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
"${FASTQC_BIN}" --threads "${THREADS}" -o "${FASTQC_TRIM_DIR}" "${TRIMMED_R1}" "${TRIMMED_R2}"

# ===================================================================
# 2) ALIGN + DEDUP
# ===================================================================
echo -e "\n---=== BƯỚC 2: BWA MEM + sort + markdup ===---"
"${BWA_BIN}" mem -Y -K 100000000 -t "${THREADS}" -R "@RG\tID:${SAMPLE_NAME}_RG\tSM:${SAMPLE_NAME}\tPL:Illumina\tLB:Lib1\tPU:${SAMPLE_NAME}_Illumina_Lib1" \
  "${REF}" "${TRIMMED_R1}" "${TRIMMED_R2}" | "${SAMTOOLS_BIN}" view -b -o "${BAM}" -
java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" SortSam I="${BAM}" O="${SORTED_BAM}" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
java ${JAVA_OPTS_PICARD} -jar "${PICARD_JAR}" MarkDuplicates I="${SORTED_BAM}" O="${BAM_DEDUP}" M="${DEDUP_METRICS}" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=2000000
[[ -f "${BAM_DEDUP_BAI}" ]] || "${SAMTOOLS_BIN}" index "${BAM_DEDUP}"

# ===================================================================
# 3) DEEPVARIANT (từ BAM dedup)
# ===================================================================
echo -e "\n---=== BƯỚC 3: DeepVariant ===---"
docker run --rm -u "$(id -u)":"$(id -g)" -v "${PROJECT_DIR}:${PROJECT_DIR}" "${DV_DOCKER_IMAGE}" \
  /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref="${REF}" \
    --reads="${BAM_DEDUP}" \
    --regions="${TARGET_BED}" \
    --output_vcf="${DV_VCF}" \
    --output_gvcf="${DV_GVCF}" \
    --num_shards="${THREADS}" \
    --vcf_stats_report=true
"${TABIX_BIN}" -f -p vcf "${DV_VCF}" || true
"${TABIX_BIN}" -f -p vcf "${DV_GVCF}" || true

# ===================================================================
# 4) BQSR (cho nhánh GATK)
# ===================================================================
echo -e "\n---=== BƯỚC 4: BQSR (GATK) ===---"
conda activate "${ENV_GATK}"
"${GATK_BIN}" BaseRecalibrator -I "${BAM_DEDUP}" -R "${REF}" \
  --known-sites "${KNOWN_SNP}" --known-sites "${KNOWN_INDEL}" --known-sites "${MILLS_1000G_INDEL}" \
  -O "${RECAL_DATA_TABLE}"
"${GATK_BIN}" ApplyBQSR -R "${REF}" -I "${BAM_DEDUP}" --bqsr-recal-file "${RECAL_DATA_TABLE}" -O "${RECAL_BAM}"
"${SAMTOOLS_BIN}" index "${RECAL_BAM}"

# ===================================================================
# 5) HaplotypeCaller + GenotypeGVCFs (GATK)
# ===================================================================
echo -e "\n---=== BƯỚC 5: HaplotypeCaller + GenotypeGVCFs ===---"
"${GATK_BIN}" HaplotypeCaller -R "${REF}" -I "${RECAL_BAM}" -O "${HAPLO_GVCF}" -L "${TARGET_BED}" -ERC GVCF
"${TABIX_BIN}" -f -p vcf "${HAPLO_GVCF}"
"${GATK_BIN}" GenotypeGVCFs \
  -R "${REF}" \
  -V "${HAPLO_GVCF}" \
  -O "${HAPLO_VCF}" \
  -L "${TARGET_BED}" \
  --only-output-calls-starting-in-intervals true


# ===================================================================
# 6) HARD FILTER (GATK) + ANNOTATION (DV & GATK)
# ===================================================================
conda activate "${ENV_GATK}"
echo -e "\n---=== BƯỚC 6A: Hard filter GATK (Best Practices) ===---"
SNP_VCF="${HAPLO_DIR}/${SAMPLE_NAME}.snps.vcf.gz"
INDEL_VCF="${HAPLO_DIR}/${SAMPLE_NAME}.indels.vcf.gz"
SNP_FILT="${HAPLO_DIR}/${SAMPLE_NAME}.snps.filtered.vcf.gz"
INDEL_FILT="${HAPLO_DIR}/${SAMPLE_NAME}.indels.filtered.vcf.gz"

"${GATK_BIN}" SelectVariants -V "${HAPLO_VCF}" --select-type-to-include SNP   -O "${SNP_VCF}"
"${GATK_BIN}" SelectVariants -V "${HAPLO_VCF}" --select-type-to-include INDEL -O "${INDEL_VCF}"

# --- SNP filters ---
"${GATK_BIN}" VariantFiltration \
  -V "${SNP_VCF}" \
  --missing-values-evaluate-as-failing false \
  -filter "QD < 2.0"                                --filter-name "QD2" \
  -filter "QUAL < 30.0"                             --filter-name "QUAL30" \
  -filter "SOR > 3.0"                               --filter-name "SOR3" \
  -filter "FS > 60.0"                               --filter-name "FS60" \
  -filter "vc.hasAttribute('MQ') && MQ < 40.0"      --filter-name "MQ40" \
  -filter "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
                                                    --filter-name "MQRankSum-12.5" \
  -filter "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
                                                    --filter-name "ReadPosRankSum-8" \
  -O "${SNP_FILT}"

# --- INDEL filters ---
"${GATK_BIN}" VariantFiltration \
  -V "${INDEL_VCF}" \
  --missing-values-evaluate-as-failing false \
  -filter "QD < 2.0"                                  --filter-name "QD2" \
  -filter "QUAL < 30.0"                               --filter-name "QUAL30" \
  -filter "FS > 200.0"                                --filter-name "FS200" \
  -filter "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
                                                      --filter-name "ReadPosRankSum-20" \
  -O "${INDEL_FILT}"

# Gộp lại
"${GATK_BIN}" MergeVcfs -I "${SNP_FILT}" -I "${INDEL_FILT}" -O "${HAPLO_VCF_FILTERED}"
"${TABIX_BIN}" -f -p vcf "${HAPLO_VCF_FILTERED}"

# ---- Annotation cho DV & GATK ----
conda activate "${ENV_ANN}"
for c in "${BCFTOOLS_BIN}" "${TABIX_BIN}" "${BGZIP_BIN}" java; do need_cmd "$c"; done
need_file "${GNOMAD_VCF}"; need_file "${CLINVAR_VCF}"; need_file "${THOUSANDG_VCF}"

# == DV ==
echo -e "\n--- Annotate DeepVariant ---"
"${BCFTOOLS_BIN}" norm -m -both -f "${REF}" -O z -o "${NORM_DV}" "${DV_VCF}"; "${TABIX_BIN}" -f -p vcf "${NORM_DV}"
ATOMIZE="${ATOMIZE:-true}"
if [[ "${ATOMIZE}" == "true" ]] && supports_bcftools_atomize; then
  "${BCFTOOLS_BIN}" norm --atomize -f "${REF}" -O z -o "${ATOM_DV}" "${NORM_DV}"; "${TABIX_BIN}" -f -p vcf "${ATOM_DV}"
else
  cp -f "${NORM_DV}" "${ATOM_DV}"; "${TABIX_BIN}" -f -p vcf "${ATOM_DV}" || true
fi
harmonize_to_resource "${ATOM_DV}" "${GNOMAD_VCF}" "${HARM_DV}"

java ${JAVA_OPTS_SNPEFF}  -jar "${SNPEFF_JAR}"  ann -c "${SNPEFF_CONFIG}" -v "${SNPEFF_DB}" -stats "${ANN_DIR}/${SAMPLE_NAME}_dv" "${HARM_DV}" > "${ANN_DV_SNPEFF}"
java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -info AF "${GNOMAD_VCF}"    "${ANN_DV_SNPEFF}"   > "${ANN_DV_GNOMAD}"
rename_AF_to_new_tag "${ANN_DV_GNOMAD}" "GNOMAD_AF" "${ANN_DV_GNOMAD_R}"
java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -info CLNSIG,CLNDN "${CLINVAR_VCF}" "${ANN_DV_GNOMAD_R}" > "${ANN_DV_CLINVAR}"
java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -info AF "${THOUSANDG_VCF}"  "${ANN_DV_CLINVAR}" > "${ANN_DV_1KG}"
rename_AF_to_new_tag "${ANN_DV_1KG}" "KG_AF" "${ANN_DV_FINAL}"
"${BGZIP_BIN}" -f "${ANN_DV_FINAL}"; "${TABIX_BIN}" -f -p vcf "${ANN_DV_FINAL}.gz"

# == GATK (sau hard-filter) ==
echo -e "\n--- Annotate GATK (sau hard-filter) ---"
"${BCFTOOLS_BIN}" norm -m -both -f "${REF}" -O z -o "${NORM_GATK}" "${HAPLO_VCF_FILTERED}"; "${TABIX_BIN}" -f -p vcf "${NORM_GATK}"
if [[ "${ATOMIZE}" == "true" ]] && supports_bcftools_atomize; then
  "${BCFTOOLS_BIN}" norm --atomize -f "${REF}" -O z -o "${ATOM_GATK}" "${NORM_GATK}"; "${TABIX_BIN}" -f -p vcf "${ATOM_GATK}"
else
  cp -f "${NORM_GATK}" "${ATOM_GATK}"; "${TABIX_BIN}" -f -p vcf "${ATOM_GATK}" || true
fi
harmonize_to_resource "${ATOM_GATK}" "${GNOMAD_VCF}" "${HARM_GATK}"

java ${JAVA_OPTS_SNPEFF}  -jar "${SNPEFF_JAR}"  ann -c "${SNPEFF_CONFIG}" -v "${SNPEFF_DB}" -stats "${ANN_DIR}/${SAMPLE_NAME}" "${HARM_GATK}" > "${ANN_GATK_SNPEFF}"
java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -info AF "${GNOMAD_VCF}"    "${ANN_GATK_SNPEFF}"   > "${ANN_GATK_GNOMAD}"
rename_AF_to_new_tag "${ANN_GATK_GNOMAD}" "GNOMAD_AF" "${ANN_GATK_GNOMAD_R}"
java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -info CLNSIG,CLNDN "${CLINVAR_VCF}" "${ANN_GATK_GNOMAD_R}" > "${ANN_GATK_CLINVAR}"
java ${JAVA_OPTS_SNPSIFT} -jar "${SNPSIFT_JAR}" annotate -info AF "${THOUSANDG_VCF}"  "${ANN_GATK_CLINVAR}" > "${ANN_GATK_1KG}"
rename_AF_to_new_tag "${ANN_GATK_1KG}" "KG_AF" "${ANN_GATK_FINAL}"
"${BGZIP_BIN}" -f "${ANN_GATK_FINAL}"; "${TABIX_BIN}" -f -p vcf "${ANN_GATK_FINAL}.gz"

# ===================================================================
# 7) MULTIQC
# ===================================================================
echo -e "\n---=== BƯỚC 7: MultiQC ===---"
conda activate "${ENV_MQC}"
"${MULTIQC_BIN}" "${SAMPLE_DIR}" --outdir "${MULTIQC_DIR}" \
  --title "Báo cáo QC cho mẫu ${SAMPLE_NAME}" --filename "${SAMPLE_NAME}_report.html" --force

# ===================================================================
# 8) CLEANUP
# ===================================================================
echo -e "\n---=== BƯỚC 8: Cleanup ===---"
if [[ "${CLEANUP}" == "true" ]]; then
  rm -f \
    "${TRIMMED_R1}" "${UNPAIRED_R1}" "${TRIMMED_R2}" "${UNPAIRED_R2}" \
    "${BAM}" "${SORTED_BAM}" \
    "${RECAL_DATA_TABLE}" \
    "${HAPLO_GVCF}" "${HAPLO_GVCF}.tbi" \
    "${HAPLO_VCF}" "${HAPLO_VCF}.tbi" \
    "${NORM_DV}" "${NORM_DV}.tbi" "${ATOM_DV}" "${ATOM_DV}.tbi" "${HARM_DV}" "${HARM_DV}.tbi" \
    "${NORM_GATK}" "${NORM_GATK}.tbi" "${ATOM_GATK}" "${ATOM_GATK}.tbi" "${HARM_GATK}" "${HARM_GATK}.tbi" \
    "${ANN_DV_SNPEFF}" "${ANN_DV_GNOMAD}" "${ANN_DV_GNOMAD_R}" "${ANN_DV_CLINVAR}" "${ANN_DV_1KG}" \
    "${ANN_GATK_SNPEFF}" "${ANN_GATK_GNOMAD}" "${ANN_GATK_GNOMAD_R}" "${ANN_GATK_CLINVAR}" "${ANN_GATK_1KG}" \
    "$MKMAP_ADDCHR" "$MKMAP_RMCHR"
else
  echo "Bỏ qua dọn dẹp (CLEANUP=false)"
fi

echo -e "\n==================================================================="
echo "✅ HOÀN TẤT — ${SAMPLE_NAME}"
echo "DV final:     ${ANN_DV_FINAL}.gz"
echo "GATK filtered:${HAPLO_VCF_FILTERED}"
echo "GATK final:   ${ANN_GATK_FINAL}.gz"
echo "MultiQC:      ${MULTIQC_DIR}/${SAMPLE_NAME}_report.html"
echo "==================================================================="


