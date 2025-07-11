#!/bin/bash

#=== Step 0: Logging utilities ===#
function log_step {
  STEP_NAME="$1"
  echo -e "\n\033[1;34m>>> ${STEP_NAME}...\033[0m"
}

function log_done {
  echo -e "\033[1;32m[Done]\033[0m"
}

function log_error {
  ERROR_MESSAGE="$1"
  echo -e "\033[1;31mERROR: ${ERROR_MESSAGE}\033[0m"
  exit 1
}

# Ngắt script ngay lập tức nếu bất kỳ lệnh nào thất bại
set -e

#=== Step 1: Check tools in Conda environments ===#
log_step "Step 1: Checking required tools in Conda environments"

function check_tool {
  TOOL_NAME="$1"
  ENV_NAME="$2"

  if ! conda run -n "$ENV_NAME" which "$TOOL_NAME" &>/dev/null; then
    log_error "$TOOL_NAME NOT FOUND in environment $ENV_NAME. Please install it."
  else
    echo -e "\033[1;32m✔ $TOOL_NAME found in environment $ENV_NAME\033[0m"
  fi
}

# Tools in BRCA env
check_tool trimmomatic BRCA
check_tool bwa BRCA
check_tool samtools BRCA
check_tool picard BRCA
check_tool mosdepth BRCA
check_tool depvariance BRCA
check_tool fastqc BRCA # Thêm kiểm tra FastQC
check_tool multiqc BRCA # Thêm kiểm tra MultiQC

# Tools in GATK env
check_tool gatk GATK
check_tool snpEff GATK

log_done

#=== Step 2: Input Configuration ===#
SAMPLE_NAME="$1"
if [[ -z "$SAMPLE_NAME" ]]; then
  log_error "Missing sample name. Usage: bash pipeline_default.sh <SAMPLE_NAME>"
fi

READ1="raw_data/${SAMPLE_NAME}_R1.fq.gz"
READ2="raw_data/${SAMPLE_NAME}_R2.fq.gz"
REF="/media/shmily/writable/BRCA_project/reference/reference.fa"
KNOWN_SITES_DIR="/media/shmily/writable/BRCA_project/reference/known_sites"
OUTDIR="results/${SAMPLE_NAME}"

# Kiểm tra sự tồn tại của các file đầu vào
if [[ ! -f "$READ1" ]]; then
  log_error "Input Read1 file not found: ${READ1}"
fi
if [[ ! -f "$READ2" ]]; then
  log_error "Input Read2 file not found: ${READ2}"
fi
if [[ ! -f "$REF" ]]; then
  log_error "Reference file not found: ${REF}"
fi
if [[ ! -f "${KNOWN_SITES_DIR}/Mills_and_1000G_gold_standard.indels.hg38.chr_FIXED.vcf.gz" ]]; then
  log_error "Known sites VCF (Mills_and_1000G) not found."
fi
if [[ ! -f "${KNOWN_SITES_DIR}/Homo_sapiens_assembly38.dbsnp138.chr_FIXED.vcf.gz" ]]; then
  log_error "Known sites VCF (dbsnp) not found."
fi

mkdir -p ${OUTDIR}

#=== Bước tiền xử lý Conda hook (chỉ chạy 1 lần) ===#
eval "$(conda shell.bash hook)"

#=== Step 3: Quality trimming ===#
log_step "Step 3: Quality trimming with Trimmomatic"
conda activate BRCA
trimmomatic PE -threads 8 \
  ${READ1} ${READ2} \
  ${OUTDIR}/${SAMPLE_NAME}_R1_paired.fq.gz ${OUTDIR}/${SAMPLE_NAME}_R1_unpaired.fq.gz \
  ${OUTDIR}/${SAMPLE_NAME}_R2_paired.fq.gz ${OUTDIR}/${SAMPLE_NAME}_R2_unpaired.fq.gz \
  SLIDINGWINDOW:4:20 MINLEN:36 \
  &>> ${OUTDIR}/trimmomatic.log # Ghi log vào file
log_done

#=== Step 4: Quality Control on Trimmed Reads (FastQC) ===#
log_step "Step 4: Running FastQC on trimmed reads"
fastqc ${OUTDIR}/${SAMPLE_NAME}_R1_paired.fq.gz -o ${OUTDIR} &>> ${OUTDIR}/fastqc_trimmed.log
fastqc ${OUTDIR}/${SAMPLE_NAME}_R2_paired.fq.gz -o ${OUTDIR} &>> ${OUTDIR}/fastqc_trimmed.log
log_done

#=== Step 5: Alignment ===#
log_step "Step 5: BWA MEM alignment"
bwa mem -t 8 -M ${REF} \
  ${OUTDIR}/${SAMPLE_NAME}_R1_paired.fq.gz ${OUTDIR}/${SAMPLE_NAME}_R2_paired.fq.gz \
  > ${OUTDIR}/${SAMPLE_NAME}.sam \
  2>> ${OUTDIR}/bwa_mem.log # Ghi log lỗi vào file

# Dọn dẹp reads unpaired (nếu không cần thiết cho pipeline của bạn)
rm ${OUTDIR}/${SAMPLE_NAME}_R1_unpaired.fq.gz || true
rm ${OUTDIR}/${SAMPLE_NAME}_R2_unpaired.fq.gz || true
log_done

#=== Step 6: Convert SAM to BAM + Sort ===#
log_step "Step 6: Convert SAM to sorted BAM"
samtools view -bS ${OUTDIR}/${SAMPLE_NAME}.sam | samtools sort -o ${OUTDIR}/${SAMPLE_NAME}_sorted.bam \
  &>> ${OUTDIR}/samtools_sort.log
rm ${OUTDIR}/${SAMPLE_NAME}.sam # Xóa file .sam trung gian
log_done

#=== Step 7: Mark Duplicates ===#
log_step "Step 7: Mark duplicates with Picard"
# Tối ưu: cấp phát thêm RAM cho Picard nếu cần (ví dụ: -Xmx10g)
picard MarkDuplicates \
  I=${OUTDIR}/${SAMPLE_NAME}_sorted.bam \
  O=${OUTDIR}/${SAMPLE_NAME}_marked.bam \
  M=${OUTDIR}/${SAMPLE_NAME}_metrics.txt \
  &>> ${OUTDIR}/picard_markduplicates.log
samtools index ${OUTDIR}/${SAMPLE_NAME}_marked.bam &>> ${OUTDIR}/samtools_index_marked.log
rm ${OUTDIR}/${SAMPLE_NAME}_sorted.bam # Xóa file .bam đã sắp xếp chưa mark
log_done

#=== Step 8: Base Recalibration ===#
log_step "Step 8: Base Quality Score Recalibration"
conda activate GATK # Kích hoạt lại môi trường GATK
gatk BaseRecalibrator \
  -I ${OUTDIR}/${SAMPLE_NAME}_marked.bam \
  -R ${REF} \
  --known-sites ${KNOWN_SITES_DIR}/Mills_and_1000G_gold_standard.indels.hg38.chr_FIXED.vcf.gz \
  --known-sites ${KNOWN_SITES_DIR}/Homo_sapiens_assembly38.dbsnp138.chr_FIXED.vcf.gz \
  -O ${OUTDIR}/recal_data.table \
  &>> ${OUTDIR}/gatk_baserecalibrator.log

gatk ApplyBQSR \
  -R ${REF} \
  -I ${OUTDIR}/${SAMPLE_NAME}_marked.bam \
  --bqsr-recal-file ${OUTDIR}/recal_data.table \
  -O ${OUTDIR}/${SAMPLE_NAME}_recal.bam \
  &>> ${OUTDIR}/gatk_applybqsr.log

# Xóa các file trung gian sau BQSR
rm ${OUTDIR}/recal_data.table || true
rm ${OUTDIR}/${SAMPLE_NAME}_marked.bam || true
rm ${OUTDIR}/${SAMPLE_NAME}_marked.bam.bai || true # Hoặc .bam.csi
log_done

#=== Step 9: Variant Calling ===#
log_step "Step 9: Variant Calling with HaplotypeCaller"
gatk HaplotypeCaller \
  -R ${REF} \
  -I ${OUTDIR}/${SAMPLE_NAME}_recal.bam \
  -O ${OUTDIR}/${SAMPLE_NAME}.g.vcf.gz \
  -ERC GVCF \
  &>> ${OUTDIR}/gatk_haplotypecaller.log
log_done

#=== Step 10: Genotype GVCF ===#
log_step "Step 10: Genotyping GVCF to final VCF"
gatk GenotypeGVCFs \
  -R ${REF} \
  -V ${OUTDIR}/${SAMPLE_NAME}.g.vcf.gz \
  -O ${OUTDIR}/${SAMPLE_NAME}.vcf.gz \
  &>> ${OUTDIR}/gatk_genotypegvcfs.log
log_done

#=== Step 11: Index final VCF ===#
log_step "Step 11: Indexing final VCF"
tabix -p vcf ${OUTDIR}/${SAMPLE_NAME}.vcf.gz &>> ${OUTDIR}/tabix_index.log
log_done

#=== Step 12: Annotation with SnpEff ===#
log_step "Step 12: Annotating variants with SnpEff"
SNPEFF_CONFIG="/home/shmily/miniconda/envs/GATK/share/snpeff-5.2-1/snpEff.config"
SNPEFF_DB="GRCh38.86" # Đảm bảo DB này đã được cài đặt trong snpEff
snpEff -c $SNPEFF_CONFIG -v $SNPEFF_DB \
  ${OUTDIR}/${SAMPLE_NAME}.vcf.gz > ${OUTDIR}/${SAMPLE_NAME}_annotated.vcf \
  2>> ${OUTDIR}/snpeff_annotation.log # Ghi lỗi SnpEff
log_done

#=== Step 13: Depth Analysis ===#
log_step "Step 13: Depth analysis with mosdepth"
conda activate BRCA # Kích hoạt lại môi trường BRCA
mosdepth ${OUTDIR}/${SAMPLE_NAME} ${OUTDIR}/${SAMPLE_NAME}_recal.bam \
  &>> ${OUTDIR}/mosdepth.log
depvariance ${OUTDIR}/${SAMPLE_NAME}.mosdepth.global.dist.txt > ${OUTDIR}/${SAMPLE_NAME}_coverage_stats.txt \
  2>> ${OUTDIR}/depvariance.log
log_done

#=== Step 14: MultiQC Report ===#
log_step "Step 14: Generating MultiQC report"
# MultiQC sẽ tìm kiếm các tệp log và báo cáo từ các công cụ trong thư mục OUTDIR
# Bạn có thể chỉ định các thư mục cụ thể nếu muốn, ví dụ: multiqc ${OUTDIR} -o ${OUTDIR}/multiqc_report
multiqc ${OUTDIR} -o ${OUTDIR} --filename multiqc_report.html &>> ${OUTDIR}/multiqc.log
log_done

#=== Finished ===#
echo -e "\n\033[1;35m>>> Pipeline completed for ${SAMPLE_NAME} <<<\033[0m"
echo -e "Check the MultiQC report: ${OUTDIR}/multiqc_report.html"
