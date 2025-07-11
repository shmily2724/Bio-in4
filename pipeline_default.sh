#!/bin/bash

#=== Step 0: Logging utilities ===#
function log_step {
  STEP_NAME="$1"
  echo -e "\n\033[1;34m>>> ${STEP_NAME}...\033[0m"
}

function log_done {
  echo -e "\033[1;32m[Done]\033[0m"
}

#=== Step 1: Check tools in Conda environments ===#
log_step "Step 1: Checking required tools in Conda environments"

function check_tool {
  TOOL_NAME="$1"
  ENV_NAME="$2"

  if conda run -n "$ENV_NAME" which "$TOOL_NAME" &>/dev/null; then
    echo -e "\033[1;32m✔ $TOOL_NAME found in environment $ENV_NAME\033[0m"
  else
    echo -e "\033[1;31m✘ $TOOL_NAME NOT FOUND in environment $ENV_NAME\033[0m"
    MISSING_TOOLS=true
  fi
}

MISSING_TOOLS=false

# Tools in BRCA env
check_tool trimmomatic BRCA
check_tool bwa BRCA
check_tool samtools BRCA
check_tool picard BRCA
check_tool mosdepth BRCA
check_tool depvariance BRCA

# Tools in GATK env
check_tool gatk GATK
check_tool snpEff GATK

if $MISSING_TOOLS; then
  echo -e "\033[1;31mERROR: One or more tools are missing. Please install them before running the pipeline.\033[0m"
  exit 1
fi
log_done

#=== Step 2: Input Configuration ===#
SAMPLE_NAME="$1"
if [[ -z "$SAMPLE_NAME" ]]; then
  echo -e "\033[1;31mERROR: Missing sample name. Usage: bash pipeline_default.sh\033[0m"
  exit 1
fi

READ1="raw_data/${SAMPLE_NAME}_R1.fq.gz"
READ2="raw_data/${SAMPLE_NAME}_R2.fq.gz"
REF="/media/shmily/writable/BRCA_project/reference/reference.fa"
KNOWN_SITES_DIR="/media/shmily/writable/BRCA_project/reference/known_sites"
OUTDIR="results/${SAMPLE_NAME}"
mkdir -p ${OUTDIR}

#=== Step 3: Quality trimming ===#
log_step "Step 3: Quality trimming with Trimmomatic"
eval "$(conda shell.bash hook)"
conda activate BRCA
trimmomatic PE -threads 8 \
  ${READ1} ${READ2} \
  ${OUTDIR}/${SAMPLE_NAME}_R1_paired.fq.gz ${OUTDIR}/${SAMPLE_NAME}_R1_unpaired.fq.gz \
  ${OUTDIR}/${SAMPLE_NAME}_R2_paired.fq.gz ${OUTDIR}/${SAMPLE_NAME}_R2_unpaired.fq.gz \
  SLIDINGWINDOW:4:20 MINLEN:36
log_done

#=== Step 4: Alignment ===#
log_step "Step 4: BWA MEM alignment"
bwa mem -t 8 -M ${REF} \
  ${OUTDIR}/${SAMPLE_NAME}_R1_paired.fq.gz ${OUTDIR}/${SAMPLE_NAME}_R2_paired.fq.gz \
  > ${OUTDIR}/${SAMPLE_NAME}.sam
log_done

#=== Step 5: Convert SAM to BAM + Sort ===#
log_step "Step 5: Convert SAM to sorted BAM"
samtools view -bS ${OUTDIR}/${SAMPLE_NAME}.sam | samtools sort -o ${OUTDIR}/${SAMPLE_NAME}_sorted.bam
rm ${OUTDIR}/${SAMPLE_NAME}.sam
log_done

#=== Step 6: Mark Duplicates ===#
log_step "Step 6: Mark duplicates with Picard"
picard MarkDuplicates \
  I=${OUTDIR}/${SAMPLE_NAME}_sorted.bam \
  O=${OUTDIR}/${SAMPLE_NAME}_marked.bam \
  M=${OUTDIR}/${SAMPLE_NAME}_metrics.txt
samtools index ${OUTDIR}/${SAMPLE_NAME}_marked.bam
log_done

#=== Step 7: Base Recalibration ===#
log_step "Step 7: Base Quality Score Recalibration"
conda activate GATK
gatk BaseRecalibrator \
  -I ${OUTDIR}/${SAMPLE_NAME}_marked.bam \
  -R ${REF} \
  --known-sites ${KNOWN_SITES_DIR}/Mills_and_1000G_gold_standard.indels.hg38.chr_FIXED.vcf.gz \
  --known-sites ${KNOWN_SITES_DIR}/Homo_sapiens_assembly38.dbsnp138.chr_FIXED.vcf.gz \
  -O ${OUTDIR}/recal_data.table

gatk ApplyBQSR \
  -R ${REF} \
  -I ${OUTDIR}/${SAMPLE_NAME}_marked.bam \
  --bqsr-recal-file ${OUTDIR}/recal_data.table \
  -O ${OUTDIR}/${SAMPLE_NAME}_recal.bam
log_done

#=== Step 8: Variant Calling ===#
log_step "Step 8: Variant Calling with HaplotypeCaller"
gatk HaplotypeCaller \
  -R ${REF} \
  -I ${OUTDIR}/${SAMPLE_NAME}_recal.bam \
  -O ${OUTDIR}/${SAMPLE_NAME}.g.vcf.gz \
  -ERC GVCF
log_done

#=== Step 9: Genotype GVCF ===#
log_step "Step 9: Genotyping GVCF to final VCF"
gatk GenotypeGVCFs \
  -R ${REF} \
  -V ${OUTDIR}/${SAMPLE_NAME}.g.vcf.gz \
  -O ${OUTDIR}/${SAMPLE_NAME}.vcf.gz
log_done

#=== Step 10: Annotation with SnpEff ===#
log_step "Step 10: Annotating variants with SnpEff"
SNPEFF_CONFIG="/home/shmily/miniconda/envs/GATK/share/snpeff-5.2-1/snpEff.config"
SNPEFF_DB="GRCh38.86"
snpEff -c $SNPEFF_CONFIG -v $SNPEFF_DB \
  ${OUTDIR}/${SAMPLE_NAME}.vcf.gz > ${OUTDIR}/${SAMPLE_NAME}_annotated.vcf
log_done

#=== Step 11: Depth Analysis ===#
log_step "Step 11: Depth analysis with mosdepth"
conda activate BRCA
mosdepth ${OUTDIR}/${SAMPLE_NAME} ${OUTDIR}/${SAMPLE_NAME}_recal.bam
depvariance ${OUTDIR}/${SAMPLE_NAME}.mosdepth.global.dist.txt > ${OUTDIR}/${SAMPLE_NAME}_coverage_stats.txt
log_done

#=== Finished ===#
echo -e "\n\033[1;35m>>> Pipeline completed for ${SAMPLE_NAME} <<<\033[0m"
