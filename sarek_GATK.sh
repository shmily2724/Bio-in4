#!/usr/bin/env bash
set -euo pipefail

# === config ===
PROJECT_DIR="/media/shmily/writable/BRCA_project"
BED="$PROJECT_DIR/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"
GENOME="GATK.GRCh38"
INPUT="$PROJECT_DIR/samplesheet_germline.final.csv"          # sample=...,status=normal
OUTDIR="$PROJECT_DIR/sarek_results/germline_hc"

# Thu gọn heap của Nextflow (tùy chọn)
export NXF_OPTS='-Xms512m -Xmx1500m'

nextflow run nf-core/sarek -r 3.5.1 \
  -profile docker \
  --custom_config_base false \
  --input "$INPUT" \
  --outdir "$OUTDIR" \
  --genome "$GENOME" \
  --wes \
  --intervals "$BED" \
  -c sarek_lowres.config \
  --tools haplotypecaller \
  -process.maxMemory '8 GB' \
  -process.maxCpus 4 \
  -process.maxTime '48 h' \
  -resume

#!/usr/bin/env bash
set -euo pipefail

# === config ===
PROJECT_DIR="/home/shmily/BRCA_project"
BED="$PROJECT_DIR/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"
GENOME="GATK.GRCh38"
INPUT="$PROJECT_DIR/samplesheet_germline.nostatus.csv"          # sample=...,status=normal
OUTDIR="$PROJECT_DIR/sarek_results/germline_hc"

# Thu gọn heap của Nextflow (tùy chọn)
export NXF_OPTS='-Xms512m -Xmx1500m'

nextflow run nf-core/sarek -r 3.5.1 \
  -profile docker \
  --custom_config_base false \
  --input "$INPUT" \
  --outdir "$OUTDIR" \
  --genome "$GENOME" \
  --wes \
  --intervals "$BED" \
  -c sarek_lowres.config \
  --tools haplotypecaller \
  -process.maxMemory '8 GB' \
  -process.maxCpus 4 \
  -process.maxTime '48 h' \
  -resume
