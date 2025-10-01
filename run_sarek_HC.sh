#!/usr/bin/env bash
# run_sarek_all.sh
# Chạy HaplotypeCaller → DeepVariant → Mutect2 (tumor-only), mỗi cái đều annotate bằng SnpEff
# Sau khi HOÀN TẤT cả 3 lượt, tự động dọn work/ để giải phóng dung lượng.
set -euo pipefail

# ===== Đường dẫn dự án =====
PROJECT_DIR="/home/shmily/BRCA_project"
CFG="$PROJECT_DIR/sarek_local.config"
PARAMS="$PROJECT_DIR/sarek_local_params.yaml"
SNPEFF_CACHE="$PROJECT_DIR/snpeff_cache/data"
INTERVALS="$PROJECT_DIR/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Samplesheets
GERM_SHEET="$PROJECT_DIR/samplesheet_germline.csv"

# Thư mục dữ liệu trên phân vùng ext4 mới
OUT_BASE="/mnt/sarek_storage/sarek_results"
WORK_DIR="/mnt/sarek_storage/work"

# ===== Helper =====
ts() { date +"%Y%m%d_%H%M%S"; }
die() { echo "❌ $*" >&2; exit 1; }

# ===== Kiểm tra đầu vào =====
for p in "$CFG" "$PARAMS" "$INTERVALS" "$SNPEFF_CACHE"; do
  [[ -e "$p" ]] || die "Thiếu resource: $p"
done
[[ -s "$GERM_SHEET" ]] || die "CSV thiếu hoặc rỗng: $GERM_SHEET"
# SOMA có thể rỗng (nếu bạn chưa có somatic), vẫn tiếp tục chạy germline.

mkdir -p "$OUT_BASE" "$WORK_DIR"

echo "===> 1/3 HaplotypeCaller + SnpEff (germline)"
nextflow run nf-core/sarek -r 3.5.1 -profile docker \
  -c "$CFG" -params-file "$PARAMS" \
  --input "$GERM_SHEET" \
  --intervals "$INTERVALS" \
  --snpeff_cache "$SNPEFF_CACHE" \
  --tools haplotypecaller,snpeff \
  --germline_filtering_tool hard \
  --outdir "$OUT_BASE/haplotypecaller_$(ts)" \
  -work-dir "$WORK_DIR" \
  -resume
echo "✅ Xong! kết quả ở: $OUT_BASE"



#!/usr/bin/env bash
# run_sarek_all.sh
# Chạy HaplotypeCaller → DeepVariant → Mutect2 (tumor-only), mỗi cái đều annotate bằng SnpEff
# Sau khi HOÀN TẤT cả 3 lượt, tự động dọn work/ để giải phóng dung lượng.
set -euo pipefail

# ===== Đường dẫn dự án =====
PROJECT_DIR="/home/shmily/BRCA_project"
CFG="$PROJECT_DIR/sarek_local.config"
PARAMS="$PROJECT_DIR/sarek_local_params.yaml"
SNPEFF_CACHE="$PROJECT_DIR/snpeff_cache/data"
INTERVALS="$PROJECT_DIR/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"

# Samplesheets
GERM_SHEET="$PROJECT_DIR/samplesheet_germline.csv"

# Thư mục dữ liệu trên phân vùng ext4 mới
OUT_BASE="/mnt/sarek_storage/sarek_results"
WORK_DIR="/mnt/sarek_storage/work"

# ===== Helper =====
ts() { date +"%Y%m%d_%H%M%S"; }
die() { echo "❌ $*" >&2; exit 1; }

# ===== Kiểm tra đầu vào =====
for p in "$CFG" "$PARAMS" "$INTERVALS" "$SNPEFF_CACHE"; do
  [[ -e "$p" ]] || die "Thiếu resource: $p"
done
[[ -s "$GERM_SHEET" ]] || die "CSV thiếu hoặc rỗng: $GERM_SHEET"
# SOMA có thể rỗng (nếu bạn chưa có somatic), vẫn tiếp tục chạy germline.

mkdir -p "$OUT_BASE" "$WORK_DIR"

echo "===> 2/3 DeepVariant + SnpEff (germline)"
nextflow run nf-core/sarek -r 3.5.1 -profile docker \
  -c "$CFG" -params-file "$PARAMS" \
  --input "$GERM_SHEET" \
  --intervals "$INTERVALS" \
  --snpeff_cache "$SNPEFF_CACHE" \
  --tools deepvariant,snpeff \
  --outdir "$OUT_BASE/deepvariant_$(ts)" \
  -work-dir "$WORK_DIR" \
  -resume




#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/home/shmily/BRCA_project"
CFG="$PROJECT_DIR/sarek_local.config"
PARAMS="$PROJECT_DIR/sarek_local_params.yaml"
SNPEFF_CACHE="$PROJECT_DIR/snpeff_cache/data"
INTERVALS="$PROJECT_DIR/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"
SOMA_SHEET="$PROJECT_DIR/samplesheet_somatic.csv"

OUT_BASE="/mnt/sarek_storage/sarek_results"
WORK_DIR="/mnt/sarek_storage/work"

ts() { date +"%Y%m%d_%H%M%S"; }
die() { echo "❌ $*" >&2; exit 1; }

for p in "$CFG" "$PARAMS" "$INTERVALS" "$SNPEFF_CACHE"; do
  [[ -e "$p" ]] || die "Thiếu resource: $p"
done
[[ -s "$SOMA_SHEET" ]] || die "CSV somatic thiếu hoặc rỗng: $SOMA_SHEET"

mkdir -p "$OUT_BASE" "$WORK_DIR"

if [[ "$(wc -l < "$SOMA_SHEET")" -ge 2 ]]; then
  echo "===> 3/3 Mutect2 (tumor-only) + SnpEff (somatic)"
  nextflow run nf-core/sarek -r 3.5.1 -profile docker \
    -c "$CFG" -params-file "$PARAMS" \
    --input "$SOMA_SHEET" \
    --intervals "$INTERVALS" \
    --snpeff_cache "$SNPEFF_CACHE" \
    --tools mutect2,snpeff \
    --outdir "$OUT_BASE/mutect2_$(ts)" \
    -work-dir "$WORK_DIR" \
    -resume
else
  echo "⏩ Bỏ qua Mutect2: $SOMA_SHEET không có mẫu (chỉ header hoặc trống)."
fi

echo "🧹 Dọn work/ sau khi hoàn tất..."
nextflow clean -f || true




