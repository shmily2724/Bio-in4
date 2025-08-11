#!/usr/bin/env bash
set -euo pipefail
PROJECT_DIR="/media/shmily/writable/BRCA_project"
BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"
URLS_TSV="/media/shmily/writable/BRCA_project/urls.tsv"
OUTDIR="${PROJECT_DIR}/reference/gnomad_v4.1_panel"
THREADS="${THREADS:-4}"

mkdir -p "$OUTDIR/parts" "$OUTDIR/logs"

# kiểm tra
command -v bcftools >/dev/null 2>&1 || { echo "Cần bcftools trong PATH"; exit 1; }
command -v tabix     >/dev/null 2>&1 || { echo "Cần tabix trong PATH"; exit 1; }
[[ -s "$BED" ]]      || { echo "Không thấy BED: $BED"; exit 1; }
[[ -s "$URLS_TSV" ]] || { echo "Không thấy URLs: $URLS_TSV"; exit 1; }

echo "[INFO] Bắt đầu trích theo BED từ các URL trong $URLS_TSV"

# vòng lặp từng chr
while read -r CHR URL; do
  [[ -z "${CHR:-}" || -z "${URL:-}" ]] && continue
  OUT="$OUTDIR/parts/gnomad.v4.1.sites.${CHR}.panel.vcf.gz"

  if [[ -s "$OUT" ]]; then
    echo "[SKIP] Đã có $OUT"
    continue
  fi

  echo "[RUN ] $CHR"
  # Trích vùng theo BED, nén bgzip
  bcftools view --threads "$THREADS" -R "$BED" -Oz \
    -o "$OUT" "$URL" 2> "$OUTDIR/logs/${CHR}.log"

  tabix -p vcf "$OUT"
done < "$URLS_TSV"

echo "[INFO] Ghép các chr lại…"
# Sắp xếp tên file theo thứ tự tự nhiên
PARTS=$(ls "$OUTDIR"/parts/gnomad.v4.1.sites.chr{1..22}.panel.vcf.gz 2>/dev/null || true)
PARTS_XY=$(ls "$OUTDIR"/parts/gnomad.v4.1.sites.chrX.panel.vcf.gz "$OUTDIR"/parts/gnomad.v4.1.sites.chrY.panel.vcf.gz 2>/dev/null || true)

bcftools concat --threads "$THREADS" -Oz \
  -o "$OUTDIR/gnomad.v4.1.panel.merged.vcf.gz" \
  $PARTS $PARTS_XY

tabix -p vcf "$OUTDIR/gnomad.v4.1.panel.merged.vcf.gz"

echo "[DONE] File cuối: $OUTDIR/gnomad.v4.1.panel.merged.vcf.gz"




URL=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20250810.vcf.gz
PROJECT_DIR="/media/shmily/writable/BRCA_project"
BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"
bcftools view -R "$BED" "$URL" -Oz -o clinvar_20250810.panel.vcf.gz
tabix -p vcf clinvar_20250810.panel.vcf.gz


#!/usr/bin/env bash
set -euo pipefail

# === CẤU HÌNH ===
PROJECT_DIR="/media/shmily/writable/BRCA_project"
BED="${PROJECT_DIR}/reference/TruSight_Cancer_TargetedRegions_v1.0.hg38.bed"
URLS_TSV="${PROJECT_DIR}/urls_1000g.tsv"
OUTDIR="${PROJECT_DIR}/reference/1000g_vcf_panel"
THREADS="${THREADS:-4}"

# Tạo thư mục output nếu chưa có
mkdir -p "$OUTDIR/parts" "$OUTDIR/logs"

# Kiểm tra điều kiện
command -v bcftools >/dev/null || { echo "Cần bcftools trong PATH"; exit 1; }
command -v tabix >/dev/null   || { echo "Cần tabix trong PATH"; exit 1; }
[[ -s "$BED" ]]      || { echo "Không thấy BED: $BED"; exit 1; }
[[ -s "$URLS_TSV" ]] || { echo "Không thấy URLs: $URLS_TSV"; exit 1; }

echo "[INFO] Bắt đầu trích xuất biến thể từ 1000G theo BED"

# === VÒNG LẶP TỪNG CHR ===
while read -r CHR URL; do
  # Bỏ qua dòng trống hoặc dòng chrY
  [[ -z "${CHR:-}" || -z "${URL:-}" ]] && continue
  if [[ "$CHR" == "chrY" ]]; then
    echo "[SKIP] Bỏ qua $CHR vì không có dữ liệu"
    continue
  fi

  OUT="$OUTDIR/parts/1000g.${CHR}.panel.vcf.gz"

  # Nếu đã có output thì bỏ qua
  if [[ -s "$OUT" ]]; then
    echo "[SKIP] Đã tồn tại $OUT"
    continue
  fi

  echo "[RUN ] Trích $CHR từ $URL"
  bcftools view --threads "$THREADS" -R "$BED" -Oz -o "$OUT" "$URL" \
    2> "$OUTDIR/logs/${CHR}.log"

  tabix -p vcf "$OUT"
done < "$URLS_TSV"

# === GHÉP FILES ===
echo "[INFO] Ghép các file chr thành 1 VCF duy nhất..."
PARTS=$(ls "$OUTDIR"/parts/1000g.chr{1..22}.panel.vcf.gz 2>/dev/null || true)

bcftools concat --threads "$THREADS" -Oz \
  -o "$OUTDIR/1000g.panel.merged.vcf.gz" \
  $PARTS

tabix -p vcf "$OUTDIR/1000g.panel.merged.vcf.gz"

echo "[DONE] File cuối: $OUTDIR/1000g.panel.merged.vcf.gz"


