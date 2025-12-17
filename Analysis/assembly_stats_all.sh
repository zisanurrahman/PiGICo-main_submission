#!/usr/bin/env bash
set -euo pipefail

DIR="${1:?Usage: bash assembly_stats_all.sh /path/to/fasta_dir}"

command -v seqkit >/dev/null 2>&1 || { echo "ERROR: seqkit not found. Install: conda install -c bioconda seqkit"; exit 1; }

OUTDIR="/Users/zisan/SGI_Paper/Full_paper_clean/Data/genome_stats_out"
mkdir -p "$OUTDIR"


# 1) Build a safe, NUL-delimited list of only real FASTA files
#    - Ignore AppleDouble files (name starts with ._)
#    - Ignore other hidden dotfiles
#    - Accept .fasta .fa .fna (case-insensitive)
FASTA_LIST="${OUTDIR}/fasta_files.list"
: > "$FASTA_LIST"
find "$DIR" -type f \( -iname "*.fasta" -o -iname "*.fa" -o -iname "*.fna" \) \
  ! -name "._*" ! -name ".*" -print0 >> "$FASTA_LIST"

# Count files in the list (NUL-based count)
N=$(tr -cd '\0' < "$FASTA_LIST" | wc -c | awk '{print $1}')
if [ "$N" -eq 0 ]; then
  echo "No FASTA files found (after filtering) in: $DIR"
  exit 1
fi
echo "Found $N FASTA files (filtered)…"

# 2) Run seqkit stats one file at a time and merge headers safely
STATS_ALL="${OUTDIR}/assembly_stats_all.tsv"
: > "$STATS_ALL"
first=1
# Read NUL-delimited filenames safely (works on macOS & Linux)
while IFS= read -r -d '' f; do
  if [ $first -eq 1 ]; then
    # include header for the first file
    seqkit stats -Ta "$f" >> "$STATS_ALL"
    first=0
  else
    # skip header for subsequent files
    seqkit stats -Ta "$f" | awk 'NR>1' >> "$STATS_ALL"
  fi
done < "$FASTA_LIST"

# Sanity check
if [ ! -s "$STATS_ALL" ]; then
  echo "ERROR: seqkit produced no output. Check file formats/paths."
  exit 1
fi

# 3) Add Group column based on basename (SGI_* vs PiBac_*)
STATS_WITH_GROUP="${OUTDIR}/assembly_stats_with_group.tsv"
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1 {print $0,"Group"; next}
{
  n=$1; sub(/^.*\//,"",n)   # basename
  g="Other"
  if (n ~ /^SGI_/)      g="SGI"
  else if (n ~ /^PiBac_/) g="PiBac"
  print $0,g
}' "$STATS_ALL" > "$STATS_WITH_GROUP"

# 4) Split per group (keep header)
STATS_SGI="${OUTDIR}/assembly_stats_SGI.tsv"
STATS_PIBAC="${OUTDIR}/assembly_stats_PiBac.tsv"
awk -F'\t' 'NR==1 || $NF=="SGI"'   "$STATS_WITH_GROUP" > "$STATS_SGI"
awk -F'\t' 'NR==1 || $NF=="PiBac"' "$STATS_WITH_GROUP" > "$STATS_PIBAC"

# 5) Quick summary log
LOG="${OUTDIR}/summary.log"
{
  echo "Input dir: $DIR"
  echo "FASTA files (filtered): $N"
  echo "Combined stats: $STATS_WITH_GROUP"
  echo "SGI stats:      $STATS_SGI"
  echo "PiBac stats:    $STATS_PIBAC"
  echo -n "SGI count:   "; awk 'NR>1{c++} END{print c+0}' "$STATS_SGI"
  echo -n "PiBac count: "; awk 'NR>1{c++} END{print c+0}' "$STATS_PIBAC"
} > "$LOG"

echo "Done."
echo "Combined: $STATS_WITH_GROUP"
echo "  SGI   : $STATS_SGI"
echo "  PiBac : $STATS_PIBAC"
echo "Summary: $LOG"

