#!/usr/bin/env bash
# Read-based pipeline from raw reads: MetaPhlAn 3 (pinned 3.0.14, marker database mpa_v30_CHOCOPhlAn_201901, default parameters,
# all reads of the sample given as unpaired input, as in curatedMetagenomicData 3), then checkgm import + score against a
# read-based bundle. No read QC: the cMD3 profiles were produced from the deposited reads. Reads processed (both mates counted)
# are counted from the input files and set the depth band. Optional HUMAnN 3 tables add the function layers.
# Usage: run_read.sh --sample ID --r1 R1.fastq.gz [--r2 R2.fastq.gz] --dbs DIR --bundle BUNDLE --out DIR [--threads 8]
#        [--humann-genefamilies FILE] [--humann-pathabundance FILE]
set -euo pipefail
THREADS=8; R2=""; GF=""; PA=""; INDEX=mpa_v30_CHOCOPhlAn_201901
while [ $# -gt 0 ]; do case "$1" in
  --sample) SAMPLE=$2; shift 2;; --r1) R1=$2; shift 2;; --r2) R2=$2; shift 2;; --dbs) DBS=$2; shift 2;; --bundle) BUNDLE=$2; shift 2;;
  --out) OUT=$2; shift 2;; --humann-genefamilies) GF=$2; shift 2;; --humann-pathabundance) PA=$2; shift 2;; --threads) THREADS=$2; shift 2;; *) echo "unknown option $1"; exit 2;; esac; done
: "${SAMPLE:?}" "${R1:?}" "${DBS:?}" "${BUNDLE:?}" "${OUT:?}"
# checkgm: the installed console script (pip install checkgm), or CHECKGM_CLI=/path/to/cli.py run with CHECKGM_PY
if [ -n "${CHECKGM_CLI:-}" ]; then CHECKGM="${CHECKGM_PY:-python3} $CHECKGM_CLI"; else CHECKGM=${CHECKGM_BIN:-checkgm}; fi
mkdir -p "$OUT"; LOG="$OUT/read.log"; step() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG"; }
metaphlan --version | tee -a "$LOG"; metaphlan --version | grep -q "version 3\.0\." || { echo "MetaPhlAn is not 3.0.x; the baseline was measured with MetaPhlAn 3 (mpa_v30)"; exit 3; }
[ -e "$DBS/$INDEX.pkl" ] || { echo "marker database $INDEX not found in $DBS (run container/fetch_dbs_read.sh)"; exit 3; }
PROFILE="$OUT/$SAMPLE.txt"
if [ ! -s "$PROFILE" ] || [ ! -s "$OUT/$SAMPLE.reads.tsv" ]; then step "MetaPhlAn 3"
  rm -f "$OUT/$SAMPLE.bowtie2.bz2"
  metaphlan "$R1${R2:+,$R2}" --input_type fastq --bowtie2db "$DBS" --index "$INDEX" --nproc "$THREADS" --bowtie2out "$OUT/$SAMPLE.bowtie2.bz2" --sample_id "$SAMPLE" -o "$OUT/$SAMPLE.profile.tmp" >> "$LOG" 2>&1
  N=$( { zcat -f "$R1" ${R2:+"$R2"} || true; } | awk 'END {printf "%d", NR / 4}' )
  printf '%s\t%s\n' "$SAMPLE" "$N" > "$OUT/$SAMPLE.reads.tsv"; mv "$OUT/$SAMPLE.profile.tmp" "$PROFILE"; fi
step "import and score"
$CHECKGM run-metaphlan --bundle "$BUNDLE" --input "$PROFILE" --reads-tsv "$OUT/$SAMPLE.reads.tsv" ${GF:+--humann-genefamilies "$GF"} ${PA:+--humann-pathabundance "$PA"} --out "$OUT/score" | tee -a "$LOG"
step "done"
