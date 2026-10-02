#!/usr/bin/env bash
# refmb tier A — raw reads (or contigs) -> MGnify-pipeline-v5-equivalent annotations -> refmb query directory.
# Pinned to the reference's measurement (docs.mgnify.org v5.0 + EBI-Metagenomics/pipeline-v5 job values, see container/pv5_ref):
#   metaSPAdes 3.13 | contigs >= 500 bp | Prodigal 2.6.3 (-p meta) + FragGeneScan 1.31 (illumina_5), Prodigal priority
#   DIAMOND 0.9.25 blastp --max-target-seqs 1 vs UniRef90 2019_11, taxonomy via db_uniref90_result
#   eggNOG-mapper 2.0.0 (-m diamond, --no_annot then --annotate_hits_table) with MGnify's eggnog.db / eggnog_proteins.dmnd
#   InterProScan 5.36-75.0, Pfam application only | hmmsearch 3.2.1 --cut_ga vs KOfam (KEGG 90.0) | MGnify KEGG module completeness
# Usage:
#   run_tierA.sh --sample S --out OUTDIR --dbs DBDIR (--r1 R1.fq.gz [--r2 R2.fq.gz] | --contigs contigs.fasta) [--threads 16] [--mem-gb 120] [--skip-ips]
set -euo pipefail
THREADS=16; MEMGB=120; R1=""; R2=""; CONTIGS=""; SKIP_IPS=0; QC=auto; ASM=auto; EC=auto; SCRATCH=""   # --scratch: node-local dir for assembly and tool temp files (NFS made steps 5-10x slower)
# Default profile (validated on anchors 2026-09-27): paired-end -> metaSPAdes 3.15.3 WITH read error correction; single-end -> MEGAHIT 1.2.9;
# no read trimming / host removal in either case (both reproduce EBI's native assemblies to within ~1%; QC removed 5-13% of assembly on paired data).
# Explicit --qc / --no-qc / --error-correction / --only-assembler / --assembler override the profile.   # auto: paired -> metaSPAdes, single-end -> MEGAHIT (MGnify's rule; contig names show it)
while [[ $# -gt 0 ]]; do case $1 in
  --sample) SAMPLE=$2; shift 2;; --out) OUT=$2; shift 2;; --dbs) DBS=$2; shift 2;; --r1) R1=$2; shift 2;; --r2) R2=$2; shift 2;;
  --contigs) CONTIGS=$2; shift 2;; --threads) THREADS=$2; shift 2;; --mem-gb) MEMGB=$2; shift 2;; --skip-ips) SKIP_IPS=1; shift;; --qc) QC=1; shift;; --no-qc) QC=0; shift;; --assembler) ASM=$2; shift 2;; --error-correction) EC=1; shift;; --only-assembler) EC=0; shift;; --scratch) SCRATCH=$2; shift 2;;
  *) echo "unknown arg $1"; exit 2;; esac; done
if [ "$ASM" = auto ] && [ -z "$CONTIGS" ]; then [ -n "$R2" ] && ASM=metaspades || ASM=megahit; fi
[ "$QC" = auto ] && QC=0    # 2026-09-27: MEGAHIT without QC matched native single-end assemblies exactly (8.8/8.8 Mb, GE 24.5/24.5); QC is off by default for both layouts
[ "$EC" = auto ] && { [ "$ASM" = metaspades ] && EC=1 || EC=0; }
HERE=$(cd "$(dirname "$0")" && pwd); PY=${REFMB_PY:-python3}; EMPY=${EMAPPER_PY:-python2}   # eggNOG-mapper 2.0.0 as vendored by pipeline v5 is Python 2
mkdir -p "$OUT"; cd "$OUT"; LOG="$OUT/tierA.log"; step() { echo "[$(date -u +%FT%TZ)] $*" | tee -a "$LOG"; }
# a database file missing from a node-local cache falls back to the shared copy (a job's cache list can be older than its workflow snapshot)
pick() { if [ -e "$DBS/$1" ]; then echo "$DBS/$1"; else echo "${DBS_FALLBACK:-$DBS}/$1"; fi; }
ASMDIR="$OUT/asm"; TMPW="$OUT"; if [ -n "$SCRATCH" ]; then mkdir -p "$SCRATCH"; ASMDIR="$SCRATCH/asm"; TMPW="$SCRATCH"; export TMPDIR="$SCRATCH"; fi
{ echo "tool versions: diamond $(diamond version 2>&1 | grep -oE '[0-9.]+' | head -1); hmmer $(hmmsearch -h 2>&1 | grep -oE 'HMMER [0-9.]+' | head -1); prodigal $(prodigal -v 2>&1 | grep -oE 'V[0-9.]+' | head -1); megahit $(megahit --version 2>&1 | grep -oE 'v[0-9.]+' | head -1); spades $(${SPADES_BIN:-$(dirname "$(command -v spades.py || echo /usr/bin/spades.py)")}/spades.py --version 2>&1 | grep -oE 'v[0-9.]+' | head -1)"; } >> "$LOG" 2>&1
STAMP() { touch "$OUT/.done_$1"; }; DONE() { [ -f "$OUT/.done_$1" ]; }

# 0. optional read pre-processing as MGnify's assembly pipeline did before metaSPAdes (docs: Trimmomatic 0.36 + host removal vs hg38) ------
#    Trimmomatic step values follow MGnify raw-read defaults except MINLEN:50 (assumption; see below); host removal with bwa-mem2 on
#    the hg38 index EBI ships in pipeline-5.0/ref-dbs/hg38 (bwa-mem2 index files). Unmapped read pairs are kept.
if [ "$QC" = 1 ] && [ -n "$R1" ] && ! DONE qc; then step "read QC: Trimmomatic 0.36 + bwa-mem2 hg38 host removal"
  TRIM="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"   # MINLEN:100 (raw-read analysis default) dropped 100% of a 75 bp library; 50 is the assembly-QC choice, recorded as an assumption
  if [ -n "$R2" ]; then trimmomatic PE -threads "$THREADS" -phred33 "$R1" "$R2" qc_1.fq.gz qc_1u.fq.gz qc_2.fq.gz qc_2u.fq.gz $TRIM >> "$LOG" 2>&1
    bwa-mem2 mem -t "$THREADS" "$DBS/hg38/hg38.fa" qc_1.fq.gz qc_2.fq.gz 2>> "$LOG" | samtools fastq -f 12 -1 clean_1.fq.gz -2 clean_2.fq.gz -0 /dev/null -s /dev/null -n - >> "$LOG" 2>&1
    R1="$OUT/clean_1.fq.gz"; R2="$OUT/clean_2.fq.gz"
  else trimmomatic SE -threads "$THREADS" -phred33 "$R1" qc_1.fq.gz $TRIM >> "$LOG" 2>&1
    bwa-mem2 mem -t "$THREADS" "$DBS/hg38/hg38.fa" qc_1.fq.gz 2>> "$LOG" | samtools fastq -f 4 -0 clean_1.fq.gz -n - >> "$LOG" 2>&1
    R1="$OUT/clean_1.fq.gz"; fi
  rm -f qc_1.fq.gz qc_2.fq.gz qc_1u.fq.gz qc_2u.fq.gz; STAMP qc; step "read QC done"
elif [ "$QC" = 1 ] && DONE qc && [ -z "$CONTIGS" ]; then R1="$OUT/clean_1.fq.gz"; [ -f "$OUT/clean_2.fq.gz" ] && R2="$OUT/clean_2.fq.gz"; fi
# 1. assembly ------------------------------------------------------------------------------------------------------
if [ -n "$CONTIGS" ]; then [ "$(readlink -f "$CONTIGS")" = "$(readlink -f contigs_raw.fasta)" ] || cp "$CONTIGS" contigs_raw.fasta; DONE asm || { STAMP asm; step "assembly: user-supplied contigs"; };
elif ! DONE asm; then step "assembly start"
  # --only-assembler: MGnify's assembly pipeline skips SPAdes read error correction (miassembler default spades_only_assembler=true)
  if [ "$ASM" = metaspades ]; then
    SPADES_BIN=${SPADES_BIN:-$(dirname "$(command -v metaspades.py || echo /usr/bin/metaspades.py)")}   # metaSPAdes 3.15.3: the reference pool was assembled with metaSPAdes 3.12-3.15.3 (mostly 3.14.1/3.15.3)
    [ -x "$SPADES_BIN/metaspades.py" ] || SPADES_BIN=$(dirname "$(command -v metaspades.py)")
    OA=$( [ "$EC" = 1 ] && echo "" || echo "--only-assembler" )
    if [ -n "$R2" ]; then "$SPADES_BIN/metaspades.py" $OA -1 "$R1" -2 "$R2" -o "$ASMDIR" -t "$THREADS" -m "$MEMGB" >> "$LOG" 2>&1
    else "$SPADES_BIN/spades.py" $OA -s "$R1" -o "$ASMDIR" -t "$THREADS" -m "$MEMGB" >> "$LOG" 2>&1; fi
    "$SPADES_BIN/spades.py" --version 2>&1 | head -1 >> "$LOG"
    cp "$ASMDIR/contigs.fasta" contigs_raw.fasta; mkdir -p "$OUT/asm"; cp "$ASMDIR/spades.log" "$ASMDIR/params.txt" "$OUT/asm/" 2>/dev/null || true
  else  # MEGAHIT (MGnify used it for single-end libraries: contig names carry 4-decimal `multi` coverage); headers renamed to NODE_i_length_L_cov_C
    if [ -n "$R2" ]; then megahit -1 "$R1" -2 "$R2" -o "$ASMDIR" -t "$THREADS" --min-contig-len 200 >> "$LOG" 2>&1
    else megahit -r "$R1" -o "$ASMDIR" -t "$THREADS" --min-contig-len 200 >> "$LOG" 2>&1; fi
    $PY "$HERE/pv5_to_query.py" rename-megahit --in "$ASMDIR/final.contigs.fa" --out contigs_raw.fasta >> "$LOG" 2>&1
  fi
  echo "$ASM ec=$EC qc=$QC" > assembler.txt; STAMP asm; step "assembly done ($ASM, error correction $EC, qc $QC)"; fi
# 2. contig filter (>= 500 bp, MGnify contig_min_length) --------------------------------------------------------------
DONE filt || { $PY "$HERE/pv5_to_query.py" filter-contigs --in contigs_raw.fasta --out contigs.fasta --min-length 500 >> "$LOG" 2>&1; STAMP filt; }
# 3. gene calling: Prodigal + FragGeneScan, merged with Prodigal priority --------------------------------------------
if ! DONE cds; then step "gene calling"
  prodigal -p meta -f sco -i contigs.fasta -o prodigal.sco -d prodigal.ffn -a prodigal.faa >> "$LOG" 2>&1
  FGS_DIR=$(dirname "$(command -v FragGeneScan)"); TRAIN_DIR="$FGS_DIR/train"; [ -d "$TRAIN_DIR" ] || TRAIN_DIR="$(dirname "$FGS_DIR")/share/fraggenescan*/train"
  ( cd "$(dirname "$(command -v FragGeneScan)")" >/dev/null; FragGeneScan -s "$OUT/contigs.fasta" -o "$OUT/fgs" -w 1 -t illumina_5 -p "$THREADS" ) >> "$LOG" 2>&1 || step "WARN FragGeneScan failed; continuing with Prodigal only"
  $PY "$HERE/pv5_to_query.py" merge-cds --prodigal-faa prodigal.faa --fgs-faa fgs.faa --out CDS.faa >> "$LOG" 2>&1; STAMP cds; step "CDS: $(grep -c '^>' CDS.faa)"; fi
# 4. taxonomy: DIAMOND blastp vs UniRef90, best hit, joined to the UniRef90 taxonomy table -----------------------------
if ! DONE dmnd; then step "DIAMOND blastp"
  diamond blastp --db "$DBS/uniref90_v2019_11_diamond-v0.9.25.dmnd" --query CDS.faa --max-target-seqs 1 --threads "$THREADS" \
    --outfmt 6 sseqid qseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --out diamond_raw.tsv >> "$LOG" 2>&1
  $PY "$HERE/pv5_to_query.py" join-taxonomy --diamond diamond_raw.tsv --uniref-tax "$(pick db_uniref90_v2019_11.txt.gz)" --out diamond.tsv >> "$LOG" 2>&1; STAMP dmnd; fi
# 6. Pfam: InterProScan 5.36-75.0, Pfam application only -------------------------------------------------------------
if [ "$SKIP_IPS" = 0 ] && ! DONE ips; then step "InterProScan (Pfam)"
  IPS="$DBS/interproscan-5.36-75.0/interproscan.sh"; sed 's/\*/X/g' CDS.faa > CDS.noast.faa
  [ -n "${JAVA8_BIN:-}" ] && export PATH="$JAVA8_BIN:$PATH"    # InterProScan 5.36 requires Java 1.8
  "$IPS" -i CDS.noast.faa -appl Pfam -f TSV -dp -o ips.tsv --cpu "$THREADS" -T "$TMPW/ips_tmp" >> "$LOG" 2>&1 || "$IPS" -i CDS.noast.faa -appl PfamA -f TSV -dp -o ips.tsv --cpu "$THREADS" -T "$TMPW/ips_tmp" >> "$LOG" 2>&1
  rm -rf ips_tmp "$TMPW/ips_tmp" CDS.noast.faa; STAMP ips; fi
# 7. KOfam: hmmsearch --cut_ga -------------------------------------------------------------------------------------
if ! DONE kofam; then step "hmmsearch KOfam"
  hmmsearch --cut_ga --noali --cpu "$THREADS" --domtblout kofam.domtbl -o /dev/null "$DBS/db_kofam.hmm" CDS.faa >> "$LOG" 2>&1; STAMP kofam; fi
# 8. KEGG module completeness (MGnify give_pathways on the per-contig KO union from KOfam) -----------------------------
if ! DONE kegg; then step "KEGG module completeness"
  $PY "$HERE/pv5_to_query.py" kofam-union --domtbl kofam.domtbl --out kofam_union.tsv --best-hits kofam_best.tsv >> "$LOG" 2>&1
  $PY "$HERE/pv5_ref/tools_Assembly_KEGG_analysis_KEGG_pathways_give_pathways.py" -i kofam_union.tsv -g "$(pick graphs-20200805.pkl)" -c "$DBS/all_pathways_class.txt" -n "$DBS/all_pathways_names.txt" -o kegg --summary-only >> "$LOG" 2>&1
  echo graphs-20200805.pkl > .kegg_graphs; STAMP kegg; fi   # graphs-20200805.pkl reproduces MGnify's module completeness exactly (146/146 on identical KO input); graphs.pkl only 139/153
# 8b. KO: eggNOG-mapper 2.0.0, diamond mode, two-step as in pipeline v5 (last: its 37 GB sqlite db is the slowest to provision) -----------------------------------------------
if ! DONE emapper; then step "eggNOG-mapper"
  $EMPY "$HERE/emapper2/emapper.py" -i CDS.faa -m diamond --no_annot --no_file_comments --cpu "$THREADS" --data_dir "$DBS" --dmnd_db "$DBS/eggnog_proteins.dmnd" -o emapper --override >> "$LOG" 2>&1
  $EMPY "$HERE/emapper2/emapper.py" --annotate_hits_table emapper.emapper.seed_orthologs --no_file_comments --cpu "$THREADS" --data_dir "$DBS" -o emapper --override >> "$LOG" 2>&1; STAMP emapper; fi
# 9. query directory -------------------------------------------------------------------------------------------------
$PY "$HERE/pv5_to_query.py" build --sample "$SAMPLE" --contigs contigs.fasta --cds CDS.faa --diamond diamond.tsv --emapper emapper.emapper.annotations \
   $( [ -f ips.tsv ] && echo --ips ips.tsv ) --kofam-best kofam_best.tsv --kegg kegg.summary.kegg_pathways.tsv --out "$OUT/query" >> "$LOG" 2>&1
step "done -> $OUT/query"
