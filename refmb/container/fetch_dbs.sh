#!/usr/bin/env bash
# refmb tier A — fetch the pinned MGnify pipeline v5.0 reference databases and tool data. Run once; resumable; 4 parallel streams.
# Usage: fetch_dbs.sh <dest_dir>
set -uo pipefail
DEST=${1:?dest dir}; mkdir -p "$DEST"; cd "$DEST"
EBI=https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipeline-5.0/ref-dbs
cat > urls.txt <<LIST
$EBI/uniref90_v2019_11_diamond-v0.9.25.dmnd.gz
$EBI/uniref90_v2019_11_diamond-v0.9.25.dmnd.gz.md5
$EBI/db_uniref90_result.txt.gz
$EBI/db_uniref90_result.txt.gz.md5
$EBI/db_uniref90_v2019_11.txt.gz
$EBI/db_uniref90_v2019_11.txt.gz.md5
$EBI/eggnog.db
$EBI/eggnog_proteins.dmnd
$EBI/db_kofam.hmm.h3f.gz
$EBI/db_kofam.hmm.h3i.gz
$EBI/db_kofam.hmm.h3m.gz
$EBI/db_kofam.hmm.h3p.gz
$EBI/kofam_ko_desc.tsv
$EBI/graphs.pkl.gz
$EBI/graphs-20200805.pkl.gz
$EBI/all_pathways-20200805.txt.gz
$EBI/all_pathways_class.txt.gz
$EBI/all_pathways_names.txt.gz
$EBI/list_pathways.txt.gz
$EBI/all_pathways.txt.gz
https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.36-75.0/interproscan-5.36-75.0-64-bit.tar.gz
https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.36-75.0/interproscan-5.36-75.0-64-bit.tar.gz.md5
https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/Pfam-A.hmm.gz
https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam32.0/md5_checksums
LIST
fetch_one() { f=$(basename "$1"); curl -sSL -C - --retry 8 --retry-delay 10 -o "$f" "$1" && echo "[ok] $f $(stat -c %s "$f")" || echo "[FAIL] $f"; }
export -f fetch_one
xargs -P 4 -I{} bash -c 'fetch_one "$@"' _ {} < urls.txt
echo "[fetch] done: $(du -sh . | cut -f1)"
