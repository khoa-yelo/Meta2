#!/usr/bin/env bash
# refmb tier A — verify checksums where provided and unpack the fetched databases in place. Usage: unpack_dbs.sh <dbdir>
set -euo pipefail
cd "${1:?dbdir}"
for f in *.md5; do base=${f%.md5}; [ -f "$base" ] || continue; exp=$(awk '{print $1}' "$f" | head -1); got=$(md5sum "$base" | cut -d' ' -f1); [ "$exp" = "$got" ] && echo "[md5 ok] $base" || { echo "[md5 FAIL] $base"; exit 1; }; done
if [ -f md5_checksums ] && [ -f Pfam-A.hmm.gz ]; then grep "Pfam-A.hmm.gz" md5_checksums | awk '{print $1"  Pfam-A.hmm.gz"}' | md5sum -c - ; fi
for g in uniref90_v2019_11_diamond-v0.9.25.dmnd.gz db_kofam.hmm.h3f.gz db_kofam.hmm.h3i.gz db_kofam.hmm.h3m.gz db_kofam.hmm.h3p.gz graphs.pkl.gz graphs-20200805.pkl.gz all_pathways-20200805.txt.gz all_pathways_class.txt.gz all_pathways_names.txt.gz list_pathways.txt.gz all_pathways.txt.gz; do
  [ -f "$g" ] && [ ! -f "${g%.gz}" ] && { echo "[gunzip] $g"; gunzip -k "$g"; }
done
[ -f db_kofam.hmm ] || touch db_kofam.hmm     # hmmsearch needs the base name to exist next to the pressed .h3* files
if [ -f interproscan-5.36-75.0-64-bit.tar.gz ] && [ ! -d interproscan-5.36-75.0 ]; then echo "[untar] interproscan"; tar -xzf interproscan-5.36-75.0-64-bit.tar.gz; fi
if [ -d interproscan-5.36-75.0 ] && [ -f interproscan-5.36-75.0/initial_setup.py ] && [ ! -f interproscan-5.36-75.0/.setup_done ]; then ( cd interproscan-5.36-75.0 && python3 initial_setup.py && touch .setup_done ); fi
echo "[unpack] done: $(du -sh . | cut -f1)"; ls -la
