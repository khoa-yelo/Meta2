#!/usr/bin/env bash
# Fetch the MetaPhlAn 3 marker database used by curatedMetagenomicData 3 (mpa_v30_CHOCOPhlAn_201901, 0.4 GB download,
# about 3 GB with the Bowtie2 index) and build the index. Usage: fetch_dbs_B.sh DIR [threads]
set -euo pipefail
D=${1:?target directory}; T=${2:-8}; U=http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases; I=mpa_v30_CHOCOPhlAn_201901
mkdir -p "$D"; cd "$D"
[ -s $I.tar ] || curl -fsSL -O $U/$I.tar; curl -fsSL -O $U/$I.md5
[ "$(md5sum $I.tar | cut -d' ' -f1)" = "$(cut -d' ' -f1 $I.md5)" ] || { echo "checksum mismatch for $I.tar"; exit 1; }
[ -s $I.pkl ] || tar xf $I.tar
if [ ! -s $I.1.bt2 ]; then [ -s $I.fna ] || bunzip2 -k $I.fna.bz2; bowtie2-build --threads "$T" $I.fna $I > bowtie2-build.log 2>&1; rm -f $I.fna; fi
ls -la "$D"; echo "marker database ready in $D"
