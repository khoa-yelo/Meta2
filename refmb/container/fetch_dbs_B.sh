#!/usr/bin/env bash
# Fetch the MetaPhlAn 3 marker database used by curatedMetagenomicData 3 (mpa_v30_CHOCOPhlAn_201901, 0.4 GB download,
# about 3 GB with the Bowtie2 index) and build the index. Usage: fetch_dbs_B.sh DIR [threads]
# The expected md5 of the archive is pinned below (taken from the copy the baseline was built with), so the check does not
# depend on the mirror's own .md5 file. MPA_DB_URL overrides the mirror; the cibio host answered http only when this
# script was written (https could not be verified from the build host), so http is the default.
set -euo pipefail
D=${1:?target directory}; T=${2:-8}; U=${MPA_DB_URL:-http://cmprod1.cibio.unitn.it/biobakery3/metaphlan_databases}; I=mpa_v30_CHOCOPhlAn_201901
EXPECTED_MD5=eacca5bc501ef36a3314db3ba561d09a   # mpa_v30_CHOCOPhlAn_201901.tar
mkdir -p "$D"; cd "$D"
[ -s $I.tar ] || curl -fsSL -O $U/$I.tar
[ "$(md5sum $I.tar | cut -d' ' -f1)" = "$EXPECTED_MD5" ] || { echo "checksum mismatch for $I.tar (expected $EXPECTED_MD5)"; exit 1; }
[ -s $I.pkl ] || tar xf $I.tar
if [ ! -s $I.1.bt2 ]; then [ -s $I.fna ] || bunzip2 -k $I.fna.bz2; bowtie2-build --threads "$T" $I.fna $I > bowtie2-build.log 2>&1; rm -f $I.fna; fi
ls -la "$D"; echo "marker database ready in $D"
