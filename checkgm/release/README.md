# Release artifacts

The baseline tarballs and the container image are not kept in git; they are built by the release step in the
project workspace and deposited on Zenodo. `zenodo.json` records the metadata and the exact file list of that
deposit, and `RELEASE.json` (alongside the artifacts) records their SHA-256 sums and the round-trip check.

To deposit: see `scripts/zenodo_deposit.py`.
