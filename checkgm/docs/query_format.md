# Query directory format (assembly-based pipeline)

A *query* is one directory per sample holding the output of an assembly pipeline in four small tab-separated tables.
It is the input of `checkgm normalize`, which turns it into normalized tables that `checkgm score` compares with an
assembly-based bundle. The shipped container (`container/checkgm_assembly.def`, run script `run_assembly.sh`) writes this
directory itself; the format is documented so that any other assembly pipeline can produce it. Tab-separated, header row
required, UTF-8. Identifiers are the user's own; nothing depends on MGnify naming.

| file | required | columns | notes |
|---|---|---|---|
| `sample.json` | no | `sample_id`, `body_site`, `pipeline` {name, version, databases}, optional `age_years`, `country`, `condition` | `body_site` must equal the bundle's (Gut); the scorer refuses others with `NO_MATCHING_REFERENCE`. The container's `pv5_to_query.py` fills `pipeline.version` with the checkgm version only when the checkGM package is importable by the Python it runs under; inside the assembly image (Python 3.8) it is not, and the field reads `unknown` |
| `contigs.tsv` | yes | `contig_id`, `length` (bp), `coverage` | ALL contigs of the assembly (used for total assembled length and so for the quality band). `coverage` is the per-contig read depth (metaSPAdes k-mer coverage, MEGAHIT `multi`, or read-mapping depth). Empty coverage: the CDS on that contig count with weight 0 |
| `cds.tsv` | yes | `cds_id`, `contig_id`, `ko`, `pfam`, optional `ko_kofam` | one row per predicted CDS; `ko` and `pfam` are `;`- or `,`-separated lists (may be empty). `ko` is the KO source the reference uses (eggNOG-mapper); `ko_kofam` is a single best KOfam KO per CDS if available |
| `cds_taxonomy.tsv` | yes | `cds_id` or `contig_id`, `taxid` | one NCBI taxid per CDS (best hit or lowest common ancestor of the protein search). Any rank works: the normalizer walks the lineage. If `contig_id` is given it is used directly |
| `modules.tsv` | no | `module_id`, `completeness` (0-1) | KEGG module completeness (for example from MGnify's kegg-pathways-completeness tool). Without it the module layer is `not_assessable` |

## What `checkgm normalize` does with it

The rules come from the bundle (`normalization/rules.json`), so a query is normalized exactly as the reference pool was.

1. **Taxonomy.** Each CDS taxid is placed on the NCBI taxonomy shipped in the bundle. Per contig and rank, a majority vote
   (at least 50 % of the resolvable CDS calls) assigns the contig; contig weight is length x coverage. Contigs with no
   majority go to an `unmapped` bucket. Proportions are taken among mapped weight and transformed to centred log-ratios
   (CLR) over the bundle's feature basis (`clr_basis_family.txt`, 584 families; genus 1,687) with multiplicative replacement
   of zeros.
2. **Gene layers.** Per KO and Pfam family: CDS count and coverage sum. Genome equivalents = median coverage sum over the
   52 single-copy marker KOs in `marker_panel.tsv`. The value scored is copies per genome; at scoring, values under 0.05
   copies per genome count as not detected.
3. **Quality.** Total assembled length over all contigs gives the band low / medium / high (5-90 / 90-180 / >= 180 Mb).
   Floors: 5 Mb assembled, N50 >= 1 kb over CDS-bearing contigs, >= 5,000 CDS, >= 20 marker KOs detected, genome
   equivalents >= 10, family mapped fraction >= 0.5. A sample failing a floor is not scored and carries the reason code in
   `rejections.tsv`.

## Does my pipeline need a calibration?

- If the tables come from the shipped container, no: it pins the tools and databases the reference was measured with, and
  its agreement with the original analyses was checked on 61 public samples (`docs/validation.md`, section 1). Scores
  carry `calibration_applied = none`.
- If the tables come from another assembly pipeline, its measurement space may differ (different gene caller, annotation
  database or assembler). The `checkgm.calibrate` module (library only, there is no command for it) can fit a calibration
  from a set of samples measured both by the user's pipeline and by the reference pipeline, deciding layer by layer what
  is scorable. `checkgm score --calibration DIR` then applies it and records its id in every score. No calibration is
  shipped in this release. Without one, scores from a different pipeline are not comparable with the reference and should
  not be reported as such.

## Minimal example

```
mysample/
  sample.json          {"sample_id": "S1", "body_site": "Gut", "pipeline": {"name": "megahit+prodigal+eggnog-mapper", "version": "..."}}
  contigs.tsv          contig_id  length  coverage
  cds.tsv              cds_id  contig_id  ko  pfam
  cds_taxonomy.tsv     cds_id  taxid
  modules.tsv          module_id  completeness
```

```bash
checkgm normalize --bundle gut-assembly-adult-global-v0.8-lenient --query mysample --out mysample_norm
checkgm score     --bundle gut-assembly-adult-global-v0.8-lenient --normalized mysample_norm --out mysample_report
```

This is the first version of the format (it has not changed since the container was validated); a later change would be
announced in `CHANGELOG.md`.
