# refmb query format (v0.1)

A query is one directory per sample. It is what any assembly pipeline has to produce for `refmb score`. Tab-separated,
header row required, UTF-8. Identifiers are the user's; nothing here depends on MGnify naming.

| file | required | columns | notes |
|---|---|---|---|
| `sample.json` | no | `sample_id`, `body_site`, `pipeline` {name, version, databases}, optional `age_years`, `country`, `condition` | `body_site` must equal the bundle's (Gut); the scorer refuses others with NO_MATCHING_REFERENCE |
| `contigs.tsv` | yes | `contig_id`, `length` (bp), `coverage` | ALL contigs of the assembly (used for total assembled length → quality band). `coverage` = per-contig read depth (metaSPAdes k-mer coverage, MEGAHIT `multi`, or read-mapping depth). Empty coverage → CDS on that contig count with weight 0 |
| `cds.tsv` | yes | `cds_id`, `contig_id`, `ko`, `pfam`, optional `ko_kofam` | one row per predicted CDS; `ko` and `pfam` are `;`- or `,`-separated lists (may be empty). `ko` is the KO source the reference uses (eggNOG-mapper); `ko_kofam` is a single best KOfam KO per CDS if available |
| `cds_taxonomy.tsv` | yes | `cds_id` or `contig_id`, `taxid` | one NCBI taxid per CDS (best-hit or LCA of the protein search). Any rank is fine: the normalizer walks the lineage. If `contig_id` is given it is used directly |
| `modules.tsv` | no | `module_id`, `completeness` (0–1) | KEGG module completeness (e.g. MGnify's kegg-pathways-completeness tool). Without it the module layer is not_assessable |

## What the normalizer does with it (bundle/normalization/rules.json)

1. **Taxonomy.** CDS taxids → lineage on the NCBI backbone shipped in the bundle. Per contig and rank, majority vote (≥ 50 % of
   resolvable CDS calls). Contig weight = length × coverage. Contigs with no majority go to `unmapped`. Proportions among mapped
   weight; CLR over the bundle's basis (`clr_basis_family.txt`, 584 families; genus 1,690) with multiplicative replacement.
2. **Gene layers.** Per KO / Pfam: CDS count and coverage sum. Genome equivalents = median coverage sum over the 52-KO validated
   marker panel. Value = copies per genome. Detection floor 0.05 copies per genome at scoring.
3. **Quality.** Total assembled length (all contigs) → band low / medium / high (5–90 / 90–180 / ≥ 180 Mb). Floors: 5 Mb, N50 ≥ 1 kb
   over CDS-bearing contigs, ≥ 5,000 CDS, ≥ 20 marker KOs detected, genome equivalents ≥ 10, family mapped fraction ≥ 0.5.
   A sample failing a floor is not scored and carries the reason code.

## Contracts
- **Tier A (calibration = none).** The tables come from a pipeline pinned to the reference's tools and databases (the `refmb run`
  container, validated on the anchor panel). Full reference.
- **Tier B (calibration = <kit>).** Any pipeline. The user first runs the anchor panel through their pipeline; `refmb calibrate`
  fits the mapping and decides, layer by layer, what is scorable. Scores carry the calibration id.

## Minimal example
```
mysample/
  sample.json          {"sample_id": "S1", "body_site": "Gut", "pipeline": {"name": "megahit+prodigal+eggnog-mapper", "version": "..."}}
  contigs.tsv          contig_id  length  coverage
  cds.tsv              cds_id  contig_id  ko  pfam
  cds_taxonomy.tsv     cds_id  taxid
  modules.tsv          module_id  completeness
```
