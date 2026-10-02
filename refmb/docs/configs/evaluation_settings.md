# Settings of the cross-study case/control evaluation

Transcribed from the header and constants of the evaluation script in the project workspace (`s8_evaluate.py`), which
produced `docs/reports/assembly_case_control_evaluation.md` and the per-study tables.

- **Folds.** Every study's samples are scored against a baseline rebuilt **without** that study, so no sample ever meets a
  reference containing its own study. The inclusion rule differs between the two pipelines and is stated here as applied:
  - assembly-based pipeline: studies with at least 20 cases and at least 10 controls (`splits_rules.yaml`, symptom-only
    labels excluded); 11 studies, of which PRJEB62821 (hemodialysis) has 24 cases and 16 controls;
  - read-based pipeline: studies with at least 20 cases and at least 20 controls; 14 studies.
  One of the 14 read-pipeline studies, `HMP_2019_ibdmdb` (inflammatory bowel disease; `HMP_2019` in the report tables), is a largely adolescent cohort:
  the median age of both its cases and its controls is 16.5 years (`docs/reports/reads_disease_by_stratum.md`), although
  the baseline covers adults (18-65). It is kept in the read-pipeline evaluation and in the disease atlas; its AUROC is
  reported per study so that a reader can judge the comparison without it.
- **Task.** Case (diseased) versus control (healthy, control, or curatedMetagenomicData control), adults in the stratum
  only. Train on all other folds (cases and controls), test on the held-out study. Reported: AUROC per study, the mean
  over studies, and per condition where at least two studies carry it. The within-study 5-fold cross-validated AUROC is
  reported as a ceiling.
- **Classifier.** Random forest, 500 trees, seed 20261020, the same for every feature set.
- **Feature sets.**
  - reference-relative: per-feature percentile (undetected and expected-but-missing set to 0) on the family, genus,
    eggNOG KO and KEGG module layers, plus per-layer counts of deviating features;
  - raw abundances: family and genus CLR, log10(copies per genome + 0.001) for eggNOG KOs, module completeness;
  - alpha diversity: Shannon index and richness on family proportions, Shannon index on KO copies;
  - health index: a genus-level approximation of the Gut Microbiome Health Index (Gupta et al. 2020) from family and
    genus proportions, with 3 health-prevalent genera (Alistipes, Bifidobacterium, Sutterella) and 18 health-scarce genera
    (Anaerotruncus, Atopobium, Blautia, Clostridium, Eggerthella, Flavonifractor, Fusobacterium, Gemella, Granulicatella,
    Holdemania, Klebsiella, Lactobacillus, Peptostreptococcus, Mediterraneibacter, Solobacterium, Streptococcus,
    Subdoligranulum, Veillonella). GMHI is defined on 50 MetaPhlAn species, so this is stated as an approximation.
- **Feature filter.** Features detected in at least 10 % of the evaluation samples are kept for modelling.
- **Paired comparison.** One-sided Wilcoxon signed-rank test over studies that reference-relative AUROC exceeds the
  comparison set, as reported in `assembly_case_control_evaluation.md`; p-values are not corrected for the number of
  comparisons.
