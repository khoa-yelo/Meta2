#!/usr/bin/env python3
"""S12.1 — out-of-sample scores for every healthy sample and every case/control sample, both pipelines, joined to metadata.

  s12_oos_scores.py A   # pipeline A: pool studies vs their LOSO bundle (refs/...-v0.7-lenient-loso-<study>, built by s5), held-out healthy
                        #             studies vs base (from results/s7/scores_gut, already out-of-sample), case/control studies from work/s8/scores
  s12_oos_scores.py B   # pipeline B: pool studies vs their LOSO bundle, held-out healthy + case/control studies vs base (refmb.score);
                        #             a study that shares specimens with a pool study is scored, in full, against the LOSO bundle of that
                        #             study as well, so no study is ever split across two baselines (see the comment before `excl_for`)
Outputs work/s12/<A|B>/oos_scores.parquet (sample, study, role, layer, feature_id, percentile, call), oos_summaries.parquet (per sample x layer
with fraction outside, excess test, landscape distance) and meta.parquet (sample, study, role, condition, is_case, age, age_bin, sex, bmi,
bmi_class, country, region, westernized)."""
import sys as _s; _s.path.insert(0, "/home/classes/bios/270/khoa/meta2/project"); from refmb.paths import VER as _VER, VTAG as _VTAG, bundle_A as _bundle_A
import os, subprocess, sys, time
import numpy as np, pandas as pd, polars as pl
sys.path.insert(0, "/home/classes/bios/270/khoa/meta2/project")
P = "/home/classes/bios/270/khoa/meta2/project"; PY = "/home/classes/bios/270/khoa/envs/refmb/bin/python"
which = sys.argv[1]; OUT = f"{P}/work/s12/{which}"; os.makedirs(OUT, exist_ok=True); t0 = time.time()
REGION = {"United States": "N. America", "Canada": "N. America", "USA": "N. America", "CAN": "N. America", "China": "E. Asia", "CHN": "E. Asia", "Japan": "E. Asia", "JPN": "E. Asia", "South Korea": "E. Asia", "KOR": "E. Asia",
          "Russia": "Europe", "RUS": "Europe", "Finland": "Europe", "FIN": "Europe", "Denmark": "Europe", "DNK": "Europe", "Estonia": "Europe", "EST": "Europe", "Netherlands": "Europe", "NLD": "Europe", "United Kingdom": "Europe", "GBR": "Europe",
          "Germany": "Europe", "DEU": "Europe", "France": "Europe", "FRA": "Europe", "Spain": "Europe", "ESP": "Europe", "Italy": "Europe", "ITA": "Europe", "Sweden": "Europe", "SWE": "Europe", "Austria": "Europe", "AUT": "Europe", "Luxembourg": "Europe", "LUX": "Europe", "Ireland": "Europe", "IRL": "Europe",
          "Israel": "W. Asia", "ISR": "W. Asia", "India": "S. Asia", "IND": "S. Asia", "Bangladesh": "S. Asia", "BGD": "S. Asia", "Kenya": "Africa", "KEN": "Africa", "Tanzania": "Africa", "TZA": "Africa", "Madagascar": "Africa", "MDG": "Africa", "Ghana": "Africa", "GHA": "Africa", "Cameroon": "Africa", "CMR": "Africa", "Ethiopia": "Africa", "ETH": "Africa",
          "Fiji": "Oceania", "FJI": "Oceania", "Australia": "Oceania", "AUS": "Oceania", "Peru": "S. America", "PER": "S. America", "El Salvador": "C. America", "SLV": "C. America", "Mongolia": "E. Asia", "MNG": "E. Asia", "Kazakhstan": "C. Asia", "KAZ": "C. Asia", "Indonesia": "SE Asia", "IDN": "SE Asia", "Malaysia": "SE Asia", "MYS": "SE Asia", "Vietnam": "SE Asia", "VNM": "SE Asia", "Thailand": "SE Asia", "THA": "SE Asia", "Philippines": "SE Asia", "PHL": "SE Asia"}


def bins(age, bmi):
    ab = pd.cut(age, [17, 39, 64, 120], labels=["18-39", "40-64", "65+"]).astype(str).replace("nan", "")
    bb = pd.cut(bmi, [0, 18.5, 25, 30, 100], labels=["<18.5", "18.5-25", "25-30", ">=30"], right=False).astype(str).replace("nan", "")
    return ab, bb


if which == "A":
    inv = pd.read_parquet(f"{P}/work/s0/inventory.parquet"); sp = pd.read_csv(f"{P}/work/s0/splits.tsv", sep="\t")
    qc = pd.read_parquet(f"{P}/work/s2/qc_per_analysis.parquet", columns=["analysis_id", "s2_status"]); ok = set(qc.loc[qc.s2_status == "retained", "analysis_id"])
    g = inv[inv.is_primary_analysis & (inv.body_site_core == "Gut") & inv.analysis_id.isin(ok)].merge(sp[["study_bioproject", "role"]], on="study_bioproject", how="left")
    cmd = pd.read_csv("/home/classes/bios/270/khoa/meta2/cmd3/cmd3_release_sampleMetadata.csv", dtype=str, keep_default_na=False)[["sample_id", "study_name", "age", "gender", "BMI", "country", "non_westernized"]]
    g = g.merge(cmd, left_on=["cmd3_sample_id", "cmd3_study_name"], right_on=["sample_id", "study_name"], how="left")
    age = pd.to_numeric(g.cur_Age_years, errors="coerce").fillna(pd.to_numeric(g.age, errors="coerce")); age = age.where(age <= 110)
    sex = g.cur_Sex.replace("", np.nan).fillna(g.gender.replace("", np.nan)).str.lower().str[0].map({"f": "female", "m": "male"})
    country = g.cur_Country.replace("", np.nan).fillna(g.country.replace("", np.nan)); bmi = pd.to_numeric(g.BMI, errors="coerce")
    healthy = g.cur_Health_status_group.isin(["Healthy", "Control"]) | (g.cmd3_study_condition == "control")
    cond = g.cur_Condition_group.replace("", np.nan).fillna(g.cmd3_study_condition.replace("", np.nan)).fillna("")
    ab, bb = bins(age, bmi)
    B = _bundle_A("lenient"); pool_st = set(pd.read_csv(f"{B}/provenance/studies.tsv", sep="\t").study_bioproject)
    excl = set(pd.read_csv(f"{B}/provenance/pool_exclusions.tsv", sep="\t").analysis_id)
    role = np.where(g.study_bioproject.isin(pool_st) & ~g.analysis_id.isin(excl) & healthy, "reference_pool", g.role.fillna("other"))
    meta = pd.DataFrame({"sample": g.analysis_id, "study": g.study_bioproject, "role": role, "is_healthy": healthy.to_numpy(), "condition": cond.to_numpy(), "age": age.to_numpy(), "age_bin": ab.to_numpy(),
                         "sex": sex.fillna("").to_numpy(), "bmi": bmi.to_numpy(), "bmi_class": bb.to_numpy(), "country": country.fillna("").to_numpy(),
                         "region": country.map(REGION).fillna(pd.Series(np.where(country.notna(), "other", ""), index=country.index)).astype(str).to_numpy(), "westernized": np.where(g.non_westernized.fillna("") == "yes", "no", np.where(g.non_westernized.fillna("") == "no", "yes", ""))})
    meta.to_parquet(f"{OUT}/meta.parquet", index=False)
    # scores: pool studies -> LOSO bundles (s7_score.py, ids per study); held-out healthy + others -> base scores from results/s7
    parts = []; summs = []
    for st in sorted(pool_st):
        bundle = _bundle_A("lenient", loso=st)
        if not os.path.exists(f"{bundle}/manifest.json"):
            print("missing LOSO bundle", st); continue
        ids = meta.loc[(meta.study == st) & (meta.role == "reference_pool"), "sample"]
        idf = f"{OUT}/ids_{st}.txt"; open(idf, "w").write("\n".join(ids) + "\n")
        subprocess.run([PY, f"{P}/scripts/s7_score.py", "--bundle", bundle, "--ids-file", idf, "--name", f"loso_{st}", "--out", f"{OUT}/parts"], check=True, stdout=subprocess.DEVNULL)
        parts.append(pd.read_parquet(f"{OUT}/parts/loso_{st}.parquet")); summs.append(pd.read_parquet(f"{OUT}/parts/loso_{st}_summaries.parquet")); print("scored", st, f"{time.time()-t0:.0f}s", flush=True)
    base_ids = meta.loc[(meta.role != "reference_pool"), "sample"]
    sc = pl.scan_parquet(f"{P}/results/s7/scores_gut.parquet").filter(pl.col("analysis_id").is_in(list(base_ids))).collect().to_pandas(); parts.append(sc)
    sm = pd.read_parquet(f"{P}/results/s7/summaries_gut.parquet"); summs.append(sm[sm.analysis_id.isin(base_ids)])
    S = pd.concat(parts, ignore_index=True).rename(columns={"analysis_id": "sample"}); M = pd.concat(summs, ignore_index=True).rename(columns={"analysis_id": "sample"})
    S[["sample", "layer", "feature_id", "band", "value", "percentile", "call"]].to_parquet(f"{OUT}/oos_scores.parquet", index=False)
    keep = [c for c in ["sample", "layer", "n_assessed", "n_low_raw", "n_high_raw", "n_missing", "frac_outside_raw", "excess_outside_p", "excess_outside_q", "core_distance_pct", "quality_band"] if c in M.columns]
    M[keep].to_parquet(f"{OUT}/oos_summaries.parquet", index=False); print("A done", len(S), len(M), f"{time.time()-t0:.0f}s")

else:
    from refmb.score import Bundle, score
    D = f"{P}/work/s11/pipeB"; D_ = D; S = pd.read_csv(f"{D}/samples.tsv", sep="\t"); Mm = pd.read_parquet(f"{P}/work/s11/metaphlan/sample_meta.parquet").drop_duplicates("sample_key").set_index("sample_key")
    BVER = os.environ.get("REFMB_B_VER", "0.2"); BUILDER = "s11_build_pipelineB.py" if BVER == "0.1" else "s13_build_pipelineB_v02.py"
    tabs = {f"taxonomy_{r}": pd.read_parquet(f"{D}/taxonomy_{r}_reads.parquet") for r in ["family", "genus", "species"]}; band_of = dict(zip(S.sample_key, S.depth_band))
    if BVER != "0.1":
        tabs.update({l: pd.read_parquet(f"{D}/function_{l}.parquet") for l in ["ko_humann", "pathway_humann"]})
    per = {l: {k: g for k, g in t[t.proportion_mapped > 0].groupby("sample_key")} for l, t in tabs.items()}; has = {l: set(t.sample_key.unique()) for l, t in tabs.items()}
    m = Mm.loc[S.sample_key]; age = pd.to_numeric(m.age, errors="coerce"); bmi = pd.to_numeric(m.BMI, errors="coerce"); ab, bb = bins(age, bmi)
    sex = m.gender.map({"female": "female", "male": "male"}).fillna("")
    meta = pd.DataFrame({"sample": S.sample_key.to_numpy(), "study": S.study.to_numpy(), "role": S.role.to_numpy(), "is_healthy": (S.condition == "control").to_numpy(), "condition": S.condition.to_numpy(), "age": age.to_numpy(), "age_bin": ab.to_numpy(),
                         "sex": sex.to_numpy(), "bmi": bmi.to_numpy(), "bmi_class": bb.to_numpy(), "country": m.country.fillna("").to_numpy(), "region": m.country.map(REGION).fillna("other").to_numpy(),
                         "westernized": m.non_westernized.map({"yes": "no", "no": "yes"}).fillna("").to_numpy()})
    meta = meta[meta.role.isin(["reference_pool", "heldout_healthy", "case_control"])]; meta.to_parquet(f"{OUT}/meta.parquet", index=False)

    def as_norm(ids):
        out = {}
        for sid in ids:
            d = {"qc": {"status": "retained", "rejections": [], "quality_band": band_of.get(sid, "unbanded"), "mapped_fraction_family": np.nan, "genome_equivalents_cov": np.nan}}
            for layer in tabs:
                if sid not in has[layer]:
                    continue   # layer not measured for this sample (function tables are missing for some studies)
                gg = per[layer].get(sid, tabs[layer].iloc[:0])
                d[layer] = pd.DataFrame({"feature_id": gg.feature_id.to_numpy(), "clr": gg.clr.to_numpy(), "proportion_mapped": gg.proportion_mapped.to_numpy()})
            out[sid] = d
        return out
    # Exclusion is a property of the STUDY, not of the sample (2026-10-06; before that date the duplicate exclusion was
    # applied sample by sample, which split GuptaA_2019 across two baselines: see below). Two rules, unioned:
    #   own        a study whose samples sit in the pool is scored against the baseline rebuilt without it;
    #   duplicate  a study that shares specimens with a pool study (work/s11/pipeB/cross_study_duplicates.tsv, built by
    #              s11_cross_study_duplicates.py) is scored against the baseline rebuilt without THAT study too.
    # Applying the duplicate rule to the whole study, not to the duplicated samples alone, is what keeps the cases and the
    # controls of one study on one baseline. GuptaA_2019 is the case that forced it: all 30 of its controls duplicate
    # DhakanDB_2019 and none of its 30 CRC cases does, so a per-sample rule scored its controls against the baseline
    # without DhakanDB_2019 and its cases against the full baseline, confounding the baseline with the case/control label.
    # Unioning the two rules is what makes a leave-one-study-out baseline hold out the whole donor cohort rather than one
    # deposit of it: LeChatelierE_2013 and NielsenHB_2014 are two cMD3 deposits of the same MetaHIT Danish cohort (171
    # specimens in the pool twice), so holding out one deposit alone would leave the other deposit's profile of the same
    # specimen in the baseline that scores it. Both now map to one -loso-LeChatelierE_2013_NielsenHB_2014 bundle.
    pool_st = set(meta.loc[meta.role == "reference_pool", "study"])
    D = pd.read_csv(f"{D_}/cross_study_duplicates.tsv", sep="\t") if os.path.exists(f"{D_}/cross_study_duplicates.tsv") else pd.DataFrame(columns=["study", "duplicate_study"])
    dup_of = D[D.duplicate_study.isin(pool_st)].groupby("study").duplicate_study.agg(lambda s: set(s)).to_dict()

    def excl_for(st):   # the studies a sample of study `st` must not see in its baseline, as a sorted tuple
        return tuple(sorted(({st} if st in pool_st else set()) | dup_of.get(st, set())))

    def bundle_for(ex):
        if not ex:
            return base
        b = f"{P}/refs/gut-reads-adult-global-v{BVER}-lenient-loso-" + "_".join(ex)
        if not os.path.exists(f"{b}/manifest.json"):
            subprocess.run([PY, f"{P}/scripts/{BUILDER}", "--exclude-study", *ex, "--bootstrap", "0"], check=True, stdout=subprocess.DEVNULL)
        return Bundle(b)

    base = Bundle(f"{P}/refs/gut-reads-adult-global-v{BVER}-lenient"); parts = []; summs = []
    for st, g in meta[meta.role == "reference_pool"].groupby("study"):
        ex = excl_for(st); assert st in ex, st
        r = score(as_norm(list(g["sample"])), bundle_for(ex)); parts.append(r["scores"]); summs.append(r["summaries"])
        print("scored", st, "against -loso-" + "_".join(ex), f"{time.time()-t0:.0f}s", flush=True)
    rest = meta[meta.role != "reference_pool"].copy()
    tup = {"_".join(excl_for(st)): excl_for(st) for st in rest.study.unique()}   # join the tuple so groupby sees a scalar key
    rest["excl"] = ["_".join(excl_for(st)) for st in rest.study]
    print("scored against a baseline without these studies:", {k or "(none)": int(v) for k, v in rest.groupby("excl").size().items()}, flush=True)
    for ex, grp in rest.groupby("excl"):
        bund = bundle_for(tup[ex])
        for i in range(0, len(grp), 400):
            r = score(as_norm(list(grp["sample"].iloc[i:i + 400])), bund); parts.append(r["scores"]); summs.append(r["summaries"])
    parts = [x for x in parts if len(x)]; summs = [x for x in summs if len(x)]
    Sc = pd.concat(parts, ignore_index=True).rename(columns={"analysis_id": "sample"}); Mu = pd.concat(summs, ignore_index=True).rename(columns={"analysis_id": "sample"})
    Sc[["sample", "layer", "feature_id", "band", "value", "percentile", "call"]].to_parquet(f"{OUT}/oos_scores.parquet", index=False)
    keep = [c for c in ["sample", "layer", "n_assessed", "n_low_raw", "n_high_raw", "n_missing", "frac_outside_raw", "excess_outside_p", "excess_outside_q", "core_distance_pct", "quality_band"] if c in Mu.columns]
    Mu[keep].to_parquet(f"{OUT}/oos_summaries.parquet", index=False); print("B done", len(Sc), len(Mu), f"{time.time()-t0:.0f}s")
