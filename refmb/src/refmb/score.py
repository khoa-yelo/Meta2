"""refmb.score — score normalized samples against a reference bundle. No dependency on project-internal tables.

Input: {sample_id: normalized dict} as produced by refmb.normalize (layers + qc). Output: long score table (one row per
sample x layer x feature), per-sample x layer summaries (with the sample-level excess-outside test), landscape coordinates and
rejections. The logic is the one the reference pool was scored with during validation (identical output verified on MGnify analyses).

Optional calibration: a Calibration object (refmb.calibrate) mapping a different pipeline's measurement space into the
reference's; it decides which layers are scorable. Without it, calibration_applied = "none", which is the contract for the
shipped containers (they reproduce the reference measurement).
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import binom

LAYER_SPEC = {  # value column, detection column (mirrors configs/s5.yaml layers)
    "taxonomy_family": ("clr", "proportion_mapped"), "taxonomy_genus": ("clr", "proportion_mapped"), "taxonomy_species": ("clr", "proportion_mapped"),
    "ko_eggnog": ("copies_per_genome", "copies_per_genome"), "ko_kofam": ("copies_per_genome", "copies_per_genome"),
    "pfam": ("copies_per_genome", "copies_per_genome"), "module": ("completeness", "completeness"),
    "ko_humann": ("clr", "proportion_mapped"), "pathway_humann": ("clr", "proportion_mapped"),   # read-based function layers (HUMAnN 3)
}
ALPHA = 0.05
MISSING_PREV = 0.70
OUTSIDE_EXPECTED = 0.05   # fraction of assessed features outside p2.5–p97.5 under the reference


def interp_percentile(vals, grid_vals, grid):
    """Percentile of each value within its feature's reference distribution.

    vals: array of n values; grid_vals: n x k array, row i holding feature i's reference values at the k grid percentiles
    (e.g. p1, p2.5, ..., p99); grid: the k percentiles. Linear interpolation on the grid; a value below the first grid
    value is placed at half the first percentile (0.5 for a grid starting at 1), one above the last at the symmetric
    point (99.5). A row with a missing grid value gives NaN (the feature is not assessable)."""
    out = np.empty(len(vals))
    for i in range(len(vals)):
        g = grid_vals[i]
        if np.isnan(g).any():
            out[i] = np.nan; continue
        out[i] = np.interp(vals[i], g, grid, left=grid[0] * 0.5, right=100 - (100 - grid[-1]) * 0.5)
    return out


def bh(p):
    """Benjamini-Hochberg adjusted p-values (q-values) for a 1-d array of p-values, capped at 1."""
    n = len(p); order = np.argsort(p); ranked = np.empty(n); ranked[order] = np.arange(1, n + 1)
    q = p * n / ranked; qs = q[order]; qs = np.minimum.accumulate(qs[::-1])[::-1]; q[order] = qs
    return np.minimum(q, 1.0)


class Bundle:
    """A reference bundle opened for scoring.

    Reads manifest.json (percentile grid, detection floor, which layers take percentiles within the sample's quality band),
    finds the available `features/<layer>.parquet` tables (loaded lazily by features()), and the optional landscape
    (PCA loadings, reference coordinates, core-distance tables). `path` is the bundle directory."""

    def __init__(self, path: str):
        self.path = path
        self.manifest = json.load(open(os.path.join(path, "manifest.json")))
        self.grid = np.array(self.manifest["percentile_grid"], dtype=float)
        self.pcols = [f"p{str(p).replace('.', '_')}" for p in self.manifest["percentile_grid"]]
        self.floor = self.manifest.get("scoring", {}).get("detection_floor_cpg", 0.0)
        self.band_layers = set(self.manifest.get("scoring", {}).get("percentiles_within_band_layers", []))
        self.layers = [l for l in self.manifest["layers"] if os.path.exists(os.path.join(path, "features", f"{l}.parquet"))]
        self._feat = {}
        lp = os.path.join(path, "landscape", "pca_loadings.parquet")
        self.landscape = None
        if os.path.exists(lp):
            load = pd.read_parquet(lp).set_index("feature_id")
            self.landscape = {"loadings": load, "pcs": [c for c in load.columns if c.startswith("PC")],
                              "coords": pd.read_parquet(os.path.join(path, "landscape", "reference_coords.parquet")),
                              "meta": json.load(open(os.path.join(path, "landscape", "core_distance.json"))),
                              "layer": self.manifest.get("landscape", {}).get("layer", "taxonomy_family") if isinstance(self.manifest.get("landscape"), dict) else "taxonomy_family",
                              "k": int(self.manifest.get("landscape", {}).get("k_neighbors", 10)) if isinstance(self.manifest.get("landscape"), dict) else 10}

    def features(self, layer: str) -> pd.DataFrame:
        if layer not in self._feat:
            self._feat[layer] = pd.read_parquet(os.path.join(self.path, "features", f"{layer}.parquet")).set_index("feature_id")
        return self._feat[layer]


def _long(norm: dict, layer: str) -> pd.DataFrame:
    rows = []
    for sid, d in norm.items():
        if layer in d and len(d[layer]):
            t = d[layer].copy(); t.insert(0, "analysis_id", sid); rows.append(t)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=["analysis_id", "feature_id"])


def score(norm: dict, bundle: Bundle, calibration=None, body_site_of=None, ignore_layer_gate: bool = False) -> dict:
    """norm: {sample_id: normalized dict}. body_site_of: optional {sample_id: body_site} to enforce the bundle's body site."""
    man = bundle.manifest; n_ref = man["n_samples"]
    rej = []; keep = []
    for sid, d in norm.items():
        q = d["qc"]
        if q["status"] != "retained":
            rej.append((sid, "NOT_NORMALIZED", "; ".join(f"{c}: {m}" for c, m in q["rejections"])))
        elif body_site_of and body_site_of.get(sid) not in (None, man["reference_key"]["body_site"]):
            rej.append((sid, "NO_MATCHING_REFERENCE", f"body site {body_site_of.get(sid)}"))
        elif q["quality_band"] in (None, "unbanded") or (q["quality_band"] == "unknown" and bundle.band_layers):
            rej.append((sid, "NO_MATCHING_REFERENCE", "no quality band" if q["quality_band"] != "unknown" else "depth unknown but the bundle scores within bands"))
        else:
            keep.append(sid)
    band_of = {sid: norm[sid]["qc"]["quality_band"] for sid in keep}
    cal_layers = None if calibration is None else (set(calibration.layers) if ignore_layer_gate else set(calibration.scorable_layers))
    scores = []; summ = []
    for name in bundle.layers:
        if name not in LAYER_SPEC or not keep:
            continue
        if cal_layers is not None and name not in cal_layers:
            for sid in keep:
                summ.append(dict(analysis_id=sid, layer=name, n_detected=0, n_assessed=0, n_low=0, n_high=0, n_low_raw=0, n_high_raw=0, n_missing=0,
                                 n_not_assessable=0, frac_outside_raw=np.nan, expected_false_positives_at_alpha=np.nan, weighted_deviation_score=np.nan,
                                 layer_status="PIPELINE_NOT_CALIBRATED"))
            continue
        vcol, dcol = LAYER_SPEC[name]
        ref = bundle.features(name)
        keep_all = keep; keep = [s for s in keep_all if name in norm[s]]   # a layer that was not measured for a sample is not scored for it (no expected-but-missing calls)
        t = _long({s: norm[s] for s in keep}, name)
        if t.empty:
            keep = keep_all; continue
        if calibration is not None:
            t = calibration.apply(name, t, vcol)
        gene = vcol == "copies_per_genome"
        t = t[t[dcol].astype(float) > 0]
        t = t[t["feature_id"].isin(ref.index)].copy()
        t["band"] = t["analysis_id"].map(band_of); t["value"] = t[vcol].astype(float)
        assessable = ~np.isnan(t["value"].to_numpy())
        if gene and bundle.floor > 0:
            assessable &= t["value"].to_numpy() >= bundle.floor
        grid = ref[bundle.pcols].reindex(t["feature_id"]).to_numpy().copy()
        if name in bundle.band_layers:
            for b in ("low", "medium", "high"):
                cols_b = [f"{c}_band_{b}" for c in bundle.pcols]
                if all(c in ref.columns for c in cols_b):
                    m = (t["band"] == b).to_numpy()
                    gb = ref[cols_b].reindex(t.loc[m, "feature_id"]).to_numpy(); ok = ~np.isnan(gb).any(axis=1)
                    sub = grid[m]; sub[ok] = gb[ok]; grid[m] = sub
        assessable &= ~np.isnan(grid).any(axis=1)
        pct = np.full(len(t), np.nan)
        pct[assessable] = interp_percentile(t["value"].to_numpy()[assessable], grid[assessable], bundle.grid)
        lo = grid[:, bundle.pcols.index("p2_5")]; hi = grid[:, bundle.pcols.index("p97_5")]; v = t["value"].to_numpy()
        t["percentile"] = pct
        t["call"] = np.where(~assessable, "not_assessable", np.where(v < lo, "low", np.where(v > hi, "high", "within")))
        if "cal_sigma" in t.columns:   # calibrated pipeline: conservative call needs the 95% mapping interval to clear the band edge
            from refmb.calibrate import KIND, fwd, inv as _inv
            k = KIND.get(name, "clr"); yv = fwd(k, v); sg = t["cal_sigma"].to_numpy()
            v_lo = _inv(k, yv - 1.96 * sg); v_hi = _inv(k, yv + 1.96 * sg)
            t["call_conservative"] = np.where(~assessable, "not_assessable", np.where(v_hi < lo, "low", np.where(v_lo > hi, "high", "within")))
        t["p_two_sided"] = np.nan; t["fdr_q"] = np.nan; pfloor = 1.0 / (n_ref + 1)
        for sid, g in t[assessable].groupby("analysis_id"):
            p = np.maximum(2 * np.minimum(g["percentile"], 100 - g["percentile"]) / 100.0, pfloor)
            t.loc[g.index, "p_two_sided"] = p.to_numpy(); t.loc[g.index, "fdr_q"] = bh(p.to_numpy())
        dev = t["call"].isin(["low", "high"]) & (t["fdr_q"] <= ALPHA)
        t["call_fdr"] = np.where(t["call"].isin(["low", "high"]) & ~dev, "within_after_fdr", t["call"])
        detected = t.groupby("analysis_id")["feature_id"].apply(set).to_dict()
        ebm = []
        for sid in keep:
            col = f"prevalence_band_{band_of[sid]}"
            if col not in ref.columns:
                continue
            cand = set(ref.index[(ref[col] >= MISSING_PREV) & ref.get("percentiles_available", True)])
            if calibration is not None:
                cand &= calibration.ebm_eligible(name)
            ebm.extend((sid, f, band_of[sid]) for f in sorted(cand - detected.get(sid, set())))
        e = pd.DataFrame(ebm, columns=["analysis_id", "feature_id", "band"])
        for c, val in (("value", np.nan), ("percentile", np.nan), ("call", "expected_but_missing"), ("call_fdr", "expected_but_missing"), ("p_two_sided", np.nan), ("fdr_q", np.nan)):
            e[c] = val
        keep_cols = ["analysis_id", "feature_id", "band", "value", "percentile", "call", "call_fdr", "p_two_sided", "fdr_q"] + [c for c in ("call_conservative", "cal_status") if c in t.columns]
        out = pd.concat([t[keep_cols], e], ignore_index=True)
        out["layer"] = name; scores.append(out)
        groups = dict(tuple(out.groupby("analysis_id")))
        for sid in keep:
            g = groups.get(sid, out.iloc[0:0])
            n_ass = int((g["call"] != "not_assessable").sum() - (g["call"] == "expected_but_missing").sum())
            n_out = int(((g["call"] == "low") | (g["call"] == "high")).sum())
            summ.append(dict(analysis_id=sid, layer=name, n_detected=int((g["call"] != "expected_but_missing").sum()), n_assessed=n_ass,
                             n_low=int((g["call_fdr"] == "low").sum()), n_high=int((g["call_fdr"] == "high").sum()),
                             n_low_raw=int((g["call"] == "low").sum()), n_high_raw=int((g["call"] == "high").sum()),
                             n_missing=int((g["call"] == "expected_but_missing").sum()), n_not_assessable=int((g["call"] == "not_assessable").sum()),
                             frac_outside_raw=n_out / n_ass if n_ass else np.nan, expected_false_positives_at_alpha=ALPHA * n_ass,
                             weighted_deviation_score=float(np.nansum(np.abs(g["percentile"] - 50) / 50 * (g["fdr_q"] <= ALPHA))) if n_ass else np.nan,
                             excess_outside_p=float(binom.sf(n_out - 1, n_ass, OUTSIDE_EXPECTED)) if n_ass else np.nan,
                             excess_outside_ratio=(n_out / n_ass) / OUTSIDE_EXPECTED if n_ass else np.nan, layer_status="scored"))
        keep = keep_all
    # landscape
    land = pd.DataFrame(); centring = None
    L = bundle.landscape
    if L is not None and keep and (cal_layers is None or L["layer"] in cal_layers):
        fam = _long({s: norm[s] for s in keep}, L["layer"])
        if calibration is not None:
            fam = calibration.apply(L["layer"], fam, "clr")
        if not fam.empty:
            meta = L["meta"]; load = L["loadings"]; pcs = L["pcs"]
            X = fam.pivot(index="analysis_id", columns="feature_id", values="clr").reindex(columns=load.index).fillna(meta["feature_fill_value"])
            if meta.get("column_means"):
                X = X - pd.Series(meta["column_means"], index=meta["column_features"]).reindex(load.index).to_numpy(); centring = "reference"
            else:
                X = X - X.mean(axis=0); centring = "query_batch"
            C = X.to_numpy() @ load[pcs].to_numpy(); cent = np.array(meta["centroid"]); core = np.linalg.norm(C - cent, axis=1)
            ref_core = np.sort(L["coords"]["core_distance"].to_numpy())
            dd, ii = cKDTree(L["coords"][pcs].to_numpy()).query(C, k=L["k"])
            rows = []
            for i, sid in enumerate(X.index):
                nb = L["coords"].iloc[ii[i]]
                rows.append(dict(analysis_id=sid, core_distance=float(core[i]), core_distance_pct=float(np.searchsorted(ref_core, core[i]) / len(ref_core) * 100),
                                 neighbor_studies=json.dumps(nb["study_bioproject"].value_counts().to_dict()), neighbor_bands=json.dumps(nb["quality_band"].value_counts().to_dict()),
                                 **{f"PC{k+1}": float(C[i, k]) for k in range(min(3, len(pcs)))}))
            land = pd.DataFrame(rows)
    S = pd.DataFrame(summ)
    if not S.empty:
        qcs = pd.DataFrame([{"analysis_id": s, "quality_band": norm[s]["qc"]["quality_band"], "mapped_fraction_family": norm[s]["qc"]["mapped_fraction_family"],
                             "genome_equivalents_cov": norm[s]["qc"]["genome_equivalents_cov"]} for s in keep])
        S = S.merge(qcs, on="analysis_id", how="left")
        S["bundle_id"] = man["bundle_id"]; S["calibration_applied"] = "none" if calibration is None else calibration.id
        S["fallback_used"] = False; S["landscape_centring"] = centring
        if not land.empty:
            S = S.merge(land, on="analysis_id", how="left")
        # BH over samples within layer for the excess test
        for c in ("excess_outside_p", "excess_outside_ratio"):
            if c not in S.columns: S[c] = np.nan
        S["excess_outside_q"] = np.nan
        for name, g in S[S["excess_outside_p"].notna()].groupby("layer"):
            S.loc[g.index, "excess_outside_q"] = bh(g["excess_outside_p"].to_numpy())
    return {"scores": pd.concat(scores, ignore_index=True) if scores else pd.DataFrame(), "summaries": S, "landscape": land,
            "rejections": pd.DataFrame(rej, columns=["analysis_id", "reason_code", "detail"])}
