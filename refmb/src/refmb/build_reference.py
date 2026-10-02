"""refmb.build_reference — generic reference-bundle builder (same bundle format as S5, readable by refmb.score).

Inputs
  samples : DataFrame sample_id, study, band            (the reference pool; equal per-study weights)
  layers  : {layer_name: DataFrame sample_id, feature_id, value, detect}   (long tables; detect > 0 = detected)
  spec    : dict with percentiles, min_detected_for_percentiles, bootstrap_n, band_layers (percentiles within band), landscape
            {layer, n_components, min_prevalence, k_neighbors}, detection_floor (gene layers), manifest fields;
            layer_samples {layer_name: iterable of sample_id} restricts a layer to the pool samples that were measured for it
            (prevalence denominators and study weights then come from that subset; the manifest records n_samples per layer)
Writes bundle_dir/{features/*.parquet, landscape/*, provenance/*, manifest.json}. Statistics only: no sample ids are stored.
Logic copied from scripts/s5_build_reference.py (layer_stats, landscape) so the two builders stay numerically identical.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import time

import numpy as np
import pandas as pd


def sha256(p):
    """Hex sha256 of a file."""
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def wilson(k, n, z=1.96):
    """Wilson score interval (low, high) for k successes out of n."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def weighted_quantiles(vals, w, probs):
    """Weighted quantiles of vals at the given probabilities (0-1), interpolating on cumulative normalized weights."""
    n, F = vals.shape
    order = np.argsort(vals, axis=0, kind="stable"); sv = np.take_along_axis(vals, order, axis=0)
    sw = w[order] * ~np.isnan(sv); cw = np.cumsum(sw, axis=0); tot = cw[-1]
    out = np.full((len(probs), F), np.nan); ok = tot > 0
    for i, p in enumerate(probs):
        idx = np.minimum((cw < p * tot).sum(axis=0), n - 1)
        out[i, ok] = sv[idx[ok], np.arange(F)[ok]]
    return out


def layer_stats(name, t, samples, spec, rng, gene_layer=False):
    ids = samples["sample_id"].to_numpy(); id_ix = {a: i for i, a in enumerate(ids)}
    t = t[t["sample_id"].isin(id_ix) & (t["detect"] > 0)]
    feats = np.array(sorted(t["feature_id"].astype(str).unique())); f_ix = {f: j for j, f in enumerate(feats)}
    V = np.full((len(ids), len(feats)), np.nan)
    V[t["sample_id"].map(id_ix).to_numpy(), t["feature_id"].astype(str).map(f_ix).to_numpy()] = t["value"].to_numpy(dtype=float)
    floor = spec.get("detection_floor", 0.0)
    if gene_layer and floor > 0:
        V[V < floor] = np.nan
    det = ~np.isnan(V)
    studies = samples["study"].to_numpy(); us, inv_s = np.unique(studies, return_inverse=True)
    w = 1.0 / np.bincount(inv_s)[inv_s]; w = w / w.sum()
    out = pd.DataFrame({"feature_id": feats})
    out["n_detected"] = det.sum(axis=0)
    out["n_studies"] = pd.DataFrame({"f": t["feature_id"].astype(str), "s": t["sample_id"].map(dict(zip(ids, studies)))}).groupby("f")["s"].nunique().reindex(feats).to_numpy()
    out["prevalence_weighted"] = (det * w[:, None]).sum(axis=0)
    bands = samples["band"].to_numpy()
    for b in sorted(set(bands)):
        m = bands == b; wb = w[m] / w[m].sum() if w[m].sum() > 0 else w[m]
        out[f"prevalence_band_{b}"] = (det[m] * wb[:, None]).sum(axis=0)
        k = det[m].sum(axis=0); n = int(m.sum()); ci = np.array([wilson(kk, n) for kk in k])
        out[f"prevalence_ci_low_{b}"] = ci[:, 0]; out[f"prevalence_ci_high_{b}"] = ci[:, 1]; out[f"n_band_{b}"] = n
    P = spec["percentiles"]; probs = np.array(P) / 100.0
    Q = weighted_quantiles(V, w, probs)
    for p, row in zip(P, Q):
        out[f"p{str(p).replace('.', '_')}"] = row
    if name in spec.get("band_layers", []):
        for b in sorted(set(bands)):
            m = bands == b; wb = w[m] / w[m].sum() if w[m].sum() > 0 else w[m]
            Qb = weighted_quantiles(V[m], wb, probs); nb = det[m].sum(axis=0)
            for p, row in zip(P, Qb):
                out[f"p{str(p).replace('.', '_')}_band_{b}"] = np.where(nb >= spec["min_detected_for_percentiles"], row, np.nan)
            out[f"n_detected_band_{b}"] = nb
    n_boot = spec.get("bootstrap_n", 0)
    if n_boot > 0:
        lo_i = P.index(2.5); hi_i = P.index(97.5)
        bl = np.empty((n_boot, len(feats))); bh_ = np.empty((n_boot, len(feats)))
        for b in range(n_boot):
            pick = rng.integers(0, len(us), len(us)); mult = np.bincount(pick, minlength=len(us))[inv_s]
            wb = w * mult; wb = wb / wb.sum() if wb.sum() > 0 else wb
            qb = weighted_quantiles(V, wb, probs[[lo_i, hi_i]]); bl[b] = qb[0]; bh_[b] = qb[1]
        out["p2_5_ci_low"] = np.nanpercentile(bl, 2.5, axis=0); out["p2_5_ci_high"] = np.nanpercentile(bl, 97.5, axis=0)
        out["p97_5_ci_low"] = np.nanpercentile(bh_, 2.5, axis=0); out["p97_5_ci_high"] = np.nanpercentile(bh_, 97.5, axis=0)
    few = out["n_detected"] < spec["min_detected_for_percentiles"]
    out.loc[few, [c for c in out.columns if re.match(r"^p\d", c)]] = np.nan
    out["percentiles_available"] = ~few; out["layer"] = name; out["reliability_tier"] = "NA"
    return out, V, feats, w


def landscape(V, feats, w, samples, spec, rng):
    L = spec["landscape"]; prev = (~np.isnan(V)).mean(axis=0); keep = prev >= L["min_prevalence"]
    fill = float(np.nanmin(V[:, keep])) if keep.any() else 0.0
    X = np.nan_to_num(V[:, keep], nan=fill); col_means = X.mean(axis=0); X = X - col_means
    k = min(L["n_components"], X.shape[1], X.shape[0] - 1)
    U, S, Vt = np.linalg.svd(X, full_matrices=False); coords = U[:, :k] * S[:k]; var = (S ** 2) / (S ** 2).sum()
    centroid = coords.mean(axis=0); core = np.linalg.norm(coords - centroid, axis=1)
    from scipy.spatial import cKDTree
    dd, _ = cKDTree(coords).query(coords, k=L["k_neighbors"] + 1); knn = dd[:, 1:].mean(axis=1)
    load = pd.DataFrame(Vt[:k].T, index=feats[keep], columns=[f"PC{i+1}" for i in range(k)]).reset_index().rename(columns={"index": "feature_id"})
    ref = pd.DataFrame(coords, columns=[f"PC{i+1}" for i in range(k)]); ref.insert(0, "study_bioproject", samples["study"].to_numpy())
    ref["core_distance"] = core; ref["knn_distance"] = knn; ref["quality_band"] = samples["band"].to_numpy()
    P = spec["percentiles"]
    return load, ref, {"variance_explained": var[:k].round(4).tolist(), "core_distance_percentiles": {str(p): float(np.percentile(core, p)) for p in P},
                       "knn_distance_percentiles": {str(p): float(np.percentile(knn, p)) for p in P}, "centroid": centroid.round(6).tolist(),
                       "feature_fill_value": fill, "column_means": col_means.round(6).tolist(), "column_features": feats[keep].tolist()}


def build(out: str, samples: pd.DataFrame, layers: dict, spec: dict, manifest_extra: dict, exclusions: pd.DataFrame | None = None, seed: int = 0) -> dict:
    """Build a reference bundle directory from long per-layer tables and return its manifest.

    out: bundle directory to create. samples: DataFrame(sample_id, study, band) of the reference pool. layers: {layer name:
    DataFrame(sample_id, feature_id, value, detect)} with detect > 0 meaning detected. spec: percentiles, min_detected_for_percentiles,
    bootstrap_n, band_layers (layers whose percentiles are also taken within each band), gene_layers, detection_floor,
    landscape settings. manifest_extra is merged into manifest.json; exclusions (sample_id, reason_code, study) is written to
    provenance/pool_exclusions.tsv. Samples are weighted so every study carries equal total weight; percentiles are weighted
    quantiles among detected samples; a study-level bootstrap gives intervals on the 2.5th and 97.5th percentiles. The caller
    adds the normalization/ directory beside the output."""
    t0 = time.time(); rng = np.random.default_rng(seed)
    for sub in ["features", "landscape", "provenance"]:
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    samples = samples.sort_values("sample_id").reset_index(drop=True)
    checks = {}; layer_summ = {}; land_meta = None
    for name, t in layers.items():
        gene = spec.get("gene_layers", []) and name in spec["gene_layers"]
        ls = samples
        if name in spec.get("layer_samples", {}):
            ls = samples[samples["sample_id"].isin(set(spec["layer_samples"][name]))].reset_index(drop=True)
        st, V, feats, w = layer_stats(name, t, ls, spec, rng, gene_layer=bool(gene))
        p = os.path.join(out, "features", f"{name}.parquet"); st.to_parquet(p, index=False); checks[f"features/{name}.parquet"] = "sha256:" + sha256(p)
        pc = [f"p{str(q).replace('.', '_')}" for q in spec["percentiles"]]
        mono = np.all(np.diff(st.loc[st["percentiles_available"], pc].to_numpy(), axis=1) >= -1e-12, axis=1)
        layer_summ[name] = {"n_samples": int(len(ls)), "n_studies": int(ls["study"].nunique()), "n_features": int(len(st)), "n_with_percentiles": int(st["percentiles_available"].sum()), "non_monotone_percentile_vectors": int((~mono).sum()),
                            "prevalence_bands": sorted(c for c in st.columns if c.startswith("prevalence_band_"))}
        if name == spec["landscape"]["layer"]:
            load, ref, land_meta = landscape(V, feats, w, ls, spec, rng)
            load.to_parquet(os.path.join(out, "landscape", "pca_loadings.parquet"), index=False); ref.to_parquet(os.path.join(out, "landscape", "reference_coords.parquet"), index=False)
            json.dump(land_meta, open(os.path.join(out, "landscape", "core_distance.json"), "w"), indent=1)
    studies = samples.groupby("study").agg(n=("sample_id", "size"), bands=("band", lambda s: dict(s.value_counts()))).reset_index().rename(columns={"study": "study_bioproject"})
    studies.to_csv(os.path.join(out, "provenance", "studies.tsv"), sep="\t", index=False)
    (exclusions if exclusions is not None else pd.DataFrame(columns=["analysis_id", "reason_code", "study_bioproject"])).to_csv(os.path.join(out, "provenance", "pool_exclusions.tsv"), sep="\t", index=False)
    manifest = {"bundle_id": os.path.basename(out), "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "builder_version": "refmb 0.2.0 (refmb/build_reference.py)",
                "n_samples": int(len(samples)), "n_studies": int(samples["study"].nunique()), "band_counts": samples["band"].value_counts().to_dict(),
                "study_weighting": "equal_per_study", "percentile_grid": spec["percentiles"], "bootstrap": {"n": spec.get("bootstrap_n", 0), "unit": "study"},
                "scoring": {"detection_floor_cpg": spec.get("detection_floor", 0.0), "percentiles_within_band": bool(spec.get("band_layers")), "percentiles_within_band_layers": spec.get("band_layers", [])},
                "layers": layer_summ, "landscape": {"variance_explained": (land_meta or {}).get("variance_explained"), "layer": spec["landscape"]["layer"], "k_neighbors": spec["landscape"]["k_neighbors"]},
                "checksums": checks, "wall_seconds": round(time.time() - t0, 1), **manifest_extra}
    json.dump(manifest, open(os.path.join(out, "manifest.json"), "w"), indent=1, sort_keys=True, default=float)
    return manifest
