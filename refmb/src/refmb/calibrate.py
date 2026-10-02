"""refmb.calibrate — tier B: map a user pipeline's measurement space into the reference's using the anchor panel.

Model (configs/s11.yaml): per layer, transformed values (log10 copies per genome; CLR as is; completeness as is) are related by a
per-feature affine map y_ref = a_f + b_f x_user, with a_f, b_f and the residual sd shrunk toward a layer-global fit. Anchor-level
K-fold cross-validation gives out-of-sample residuals; those drive the feature gate (residual sd vs reference band half-width) and
the layer gate (fraction outside, resolution Spearman, surviving fraction). A detection model gives per-feature sensitivity
P(user detects | reference detects), which decides whether expected_but_missing may be called.

A Calibration is saved as a directory: calibration.json (id, config, per-layer gate results) + features.parquet (per-feature params).
"""
from __future__ import annotations

import hashlib
import json
import os
import time

import numpy as np
import pandas as pd
import yaml

from refmb.score import LAYER_SPEC, interp_percentile

KIND = {"taxonomy_family": "clr", "taxonomy_genus": "clr", "taxonomy_species": "clr", "ko_eggnog": "gene", "ko_kofam": "gene", "pfam": "gene", "module": "completeness"}


def fwd(kind, v):
    v = np.asarray(v, dtype=float)
    return np.log10(np.where(v > 0, v, np.nan)) if kind == "gene" else v


def inv(kind, y):
    return np.power(10.0, y) if kind == "gene" else y


def _ols(x, y):
    if len(x) < 2 or np.nanvar(x) == 0:
        return np.nan, np.nan
    b = np.cov(x, y, bias=True)[0, 1] / np.var(x); a = y.mean() - b * x.mean()
    return a, b


def _fit_params(pairs: pd.DataFrame, cfg: dict) -> tuple[dict, pd.DataFrame]:
    """pairs: feature_id, x (user), y (reference), both detected. Returns global params and per-feature shrunk params for the
    prediction y = a + b x.

    kind = inverse   : regress y on x (conditional mean; shrinks predictions toward the mean, under-calls extremes)
    kind = classical : regress x on y and invert (classical instrument calibration; preserves the spread of y, noisier)
    Shrinkage toward the layer-global fit is applied to the fitted slope/offset in the regression's own parameterization."""
    m = cfg["model"]; kind = m.get("kind", "inverse"); classical = kind.startswith("classical")
    X = pairs["x"].to_numpy(); Y = pairs["y"].to_numpy()
    if kind == "matched":
        return _fit_matched(pairs, cfg)
    if classical:
        al_g, be_g = _ols(Y, X)                     # x = alpha + beta y
        s_g_pred = None
    else:
        al_g, be_g = _ols(X, Y)                     # y = alpha + beta x
    def to_pred(al, be):
        return (-al / be, 1.0 / be) if classical else (al, be)
    a_g, b_g = to_pred(al_g, be_g)
    res_g = Y - (a_g + b_g * X); s_g = float(np.std(res_g, ddof=1)) if len(res_g) > 2 else np.nan
    rows = []
    for f, g in pairs.groupby("feature_id", sort=True):
        n = len(g); x = g["x"].to_numpy(); y = g["y"].to_numpy()
        al_f, be_f = (_ols(y, x) if classical else _ols(x, y)) if n >= 3 else (np.nan, np.nan)
        be_s = be_g if (np.isnan(be_f) or (classical and be_f <= 0.05 * be_g)) else (n * be_f + m["shrink_k_slope"] * be_g) / (n + m["shrink_k_slope"])
        if classical:
            al_raw = float(np.mean(x - be_s * y)) if n else al_g
        else:
            al_raw = float(np.mean(y - be_s * x)) if n else al_g
        al_s = (n * al_raw + m["shrink_k_offset"] * al_g) / (n + m["shrink_k_offset"])
        a_s, b_s = to_pred(al_s, be_s)
        r = y - (a_s + b_s * x); s_f = float(np.sqrt((n * r.var(ddof=0) + m["shrink_k_sigma"] * s_g ** 2) / (n + m["shrink_k_sigma"]))) if n else s_g
        rows.append((f, n, a_s, b_s, s_f, n < m["min_pairs_per_feature"]))
    F = pd.DataFrame(rows, columns=["feature_id", "n_pairs", "a", "b", "sigma", "inherited"])
    return {"a": float(a_g), "b": float(b_g), "sigma": s_g, "n_pairs": int(len(pairs)), "kind": "classical" if classical else "inverse"}, F


def _fit_matched(pairs: pd.DataFrame, cfg: dict) -> tuple[dict, pd.DataFrame]:
    """kind = matched: marginal-variance matching (geometric-mean / Deming-type map). y = a_f + b_f x with b_f = sd(y)/sd(x) and
    a_f = mean(y) - b_f mean(x), so the calibrated values reproduce the reference spread by construction; resolution is then set
    by the anchor correlation alone. log(b_f) is shrunk toward the layer-global log ratio, a_f toward the global offset."""
    m = cfg["model"]; X = pairs["x"].to_numpy(); Y = pairs["y"].to_numpy()
    b_g = float(np.std(Y, ddof=1) / np.std(X, ddof=1)); a_g = float(Y.mean() - b_g * X.mean())
    res_g = Y - (a_g + b_g * X); s_g = float(np.std(res_g, ddof=1)) if len(res_g) > 2 else np.nan
    rows = []
    for f, g in pairs.groupby("feature_id", sort=True):
        n = len(g); x = g["x"].to_numpy(); y = g["y"].to_numpy()
        sx = np.std(x, ddof=1) if n >= 3 else 0.0; sy = np.std(y, ddof=1) if n >= 3 else 0.0
        lb_f = np.log(sy / sx) if (sx > 0 and sy > 0) else np.log(b_g)
        lb_s = (n * lb_f + m["shrink_k_slope"] * np.log(b_g)) / (n + m["shrink_k_slope"]); b_s = float(np.exp(lb_s))
        a_raw = float(np.mean(y - b_s * x)) if n else a_g
        a_s = (n * a_raw + m["shrink_k_offset"] * a_g) / (n + m["shrink_k_offset"])
        r = y - (a_s + b_s * x); s_f = float(np.sqrt((n * r.var(ddof=0) + m["shrink_k_sigma"] * s_g ** 2) / (n + m["shrink_k_sigma"]))) if n else s_g
        rows.append((f, n, a_s, b_s, s_f, n < m["min_pairs_per_feature"]))
    F = pd.DataFrame(rows, columns=["feature_id", "n_pairs", "a", "b", "sigma", "inherited"])
    return {"a": a_g, "b": b_g, "sigma": s_g, "n_pairs": int(len(pairs)), "kind": "matched"}, F


class Calibration:
    def __init__(self, cid: str, cfg: dict, layers: dict, features: pd.DataFrame, meta: dict):
        self.id = cid; self.cfg = cfg; self.layers = layers; self.features = features; self.meta = meta
        self._F = {l: g.set_index("feature_id") for l, g in features.groupby("layer")}

    @property
    def scorable_layers(self):
        return [l for l, r in self.layers.items() if r["status"] == "scorable"]

    def ebm_eligible(self, layer: str) -> set:
        F = self._F.get(layer)
        if F is None:
            return set()
        return set(F.index[(F["sensitivity"] >= self.cfg["detection"]["min_sensitivity_for_missing_call"]) & (F["status"] == "calibrated")])

    def apply(self, layer: str, t: pd.DataFrame, vcol: str) -> pd.DataFrame:
        """t: analysis_id, feature_id, vcol in user space -> vcol replaced by the calibrated reference-space value; adds
        cal_sigma (transformed scale) and cal_status. Uncalibrated features get NaN (scored not_assessable)."""
        F = self._F.get(layer); kind = KIND[layer]; t = t.copy()
        if F is None:
            t[vcol] = np.nan; t["cal_sigma"] = np.nan; t["cal_status"] = "PIPELINE_NOT_CALIBRATED"; return t
        P = F.reindex(t["feature_id"])
        x = fwd(kind, t[vcol].to_numpy())
        y = P["a"].to_numpy() + P["b"].to_numpy() * x
        ok = (P["status"].to_numpy() == "calibrated")
        val = np.where(ok, inv(kind, y), np.nan)
        if kind == "completeness":
            val = np.clip(val, 0, 1)
        t[vcol] = val
        if vcol != LAYER_SPEC[layer][1]:      # keep detection column consistent with the calibrated value where they differ
            pass
        t["cal_sigma"] = np.where(ok, P["sigma"].to_numpy(), np.nan)
        t["cal_status"] = np.where(P["status"].isna(), "FEATURE_NOT_IN_ANCHORS", P["status"].fillna("").to_numpy())
        return t

    def save(self, out: str):
        os.makedirs(out, exist_ok=True)
        json.dump({"id": self.id, "config": self.cfg, "layers": self.layers, "meta": self.meta}, open(os.path.join(out, "calibration.json"), "w"), indent=1, default=float)
        self.features.to_parquet(os.path.join(out, "features.parquet"), index=False)

    @classmethod
    def load(cls, path: str):
        j = json.load(open(os.path.join(path, "calibration.json")))
        return cls(j["id"], j["config"], j["layers"], pd.read_parquet(os.path.join(path, "features.parquet")), j["meta"])


def fit(user: dict, ref: dict, anchors: list, band_of: dict, bundle, cfg: dict, pipeline_label: str, rng_seed: int | None = None) -> Calibration:
    """user / ref: {layer: DataFrame(analysis_id, feature_id, value)} restricted to anchors (any extra ids ignored).
    band_of: analysis_id -> quality band (for band-specific percentile grids). bundle: refmb.score.Bundle."""
    t0 = time.time(); rng = np.random.default_rng(cfg.get("seed", 0) if rng_seed is None else rng_seed)
    anchors = sorted(anchors); folds = np.array_split(rng.permutation(anchors), cfg["model"]["cv_folds"])
    fold_of = {a: i for i, f in enumerate(folds) for a in f}
    layers = {}; feats = []
    for layer in ref:
        if layer not in user or layer not in bundle.layers:
            layers[layer] = {"status": "PIPELINE_NOT_CALIBRATED", "reason": "layer not produced by the user pipeline"}; continue
        kind = KIND[layer]; refF = bundle.features(layer)
        u = user[layer][user[layer]["analysis_id"].isin(anchors)]; r = ref[layer][ref[layer]["analysis_id"].isin(anchors)]
        u = u[u["feature_id"].isin(refF.index)]; r = r[r["feature_id"].isin(refF.index)]
        floor = bundle.floor if kind == "gene" else 0.0
        u_det = u[u["value"] > 0]; r_det = r[r["value"] >= max(floor, np.finfo(float).tiny)]
        j = r_det.merge(u_det, on=["analysis_id", "feature_id"], how="outer", suffixes=("_ref", "_user"), indicator=True)
        # detection model: sensitivity per feature among anchors where the reference detects
        det = j.assign(ref_det=j["_merge"].isin(["both", "left_only"]), user_det=j["_merge"].isin(["both", "right_only"]))
        sens = det[det["ref_det"]].groupby("feature_id")["user_det"].agg(["mean", "size"]).rename(columns={"mean": "sensitivity", "size": "n_ref_detected"})
        pairs = j[j["_merge"] == "both"].copy()
        pairs["x"] = fwd(kind, pairs["value_user"].to_numpy()); pairs["y"] = fwd(kind, pairs["value_ref"].to_numpy())
        pairs = pairs.dropna(subset=["x", "y"])
        if len(pairs) < 50:
            layers[layer] = {"status": "PIPELINE_NOT_CALIBRATED", "reason": f"only {len(pairs)} anchor pairs"}; continue
        G, F = _fit_params(pairs, cfg)
        # cross-validated predictions
        pairs["fold"] = pairs["analysis_id"].map(fold_of); yhat = np.full(len(pairs), np.nan)
        for k in range(len(folds)):
            tr = pairs[pairs["fold"] != k]; te = pairs["fold"] == k
            if not te.any() or len(tr) < 50:
                continue
            Gk, Fk = _fit_params(tr, cfg); Pk = Fk.set_index("feature_id").reindex(pairs.loc[te, "feature_id"])
            a = Pk["a"].fillna(Gk["a"]).to_numpy(); b = Pk["b"].fillna(Gk["b"]).to_numpy()
            yhat[te.to_numpy()] = a + b * pairs.loc[te, "x"].to_numpy()
        pairs["yhat_cv"] = yhat; pairs["res_cv"] = pairs["y"] - pairs["yhat_cv"]
        cvs = pairs.groupby("feature_id")["res_cv"].agg(lambda s: float(np.sqrt(np.nanmean(np.square(s)))) if s.notna().sum() >= 3 else np.nan).rename("cv_sigma")
        F = F.set_index("feature_id").join(cvs).join(sens, how="outer")
        # reference band half-width on the transformed scale (pooled grid)
        lo = fwd(kind, refF["p2_5"].reindex(F.index).to_numpy()); hi = fwd(kind, refF["p97_5"].reindex(F.index).to_numpy())
        F["halfwidth"] = (hi - lo) / 2
        F["cv_sigma_filled"] = F["cv_sigma"].fillna(F["sigma"]).fillna(G["sigma"])
        has_pct = refF["percentiles_available"].reindex(F.index).fillna(False).astype(bool) if "percentiles_available" in refF.columns else pd.Series(True, index=F.index)
        F["status"] = np.where(F["a"].isna(), "FEATURE_NOT_PRODUCED",
                       np.where(F["cv_sigma_filled"] > cfg["feature_gate"]["max_cv_sigma_over_halfwidth"] * F["halfwidth"], "FEATURE_NOT_CALIBRATED", "calibrated"))
        F["status"] = np.where(~has_pct.to_numpy() & (F["status"] == "calibrated"), "NO_REFERENCE_PERCENTILES", F["status"])
        # layer gate on cross-validated anchors: percentiles of calibrated vs native values, band-aware
        cal_ok = pairs["feature_id"].map(F["status"]).eq("calibrated").to_numpy() & ~np.isnan(pairs["yhat_cv"].to_numpy())
        pc = pairs[cal_ok].copy(); pc["band"] = pc["analysis_id"].map(band_of)
        grid = refF[bundle.pcols].reindex(pc["feature_id"]).to_numpy().copy()
        if layer in bundle.band_layers:
            for b in ("low", "medium", "high"):
                cols_b = [f"{c}_band_{b}" for c in bundle.pcols]
                if all(c in refF.columns for c in cols_b):
                    m = (pc["band"] == b).to_numpy(); gb = refF[cols_b].reindex(pc.loc[m, "feature_id"]).to_numpy(); ok = ~np.isnan(gb).any(axis=1)
                    sub = grid[m]; sub[ok] = gb[ok]; grid[m] = sub
        okg = ~np.isnan(grid).any(axis=1); pc = pc[okg]; grid = grid[okg]
        v_cal = inv(kind, pc["yhat_cv"].to_numpy()); v_nat = pc["value_ref"].to_numpy()
        p_cal = interp_percentile(v_cal, grid, bundle.grid); p_nat = interp_percentile(v_nat, grid, bundle.grid)
        lo_b = grid[:, bundle.pcols.index("p2_5")]; hi_b = grid[:, bundle.pcols.index("p97_5")]
        out_cal = ((v_cal < lo_b) | (v_cal > hi_b)); out_nat = ((v_nat < lo_b) | (v_nat > hi_b))
        per_sample_cal = pd.Series(out_cal).groupby(pc["analysis_id"].to_numpy()).mean(); per_sample_nat = pd.Series(out_nat).groupby(pc["analysis_id"].to_numpy()).mean()
        from scipy.stats import spearmanr
        rho = float(spearmanr(p_cal, p_nat).statistic) if len(p_cal) > 10 else np.nan
        n_pct = int(has_pct.sum()) if len(has_pct) else int((refF.get("percentiles_available", pd.Series(True, index=refF.index))).sum())
        surviving = int((F["status"] == "calibrated").sum()) / max(1, int(refF["percentiles_available"].sum() if "percentiles_available" in refF.columns else len(refF)))
        g = cfg["layer_gate"]; fo = float(per_sample_cal.median()) if len(per_sample_cal) else np.nan
        fn = float(per_sample_nat.median()) if len(per_sample_nat) else np.nan
        ratio = fo / fn if fn and fn > 0 else np.nan
        passed = (fo <= g["max_frac_outside"]) and (g["min_ratio_to_native"] <= ratio <= g["max_ratio_to_native"]) and (rho >= g["min_resolution_spearman"]) and (surviving >= g["min_surviving_fraction"])
        layers[layer] = {"status": "scorable" if passed else "PIPELINE_NOT_CALIBRATED", "kind": kind, "global": G,
                         "cv": {"n_pairs": int(len(pc)), "median_frac_outside_calibrated": fo, "median_frac_outside_native": fn, "ratio_to_native": ratio,
                                "resolution_spearman": rho, "surviving_fraction": surviving, "n_calibrated": int((F["status"] == "calibrated").sum()),
                                "n_not_calibrated": int((F["status"] == "FEATURE_NOT_CALIBRATED").sum()), "n_not_produced": int((F["status"] == "FEATURE_NOT_PRODUCED").sum()),
                                "n_reference_features_with_percentiles": int(refF["percentiles_available"].sum()) if "percentiles_available" in refF.columns else len(refF)},
                         "gate": g, "reason": None if passed else "layer gate failed"}
        F["layer"] = layer; feats.append(F.reset_index())
    features = pd.concat(feats, ignore_index=True) if feats else pd.DataFrame(columns=["layer", "feature_id", "status"])
    cid = f"cal-{pipeline_label}-{hashlib.sha256((bundle.manifest['bundle_id'] + json.dumps(cfg, sort_keys=True) + ','.join(anchors)).encode()).hexdigest()[:10]}"
    meta = {"bundle_id": bundle.manifest["bundle_id"], "pipeline": pipeline_label, "n_anchors": len(anchors), "anchors": anchors, "fitted": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "wall_s": round(time.time() - t0, 1)}
    return Calibration(cid, cfg, layers, features, meta)


def load_config(path: str | None = None) -> dict:
    return yaml.safe_load(open(path or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "configs", "s11.yaml")))
