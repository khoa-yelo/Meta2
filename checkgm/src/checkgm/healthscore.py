"""checkgm.healthscore — the shipped supervised health scores, applied to what `checkgm assess` wrote.

Each model is an L1-penalised logistic regression on standardized checkgm percentiles. They were fitted in the analysis
workspace (read pipeline: scripts/s16_health_score_B.py; assembly pipeline: scripts/s8_export_assembly_model.py, the
leave-one-study-out pipeline of scripts/s8_evaluate.py refitted on all 11 case/control studies) and exported whole — intercept, and per feature its coefficient together with
the StandardScaler centre and scale — so that applying it needs neither sklearn nor that workspace. The score is the
decision function of the disease class, `intercept + sum coef * (value - mean) / scale`: higher is more case-like. It is
not a probability, and the export carries no calibration that would make it one; only its ordering was validated
(leave-one-study-out AUROC, results/s16/health_score_B.md and results/s8/report_loso_l1.md).

One model per pipeline, and the match is enforced rather than documented. The read model's features are percentiles of
MetaPhlAn 3 taxonomy against a gut-reads baseline; the assembly model's are percentiles of family, genus, eggNOG KO and
KEGG module against a gut-assembly baseline, plus per-layer deviation counts. The two measure different spaces, in
which each other's coefficients mean nothing, so the model is chosen from the assessment's baseline, an assessment that
identifies as neither (or both) is refused, and an explicitly supplied model fitted on the other pipeline is refused:
a wrong number here would look exactly like a right one.

The shipped model is the percentile-only variant. The pre-specified primary feature set also carried species presence
and scored marginally better overall (0.747 against 0.719 mean AUROC over the 14 case/control studies), but its single
largest coefficient was the presence of [Collinsella] massiliensis — a batch artifact rather than biology — so the
percentile-only model is the one that ships, as the paper reports. load_model() refuses a model file carrying
`present:` features both for that reason and because assess keys taxa by NCBI taxid, which no MetaPhlAn species name can
be matched against without the profiler's own naming table.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

MODEL_FILENAME = "health_score_model.json"   # the copy packaged beside this module; results/s16/health_score_B_model.json in the workspace
MODEL_FILES = {"reads": MODEL_FILENAME,
               "assembly": "health_score_model_assembly.json"}   # results/s8/health_score_A_model.json in the workspace
PIPELINE_WORD = {"reads": "read", "assembly": "assembly"}

# A percentile feature is named for its layer with the layer's generic suffix dropped, as the training script built the
# column names ("taxonomy_species" -> "species:1234", "ko_humann" -> "ko:K01810"); the assembly layers keep their full
# name ("ko_eggnog:K00001", "module:M00001"), since "ko" already means the HUMAnN layer. The per-layer summary features
# keep the full layer name ("frac_outside:taxonomy_species"), so the two namespaces cannot collide.
PCT_LAYER = {"family": "taxonomy_family", "genus": "taxonomy_genus", "species": "taxonomy_species",
             "ko": "ko_humann", "pathway": "pathway_humann", "ko_eggnog": "ko_eggnog", "module": "module"}
SUMMARY_COLUMN = {"frac_outside": "frac_outside_raw", "n_missing": "n_missing",    # prefix -> column of summaries.tsv
                  "n_low": "n_low_raw", "n_high": "n_high_raw"}

# What the training used for a feature the sample does not carry. Every feature matrix in s16_health_score_B.py and in
# s8_evaluate.py is built with .fillna(0.0) and the percentile column with percentile.fillna(0.0), so an undetected
# feature, an expected-but-missing one and a detected-but-not-assessable one all entered the fit as 0, which standardizes to
# (0 - mean) / scale and is a real contribution rather than a neutral one. Applying the model has to agree.
ABSENT = 0.0

ASSEMBLY_LAYERS = {"ko_eggnog", "ko_kofam", "pfam", "module"}   # measured only by the assembly pipeline
READS_LAYERS = {"ko_humann", "pathway_humann", "taxonomy_species"}   # measured only by the read pipeline


# --------------------------------------------------------------------------------------------- the model

def _feature_kind(name: str) -> tuple[str, str, str]:
    """(kind, layer, key) for one model feature name; ValueError for a name this module cannot map onto assess output.

    kind is "pct" (key = feature_id within layer) or "summary" (key = the summaries.tsv column). The error is raised at
    load time rather than at scoring time so that a model file is rejected before any sample is read."""
    prefix, _, rest = name.partition(":")
    if not rest:
        raise ValueError(f"model feature '{name}' is not '<kind>:<id>'")
    if prefix in PCT_LAYER:
        return "pct", PCT_LAYER[prefix], rest
    if prefix in SUMMARY_COLUMN:
        return "summary", rest, SUMMARY_COLUMN[prefix]
    if prefix == "present":
        raise ValueError(f"model feature '{name}' is a species-presence feature: this is the feature set whose strongest "
                         f"coefficient is a batch artifact, and presence is keyed by MetaPhlAn species name while checkgm "
                         f"keys taxa by NCBI taxid, so it cannot be applied to an assessment at all; ship the "
                         f"percentile-only model (feature_set 'pct')")
    raise ValueError(f"model feature '{name}' has an unknown kind '{prefix}' (expected one of "
                     f"{', '.join(sorted(PCT_LAYER) + sorted(SUMMARY_COLUMN))})")


def default_model_path(pipeline: str = "reads") -> str:
    """Path of the model copy packaged beside this module for one pipeline ('reads' or 'assembly')."""
    if pipeline not in MODEL_FILES:
        raise ValueError(f"no health-score model is shipped for pipeline '{pipeline}' (expected one of {', '.join(MODEL_FILES)})")
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), MODEL_FILES[pipeline])


def model_pipeline(model: dict) -> str:
    """'reads' or 'assembly': the pipeline a model was fitted on.

    Taken from the export's "pipeline" field when it carries one, otherwise from the layers its features read: a model
    reading an assembly-only layer is an assembly model; any other is a read model, since every model exported before
    the assembly one was read-based and carries no such field."""
    declared = model.get("pipeline")
    if declared is not None:
        if declared not in MODEL_FILES:
            raise ValueError(f"health-score model declares an unknown pipeline {declared!r} (expected one of {', '.join(MODEL_FILES)})")
        return declared
    layers = {lay for _, lay, _ in model["parsed"]}
    return "assembly" if layers & ASSEMBLY_LAYERS else "reads"


def load_model(path=None, pipeline: str = "reads") -> dict:
    """A health-score model (default: the copy packaged beside this module for `pipeline`), validated.

    Returns the exported object with three keys added: "path", "parsed", the per-feature (kind, layer, key) triples in
    the order of "features", and "pipeline" ('reads' or 'assembly', see model_pipeline). Raises ValueError naming the file for anything malformed — a missing field, a non-finite or
    non-positive scale (the standardization would divide by it), a duplicated feature, a presence feature — because a
    model that is wrong in any of those ways still produces a number, and a number is what the caller will believe."""
    p = path or default_model_path(pipeline)
    if not os.path.isfile(p):
        raise FileNotFoundError(f"health-score model not found: {p}" + ("" if path else
                                " (the packaged copy is missing from this installation; reinstall checkgm or pass the "
                                "model file explicitly)"))
    try:
        with open(p) as fh:   # a leaked handle raises a ResourceWarning on stderr, and the one-line contract forbids that
            m = json.load(fh)
    except json.JSONDecodeError as e:
        raise ValueError(f"health-score model is not readable JSON: {p} ({e.msg}: line {e.lineno}, column {e.colno})") from None
    if not isinstance(m, dict):
        raise ValueError(f"health-score model is not a JSON object: {p}")
    try:
        intercept = float(m["intercept"])
    except KeyError:
        raise ValueError(f"health-score model has no 'intercept': {p}") from None
    except (TypeError, ValueError):
        raise ValueError(f"health-score model has a non-numeric 'intercept' {m['intercept']!r}: {p}") from None
    if not np.isfinite(intercept):
        raise ValueError(f"health-score model has a non-finite 'intercept': {p}")
    feats = m.get("features")
    if not isinstance(feats, list) or not feats:
        raise ValueError(f"health-score model has no 'features' list: {p}")
    parsed = []; seen = set()
    for i, f in enumerate(feats):
        if not isinstance(f, dict):
            raise ValueError(f"health-score model feature {i} is not an object: {p}")
        missing = [k for k in ("name", "coef", "mean", "scale") if k not in f]
        if missing:
            raise ValueError(f"health-score model feature {i} is missing {', '.join(missing)}: {p}")
        name = str(f["name"])
        try:
            coef, mean, scale = float(f["coef"]), float(f["mean"]), float(f["scale"])
        except (TypeError, ValueError):
            raise ValueError(f"health-score model feature '{name}' has a non-numeric coef, mean or scale: {p}") from None
        if not (np.isfinite(coef) and np.isfinite(mean) and np.isfinite(scale)):
            raise ValueError(f"health-score model feature '{name}' has a non-finite coef, mean or scale: {p}")
        if scale <= 0:
            raise ValueError(f"health-score model feature '{name}' has scale {scale!r}, which cannot standardize a value: {p}")
        if name in seen:
            raise ValueError(f"health-score model names feature '{name}' twice: {p}")
        seen.add(name)
        try:
            parsed.append(_feature_kind(name))
        except ValueError as e:
            raise ValueError(f"{e} ({p})") from None
    out = dict(m); out["intercept"] = intercept; out["path"] = p; out["parsed"] = parsed
    try:
        out["pipeline"] = model_pipeline(out)
    except ValueError as e:
        raise ValueError(f"{e}: {p}") from None
    return out


def _model(model, pipeline: str, assessed_dir: str = "") -> dict:
    """The caller's model, or the shipped one for the assessment's pipeline; refused when fitted on the other pipeline.

    A model handed over as a bare exported object (plain json.load rather than load_model(), easy to do by accident) is
    re-validated on the way in."""
    if model is None:
        return load_model(pipeline=pipeline)
    if not isinstance(model, dict):
        raise ValueError("model must be the dict returned by checkgm.healthscore.load_model()")
    if "parsed" not in model:
        m = dict(model); m["parsed"] = [_feature_kind(str(f["name"])) for f in model["features"]]
        m["intercept"] = float(model["intercept"]); m["pipeline"] = model_pipeline(m)
        model = m
    have = model.get("pipeline") or model_pipeline(model)
    if have != pipeline:
        src = model.get("path") or "the supplied model"
        raise ValueError(f"{src} is a {PIPELINE_WORD[have]}-pipeline health-score model, and {assessed_dir or 'the assessment'} "
                         f"was assessed against a {PIPELINE_WORD[pipeline]}-pipeline baseline; its coefficients are defined "
                         f"on the other pipeline's percentiles and mean nothing here, so no number is produced. Omit the "
                         f"model to use the shipped {PIPELINE_WORD[pipeline]}-pipeline one, or pass a "
                         f"{PIPELINE_WORD[pipeline]}-pipeline model")
    return model


# --------------------------------------------------------------------------------------------- the assessment

def _read_assessment(assessed_dir: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(scores, summaries) of an assess output directory, with the two files' absence told apart from their contents."""
    sp = os.path.join(assessed_dir, "scores.parquet"); mp = os.path.join(assessed_dir, "summaries.tsv")
    if not os.path.isdir(assessed_dir):
        raise FileNotFoundError(f"assessed directory not found: {assessed_dir}")
    for q, what in ((sp, "scores.parquet"), (mp, "summaries.tsv")):
        if not os.path.isfile(q):
            raise FileNotFoundError(f"assessed directory has no {what}: {assessed_dir} (point it at the directory "
                                    f"'checkgm assess' wrote)")
    sc = pd.read_parquet(sp)
    if os.path.getsize(mp) == 0:
        raise ValueError(f"summaries.tsv is empty: {mp} (re-run 'checkgm assess'; a run in which every sample was "
                         f"refused still writes the header)")
    sm = pd.read_csv(mp, sep="\t")
    for t, what, cols in ((sc, "scores.parquet", ("analysis_id", "layer", "feature_id", "percentile")),
                          (sm, "summaries.tsv", ("analysis_id", "layer"))):
        miss = [c for c in cols if c not in t.columns]
        if miss:
            raise ValueError(f"{what} in {assessed_dir} has no column {', '.join(miss)}; it was not written by this "
                             f"version of 'checkgm assess'")
    return sc, sm


def assessment_pipeline(assessed_dir: str, scores: pd.DataFrame, summaries: pd.DataFrame) -> str:
    """'reads' or 'assembly' for an assessment; ValueError when it identifies as neither or as both.

    Two independent signals, because either can be absent: the bundle id stamped in summaries.tsv, and the layers that
    were scored (the assembly baseline has no species layer and the read baseline no eggNOG/KOfam/Pfam/module one). An
    assessment that matches neither is refused — scoring it with either model would be a guess, and the one thing worse
    than refusing an assessment that failed to identify itself is scoring it with the other pipeline's coefficients."""
    bids = {str(b) for b in summaries.get("bundle_id", pd.Series(dtype=str)).dropna().unique()}
    layers = set(scores["layer"].astype(str)) | set(summaries["layer"].astype(str))
    assembly = {b for b in bids if b.startswith("gut-assembly")} or (layers & ASSEMBLY_LAYERS)
    reads = {b for b in bids if b.startswith("gut-reads")} or (layers & READS_LAYERS)
    if reads and not assembly:
        return "reads"
    if assembly and not reads:
        return "assembly"
    what = (f"against {', '.join(sorted(bids))}" if bids else f"on layers {', '.join(sorted(layers)) or 'none'}")
    why = ("both a read-based and an assembly-based baseline" if assembly and reads else
           "a baseline this module cannot identify as read-based or assembly-based")
    raise ValueError(f"a health score is shipped for the read pipeline (gut-reads-* baselines) and for the assembly "
                     f"pipeline (gut-assembly-* baselines), and {assessed_dir} was assessed {what}, {why}; each model's "
                     f"features are percentiles against its own pipeline's baseline, so no number is produced rather "
                     f"than one that cannot mean what it says. Assess the sample against a gut-reads-* or gut-assembly-* "
                     f"bundle, then score that result")


def _load(assessed_dir: str, model) -> tuple[pd.DataFrame, np.ndarray, dict]:
    """(values, present, model) for an assessment: its pipeline identified, the matching model chosen or checked."""
    scores, summaries = _read_assessment(assessed_dir)
    m = _model(model, assessment_pipeline(assessed_dir, scores, summaries), assessed_dir)
    V, present = matrix(scores, summaries, m, assessed_dir)
    return V, present, m


def matrix(scores: pd.DataFrame, summaries: pd.DataFrame, model: dict, assessed_dir: str = "the assessment") -> tuple[pd.DataFrame, np.ndarray]:
    """(values, present) for the model's features: both n_samples x n_features, indexed by analysis_id.

    `scores` and `summaries` are the two tables 'checkgm assess' writes (scores.parquet, summaries.tsv); no pipeline
    check is made here, which is _load's job. `present` records which cells carried a usable value — a detected but
    not-assessable feature counts as absent, which is how the fit treated it — because an absent cell is filled with
    ABSENT and is indistinguishable afterwards from a measured 0, and absences are the thing a reader of the score needs
    counted."""
    names = [str(f["name"]) for f in model["features"]]
    samples = pd.Index(pd.unique(summaries["analysis_id"].astype(str)), name="analysis_id")
    if not len(samples):
        raise ValueError(f"no sample was placed in {assessed_dir} (summaries.tsv holds no rows; the samples are in "
                         f"rejections.tsv with the reason each could not be placed)")
    V = pd.DataFrame(np.nan, index=samples, columns=names, dtype=float)

    want = {(lay, key): n for n, (kind, lay, key) in zip(names, model["parsed"]) if kind == "pct"}
    if want:
        sc = scores[["analysis_id", "layer", "feature_id", "percentile"]].copy()
        sc["analysis_id"] = sc["analysis_id"].astype(str); sc["layer"] = sc["layer"].astype(str)
        sc["name"] = [want.get((l, str(f))) for l, f in zip(sc["layer"], sc["feature_id"])]
        sc = sc[sc["name"].notna()]
        if len(sc):
            # One row per sample x layer x feature in a well-formed assessment, but "first" keeps a duplicated input
            # from turning into a mean of two percentiles, which is what the training pivot did as well.
            W = sc.groupby(["analysis_id", "name"], sort=False)["percentile"].first().unstack()
            V.update(W.reindex(index=samples, columns=[c for c in names if c in W.columns]))

    sumwant = [(n, lay, col) for n, (kind, lay, col) in zip(names, model["parsed"]) if kind == "summary"]
    if sumwant:
        sm = summaries.copy(); sm["analysis_id"] = sm["analysis_id"].astype(str); sm["layer"] = sm["layer"].astype(str)
        for n, lay, col in sumwant:
            if col not in sm.columns:
                continue
            g = sm[sm["layer"] == lay]
            if len(g):
                V[n] = pd.to_numeric(g.groupby("analysis_id", sort=False)[col].first(), errors="coerce").reindex(samples).to_numpy()

    present = V.notna().to_numpy()
    if not present.any():
        raise ValueError(f"none of the model's {len(names)} features was measured in {assessed_dir}; the assessment and "
                         f"the model do not share a feature space, so there is nothing to score (the model was fitted on "
                         f"{model.get('feature_set', 'its own')} features of the "
                         f"{PIPELINE_WORD.get(model.get('pipeline') or model_pipeline(model), 'read')} pipeline)")
    return V.fillna(ABSENT), present


def _contributions(V: pd.DataFrame, model: dict) -> tuple[np.ndarray, np.ndarray]:
    """(contributions, scores): coef * (value - mean) / scale per cell, and the intercept plus each row's sum."""
    coef = np.array([float(f["coef"]) for f in model["features"]])
    mean = np.array([float(f["mean"]) for f in model["features"]])
    scale = np.array([float(f["scale"]) for f in model["features"]])
    with np.errstate(all="ignore"):   # scale > 0 is validated, so this only guards against a caller-built model
        C = (V.to_numpy(dtype=float) - mean) / scale * coef
    return C, float(model["intercept"]) + C.sum(axis=1)


# --------------------------------------------------------------------------------------------- public API

TOP_IN_SUMMARY = 5   # features named per sample in score_assessment's readable column; explain() gives the full ranking


def model_for(assessed_dir, model=None) -> dict:
    """The model score_assessment would apply to this assessment: the shipped one for its pipeline, or `model` (a path
    or a load_model() dict) once checked to belong to that pipeline. ValueError otherwise."""
    if isinstance(model, (str, os.PathLike)):
        model = load_model(str(model))
    scores, summaries = _read_assessment(assessed_dir)
    return _model(model, assessment_pipeline(assessed_dir, scores, summaries), assessed_dir)


def score_assessment(assessed_dir, model=None) -> pd.DataFrame:
    """Score every sample an assessment placed: one row per sample with the score and the features that moved it most.

    `model` defaults to the shipped model of the assessment's pipeline (read or assembly); a model given explicitly must
    have been fitted on that pipeline.

    Columns: analysis_id (the sample id, named as every other checkgm table names it, so that this merges with
    summaries.tsv), score (the decision function; higher = more case-like), n_features_present / n_features_absent over
    the model's features, and top_features, the TOP_IN_SUMMARY largest contributions as '<feature> <signed amount>'.
    The absence count belongs next to the score: a sample assessed against a bundle that carries few of the model's
    features still gets a number, and the number is then mostly the model's view of what the sample does not have."""
    V, present, m = _load(assessed_dir, model)
    C, s = _contributions(V, m)
    names = np.array(V.columns)
    order = np.argsort(-np.abs(C), axis=1)[:, :TOP_IN_SUMMARY]
    top = ["; ".join(f"{names[j]} {C[i, j]:+.3f}" for j in order[i]) for i in range(len(V))]
    return pd.DataFrame({"analysis_id": V.index.to_numpy(), "score": s,
                         "n_features_present": present.sum(axis=1), "n_features_absent": (~present).sum(axis=1),
                         "top_features": top})


def explain(assessed_dir, sample, model=None, top=15) -> pd.DataFrame:
    """Per-feature contributions for one sample, largest absolute contribution first.

    Columns: feature, layer, feature_id (the summaries.tsv column for a summary feature), present, value, mean, scale,
    z, coef, contribution. `top` rows, or all of them for top=None or top <= 0. The score is the intercept plus the sum
    of every contribution, not of the rows shown, so both are attached as .attrs['intercept'] and .attrs['score']."""
    V, present, m = _load(assessed_dir, model)
    sid = str(sample)
    if sid not in V.index:
        have = list(V.index[:10])
        raise ValueError(f"sample '{sid}' was not placed in {assessed_dir}; it holds "
                         f"{', '.join(have)}{', ...' if len(V.index) > 10 else ''} "
                         f"({len(V.index)} sample(s), and any sample that could not be placed is in rejections.tsv)")
    i = V.index.get_loc(sid)
    C, s = _contributions(V, m)
    mean = np.array([float(f["mean"]) for f in m["features"]]); scale = np.array([float(f["scale"]) for f in m["features"]])
    v = V.to_numpy(dtype=float)[i]
    out = pd.DataFrame({"feature": V.columns.to_numpy(), "layer": [lay for _, lay, _ in m["parsed"]],
                        "feature_id": [key for _, _, key in m["parsed"]], "present": present[i],
                        "value": v, "mean": mean, "scale": scale, "z": (v - mean) / scale,
                        "coef": [float(f["coef"]) for f in m["features"]], "contribution": C[i]})
    out = out.iloc[np.argsort(-np.abs(out["contribution"].to_numpy()), kind="stable")].reset_index(drop=True)
    if top is not None and top > 0:
        out = out.head(int(top)).copy()
    out.attrs["intercept"] = float(m["intercept"]); out.attrs["score"] = float(s[i])
    return out
