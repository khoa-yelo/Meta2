"""Figure 2 — what the baselines are made of (double column). a: assembly baseline by study (samples, stacked by
assembly-size class as three ordinal blue steps). b: by country (filled bar = assembly pipeline; countries from the curated
sample table work/s12/A/meta.parquet, reference_pool role, ISO codes mapped to names so that e.g. TZA and Tanzania are one
bar). c: metadata completeness of the gut sample set from the curated metadata vs after the cMD3 accession join (two greys).
d, e: the read baseline (cMD3 MetaPhlAn 3 profiles) by study, stacked by read-depth class, and by country (outlined bar =
read pipeline)."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T
sys.path.insert(0, T.P); from refmb.paths import bundle_A as _bundle_A

P = T.P
BUNDLE = _bundle_A("lenient")
BAND_COL = {"low": T.SEQ[0], "medium": T.SEQ[2], "high": T.SEQ[4]}
ISO = {"NLD": "Netherlands", "GBR": "United Kingdom", "ISR": "Israel", "DNK": "Denmark", "CHN": "China", "USA": "United States", "DEU": "Germany", "FRA": "France", "ESP": "Spain", "ITA": "Italy", "SWE": "Sweden",
       "AUT": "Austria", "FIN": "Finland", "IND": "India", "JPN": "Japan", "CAN": "Canada", "AUS": "Australia", "KAZ": "Kazakhstan", "FJI": "Fiji", "TZA": "Tanzania", "MDG": "Madagascar", "PER": "Peru", "SLV": "El Salvador", "HUN": "Hungary", "IRL": "Ireland", "LUX": "Luxembourg", "NOR": "Norway", "RUS": "Russia", "SGP": "Singapore", "KOR": "South Korea", "BGD": "Bangladesh", "MNG": "Mongolia", "LBR": "Liberia", "CMR": "Cameroon", "GHA": "Ghana", "ETH": "Ethiopia", "SVK": "Slovakia", "EST": "Estonia", "ISL": "Iceland", "BRN": "Brunei", "MYS": "Malaysia", "PHL": "Philippines", "IDN": "Indonesia", "THA": "Thailand", "VNM": "Vietnam", "PAK": "Pakistan", "NPL": "Nepal", "SAU": "Saudi Arabia", "ARE": "UAE", "EGY": "Egypt", "MAR": "Morocco", "ZAF": "South Africa", "NGA": "Nigeria", "KEN": "Kenya", "UGA": "Uganda", "BFA": "Burkina Faso", "MWI": "Malawi", "BRA": "Brazil", "ARG": "Argentina", "COL": "Colombia", "MEX": "Mexico", "CHL": "Chile", "VEN": "Venezuela", "ECU": "Ecuador", "BOL": "Bolivia", "GRL": "Greenland", "NZL": "New Zealand", "PNG": "Papua New Guinea", "CZE": "Czechia", "POL": "Poland", "PRT": "Portugal", "GRC": "Greece", "TUR": "Turkey", "ROU": "Romania", "BGR": "Bulgaria", "BEL": "Belgium", "CHE": "Switzerland", "LTU": "Lithuania", "LVA": "Latvia", "HRV": "Croatia", "SRB": "Serbia", "SVN": "Slovenia", "UKR": "Ukraine", "BLR": "Belarus", "GEO": "Georgia", "ARM": "Armenia", "AZE": "Azerbaijan", "UZB": "Uzbekistan", "KGZ": "Kyrgyzstan", "TJK": "Tajikistan", "TKM": "Turkmenistan", "IRN": "Iran", "IRQ": "Iraq", "JOR": "Jordan", "LBN": "Lebanon", "SYR": "Syria", "YEM": "Yemen", "OMN": "Oman", "KWT": "Kuwait", "QAT": "Qatar", "BHR": "Bahrain", "AFG": "Afghanistan", "LKA": "Sri Lanka", "MMR": "Myanmar", "KHM": "Cambodia", "LAO": "Laos", "TWN": "Taiwan", "HKG": "Hong Kong", "MAC": "Macao", "PRK": "North Korea"}


def make(out_dir):
    T.print_profile()
    inv = pd.read_parquet(f"{P}/work/s0/inventory.parquet")
    excl = set(pd.read_csv(f"{BUNDLE}/provenance/pool_exclusions.tsv", sep="\t")["analysis_id"])
    st = pd.read_csv(f"{BUNDLE}/provenance/studies.tsv", sep="\t")
    bands = pd.read_parquet(f"{P}/work/s2/assembly_quality_bands.parquet", columns=["analysis_id", "quality_band"])
    pool = inv[inv["study_bioproject"].isin(st["study_bioproject"]) & inv["is_primary_analysis"] & (inv["body_site_core"] == "Gut") & ~inv["analysis_id"].isin(excl)]
    pool = pool[pool["analysis_id"].isin(bands["analysis_id"])].merge(bands, on="analysis_id")
    fig = plt.figure(figsize=(T.DOUBLE_IN, 5.8))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.6, 1.0, 1.35], height_ratios=[1.0, 0.8], wspace=0.75, hspace=0.5, left=0.14, right=0.975, top=0.935, bottom=0.07)
    # a — by study, stacked by band
    ax = fig.add_subplot(gs[0, 0])
    tab = pool.pivot_table(index="study_bioproject", columns="quality_band", values="analysis_id", aggfunc="count", fill_value=0).reindex(columns=["low", "medium", "high"], fill_value=0)
    tab = tab.loc[tab.sum(axis=1).sort_values(ascending=False).index]
    left = np.zeros(len(tab))
    for b in ["low", "medium", "high"]:
        ax.barh(np.arange(len(tab)), tab[b], left=left, color=BAND_COL[b], height=0.72, linewidth=0.6, edgecolor=T.INK["surface"], label=f"{b} size class")
        left += tab[b].to_numpy()
    ax.set_yticks(np.arange(len(tab))); ax.set_yticklabels(tab.index); ax.invert_yaxis(); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("reference samples"); ax.set_title(f"a  assembly pipeline:\n{len(pool):,} healthy adults, {len(tab)} studies", loc="left")
    ax.legend(loc="lower right", frameon=False, handlelength=1.0, handletextpad=0.4)
    # b — by country
    ax = fig.add_subplot(gs[0, 1])
    meta = pd.read_parquet(f"{P}/work/s12/A/meta.parquet"); meta = meta[meta["role"] == "reference_pool"]
    assert len(meta) == len(pool), f"pool size differs between inventory ({len(pool)}) and work/s12/A/meta.parquet ({len(meta)})"
    c = meta["country"].fillna("").replace("", "not reported").map(lambda k: ISO.get(k, k)).value_counts()
    top = c.head(7); other = c.iloc[7:].sum()
    if other: top = pd.concat([top, pd.Series({"other": other})])
    ax.barh(np.arange(len(top)), top.to_numpy(), height=0.72, **T.BAR_ASSEMBLY)
    for i, v in enumerate(top.to_numpy()): ax.text(v + 8, i, f"{v}", va="center", fontsize=T.PT_MIN, color=T.INK["secondary"])
    ax.set_yticks(np.arange(len(top))); ax.set_yticklabels(top.index); ax.invert_yaxis(); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, top.max() * 1.25); ax.set_xlabel("reference samples"); ax.set_title("b  assembly pipeline:\nby country", loc="left")
    # c — metadata completeness: curation vs curation + cMD3 join, gut query set
    ax = fig.add_subplot(gs[0, 2])
    gut = inv[inv["is_primary_analysis"] & (inv["s0_status"] == "retained") & (inv["body_site_core"] == "Gut")]
    cmd = pd.read_csv(T.CMD3_META, dtype=str, keep_default_na=False)[["sample_id", "study_name", "age", "gender", "BMI", "antibiotics_current_use", "study_condition"]]
    g = gut.merge(cmd, left_on=["cmd3_sample_id", "cmd3_study_name"], right_on=["sample_id", "study_name"], how="left")
    fields = [("health status", g["cur_Health_status_group"] != "", (g["cur_Health_status_group"] != "") | g["study_condition"].fillna("").ne("")),
              ("age", g["cur_Age_years"] != "", (g["cur_Age_years"] != "") | g["age"].fillna("").ne("")),
              ("sex", g["cur_Sex"] != "", (g["cur_Sex"] != "") | g["gender"].fillna("").ne("")),
              ("antibiotics", g["cur_Antibiotics"] != "", (g["cur_Antibiotics"] != "") | g["antibiotics_current_use"].fillna("").ne("")),
              ("BMI", pd.Series(False, index=g.index), g["BMI"].fillna("").ne("")),
              ("country", g["cur_Country"] != "", g["cur_Country"] != "")]
    y = np.arange(len(fields))
    ax.barh(y - 0.19, [a.mean() for _, a, _ in fields], height=0.36, color=T.INK["axis"], linewidth=0, label="curated metadata")
    ax.barh(y + 0.19, [b.mean() for _, _, b in fields], height=0.36, color=T.INK["secondary"], linewidth=0, label="+ cMD3 join")
    ax.set_yticks(y); ax.set_yticklabels([f for f, _, _ in fields]); ax.invert_yaxis(); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["0", "50%", "100%"]); ax.set_xlabel(f"share of gut samples (n={len(gut):,})")
    ax.set_title("c  metadata completeness,\nall gut samples", loc="left"); ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, handlelength=1.0, handletextpad=0.4, columnspacing=0.7, borderaxespad=0)
    # d, e — pipeline B pool
    B = pd.read_csv(f"{P}/work/s11/pipeB/samples.tsv", sep="\t", low_memory=False); B = B[(B["role"] == "reference_pool") & (B["status"] == "included")]
    ax = fig.add_subplot(gs[1, 0])
    tab = B.pivot_table(index="study", columns="depth_band", values="sample_key", aggfunc="count", fill_value=0).reindex(columns=["low", "medium", "high"], fill_value=0)
    tab = tab.loc[tab.sum(axis=1).sort_values(ascending=False).index]; left = np.zeros(len(tab))
    for b in ["low", "medium", "high"]:
        ax.barh(np.arange(len(tab)), tab[b], left=left, color=BAND_COL[b], height=0.72, linewidth=0.6, edgecolor=T.INK["surface"], label=f"{b} depth class")
        left += tab[b].to_numpy()
    ax.set_yticks(np.arange(len(tab))); ax.set_yticklabels(tab.index); ax.invert_yaxis(); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("reference samples"); ax.set_title(f"d  read pipeline:\n{len(B):,} healthy adults, {len(tab)} studies", loc="left")
    ax.legend(loc="lower right", frameon=False, handlelength=1.0, handletextpad=0.4)
    ax = fig.add_subplot(gs[1, 1])
    c = B["country"].fillna("not reported").map(lambda k: ISO.get(k, k)).value_counts(); top = c.head(7); other = c.iloc[7:].sum()
    if other: top = pd.concat([top, pd.Series({"other": other})])
    ax.barh(np.arange(len(top)), top.to_numpy(), height=0.72, **T.BAR_READ)
    for i, v in enumerate(top.to_numpy()): ax.text(v + top.max() * 0.02, i, f"{v:,}", va="center", fontsize=T.PT_MIN, color=T.INK["secondary"])
    ax.set_yticks(np.arange(len(top))); ax.set_yticklabels(top.index); ax.invert_yaxis(); ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, top.max() * 1.3); ax.set_xlabel("reference samples"); ax.set_title("e  read pipeline:\nby country", loc="left")
    T.save(fig, "fig2_composition", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
