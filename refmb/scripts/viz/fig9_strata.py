"""Figure 9 — stratified baseline (double column). Out-of-sample fraction of features outside the reference 2.5–97.5
band for healthy adults, by age, sex, BMI, westernization and region. a: pipeline A (assembly), b: pipeline B (reads).
Two layers per pipeline; the 8 % gate is drawn; n samples (studies) is written beside each stratum."""
import os, sys
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import tokens as T

P = T.P
GROUPS = [("age_bin", "age", ["18-39", "40-64", "65+"]), ("sex", "sex", ["female", "male"]), ("bmi_class", "BMI", ["<18.5", "18.5-25", "25-30", ">=30"]),
          ("westernized", "westernized", ["yes", "no"]), ("region", "region", ["Europe", "N. America", "E. Asia", "W. Asia", "C. Asia", "S. Asia", "Oceania", "Africa"])]
LAYERS = {"A": [("taxonomy_family", "family"), ("ko_eggnog", "KO")], "B": [("taxonomy_family", "family"), ("ko_humann", "KO")]}
TITLE = {"A": "a  assembly pipeline (A)", "B": "b  read pipeline (B)"}


def make(out_dir):
    T.print_profile()
    fig = plt.figure(figsize=(T.DOUBLE_IN, 4.3))
    gs = fig.add_gridspec(1, 2, wspace=0.62, left=0.13, right=0.9, top=0.93, bottom=0.14)
    for k, pl in enumerate("AB"):
        ax = fig.add_subplot(gs[k]); y = 0; ticks = []; labels = []; seps = []
        for col, gname, order in GROUPS:
            d = pd.read_csv(f"{P}/results/s12/strata_{pl}_{col}.tsv", sep="\t"); d[col] = d[col].astype(str)
            for s in order:
                r = d[d[col] == s]
                if r.empty:
                    continue
                for j, (layer, _) in enumerate(LAYERS[pl]):
                    v = r[r["layer"] == layer]
                    if len(v):
                        ax.scatter(100 * v["median_frac_outside"].iloc[0], y + (-0.17 if j == 0 else 0.17), s=14 if j == 0 else 12, marker=T.LAYER_MARK["family" if j == 0 else "KO"], facecolor=T.INK["primary"] if j == 0 else T.INK["surface"], edgecolor=T.INK["primary"], linewidth=0.6, zorder=5, clip_on=False)
                r0 = r.iloc[0]
                ax.text(1.02, y, f"{int(r0['n']):,} ({int(r0['studies'])})", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.PT_MIN, color=T.INK["secondary"])
                ticks.append(y); labels.append(f"{gname}  {s.replace('>=', '≥')}"); y += 1
            seps.append(y - 0.5)
        for s in seps[:-1]:
            ax.axhline(s, color=T.INK["grid"], lw=0.4, zorder=1)
        ax.axvline(8, color=T.INK["muted"], lw=0.6, ls=(0, (3, 2)), zorder=1)
        ax.text(8.7, -0.75, "8 % limit", ha="left", va="center", fontsize=T.PT_MIN, color=T.INK["secondary"])
        ax.set_yticks(ticks); ax.set_yticklabels(labels); ax.set_ylim(y - 0.4, -1.2)
        ax.set_xlim(-0.5, 21); ax.set_xlabel("median share of features outside the healthy range (%)")
        ax.grid(axis="x", zorder=0); ax.set_axisbelow(True); ax.tick_params(length=2, width=0.4)
        ax.set_title(TITLE[pl], loc="left")
        ax.text(1.02, -0.75, "n (studies)", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.PT_MIN, color=T.INK["muted"])
        h = [plt.Line2D([], [], marker=T.LAYER_MARK[lab], ls="", ms=3.6, markerfacecolor=T.INK["primary"] if lab == "family" else T.INK["surface"], markeredgecolor=T.INK["primary"], markeredgewidth=0.6, label=lab) for _, lab in LAYERS[pl]]
        ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=2, frameon=False, handletextpad=0.2, columnspacing=1.2, borderaxespad=0)
    T.save(fig, "fig9_strata", out_dir); plt.close(fig)


if __name__ == "__main__":
    make(sys.argv[1] if len(sys.argv) > 1 else f"{P}/figures")
