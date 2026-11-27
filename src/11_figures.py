"""Figures for the exploratory graph analyses (reads outputs/net_*.csv; writes outputs/figures/*.png and *.pdf).

Fig A  Measurement DAG with Part A ratings (rating printed in each node; fill is an ordinal single-hue ramp).
Fig B  Panel similarity (Jaccard) between surveillance units, ordered by organism group, organism and network.
Fig C  Resistance-profile dendrograms (E. coli; staphylococci); MST with bootstrap support is tabulated.
Fig D  Human-minus-animal differences with Newcombe 95% intervals and the method-effect yardstick band.

Categorical palette (human, fisheries, livestock) validated with the dataviz validator (light mode: all checks pass;
orange below 3:1 contrast, so every mark is direct-labelled or sits in a labelled row).

Usage:  python -I src/11_figures.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Patch

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "outputs")
FIG = os.path.join(OUT, "figures")
INK, INK2, GRID = "#191919", "#666666", "#e6e6e6"
GROUP = {"human": "#0072B2", "fisheries": "#E69F00", "livestock": "#009E73"}
RAMP = {"aligned": "#ffffff", "partly aligned": "#fde0c5", "not reported": "#f4a582", "not aligned": "#d6604d"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.spines.top": False,
                     "axes.spines.right": False})


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def box(ax, xy, text, fc, ec=INK2, w=0.2, h=0.055, fs=7.5, bold=False):
    x, y = xy
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.004,rounding_size=0.01",
                                fc=fc, ec=ec, lw=0.8, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, zorder=4, weight="bold" if bold else "normal")


def fig_dag():
    e = pd.read_csv(os.path.join(OUT, "net_dag_edges.csv"))
    n = pd.read_csv(os.path.join(OUT, "net_dag_nodes.csv"))
    rating = dict(zip(n.node, n.rating))
    meas = [m for m in n.node if m.endswith("(human)")]
    pos = {"Antimicrobial use (human)": (0.36, 0.95), "Antimicrobial use (animal)": (0.64, 0.95),
           "Shared environment (water, soil, waste)": (0.5, 0.83), "Food chain and animal contact": (0.5, 0.71),
           "True resistance (human)": (0.38, 0.58), "True resistance (animal)": (0.62, 0.58),
           "Reported % (human)": (0.38, 0.16), "Reported % (animal)": (0.62, 0.16),
           "Reported difference (human - animal)": (0.5, 0.03)}
    for i, m in enumerate(meas):
        y = 0.86 - i * 0.092
        pos[m] = (0.1, y)
        pos[m.replace("(human)", "(animal)")] = (0.9, y)
    fig, ax = plt.subplots(figsize=(11, 8.5))
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.07, 1.0); ax.axis("off")
    for _, r in e.iterrows():
        (x1, y1), (x2, y2) = pos[r["from"]], pos[r["to"]]
        mline = r["from"] in rating
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>", color="#c0392b" if mline else INK2, lw=0.7 if mline else 1.0,
                                    alpha=0.55 if mline else 0.9, shrinkA=14, shrinkB=14), zorder=1)
    for node, (x, y) in pos.items():
        if node in rating:
            box(ax, (x, y), node.rsplit(" (", 1)[0] + f"\n[{rating[node]}]", RAMP[rating[node]], w=0.19, h=0.07, fs=7)
        elif node.startswith("True"):
            box(ax, (x, y), node, "#dbe9f6", w=0.2, bold=True)
        elif node.startswith("Reported"):
            box(ax, (x, y), node, "#f2f2f2", w=0.3 if "difference" in node else 0.18, bold=True)
        else:
            box(ax, (x, y), node, "#ffffff", w=0.3 if "environment" in node else 0.22)
    ax.text(0.1, 0.97, "Human surveillance (AMRSN)", ha="center", fontsize=9, weight="bold")
    ax.text(0.9, 0.97, "Animal surveillance (INFAAR)", ha="center", fontsize=9, weight="bold")
    ax.legend(handles=[Patch(fc=RAMP[k], ec=INK2, label=k) for k in ("partly aligned", "not reported", "not aligned")],
              title="Part A rating of the measurement step", loc="lower left", bbox_to_anchor=(0.0, 0.02),
              frameon=False, fontsize=8, title_fontsize=8)
    ax.text(0.5, -0.06, "Red arrows: measurement paths to the reported difference that bypass true resistance. "
            "Common causes link the two reservoirs without any arrow between them.", ha="center", fontsize=7.5, color=INK2)
    save(fig, "figA_measurement_dag")


def fig_jaccard():
    p = pd.read_csv(os.path.join(OUT, "net_panel_pairs.csv"))
    e = pd.read_csv(os.path.join(OUT, "net_panel_edges.csv"))
    units = sorted(set(e.a) | set(e.b) | set(p.unit_a) | set(p.unit_b))
    grp = lambda u: (0 if any(o in u for o in ("S. aureus", "CoNS")) else 1, u.split(" | ")[1], u.split(" | ")[0], u)
    units = sorted(units, key=grp)
    M = pd.DataFrame(np.nan, index=units, columns=units)
    for _, r in e.iterrows():
        M.loc[r.a, r.b] = M.loc[r.b, r.a] = r.jaccard
    for u in units:
        M.loc[u, u] = 1.0
    M = M.fillna(0)
    lab = [u.replace("all specimens except urine and faeces", "non-urine").replace("all specimens", "all")
           .replace(" | ", " · ") for u in units]
    fig, ax = plt.subplots(figsize=(9.5, 8))
    im = ax.imshow(M.values, cmap="Blues", vmin=0, vmax=1)
    ax.set_xticks(range(len(units)), lab, rotation=55, ha="right", fontsize=7.5)
    ax.set_yticks(range(len(units)), lab, fontsize=7.5)
    for i in range(len(units)):
        for j in range(len(units)):
            v = M.values[i, j]
            if i != j and v > 0:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.5, color="white" if v > 0.55 else INK)
    ax.spines[:].set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Panel similarity (Jaccard: shared drugs / all drugs)", fontsize=8)
    ax.set_title("Antimicrobial panels of surveillance units", loc="left", fontsize=10, weight="bold")
    save(fig, "figB_panel_jaccard")


def fig_mst():
    """Average-linkage dendrograms of the profile distances (the MST and its bootstrap support are tabulated in
    outputs/net_profile_mst_*.csv)."""
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), gridspec_kw={"wspace": 0.75})
    for ax, key, title in ((axes[0], "ecoli", "E. coli: amikacin, cefotaxime, imipenem, FQ"),
                           (axes[1], "staphylococci", "Staphylococci: PEN, ERY, LIN, TET, COT")):
        D = pd.read_csv(os.path.join(OUT, f"net_profile_distance_{key}.csv"), index_col=0)
        val = pd.read_csv(os.path.join(OUT, f"net_profile_values_{key}.csv")).fillna("")
        g = dict(zip(val.node, val.group))
        flag = dict(zip(val.node, val.n_flag))
        Zl = linkage(squareform(D.values, checks=False), method="average")
        dn = dendrogram(Zl, labels=[n + (" *" if flag[n] else "") for n in D.index], orientation="left", ax=ax,
                        color_threshold=0, above_threshold_color=INK2, leaf_font_size=8)
        for t in ax.get_yticklabels():
            name = t.get_text().rstrip(" *")
            t.set_color(INK)
            ax.plot([-1.5], [t.get_position()[1]], marker="o", ms=7, color=GROUP[g[name]], clip_on=False, transform=ax.transData)
        ax.tick_params(axis="y", pad=14, length=0)
        ax.set_xlabel("Average distance between resistance profiles (percentage points)")
        ax.set_title(title, loc="left", fontsize=9, weight="bold")
        ax.spines[["left"]].set_visible(False)
        ax.grid(axis="x", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
    fig.legend(handles=[Patch(color=c, label=k) for k, c in GROUP.items()], loc="lower center", bbox_to_anchor=(0.5, -0.1),
               ncol=3, frameon=False, fontsize=8)
    fig.text(0.01, -0.14, "Euclidean distance between resistance profiles; human: % non-susceptible, animal: % resistant. "
             "* isolate count not reported. Minimum spanning trees with bootstrap support are tabulated separately.",
             fontsize=7, color=INK2)
    save(fig, "figC_profile_dendrogram")


def fig_diff():
    d = pd.read_csv(os.path.join(OUT, "net_differences.csv"))
    y = pd.read_csv(os.path.join(OUT, "net_yardstick.csv"))
    ys = y.effect_pp.max()
    d = d[(d.human_specimen != "urine") & (d.organism != "Vibrio spp.")].copy()
    d["row"] = d.organism + " · " + np.where(d.human_drug == d.animal_drug, d.human_drug, d.human_drug + " vs " + d.animal_drug)
    order = list(dict.fromkeys(d.sort_values(["organism"], key=lambda s: s.map({"E. coli": 0, "Aeromonas spp.": 1, "S. aureus": 2, "CoNS": 3})).row))
    fig, ax = plt.subplots(figsize=(9, 0.32 * len(order) + 1.5))
    ax.axvspan(-ys, ys, color="#eeeeee", zorder=0)
    ax.axvline(0, color=INK2, lw=0.8, zorder=1)
    for i, r in enumerate(order):
        g = d[d.row == r].reset_index(drop=True)
        offs = np.linspace(-0.28, 0.28, len(g)) if len(g) > 1 else [0]
        for k, (_, x) in enumerate(g.iterrows()):
            yy = i + offs[k]
            if pd.notna(x.diff_ci_low):
                ax.plot([x.diff_ci_low, x.diff_ci_high], [yy, yy], color=GROUP[x.sector], lw=1, zorder=2)
            ax.scatter([x.difference_pp], [yy], s=16, color=GROUP[x.sector], zorder=3,
                       marker="o" if pd.notna(x.diff_ci_low) else "D")
    ax.set_yticks(range(len(order)), order, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("Human % non-susceptible minus animal % resistant (percentage points)")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(handles=[Patch(color=GROUP["fisheries"], label="fisheries source"), Patch(color=GROUP["livestock"], label="livestock source"),
                       Patch(color="#eeeeee", label=f"method-effect yardstick (±{ys:.1f} pp)")],
              loc="upper center", bbox_to_anchor=(0.45, -0.06), ncol=3, frameon=False, fontsize=7.5)
    ax.set_title("Human minus animal resistance, 2019-2022 (Newcombe 95% intervals; diamonds: no isolate count)",
                 loc="left", fontsize=9, weight="bold")
    save(fig, "figD_differences")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    fig_dag(); fig_jaccard(); fig_mst(); fig_diff()
    print("\n".join(sorted(os.listdir(FIG))))


if __name__ == "__main__":
    main()
