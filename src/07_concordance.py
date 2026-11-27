"""Concordance of resistance in human clinical isolates (AMRSN) and animal and fishery isolates (INFAAR), 2019-2022
(plan section 4.4; definitions in docs/DEVIATIONS.md C1-C7).

For each shared organism-antimicrobial pair: AMRSN non-susceptible % by year and pooled 2019-2022 (sum of
non-susceptible over sum tested, Wilson 95% CI); INFAAR resistant % by source (pooled 2019-2022 in the report), with an
approximate Wilson CI where an isolate count is reported; difference human minus animal in percentage points. Also the
rank order of the shared antimicrobials within each network (Spearman rho, descriptive; no test with 4-6 drugs). Descriptive only: the two networks sample
different populations (patients' clinical isolates against healthy animals, food and fish), so differences in level
are not evidence of transmission.

Usage:  python -I src/07_concordance.py      writes outputs/concordance.csv, outputs/amrsn_series.csv,
                                               outputs/concordance_ranks.csv, outputs/concordance_fig8_by_period.csv
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "outputs")
Z = 1.959964

# Isolates per INFAAR source (DEVIATIONS C7): fisheries figure legends; livestock text (p25-27); S. aureus p23.
N_INFAAR = {("Figure 2", "Freshwater fish"): 35, ("Figure 2", "Shrimps"): 193, ("Figure 2", "Marine fish"): 48,
            ("Figure 4", "Freshwater fish"): 1441, ("Figure 4", "Shrimps"): 481, ("Figure 4", "Marine fish"): 157,
            ("Figure 11", "Cattle"): 420, ("Figure 11", "Buffalo"): 125, ("Figure 11", "Goat"): 434,
            ("Figure 11", "Poultry"): 722, ("Figure 10A", "Food animals"): 452}
SECTOR = {"Freshwater fish": "fisheries", "Shrimps": "fisheries", "Marine fish": "fisheries", "Food animals": "livestock",
          "Cattle": "livestock", "Buffalo": "livestock", "Goat": "livestock", "Sheep": "livestock", "Poultry": "livestock",
          "Pigs": "livestock"}
# (organism, human antimicrobial, animal antimicrobial, comparison, INFAAR figures)
PAIRS = [
    ("E. coli", "amikacin", "amikacin", "same drug", ["Figure 4", "Figure 11"]),
    ("E. coli", "cefotaxime", "cefotaxime", "same drug", ["Figure 4", "Figure 11"]),
    ("E. coli", "ceftazidime", "ceftazidime", "same drug", ["Figure 4"]),
    ("E. coli", "imipenem", "imipenem", "same drug", ["Figure 4", "Figure 11"]),
    ("E. coli", "ciprofloxacin", "enrofloxacin", "fluoroquinolone class", ["Figure 4", "Figure 11"]),
    ("E. coli", "ceftriaxone", "ceftriaxone", "same drug (human: urine, 2022 only)", ["Figure 4", "Figure 11"]),
    ("E. coli", "co-trimoxazole", "co-trimoxazole", "same drug (human: urine)", ["Figure 4", "Figure 11"]),
    ("S. aureus", "cefoxitin", "cefoxitin", "methicillin resistance, phenotypic", ["Figure 2"]),
    ("S. aureus", "penicillin", "penicillin", "same drug", ["Figure 2", "Figure 10A"]),
    ("S. aureus", "erythromycin", "erythromycin", "same drug", ["Figure 2", "Figure 10A"]),
    ("S. aureus", "linezolid", "linezolid", "same drug", ["Figure 2", "Figure 10A"]),
    ("S. aureus", "tetracycline", "tetracycline", "same drug", ["Figure 2", "Figure 10A"]),
    ("S. aureus", "co-trimoxazole", "co-trimoxazole", "same drug", ["Figure 2", "Figure 10A"]),
    ("S. aureus", "ciprofloxacin", "ciprofloxacin", "same drug", ["Figure 2"]),
    ("S. aureus", "ciprofloxacin", "enrofloxacin", "fluoroquinolone class", ["Figure 10A"]),
]


def wilson(k, n):
    if not n:
        return (np.nan, np.nan)
    p = k / n
    den = 1 + Z ** 2 / n
    c = (p + Z ** 2 / (2 * n)) / den
    h = Z * np.sqrt(p * (1 - p) / n + Z ** 2 / (4 * n ** 2)) / den
    return (round(100 * (c - h), 1), round(100 * (c + h), 1))


def amrsn_series():
    """Rows: organism, specimen, antimicrobial, year, n_ns, n_tested, pct_ns, source, flag."""
    p = pd.read_csv(os.path.join(OUT, "amrsn_primary.csv"))
    p["n_ns"] = p.n_tested - p.n_susceptible
    p["source"] = "AMRSN " + p.report.astype(str) + " Table " + p.table.astype(str) + " p" + p.page.astype(str)
    p["flag"] = p.flag.fillna("")
    rows = [p[["organism", "specimen", "antimicrobial", "year", "n_ns", "n_tested", "pct_ns", "source", "flag"]]]
    # urine: 2020-2021 image "Overall" rows (whole-number %S), 2022 text Table 2.13 (OPD + ward + ICU)
    img = pd.read_csv(os.path.join(REPO, "docs", "amrsn_urine_image_transcription.csv"))
    img = img[img.centre == "Overall"]
    rows.append(pd.DataFrame({"organism": "E. coli", "specimen": "urine", "antimicrobial": img.antimicrobial,
                              "year": img.report, "n_ns": (img.n_tested * (100 - img.pct_s_rounded) / 100).round().astype(int),
                              "n_tested": img.n_tested, "pct_ns": 100 - img.pct_s_rounded,
                              "source": "AMRSN " + img.report.astype(str) + " Table " + img.table.astype(str) + " p" + img.page.astype(str),
                              "flag": "image-read; %S printed as whole number"}))
    lg = pd.read_csv(os.path.join(OUT, "amrsn_long.csv"))
    u = lg[(lg.specimen == "urine") & (lg.report == 2022)].copy()
    u["antimicrobial"] = u.antimicrobial_raw.str.replace(r"\s+", "", regex=True).str.lower().replace(
        {"trimethoprim-sulfamethoxazole": "co-trimoxazole"})
    g = u.groupby("antimicrobial")[["n_susceptible", "n_tested"]].sum().reset_index()
    rows.append(pd.DataFrame({"organism": "E. coli", "specimen": "urine", "antimicrobial": g.antimicrobial, "year": 2022,
                              "n_ns": g.n_tested - g.n_susceptible, "n_tested": g.n_tested,
                              "pct_ns": (100 * (g.n_tested - g.n_susceptible) / g.n_tested).round(1),
                              "source": "AMRSN 2022 Table 2.13 p62 (OPD+ward+ICU)", "flag": ""}))
    s = pd.concat(rows, ignore_index=True)
    return s[s.year.between(2019, 2022)]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    pd.set_option("display.width", 220)
    s = amrsn_series()
    s.to_csv(os.path.join(OUT, "amrsn_series.csv"), index=False)
    v = pd.read_csv(os.path.join(OUT, "infaar_values.csv")).fillna({"flag": ""})
    out = []
    for org, hdrug, adrug, comp, figs in PAIRS:
        for spec in (["urine"] if "urine" in comp else (["all specimens except urine and faeces", "urine"] if org == "E. coli" else ["all specimens"])):
            h = s[(s.organism == org) & (s.antimicrobial == hdrug) & (s.specimen == spec)]
            if h.empty:
                continue
            k, n = int(h.n_ns.sum()), int(h.n_tested.sum())
            lo, hi = wilson(k, n)
            hum = {"organism": org, "human_drug": hdrug, "animal_drug": adrug, "comparison": comp, "human_specimen": spec,
                   "human_years": f"{h.year.min()}-{h.year.max()}" if h.year.nunique() > 1 else str(h.year.iloc[0]),
                   "human_pct_ns": round(100 * k / n, 1), "human_ci": f"{lo}-{hi}", "human_n": n,
                   "human_by_year": "; ".join(f"{int(r.year)}: {r.pct_ns:.1f}" for _, r in h.sort_values("year").iterrows()),
                   "human_flag": "; ".join(sorted(set(f for f in h.flag if f)))}
            a = v[v.figure.isin(figs) & (v.organism.isin(["E. coli", "S. aureus"])) & (v.antimicrobial == adrug)]
            a = a[a.organism == org]
            for _, r in a.iterrows():
                na = N_INFAAR.get((r.figure, r.source))
                alo, ahi = wilson(round(r.value_pct * na / 100), na) if na else (np.nan, np.nan)
                out.append({**hum, "animal_source": r.source, "sector": SECTOR[r.source], "infaar_figure": r.figure,
                            "animal_pct_r": r.value_pct, "animal_ci_approx": "" if not na else f"{alo}-{ahi}",
                            "animal_n_isolates": na, "difference_pp": round(hum["human_pct_ns"] - r.value_pct, 1),
                            "animal_flag": r.flag})
    # livestock MRSA by mecA/mecC (INFAAR p24: 41 confirmed among 452 S. aureus tested)
    h = s[(s.organism == "S. aureus") & (s.antimicrobial == "cefoxitin")]
    k, n = int(h.n_ns.sum()), int(h.n_tested.sum())
    lo, hi = wilson(k, n)
    alo, ahi = wilson(41, 452)
    out.append({"organism": "S. aureus", "human_drug": "cefoxitin", "animal_drug": "mecA/mecC", "comparison":
                "methicillin resistance: human phenotypic (cefoxitin) vs animal genotypic (mecA/mecC)",
                "human_specimen": "all specimens", "human_years": "2019-2022", "human_pct_ns": round(100 * k / n, 1),
                "human_ci": f"{lo}-{hi}", "human_n": n, "human_by_year": "", "human_flag": "",
                "animal_source": "Food animals", "sector": "livestock", "infaar_figure": "text p24 (41/452)",
                "animal_pct_r": round(100 * 41 / 452, 1), "animal_ci_approx": f"{alo}-{ahi}", "animal_n_isolates": 452,
                "difference_pp": round(100 * k / n - 100 * 41 / 452, 1), "animal_flag": "genotypic definition"})
    c = pd.DataFrame(out)
    c.to_csv(os.path.join(OUT, "concordance.csv"), index=False)

    # rank order of shared same-drug antimicrobials, primary human series vs each animal source
    ranks = []
    for org in ("E. coli", "S. aureus"):
        prim = c[(c.organism == org) & (c.comparison == "same drug") & (c.human_specimen != "urine")]
        for src, g in prim.groupby("animal_source"):
            if len(g) >= 4:
                rho = spearmanr(g.human_pct_ns, g.animal_pct_r)[0]     # descriptive only: 4-6 drugs, no test
                ranks.append({"organism": org, "animal_source": src, "n_drugs": len(g), "drugs": ", ".join(g.human_drug),
                              "spearman_rho": round(rho, 2)})
    pd.DataFrame(ranks).to_csv(os.path.join(OUT, "concordance_ranks.csv"), index=False)

    # Figure 8: aquaculture E. coli by period, against AMRSN by year (descriptive, no denominators)
    f8 = v[v.figure == "Figure 8"]
    yr = {"Oct 2018-Mar 2020": "2019", "Apr 2020-Mar 2021": "2020", "Apr 2021-Mar 2022": "2021"}
    rows = []
    for d8 in ("amikacin", "cefotaxime", "ceftazidime", "imipenem"):
        for per, y in yr.items():
            x = f8[(f8.antimicrobial == d8) & (f8.source == per)]
            hh = s[(s.organism == "E. coli") & (s.antimicrobial == d8) & (s.specimen == "all specimens except urine and faeces") & (s.year == int(y))]
            if len(x) and len(hh):
                rows.append({"antimicrobial": d8, "infaar_period": per, "infaar_pct_r": float(x.value_pct.iloc[0]),
                             "amrsn_year": int(y), "amrsn_pct_ns": float(hh.pct_ns.iloc[0])})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "concordance_fig8_by_period.csv"), index=False)

    show = c[c.human_specimen != "urine"]
    print(show.pivot_table(index=["organism", "human_drug", "animal_drug", "human_pct_ns"], columns="animal_source",
                           values="animal_pct_r").round(1).to_string())
    print("\nUrine series (human):")
    print(c[c.human_specimen == "urine"].pivot_table(index=["human_drug", "human_years", "human_pct_ns"], columns="animal_source",
                                                     values="animal_pct_r").round(1).to_string())
    print("\nRank concordance:")
    print(pd.DataFrame(ranks).to_string(index=False))
    print("\nFigure 8 by period:")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
