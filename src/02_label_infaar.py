"""Assign each extracted INFAAR bar to figure, organism, source and antimicrobial (plan section 4.1).

Uses an explicit chart map (below) read from the report's figure titles and legends. Within each panel, category
centres are the x-centres of the category label words; each bar goes to the nearest category, and in grouped charts
its series is its position within the category slot (left to right in legend order). A series with no bar in a
category is recorded as 0% (bars of zero height are not drawn), flagged as inferred.

Usage:  python -I src/02_label_infaar.py      writes outputs/infaar_values.csv
"""
import os

import fitz
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PDF = os.path.join(REPO, "data_raw", "infaar_2019_22.pdf")
OUT = os.path.join(REPO, "outputs")

# Canonical antimicrobial names for every label used in the report.
DRUG = {
    "AMC": "amoxicillin-clavulanate", "Amox-clav": "amoxicillin-clavulanate", "Amoxi-Clav": "amoxicillin-clavulanate",
    "AMK": "amikacin", "Amikacin": "amikacin", "AMP": "ampicillin", "Ampicillinnam": "ampicillin", "Ampicilin": "ampicillin",
    "ATM": "aztreonam", "Aztreonam": "aztreonam", "CAX": "CAX (undefined in report)", "CHL": "chloramphenicol",
    "Chloramphenicol": "chloramphenicol", "COT": "co-trimoxazole", "SXT": "co-trimoxazole", "TMP": "co-trimoxazole",
    "Co-trimoxazole": "co-trimoxazole", "CPD": "cefpodoxime", "Cefpodoxime": "cefpodoxime", "CTR": "ceftriaxone",
    "Ceftriaxone": "ceftriaxone", "CTX": "cefotaxime", "Cefotaxime": "cefotaxime", "CTZ": "ceftazidime",
    "Ceftazidime": "ceftazidime", "ENR": "enrofloxacin", "Enrofloxacin": "enrofloxacin", "FOX": "cefoxitin",
    "Cefoxitin": "cefoxitin", "IPM": "imipenem", "Imipenem": "imipenem", "NAL": "nalidixic acid", "Nalidixic": "nalidixic acid",
    "TET": "tetracycline", "Tetracycline": "tetracycline", "ERY": "erythromycin", "Erythromycin": "erythromycin",
    "GEN": "gentamicin", "Gentamicin": "gentamicin", "LIN": "linezolid", "Linezolid": "linezolid", "PEN": "penicillin",
    "Penicillin": "penicillin", "Ciprofloxacin": "ciprofloxacin", "Cefepime": "cefepime", "Meropenem": "meropenem",
}

# (page, axis_x, axis_y0) -> figure, organism, series list (legend order) or single source name.
CHARTS = [
    (16, 86.1, 696.5, "Figure 2", "S. aureus", ["Freshwater fish", "Shrimps", "Marine fish"]),
    (17, 88.2, 193.4, "Figure 3", "Coagulase-negative staphylococci", ["Freshwater fish", "Shrimps", "Marine fish"]),
    (18, 87.5, 202.2, "Figure 4", "E. coli", ["Freshwater fish", "Shrimps", "Marine fish"]),
    (19, 83.7, 198.5, "Figure 5", "V. parahaemolyticus", ["Shrimps", "Marine fish"]),
    (19, 83.7, 638.2, "Figure 6", "Vibrio species", ["Shrimps", "Marine fish"]),
    (20, 88.5, 468.9, "Figure 7", "Aeromonas species", ["Aquaculture (all)"]),
    (21, 82.9, 271.8, "Figure 8", "E. coli", ["Oct 2018-Mar 2020", "Apr 2020-Mar 2021", "Apr 2021-Mar 2022"]),
    (24, 84.9, 524.0, "Figure 10A", "S. aureus", ["Food animals"]),
    (24, 84.9, 698.1, "Figure 10B", "MRSA", ["Food animals"]),
    (24, 349.7, 523.6, "Figure 10C", "Coagulase-negative staphylococci (caption) / MRCoNS (text)", ["Food animals"]),
    (24, 349.7, 698.1, "Figure 10D", "MR-CoNS (caption) / CoNS (text)", ["Food animals"]),
    (26, 79.3, 365.4, "Figure 11", "E. coli", ["Cattle"]),
    (26, 325.7, 365.4, "Figure 11", "E. coli", ["Buffalo"]),
    (26, 79.3, 522.0, "Figure 11", "E. coli", ["Goat"]),
    (26, 325.7, 522.0, "Figure 11", "E. coli", ["Sheep"]),
    (26, 79.3, 677.3, "Figure 11", "E. coli", ["Poultry"]),
    (26, 325.7, 677.3, "Figure 11", "E. coli", ["Pigs"]),
]


def main():
    bars = pd.read_csv(os.path.join(OUT, "infaar_extracted.csv")).fillna("")
    d = fitz.open(PDF)
    rows = []
    for page, ax, ay, fig, org, series in CHARTS:
        b = bars[(bars.page == page) & (abs(bars.axis_x - ax) < 1) & (abs(bars.axis_y0 - ay) < 1)].copy()
        assert len(b), (page, ax, ay)
        others = sorted({x for p, x, y, *_ in CHARTS if p == page and abs(y - ay) < 1 and x > ax + 1})
        x_right = others[0] - 10 if others else 10_000
        words = d[page - 1].get_text("words")
        cats = {}
        for w in words:
            if w[4] in DRUG and ay < w[1] < ay + 80 and ax < (w[0] + w[2]) / 2 < x_right:
                cats.setdefault(w[4], []).append((w[0] + w[2]) / 2)
        centres = {k: float(np.mean(v)) for k, v in cats.items()}
        b["category"] = [min(centres, key=lambda k: abs(centres[k] - x)) for x in b.bar_x]
        n = len(series)
        if n > 1:
            # Grouped chart: one fill colour per series; legend order = left-to-right order of the colours within
            # the category slots (mean offset of each colour's bars from its category centre).
            b["offset"] = b.bar_x - b.category.map(centres)
            order = b.groupby("fill").offset.mean().sort_values().index.tolist()
            assert len(order) == n, f"{fig} {org}: {len(order)} fill colours for {n} series"
            b["series_idx"] = b.fill.map({f: i for i, f in enumerate(order)})
        else:
            b["series_idx"] = 0
        for cat in sorted(centres, key=centres.get):          # every labelled category, even if no bar is drawn
            g = b[b.category == cat]
            got = {}
            for _, r in g.iterrows():
                assert r.series_idx not in got, f"two bars of one series in a category: {fig} {org} {cat}"
                got[r.series_idx] = r
            for i, s in enumerate(series):
                if i in got:
                    r = got[i]
                    rows.append({"figure": fig, "page": page, "organism": org, "source": s, "label": cat, "antimicrobial": DRUG[cat],
                                 "value_pct": r.value_pct, "flag": r.flag})
                else:
                    rows.append({"figure": fig, "page": page, "organism": org, "source": s, "label": cat, "antimicrobial": DRUG[cat],
                                 "value_pct": 0.0, "flag": "no bar drawn; recorded as 0 (inferred)"})
        dup = pd.DataFrame(rows)
        dup = dup[(dup.figure == fig) & (dup.organism == org)]
        clash = dup.groupby(["source", "label"]).size()
        assert (clash == 1).all(), f"two bars in one slot: {fig} {org}\n{clash[clash > 1]}"
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "infaar_values.csv"), index=False)
    print(out.groupby(["figure", "organism", "source"]).size().to_string())
    print(f"\n{len(out)} values; {int((out.flag.str.contains('inferred')).sum())} inferred zeros; "
          f"{int((out.flag.str.contains('out of order')).sum())} from mislabelled axes")


if __name__ == "__main__":
    main()
