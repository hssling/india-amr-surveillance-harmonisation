"""Validate the INFAAR extraction against the report's own text (plan section 4.2).

1. Pre-specified criterion: every text-stated percentage (docs/infaar_validation_set.csv) is compared with the extracted
   bar it describes, using the panel assignment printed in the report's figure captions. Accepted if >= 90% agree within
   2 percentage points. A range ("23 to 24 percent") agrees if the value is within 2 pp of the range.
2. Investigation of disagreements (added after the plan; reported as such):
   a. Figure 10 C and D: the same comparison with the two panels exchanged.
   b. Pooled statements (docs/infaar_pooled_statements.csv): the isolate-weighted mean of the extracted per-source
      values, using the denominators in the figure legends, is compared with each pooled percentage in the text; and
      where a per-source text value disagrees with its bar, the per-source value implied by the pooled statement is
      back-calculated.

Usage:  python -I src/03_validate_infaar.py      writes outputs/infaar_validation.csv, outputs/infaar_pooled_check.csv
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "outputs")
TOL = 2.0

# Isolates per source, from the figure legends (page 16-19; "G" in the PDF text layer is a broken glyph for 6).
N = {"Figure 2": {"Freshwater fish": 35, "Shrimps": 193, "Marine fish": 48},
     "Figure 3": {"Freshwater fish": 1389, "Shrimps": 444, "Marine fish": 61},
     "Figure 4": {"Freshwater fish": 1441, "Shrimps": 481, "Marine fish": 157},
     "Figure 5": {"Shrimps": 226, "Marine fish": 89},
     "Figure 6": {"Shrimps": 464, "Marine fish": 102}}


def compare(vs, vals, swap_10cd=False):
    key = vals.set_index(["figure", "source", "antimicrobial"])["value_pct"]
    rows = []
    for _, r in vs.iterrows():
        fig = r.figure
        if swap_10cd and fig in ("Figure 10C", "Figure 10D"):
            fig = "Figure 10D" if fig == "Figure 10C" else "Figure 10C"
        k = (fig, r.source, r.antimicrobial)
        if k not in key.index:
            raise KeyError(f"no extracted value for {k}")
        x = float(key[k])
        lo, hi = (r.low, r.high) if r.kind == "range" else (r.stated, r.stated)
        diff = 0.0 if lo <= x <= hi else (x - hi if x > hi else x - lo)
        rows.append({**r.to_dict(), "panel_compared": fig, "extracted": x, "diff_pp": round(diff, 1),
                     "agree_2pp": abs(diff) <= TOL})
    return pd.DataFrame(rows)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    vals = pd.read_csv(os.path.join(OUT, "infaar_values.csv"))
    vs = pd.read_csv(os.path.join(REPO, "docs", "infaar_validation_set.csv"))

    res = compare(vs, vals)
    n, k = len(res), int(res.agree_2pp.sum())
    print(f"PRE-SPECIFIED CRITERION (caption panel assignment): {k}/{n} = {100 * k / n:.1f}% within {TOL} pp "
          f"-> {'ACCEPTED' if k / n >= 0.90 else 'NOT ACCEPTED'}")
    print(f"  mean |diff| {res.diff_pp.abs().mean():.2f} pp, median {res.diff_pp.abs().median():.2f}, max {res.diff_pp.abs().max():.1f}")
    print(res.groupby("kind").agg(n=("agree_2pp", "size"), agree=("agree_2pp", "sum"),
                                  mean_abs_diff=("diff_pp", lambda d: d.abs().mean())).round(2).to_string())
    print("\nDisagreements:")
    print(res[~res.agree_2pp][["figure", "source", "antimicrobial", "kind", "stated", "low", "high", "extracted", "diff_pp"]].to_string(index=False))

    # a. Figure 10 C/D exchanged
    sw = compare(vs, vals, swap_10cd=True)
    cd = res.figure.isin(["Figure 10C", "Figure 10D"])
    print(f"\nFigure 10 C/D, as captioned: {int(res[cd].agree_2pp.sum())}/{int(cd.sum())} agree; "
          f"panels exchanged: {int(sw[cd].agree_2pp.sum())}/{int(cd.sum())} agree")

    # b. pooled statements
    ps = pd.read_csv(os.path.join(REPO, "docs", "infaar_pooled_statements.csv"))
    prow = []
    for _, r in ps.iterrows():
        v = vals[(vals.figure == r.figure) & (vals.antimicrobial == r.antimicrobial)].set_index("source")["value_pct"]
        w = pd.Series(N[r.figure])
        pooled = float((v[w.index] * w).sum() / w.sum())
        prow.append({**r.to_dict(), "weighted_extracted": round(pooled, 1), "diff_pp": round(pooled - r.stated, 1),
                     "agree_2pp": abs(pooled - r.stated) <= TOL})
    pc = pd.DataFrame(prow)
    print(f"\nPooled statements: {int(pc.agree_2pp.sum())}/{len(pc)} agree within {TOL} pp with the isolate-weighted extracted values")
    print(pc[["figure", "antimicrobial", "stated", "weighted_extracted", "diff_pp"]].to_string(index=False))

    # back-calculate disputed per-source values from pooled statements
    print("\nPer-source values implied by the pooled text (where the per-source text disagrees with its bar):")
    for _, r in res[~res.agree_2pp].iterrows():
        p = pc[(pc.figure == r.figure) & (pc.antimicrobial == r.antimicrobial)]
        if p.empty or r.figure not in N:
            continue
        w = pd.Series(N[r.figure])
        others = vals[(vals.figure == r.figure) & (vals.antimicrobial == r.antimicrobial) & (vals.source != r.source)].set_index("source")["value_pct"]
        implied = (p.stated.iloc[0] * w.sum() - (others * w[others.index]).sum()) / w[r.source]
        print(f"  {r.figure} {r.source} {r.antimicrobial}: text {r.stated}, bar {r.extracted}, implied by pooled text {implied:.1f}")

    res.to_csv(os.path.join(OUT, "infaar_validation.csv"), index=False)
    pc.to_csv(os.path.join(OUT, "infaar_pooled_check.csv"), index=False)


if __name__ == "__main__":
    main()
