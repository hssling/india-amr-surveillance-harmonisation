"""Build the primary AMRSN series and check it (plan sections 4.1 and 4.2).

Primary value for each organism, antimicrobial and year = the value in that year's own annual report:
  E. coli   - own-year column of the yearly trend table, all specimens except urine and faeces (the specimen category
              with the most isolates among those with reported susceptibility, every year; DEVIATIONS C1)
  S. aureus - the "S. aureus, all specimens" column of the percentage susceptibility table
Checks:
  1. every value stated in the reports' text (docs/amrsn_text_checks.csv) against the extracted table value;
  2. S. aureus: own-year column of the trend table against the all-specimens table (two tables in one report);
  3. revisions: each year's value in later reports' trend tables against its own report.

Usage:  python -I src/06_amrsn_primary.py      writes outputs/amrsn_primary.csv, outputs/amrsn_checks.csv,
                                                  outputs/amrsn_revisions.csv
"""
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "outputs")


def canon(name):
    n = re.sub(r"\s+", "", name).lower()
    for key, c in [("ciprofloxac", "ciprofloxacin"), ("cefazolin", "cefazolin"), ("mupirocin", "mupirocin high-level"),
                   ("trimethoprim", "co-trimoxazole"), ("piperacillin", "piperacillin-tazobactam")]:
        if n.startswith(key):
            return c
    return n


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    df = pd.read_csv(os.path.join(OUT, "amrsn_long.csv"))
    df["antimicrobial"] = df.antimicrobial_raw.map(canon)
    df["year"] = pd.to_numeric(df.column, errors="coerce")
    df.loc[df.layout == "groups", "year"] = df.report
    df["year"] = df.year.astype(int)
    df["pct_ns"] = (100 - df.pct_s).round(1)

    ec = df[(df.organism == "E. coli") & (df.specimen != "urine") & (df.layout == "trend") & (df.year == df.report)]
    sa = df[(df.organism == "S. aureus") & (df.layout == "groups") & (df.column == "S. aureus")]
    prim = pd.concat([ec, sa])
    prim = prim[["organism", "specimen", "antimicrobial", "year", "n_susceptible", "n_tested", "pct_s", "pct_ns",
                 "report", "table", "page", "flag"]].sort_values(["organism", "antimicrobial", "year"])
    assert not prim.duplicated(["organism", "antimicrobial", "year"]).any()
    prim.to_csv(os.path.join(OUT, "amrsn_primary.csv"), index=False)

    # 1. text statements
    tc = pd.read_csv(os.path.join(REPO, "docs", "amrsn_text_checks.csv"))
    rows = []
    for _, r in tc.iterrows():
        src = df[(df.organism == r.organism) & (df.antimicrobial == r.antimicrobial) & (df.year == r.year)
                 & (df.report == r.report) & (df.specimen != "urine")
                 & ((df.layout == "trend") if r.organism == "E. coli" else (df.column == "S. aureus"))]
        if src.empty:     # a statement about an earlier year: compare with that year's own report
            src = prim[(prim.organism == r.organism) & (prim.antimicrobial == r.antimicrobial) & (prim.year == r.year)]
        x = float(src.iloc[0]["pct_s" if r.measure == "S" else "pct_ns"])
        lo, hi = (r.low, r.high) if pd.notna(r.low) else (r.stated, r.stated)
        # integer statements are rounded: agreement = within rounding (0.5) of the stated value or range
        tol = 0.5 if float(lo).is_integer() and float(hi).is_integer() else 0.05
        diff = 0.0 if lo <= x <= hi else (x - hi if x > hi else x - lo)
        rows.append({**r.to_dict(), "extracted": x, "diff": round(diff, 2), "agree": abs(diff) <= tol + 1e-9})
    chk = pd.DataFrame(rows)
    chk.to_csv(os.path.join(OUT, "amrsn_checks.csv"), index=False)
    print(f"1. Text statements: {int(chk.agree.sum())}/{len(chk)} agree within rounding")
    if (~chk.agree).any():
        print(chk[~chk.agree][["report", "page", "antimicrobial", "year", "measure", "stated", "low", "high", "extracted", "diff"]].to_string(index=False))

    # 2. S. aureus trend own-year vs all-specimens table
    st = df[(df.organism == "S. aureus") & (df.layout == "trend") & (df.year == df.report)].set_index(["report", "antimicrobial"])
    sg = sa.set_index(["report", "antimicrobial"])
    j = st[["n_susceptible", "n_tested"]].join(sg[["n_susceptible", "n_tested"]], rsuffix="_tab", how="inner")
    mism = j[(j.n_susceptible != j.n_susceptible_tab) | (j.n_tested != j.n_tested_tab)]
    print(f"2. S. aureus trend own-year column vs all-specimens table: {len(j) - len(mism)}/{len(j)} identical counts")
    if len(mism):
        print(mism.to_string())

    # 3. revisions in later reports
    tr = df[df.layout == "trend"]
    rev = tr.merge(prim[["organism", "antimicrobial", "year", "pct_s", "n_tested"]], on=["organism", "antimicrobial", "year"],
                   suffixes=("", "_own"))
    rev = rev[rev.report > rev.year]
    rev["diff_pct_s"] = (rev.pct_s - rev.pct_s_own).round(1)
    rev["diff_n"] = rev.n_tested - rev.n_tested_own
    rev[["organism", "antimicrobial", "year", "report", "pct_s_own", "pct_s", "diff_pct_s", "n_tested_own", "n_tested", "diff_n"]] \
        .to_csv(os.path.join(OUT, "amrsn_revisions.csv"), index=False)
    print(f"3. Later-report revisions of 2019-2021 values: {len(rev)} comparisons; |diff %S| mean {rev.diff_pct_s.abs().mean():.2f}, "
          f"max {rev.diff_pct_s.abs().max():.1f} pp; {int((rev.diff_n != 0).sum())} with a changed denominator")
    print(rev.loc[rev.diff_pct_s.abs().nlargest(5).index, ["organism", "antimicrobial", "year", "report", "pct_s_own", "pct_s", "n_tested_own", "n_tested"]].to_string(index=False))

    pd.set_option("display.width", 200)
    print("\nPrimary series (% non-susceptible, n tested):")
    print(prim.pivot_table(index=["organism", "antimicrobial"], columns="year", values="pct_ns").to_string())


if __name__ == "__main__":
    main()
