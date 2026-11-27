"""Extract AMRSN susceptibility tables (human clinical isolates) from the annual report PDFs (plan section 4.1).

Only tables with a text layer are read here, by cell position (PyMuPDF find_tables); tables printed as images are
transcribed separately and flagged. Each cell holds "susceptible/tested (percent)"; line breaks inside numbers are
removed, and every parsed fraction is checked against the percentage printed beside it.

AMRSN reports percentage SUSCEPTIBLE. Non-susceptible (resistant plus intermediate) = 100 - %S.

Usage:  python -I src/05_extract_amrsn.py      writes outputs/amrsn_long.csv (every cell of every table read)
"""
import os
import re
import sys

import fitz
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
RAW = os.path.join(REPO, "data_raw")
OUT = os.path.join(REPO, "outputs")

EXCL = "all specimens except urine and faeces"
# (report year, page, table number, organism, specimen, layout, column labels for non-trend layouts)
TABLES = [
    (2019, 33, "2.4.1", "E. coli", EXCL, "trend", None),
    (2020, 75, "3.6", "E. coli", EXCL, "trend", None),
    (2021, 61, "2.6", "E. coli", EXCL, "trend", None),
    (2022, 81, "3.6", "E. coli", EXCL, "trend", None),
    (2019, 90, "6.1", "S. aureus", "all specimens", "groups", ["S. aureus", "MSSA", "MRSA", "CoNS"]),
    (2020, 140, "7.1", "S. aureus", "all specimens", "groups", ["S. aureus", "MSSA", "MRSA", "CoNS"]),
    (2021, 100, "4.1", "S. aureus", "all specimens", "groups", ["S. aureus", "MSSA", "MRSA", "CoNS"]),
    (2022, 145, "6.1", "S. aureus", "all specimens", "groups", ["S. aureus", "MSSA", "MRSA", "CoNS"]),
    (2019, 93, "6.4", "S. aureus", "all specimens", "trend", None),
    (2020, 146, "7.4", "S. aureus", "all specimens", "trend", None),
    (2021, 106, "4.4", "S. aureus", "all specimens", "trend", None),
    (2022, 151, "6.4", "S. aureus", "all specimens", "trend", None),
    (2022, 62, "2.13", "E. coli", "urine", "groups", ["OPD", "Ward", "ICU"]),
]
CELL = re.compile(r"^(\*?)(\d+)/(\d+)(?:\(([\d.]+)%?\)|\(-\))?$")


def parse_table(page, layout, groups):
    tabs = page.find_tables().tables
    assert tabs, "no table found"
    tb = max(tabs, key=lambda t: len(t.extract()))
    grid = [[re.sub(r"\s+", " ", c or "").strip() for c in row] for row in tb.extract()]
    is_val = lambda c: bool(re.match(r"^\*?\d+ ?/", c.replace(" ", "")))
    first_data = next(i for i, row in enumerate(grid) if any(is_val(c) for c in row))
    # Column headers may be split over several rows ("Year-" above "2016"): join each column's header cells.
    ncol = max(len(r) for r in grid)
    heads = [" ".join(grid[i][j] for i in range(first_data) if j < len(grid[i]) and grid[i][j]) for j in range(ncol)]
    years = [int(m.group(1)) for h in heads for m in [re.search(r"Year-?\s*(\d{4})", h)] if m]
    labels = years if layout == "trend" else groups
    assert labels, "no column labels"
    words = page.get_text("words")
    out, last = [], None
    for i, row in enumerate(grid):
        vals = [c.replace(" ", "") for c in row if is_val(c)]
        names = [c for c in row if c and not re.match(r"^\*?\d", c.replace(" ", "")) and not c.startswith(("(", "n=", "Year", "Total"))]
        if vals and not names:
            # Name cell not captured by the table finder: take the words left of the first value cell in this row.
            rb = tb.rows[i].bbox
            x_first = min(c[0] for c in tb.rows[i].cells if c is not None and is_val(re.sub(r"\s+", "", page.get_textbox(fitz.Rect(c)))))
            names = [" ".join(w[4] for w in sorted(words, key=lambda w: (round(w[1]), w[0]))
                              if rb[1] - 1 < (w[1] + w[3]) / 2 < rb[3] + 1 and w[2] < x_first - 1)]
        if not vals:
            if names and last is not None and last["name"].endswith("-") and len(names) == 1:
                last["name"] += names[0]                        # "Piperacillin-" + "tazobactam" on the next row
            continue
        name = " ".join(names).strip()
        last = {"name": name, "vals": vals}
        out.append(last)
    rows = []
    for r in out:
        flag = "" if len(r["vals"]) == len(labels) else f"{len(r['vals'])} values for {len(labels)} columns"
        for lab, v in zip(labels, r["vals"]):
            m = CELL.match(v)
            if not m:
                rows.append({"antimicrobial_raw": r["name"], "column": lab, "cell": v, "flag": "unparsed cell"})
                continue
            star, s, n, pct = m.groups()
            s, n = int(s), int(n)
            calc = round(100 * s / n, 1) if n else None
            f = flag
            if pct is not None and calc is not None and abs(float(pct) - calc) > 0.15:
                f = (f + "; " if f else "") + f"printed {pct} != computed {calc}"
            if star:
                f = (f + "; " if f else "") + "asterisk in source (small denominator)"
            rows.append({"antimicrobial_raw": r["name"], "column": lab, "cell": v, "n_susceptible": s, "n_tested": n,
                         "pct_s_printed": float(pct) if pct else None, "pct_s": calc, "flag": f})
    return rows


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    allrows = []
    for yr, pg, tno, org, spec, layout, groups in TABLES:
        d = fitz.open(os.path.join(RAW, f"amrsn_{yr}.pdf"))
        page = d[pg - 1]
        cap = re.search(r"Table " + re.escape(tno) + r"\b", page.get_text())
        assert cap, f"Table {tno} caption not on {yr} p{pg}"
        for r in parse_table(page, layout, groups):
            allrows.append({"report": yr, "page": pg, "table": tno, "organism": org, "specimen": spec, "layout": layout, **r})
    df = pd.DataFrame(allrows)
    df["antimicrobial_raw"] = df.antimicrobial_raw.str.replace(r"\s+", " ", regex=True).str.replace("- ", "-")
    os.makedirs(OUT, exist_ok=True)
    df.to_csv(os.path.join(OUT, "amrsn_long.csv"), index=False)
    pd.set_option("display.width", 220)
    print(df.groupby(["report", "table", "organism"]).agg(cells=("cell", "size"), drugs=("antimicrobial_raw", "nunique")).to_string())
    bad = df[df.flag.fillna("").str.contains("!=|unparsed|values for")]
    print(f"\n{len(df)} cells; {len(bad)} with a parsing or arithmetic flag")
    if len(bad):
        print(bad[["report", "table", "antimicrobial_raw", "column", "cell", "flag"]].to_string(index=False))
    print("\nDrugs per table:")
    for (yr, t), g in df.groupby(["report", "table"]):
        print(yr, t, sorted(g.antimicrobial_raw.unique()))


if __name__ == "__main__":
    main()
