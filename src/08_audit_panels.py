"""Audit dimension 2: antimicrobial panels for each shared organism (plan section 3).

AMRSN panel = antimicrobials tested on at least 30 isolates in at least one year (excludes stray results such as
2 urine-panel isolates in a non-urine table) in the tables used, 2019-2022 (E. coli: non-urine and urine tables
listed separately). INFAAR panel = antimicrobials on the bars of the figure(s) for that organism (labels as printed; "CAX" is undefined in the report and kept apart).
Overlap = count and Jaccard index (shared / union). Class coverage: for each organism, the antimicrobial classes
each network tests (all specimens and sectors combined), and the classes only one network covers.

Usage:  python -I src/08_audit_panels.py      writes outputs/audit_panels.csv, outputs/audit_class_coverage.csv
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "outputs")
MIN_N = 30
FQ = {"ciprofloxacin", "levofloxacin", "norfloxacin", "enrofloxacin", "nalidixic acid"}

CLASS = {
    "ampicillin": "penicillins", "penicillin": "penicillins", "amoxicillin-clavulanate": "beta-lactam/inhibitor",
    "piperacillin-tazobactam": "beta-lactam/inhibitor", "cefazolin": "1st-gen cephalosporins",
    "cefoxitin": "cephamycins (methicillin marker in staphylococci)", "oxacillin": "isoxazolyl penicillins (methicillin marker)",
    "cefotaxime": "3rd-gen cephalosporins", "ceftriaxone": "3rd-gen cephalosporins", "ceftazidime": "3rd-gen cephalosporins",
    "cefpodoxime": "3rd-gen cephalosporins", "cefixime": "3rd-gen cephalosporins", "cefepime": "4th-gen cephalosporins",
    "aztreonam": "monobactams", "imipenem": "carbapenems", "meropenem": "carbapenems", "ertapenem": "carbapenems",
    "amikacin": "aminoglycosides", "gentamicin": "aminoglycosides", "ciprofloxacin": "fluoroquinolones",
    "levofloxacin": "fluoroquinolones", "norfloxacin": "fluoroquinolones", "enrofloxacin": "fluoroquinolones",
    "nalidixic acid": "quinolones (non-fluorinated)", "co-trimoxazole": "folate-pathway inhibitors",
    "tetracycline": "tetracyclines", "minocycline": "tetracyclines", "tigecycline": "glycylcyclines",
    "chloramphenicol": "phenicols", "colistin": "polymyxins", "nitrofurantoin": "nitrofurans", "fosfomycin": "phosphonic acids",
    "erythromycin": "macrolides", "clindamycin": "lincosamides", "linezolid": "oxazolidinones", "vancomycin": "glycopeptides",
    "teicoplanin": "glycopeptides", "mupirocin high-level": "mupirocin",
}

INFAAR_FIGS = {
    ("E. coli", "fisheries"): ["Figure 4"], ("E. coli", "livestock"): ["Figure 11"],
    ("S. aureus", "fisheries"): ["Figure 2"], ("S. aureus", "livestock"): ["Figure 10A"],
    ("CoNS", "fisheries"): ["Figure 3"], ("CoNS", "livestock"): ["Figure 10D"],
    ("Aeromonas spp.", "fisheries"): ["Figure 7"], ("Vibrio spp.", "fisheries"): ["Figure 6"],
}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    s = pd.read_csv(os.path.join(OUT, "amrsn_series.csv")).fillna({"flag": ""})
    v = pd.read_csv(os.path.join(OUT, "infaar_values.csv"))
    rows = []
    for (org, sector), figs in INFAAR_FIGS.items():
        a = set(v[v.figure.isin(figs)].antimicrobial)
        undefined = {x for x in a if "undefined" in x}
        a -= undefined
        specs = s[s.organism == org].specimen.unique()
        for spec in specs:
            h = s[(s.organism == org) & (s.specimen == spec) & (s.n_tested >= MIN_N)]
            hset = set(h.antimicrobial)
            shared = sorted(hset & a)
            rows.append({"organism": org, "sector": sector, "human_specimen": spec, "n_human": len(hset), "n_animal": len(a),
                         "n_shared": len(shared), "jaccard": round(len(shared) / len(hset | a), 2),
                         "shared": ", ".join(shared), "human_only": ", ".join(sorted(hset - a)),
                         "animal_only": ", ".join(sorted(a - hset)) + (f" (+ {', '.join(sorted(undefined))})" if undefined else ""),
                         "fluoroquinolone_class_both": bool(hset & FQ) and bool(a & FQ)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "audit_panels.csv"), index=False)

    cov = []
    for org in dict.fromkeys(o for o, _ in INFAAR_FIGS):
        figs = [f for (o, _), fs in INFAAR_FIGS.items() if o == org for f in fs]
        a = {x for x in v[v.figure.isin(figs)].antimicrobial if "undefined" not in x}
        h = set(s[(s.organism == org) & (s.n_tested >= MIN_N)].antimicrobial)
        unknown = (a | h) - set(CLASS)
        assert not unknown, f"no class for {unknown}"
        ca, ch = {CLASS[x] for x in a}, {CLASS[x] for x in h}
        for c in sorted(ca | ch):
            cov.append({"organism": org, "class": c, "AMRSN": ", ".join(sorted(x for x in h if CLASS[x] == c)),
                        "INFAAR": ", ".join(sorted(x for x in a if CLASS[x] == c)),
                        "coverage": "both" if c in ca and c in ch else ("AMRSN only" if c in ch else "INFAAR only")})
    cv = pd.DataFrame(cov)
    cv.to_csv(os.path.join(OUT, "audit_class_coverage.csv"), index=False)
    print()
    print(cv.groupby(["organism", "coverage"]).size().unstack(fill_value=0).to_string())
    pd.set_option("display.width", 250, "display.max_colwidth", 90)
    print(out[["organism", "sector", "human_specimen", "n_human", "n_animal", "n_shared", "jaccard", "shared"]].to_string(index=False))


if __name__ == "__main__":
    main()
