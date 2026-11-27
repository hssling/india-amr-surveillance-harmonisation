"""Minimum common panel (Aim 3; rule fixed in docs/DEVIATIONS.md D1).

Rule: for each organism shared by both networks, every antimicrobial class tested by either network is tested by
both, with one representative drug per class: the drug already tested by the most surveillance units (network x
sector or specimen), ties broken in favour of a drug on the EU harmonised panel, then alphabetically. Organisms
for which one network has no routine panel are excluded. Methicillin
markers (cefoxitin, oxacillin) are kept as one "methicillin resistance" class for staphylococci.

E. coli result is benchmarked against the EU harmonised monitoring panel for indicator commensal E. coli
(Commission Implementing Decision (EU) 2020/1729, Annex, first panel). EU_ECOLI below was taken from the EUR-Lex
consolidated text (01.01.2025) and must be checked against OJ L 387, 19.11.2020, p. 8 before publication.

Usage:  python -I src/09_min_panel.py      writes outputs/min_panel.csv, outputs/min_panel_eu_benchmark.csv
"""
import importlib.util
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "outputs")
spec = importlib.util.spec_from_file_location("aud", os.path.join(HERE, "08_audit_panels.py"))
aud = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aud)

EU_ECOLI = {  # drug: class (EU 2020/1729 first panel, Salmonella and indicator commensal E. coli) - VERIFY against OJ
    "amikacin": "aminoglycosides", "ampicillin": "penicillins", "azithromycin": "macrolides",
    "cefotaxime": "3rd-gen cephalosporins", "ceftazidime": "3rd-gen cephalosporins", "chloramphenicol": "phenicols",
    "ciprofloxacin": "fluoroquinolones", "colistin": "polymyxins", "gentamicin": "aminoglycosides",
    "meropenem": "carbapenems", "nalidixic acid": "quinolones (non-fluorinated)",
    "sulfamethoxazole": "folate-pathway inhibitors", "tetracycline": "tetracyclines", "tigecycline": "glycylcyclines",
    "trimethoprim": "folate-pathway inhibitors",
}
STAPH = {"S. aureus", "CoNS"}
METH = {"cefoxitin", "oxacillin"}


def klass(org, drug):
    if org in STAPH and drug in METH:
        return "methicillin resistance (cefoxitin/oxacillin)"
    return aud.CLASS[drug]


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    pd.set_option("display.width", 220, "display.max_colwidth", 60)
    s = pd.read_csv(os.path.join(OUT, "amrsn_series.csv"))
    v = pd.read_csv(os.path.join(OUT, "infaar_values.csv"))
    units = []   # (organism, network, unit, drug)
    for org in dict.fromkeys(o for o, _ in aud.INFAAR_FIGS):
        for spec_ in s[s.organism == org].specimen.unique():
            for d in set(s[(s.organism == org) & (s.specimen == spec_) & (s.n_tested >= aud.MIN_N)].antimicrobial):
                units.append((org, "AMRSN", spec_, d))
        for (o, sector), figs in aud.INFAAR_FIGS.items():
            if o == org:
                for d in {x for x in v[v.figure.isin(figs)].antimicrobial if "undefined" not in x}:
                    units.append((org, "INFAAR", sector, d))
    u = pd.DataFrame(units, columns=["organism", "network", "unit", "drug"])
    u["class"] = [klass(o, d) for o, d in zip(u.organism, u.drug)]
    # only organisms for which both networks have a routine panel (AMRSN tests non-cholera Vibrio on <30 isolates/year)
    both = [o for o in u.organism.unique() if {"AMRSN", "INFAAR"} <= set(u[u.organism == o].network)]
    u = u[u.organism.isin(both)]
    eu_cl = set(EU_ECOLI.values())
    rows = []
    for (org, cl), g in u.groupby(["organism", "class"]):
        counts = g.groupby("drug").unit.nunique()
        rep = sorted(counts.index, key=lambda d: (-counts[d], d not in EU_ECOLI, d))[0]
        nets = set(g.network)
        rep_nets = set(g[g.drug == rep].network)
        rows.append({"organism": org, "class": cl, "representative": rep,
                     "units_testing_rep": int(counts[rep]), "class_tested_by": " + ".join(sorted(nets)),
                     "rep_tested_by": " + ".join(sorted(rep_nets)),
                     "to_add": ", ".join(f"{n} adds {rep}" for n in ("AMRSN", "INFAAR") if n not in rep_nets),
                     "other_drugs_in_class": ", ".join(sorted(set(g.drug) - {rep})),
                     # tiering added after derivation (DEVIATIONS D1 addendum): E. coli core = class also on EU panel
                     "tier": ("core" if cl in eu_cl else "extended") if org == "E. coli" else "not benchmarked"})
    mp = pd.DataFrame(rows)
    mp.to_csv(os.path.join(OUT, "min_panel.csv"), index=False)
    for org, g in mp.groupby("organism", sort=False):
        both_cls = int((g.class_tested_by == "AMRSN + INFAAR").sum())
        print(f"\n{org}: {len(g)} classes; class tested by both {both_cls}; same drug in both {int((g.to_add == '').sum())}; "
              f"additions AMRSN {int(g.to_add.str.contains('AMRSN').sum())}, INFAAR {int(g.to_add.str.contains('INFAAR').sum())}")
        print(g[["class", "representative", "rep_tested_by", "to_add", "tier"]].to_string(index=False))

    # EU benchmark for E. coli: class by class
    ec = mp[mp.organism == "E. coli"]
    eu_classes = sorted(set(EU_ECOLI.values()))
    bench = []
    for cl in sorted(set(eu_classes) | set(ec["class"])):
        r = ec[ec["class"] == cl]
        bench.append({"class": cl, "EU_drugs": ", ".join(sorted(d for d, c in EU_ECOLI.items() if c == cl)),
                      "derived_representative": r.representative.iloc[0] if len(r) else "",
                      "AMRSN_tests_class": bool(len(r)) and "AMRSN" in r.class_tested_by.iloc[0],
                      "INFAAR_tests_class": bool(len(r)) and "INFAAR" in r.class_tested_by.iloc[0],
                      "in_EU": cl in eu_classes, "in_derived": bool(len(r))})
    b = pd.DataFrame(bench)
    b.to_csv(os.path.join(OUT, "min_panel_eu_benchmark.csv"), index=False)
    print("\nE. coli derived panel vs EU 2020/1729 first panel (by class):")
    print(b.to_string(index=False))


if __name__ == "__main__":
    main()
