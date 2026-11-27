"""Exploratory graph and network analyses (docs/DEVIATIONS.md D2-D5).

D2  Measurement DAG: how true resistance in each reservoir becomes the reported percentage; measurement nodes carry
    their Part A rating; common causes of resistance in both reservoirs. Checks acyclicity and lists, for the
    human-minus-animal contrast, every path that does not pass through true resistance (measurement paths).
D3  Panel network and silo index: bipartite units x antimicrobials, projected onto units (Jaccard weights). Weighted
    modularity of the partition by network vs by organism and organism group (Gram-negative, staphylococci);
    label-permutation test; within-organism panel similarity
    within vs across networks; Louvain communities.
D4  Resistance-profile network: Euclidean distance between source profiles; minimum spanning tree; nearest animal
    source to each human series; parametric bootstrap for stability.
D5  Newcombe hybrid-score intervals for human-minus-animal differences; method-effect yardstick from measurement
    effects documented inside the networks.

Usage:  python -I src/10_network_analyses.py      writes outputs/net_*.csv
"""
import importlib.util
import itertools
import os
import sys

import networkx as nx
import numpy as np
import pandas as pd
from networkx.algorithms.community import louvain_communities, modularity

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "outputs")
SEED = 2026
Z = 1.959964


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + ".py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


aud = load("08_audit_panels")
con = load("07_concordance")

# ---------------------------------------------------------------- D2 measurement DAG
RATING = {"Sampling frame": "not aligned", "Isolate selection / de-duplication": "not reported",
          "AST method": "partly aligned",
          # partly aligned on paper: ICMR human SOP 2019 (CLSI M100 29th ed) and ICMR-ICAR-FAO veterinary SOP 2019
          # (CLSI 2018 adapted) give identical surveillance breakpoints for the shared E. coli drugs; INFAAR does
          # not cite the SOP, so adherence is unverified
          "Interpretive standard and edition": "partly aligned",
          "Intermediate category handling": "not aligned", "Antimicrobial panel": "partly aligned",
          "Reporting and aggregation": "not aligned", "Quality assurance": "partly aligned"}
MEAS = list(RATING)


def build_dag():
    g = nx.DiGraph()
    for side in ("human", "animal"):
        true, rep = f"True resistance ({side})", f"Reported % ({side})"
        g.add_edge(true, rep)
        for m in MEAS:
            g.add_node(f"{m} ({side})", kind="measurement", rating=RATING[m], side=side)
            g.add_edge(f"{m} ({side})", rep)
        g.add_edge(f"Antimicrobial use ({side})", true)
        g.add_edge(rep, "Reported difference (human - animal)")
    for common in ("Shared environment (water, soil, waste)", "Food chain and animal contact"):
        g.add_edge(common, "True resistance (human)")
        g.add_edge(common, "True resistance (animal)")
    g.add_edge("Antimicrobial use (animal)", "Shared environment (water, soil, waste)")
    g.add_edge("Antimicrobial use (human)", "Shared environment (water, soil, waste)")
    assert nx.is_directed_acyclic_graph(g)
    return g


def dag_summary(g):
    target = "Reported difference (human - animal)"
    rows = []
    for n, d in g.nodes(data=True):
        if d.get("kind") == "measurement":
            rows.append({"node": n, "side": d["side"], "rating": d["rating"],
                         "paths_to_contrast": len(list(nx.all_simple_paths(g, n, target))),
                         "bypasses_true_resistance": True})
    out = pd.DataFrame(rows)
    # Can the contrast identify a difference in TRUE resistance? Only if every measurement node is aligned (equal on
    # both arms) so that its effect cancels in the difference.
    misaligned = out[out.rating != "aligned"].node.str.replace(r" \((human|animal)\)$", "", regex=True).unique()
    # Direction of transmission: true resistance nodes share common causes -> d-connected without any arrow between
    # them, so association cannot establish direction.
    th, ta = "True resistance (human)", "True resistance (animal)"
    d_sep = nx.is_d_separator(g, {th}, {ta}, set()) if hasattr(nx, "is_d_separator") else nx.d_separated(g, {th}, {ta}, set())
    return out, sorted(misaligned), d_sep


# ---------------------------------------------------------------- D3 panel network and silo index
def units_table():
    s = pd.read_csv(os.path.join(OUT, "amrsn_series.csv"))
    v = pd.read_csv(os.path.join(OUT, "infaar_values.csv"))
    rows = []
    for org in dict.fromkeys(o for o, _ in aud.INFAAR_FIGS):
        for sp in s[s.organism == org].specimen.unique():
            drugs = set(s[(s.organism == org) & (s.specimen == sp) & (s.n_tested >= aud.MIN_N)].antimicrobial)
            if drugs:
                rows.append({"unit": f"AMRSN | {org} | {sp}", "network": "AMRSN", "organism": org, "drugs": drugs})
        for (o, sector), figs in aud.INFAAR_FIGS.items():
            if o == org:
                rows.append({"unit": f"INFAAR | {org} | {sector}", "network": "INFAAR", "organism": org,
                             "drugs": {x for x in v[v.figure.isin(figs)].antimicrobial if "undefined" not in x}})
    return pd.DataFrame(rows)


def silo(units, n_perm=10000):
    G = nx.Graph()
    for _, r in units.iterrows():
        G.add_node(r.unit, network=r.network, organism=r.organism)
    for (_, a), (_, b) in itertools.combinations(units.iterrows(), 2):
        j = len(a.drugs & b.drugs) / len(a.drugs | b.drugs)
        if j > 0:
            G.add_edge(a.unit, b.unit, weight=j)
    part = lambda key: [set(units[units[key] == k].unit) for k in units[key].unique()]
    units = units.assign(org_group=units.organism.map(lambda o: "staphylococci" if o in ("S. aureus", "CoNS") else "Gram-negative"))
    q_net, q_org = modularity(G, part("network"), weight="weight"), modularity(G, part("organism"), weight="weight")
    q_grp = modularity(G, part("org_group"), weight="weight")
    rng = np.random.default_rng(SEED)
    labels = units.network.to_numpy()
    null = []
    for _ in range(n_perm):
        perm = rng.permutation(labels)
        null.append(modularity(G, [set(units.unit[perm == k]) for k in np.unique(perm)], weight="weight"))
    null = np.array(null)
    p = (np.sum(null >= q_net) + 1) / (n_perm + 1)
    # within-organism similarity, within vs across networks
    pairs = []
    for (_, a), (_, b) in itertools.combinations(units.iterrows(), 2):
        if a.organism == b.organism:
            pairs.append({"organism": a.organism, "unit_a": a.unit, "unit_b": b.unit,
                          "same_network": a.network == b.network,
                          "jaccard": round(len(a.drugs & b.drugs) / len(a.drugs | b.drugs), 3)})
    pairs = pd.DataFrame(pairs)
    comms = louvain_communities(G, weight="weight", seed=SEED)
    summary = {"Q_by_network": round(float(q_net), 3), "Q_by_organism": round(float(q_org), 3),
               "Q_by_organism_group": round(float(q_grp), 3),
               "perm_null_mean": round(float(null.mean()), 3), "perm_null_95pct": round(float(np.quantile(null, .95)), 3),
               "perm_p_network": round(float(p), 4), "n_units": len(units), "n_edges": G.number_of_edges(),
               "mean_J_same_org_same_network": round(float(pairs[pairs.same_network].jaccard.mean()), 3),
               "mean_J_same_org_cross_network": round(float(pairs[~pairs.same_network].jaccard.mean()), 3)}
    edges = pd.DataFrame([{"a": a, "b": b, "jaccard": round(d["weight"], 3)} for a, b, d in G.edges(data=True)])
    comm = pd.DataFrame([{"community": i + 1, "unit": u} for i, c in enumerate(comms) for u in sorted(c)])
    return summary, pairs, edges, comm


# ---------------------------------------------------------------- D4 resistance-profile network
def profiles():
    s = pd.read_csv(os.path.join(OUT, "amrsn_series.csv"))
    v = pd.read_csv(os.path.join(OUT, "infaar_values.csv"))
    sets = {
        "E. coli": {"features": {"amikacin": ("amikacin", "amikacin"), "cefotaxime": ("cefotaxime", "cefotaxime"),
                                 "imipenem": ("imipenem", "imipenem"), "fluoroquinolone": ("ciprofloxacin", "enrofloxacin")},
                    "human": [("E. coli", "all specimens except urine and faeces", "Human E. coli (non-urine)"),
                              ("E. coli", "urine", "Human E. coli (urine)")],
                    "animal": ["Figure 4", "Figure 11"]},
        "Staphylococci": {"features": {d: (d, d) for d in ("penicillin", "erythromycin", "linezolid", "tetracycline", "co-trimoxazole")},
                          "human": [("S. aureus", "all specimens", "Human S. aureus"), ("CoNS", "all specimens", "Human CoNS")],
                          "animal": ["Figure 2", "Figure 10A", "Figure 3", "Figure 10D"]},
    }
    FIGORG = {"Figure 2": "S. aureus", "Figure 10A": "S. aureus", "Figure 3": "CoNS", "Figure 10D": "CoNS",
              "Figure 4": "E. coli", "Figure 11": "E. coli"}
    out = {}
    for name, cfg in sets.items():
        nodes = []
        for org, sp, lab in cfg["human"]:
            vals, ns = {}, {}
            for f, (hd, _) in cfg["features"].items():
                h = s[(s.organism == org) & (s.specimen == sp) & (s.antimicrobial == hd)]
                if h.n_tested.sum():
                    vals[f], ns[f] = 100 * h.n_ns.sum() / h.n_tested.sum(), int(h.n_tested.sum())
            if len(vals) == len(cfg["features"]):
                nodes.append({"node": lab, "group": "human", "vals": vals, "ns": ns, "n_flag": ""})
        for fig in cfg["animal"]:
            for src in v[v.figure == fig].source.unique():
                vals, ns = {}, {}
                for f, (_, ad) in cfg["features"].items():
                    x = v[(v.figure == fig) & (v.source == src) & (v.antimicrobial == ad)]
                    if len(x):
                        vals[f] = float(x.value_pct.iloc[0])
                n = con.N_INFAAR.get((fig, src))
                flag = "" if n else "count not reported; 125 assumed"
                ns = {f: (n or 125) for f in vals}
                sector = con.SECTOR[src]
                lab = f"{FIGORG[fig]} {src}" if name == "Staphylococci" else src
                if len(vals) == len(cfg["features"]):
                    nodes.append({"node": f"{lab} ({sector})", "group": sector, "vals": vals, "ns": ns, "n_flag": flag})
        out[name] = (list(cfg["features"]), nodes)
    return out


def profile_network(features, nodes, n_boot=2000):
    X = np.array([[nd["vals"][f] for f in features] for nd in nodes])
    N = np.array([[nd["ns"][f] for f in features] for nd in nodes])
    labs = [nd["node"] for nd in nodes]

    def mst_and_nn(M):
        D = np.sqrt(((M[:, None, :] - M[None, :, :]) ** 2).sum(-1))
        G = nx.Graph()
        for i, j in itertools.combinations(range(len(labs)), 2):
            G.add_edge(labs[i], labs[j], weight=D[i, j])
        T = nx.minimum_spanning_tree(G)
        nn = {}
        for i, nd in enumerate(nodes):
            if nd["group"] == "human":
                cand = [(D[i, k], labs[k]) for k in range(len(labs)) if nodes[k]["group"] != "human"]
                nn[labs[i]] = min(cand)[1]
        return D, T, nn

    D, T, nn = mst_and_nn(X)
    rng = np.random.default_rng(SEED)
    edge_freq, nn_freq = {}, {h: {} for h in nn}
    for _ in range(n_boot):
        Xb = rng.binomial(N, np.clip(X / 100, 0, 1)) / N * 100
        _, Tb, nnb = mst_and_nn(Xb)
        for a, b in Tb.edges():
            k = tuple(sorted((a, b)))
            edge_freq[k] = edge_freq.get(k, 0) + 1
        for h, a in nnb.items():
            nn_freq[h][a] = nn_freq[h].get(a, 0) + 1
    mst = pd.DataFrame([{"a": a, "b": b, "distance_pp": round(d["weight"], 1),
                         "bootstrap_freq": round(edge_freq.get(tuple(sorted((a, b))), 0) / n_boot, 3)}
                        for a, b, d in T.edges(data=True)]).sort_values("distance_pp")
    near = pd.DataFrame([{"human_series": h, "nearest_animal_source": nn[h],
                          "bootstrap_prob_nearest": round(nn_freq[h].get(nn[h], 0) / n_boot, 3),
                          "runner_up": max(((a, c) for a, c in nn_freq[h].items() if a != nn[h]), default=("", 0), key=lambda t: t[1])[0]}
                         for h in nn])
    prof = pd.DataFrame([{"node": nd["node"], "group": nd["group"], "n_flag": nd["n_flag"], **{f: round(nd["vals"][f], 1) for f in features}}
                         for nd in nodes])
    dist = pd.DataFrame(np.round(D, 1), index=labs, columns=labs)
    return mst, near, prof, dist


# ---------------------------------------------------------------- D5 Newcombe intervals and yardstick
def wilson(k, n):
    p = k / n
    den = 1 + Z ** 2 / n
    c = (p + Z ** 2 / (2 * n)) / den
    h = Z * np.sqrt(p * (1 - p) / n + Z ** 2 / (4 * n ** 2)) / den
    return c - h, c + h


def newcombe(k1, n1, k2, n2):
    p1, p2 = k1 / n1, k2 / n2
    l1, u1 = wilson(k1, n1)
    l2, u2 = wilson(k2, n2)
    d = p1 - p2
    return 100 * (d - np.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2)), 100 * (d + np.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2))


def yardstick():
    lg = pd.read_csv(os.path.join(OUT, "amrsn_long.csv"))
    sa = lg[(lg.organism == "S. aureus") & (lg.layout == "groups") & (lg.column == "S. aureus")]
    rows = []
    for yr, g in sa.groupby("report"):
        fox = g[g.antimicrobial_raw == "Cefoxitin"].pct_s.iloc[0]
        oxa = g[g.antimicrobial_raw == "Oxacillin"].pct_s.iloc[0]
        rows.append({"source": f"AMRSN {yr}: MRSA by cefoxitin vs oxacillin (same network-year, overlapping isolates)",
                     "effect_pp": round(abs(fox - oxa), 1)})
    rev = pd.read_csv(os.path.join(OUT, "amrsn_revisions.csv"))
    rows.append({"source": "AMRSN: revision of a year's value in later reports (max)", "effect_pp": float(rev.diff_pct_s.abs().max())})
    val = pd.read_csv(os.path.join(OUT, "infaar_validation.csv"))
    typo = val[(~val.agree_2pp) & ~val.figure.isin(["Figure 10C", "Figure 10D"])]
    rows.append({"source": "INFAAR: text value vs its own bar, reporting errors (max)", "effect_pp": float(typo.diff_pp.abs().max())})
    pix = pd.read_csv(os.path.join(OUT, "infaar_pixel_check.csv"))
    rows.append({"source": "This study: extraction, vector vs pixel reading (max)", "effect_pp": float(pix.abs_diff.max())})
    y = pd.DataFrame(rows)
    return y, float(y.effect_pp.max())


def differences(ystick):
    c = pd.read_csv(os.path.join(OUT, "concordance.csv"))
    out = []
    for _, r in c.iterrows():
        k1, n1 = round(r.human_pct_ns * r.human_n / 100), int(r.human_n)
        lo = hi = np.nan
        if pd.notna(r.animal_n_isolates):
            n2 = int(r.animal_n_isolates)
            lo, hi = newcombe(k1, n1, round(r.animal_pct_r * n2 / 100), n2)
        d = r.difference_pp
        if np.isnan(lo):
            cls = "no interval (count not reported)"
        elif lo > ystick:
            cls = "human higher, beyond yardstick"
        elif hi < -ystick:
            cls = "animal higher, beyond yardstick"
        elif lo > 0:
            cls = "human higher, within yardstick"
        elif hi < 0:
            cls = "animal higher, within yardstick"
        else:
            cls = "no clear difference"
        out.append({**r.to_dict(), "diff_ci_low": round(lo, 1), "diff_ci_high": round(hi, 1), "reading": cls})
    return pd.DataFrame(out)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    pd.set_option("display.width", 230, "display.max_colwidth", 70)

    g = build_dag()
    dag_nodes, misaligned, d_sep = dag_summary(g)
    pd.DataFrame([{"from": a, "to": b} for a, b in g.edges()]).to_csv(os.path.join(OUT, "net_dag_edges.csv"), index=False)
    dag_nodes.to_csv(os.path.join(OUT, "net_dag_nodes.csv"), index=False)
    print(f"D2 DAG: {g.number_of_nodes()} nodes, {g.number_of_edges()} edges, acyclic. Measurement dimensions not aligned or "
          f"not reported on both arms (open paths to the contrast): {len(misaligned)}/{len(MEAS)}: {', '.join(misaligned)}")
    print(f"   True resistance (human) and (animal) d-separated given nothing: {d_sep} (common causes -> association "
          f"without any arrow between them; direction not identifiable from these data)")

    units = units_table()
    summ, pairs, edges, comm = silo(units)
    pairs.to_csv(os.path.join(OUT, "net_panel_pairs.csv"), index=False)
    edges.to_csv(os.path.join(OUT, "net_panel_edges.csv"), index=False)
    comm.to_csv(os.path.join(OUT, "net_panel_communities.csv"), index=False)
    pd.DataFrame([summ]).to_csv(os.path.join(OUT, "net_silo_summary.csv"), index=False)
    print("\nD3 panel network / silo index:", summ)
    print(pairs.sort_values(["organism", "same_network"]).to_string(index=False))
    print(comm.groupby("community").unit.apply(lambda u: "; ".join(u)).to_string())

    for name, (feats, nodes) in profiles().items():
        mst, near, prof, dist = profile_network(feats, nodes)
        key = name.lower().replace(". ", "").replace(" ", "_")
        mst.to_csv(os.path.join(OUT, f"net_profile_mst_{key}.csv"), index=False)
        near.to_csv(os.path.join(OUT, f"net_profile_nearest_{key}.csv"), index=False)
        prof.to_csv(os.path.join(OUT, f"net_profile_values_{key}.csv"), index=False)
        dist.to_csv(os.path.join(OUT, f"net_profile_distance_{key}.csv"))
        print(f"\nD4 {name}: features {feats}")
        print(prof.to_string(index=False))
        print(mst.to_string(index=False))
        print(near.to_string(index=False))

    y, ys = yardstick()
    y.to_csv(os.path.join(OUT, "net_yardstick.csv"), index=False)
    print("\nD5 method-effect yardstick sources:")
    print(y.to_string(index=False))
    print(f"   yardstick = {ys} pp")
    d = differences(ys)
    d.to_csv(os.path.join(OUT, "net_differences.csv"), index=False)
    print(d.groupby(["organism", "reading"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main()
