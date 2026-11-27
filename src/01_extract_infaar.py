"""Extract resistance percentages from the vector bar charts of the INFAAR 2019-2022 report (plan section 4.1).

For each page: every y-axis is found as a vertical column of numeric tick labels; the value scale is fitted by least
squares to all its ticks. Each filled bar is assigned to the nearest axis to its left whose range contains the bar's
base, and its value is read from the scale. The bar's fill colour is kept; figure, source and antimicrobial are
assigned in 02_label_infaar.py (the category and series columns here are provisional). Pages with raster images are
listed in infaar_pages.csv.

Usage:  python -I src/01_extract_infaar.py            writes outputs/infaar_extracted.csv and outputs/infaar_pages.csv
"""
import os
import re

import fitz
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PDF = os.path.join(REPO, "data_raw", "infaar_2019_22.pdf")
OUT = os.path.join(REPO, "outputs")
NUM = re.compile(r"^\d{1,3}$")


def axes_on(words):
    """Y-axes: numeric ticks in a vertical column (right edges within 7 pt), split into separate axes wherever the
    values stop rising upwards or the spacing jumps; an axis needs at least 3 ticks including 0."""
    ticks = [w for w in words if NUM.match(w[4]) and int(w[4]) <= 100]
    cols = []
    for w in sorted(ticks, key=lambda w: w[2]):
        for c in cols:
            if abs(c["x1"] - w[2]) < 7:
                c["ticks"].append(w)
                break
        else:
            cols.append({"x1": w[2], "ticks": [w]})
    runs = []
    for c in cols:
        ts = sorted(c["ticks"], key=lambda w: w[1])              # top to bottom: values should fall
        run = [ts[0]]
        for prev, w in zip(ts, ts[1:]):
            gaps = [b[1] - a[1] for a, b in zip(run, run[1:])] or [w[1] - prev[1]]
            if int(w[4]) < int(prev[4]) and (w[1] - prev[1]) < 2.2 * np.median(gaps):
                run.append(w)
            else:
                runs.append({"x1": c["x1"], "ticks": run})
                run = [w]
        runs.append({"x1": c["x1"], "ticks": run})
    out = []
    # Fallback for axes whose intermediate tick labels are out of order in the source (Figure 10 B and D read
    # 0, 50, 25, 75, 100 from the bottom): calibrate from the 0 and 100 labels only and flag every value.
    for c in cols:
        zeros = [w for w in c["ticks"] if w[4] == "0"]
        hundreds = [w for w in c["ticks"] if w[4] == "100"]
        for z in zeros:
            h = [w for w in hundreds if 60 < z[1] - w[1] < 400]
            if not h:
                continue
            h = min(h, key=lambda w: z[1] - w[1])
            between = sorted([w for w in c["ticks"] if h[1] < w[1] < z[1]], key=lambda w: w[1])
            vals_down = [int(w[4]) for w in between]
            if vals_down and vals_down != sorted(vals_down, reverse=True):
                yz, yh = (z[1] + z[3]) / 2, (h[1] + h[3]) / 2
                slope = (100 - 0) / (yh - yz)
                out.append({"x1": c["x1"], "y_zero": yz, "y_top": yh, "slope": slope, "icpt": -slope * yz,
                            "vmax": 100, "ticks": 2, "flag": "intermediate tick labels out of order in source; "
                                                               "calibrated on 0 and 100 only"})
    for c in runs:
        vals = [int(t[4]) for t in c["ticks"]]
        if len(c["ticks"]) >= 3 and 0 in vals and len(set(vals)) == len(vals):
            ys = np.array([(t[1] + t[3]) / 2 for t in c["ticks"]])
            v = np.array(vals, float)
            slope, icpt = np.polyfit(ys, v, 1)          # value = slope * y + icpt
            if slope < 0:
                out.append({"x1": c["x1"], "y_zero": (0 - icpt) / slope, "y_top": (max(v) - icpt) / slope,
                            "slope": slope, "icpt": icpt, "vmax": max(v), "ticks": len(vals)})
    return out


def main():
    d = fitz.open(PDF)
    rows, pages = [], []
    for pno in range(d.page_count):
        page = d[pno]
        words = page.get_text("words")
        axes = axes_on(words)
        drawings = page.get_drawings()
        rects = [(it[1], dr["fill"]) for dr in drawings if dr.get("fill") is not None for it in dr["items"]
                 if it[0] == "re" and 2.5 < it[1].width < 25 and it[1].height > 0.3]
        raster = len(page.get_images(full=True))
        pages.append({"page": pno + 1, "axes": len(axes), "bar_rects": len(rects), "raster_images": raster})
        if not axes or not rects:
            continue
        # legend: small coloured squares/lines next to text, matched by colour
        legend = {}
        for r, fill in rects:
            if r.width < 12 and r.height < 12 and abs(r.width - r.height) < 6:
                txt = [w for w in words if abs((w[1] + w[3]) / 2 - (r.y0 + r.y1) / 2) < 5 and 0 < w[0] - r.x1 < 12]
                if txt:
                    line = [w for w in words if abs(w[1] - txt[0][1]) < 3 and w[0] >= txt[0][0] and w[0] - txt[0][0] < 160]
                    legend[tuple(round(c, 3) for c in fill)] = " ".join(w[4] for w in sorted(line, key=lambda w: w[0]))
        for r, fill in rects:
            cand = [a for a in axes if a["x1"] < r.x0 and abs(r.y1 - a["y_zero"]) < 3 and a["y_top"] - 5 < r.y0 < a["y_zero"] + 1]
            if not cand:
                continue
            a = min(cand, key=lambda a: r.x0 - a["x1"])
            value = a["slope"] * r.y0 + a["icpt"]
            cx = (r.x0 + r.x1) / 2
            below = [w for w in words if a["y_zero"] < w[1] < a["y_zero"] + 70 and w[0] - 6 < cx < w[2] + 6
                     and not NUM.match(w[4])]
            label = min(below, key=lambda w: w[1]) [4] if below else ""
            above = [w for w in words if w[3] < a["y_top"] and a["y_top"] - w[3] < 40 and a["x1"] - 20 < w[0] < a["x1"] + 260]
            title = " ".join(w[4] for w in sorted(above, key=lambda w: (round(w[1]), w[0])))
            rows.append({"page": pno + 1, "axis_x": round(a["x1"], 1), "panel_title": title[:120],
                         "series": legend.get(tuple(round(c, 3) for c in fill), ""), "category": label,
                         "value_pct": round(value, 1), "bar_x": round(cx, 1), "bar_w": round(r.width, 2), "bar_top": round(r.y0, 2),
                         "axis_y0": round(a["y_zero"], 1), "axis_ticks": a["ticks"],
                         "fill": "#" + "".join(f"{round(255 * c):02x}" for c in fill[:3]),
                         "flag": a.get("flag", "")})
    os.makedirs(OUT, exist_ok=True)
    df = pd.DataFrame(rows).sort_values(["page", "axis_x", "bar_x"])
    df.to_csv(os.path.join(OUT, "infaar_extracted.csv"), index=False)
    pd.DataFrame(pages).to_csv(os.path.join(OUT, "infaar_pages.csv"), index=False)
    print(pd.DataFrame(pages).query("axes > 0 or raster_images > 0").to_string(index=False))
    print(f"\n{len(df)} bars extracted")
    print(df.groupby(["page", "panel_title"]).size().to_string())


if __name__ == "__main__":
    main()
