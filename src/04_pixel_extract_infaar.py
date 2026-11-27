"""Second extraction of every INFAAR bar, reading bar heights from the rendered page image (plan section 4.2).

Each chart page is rendered at 288 dpi. For each bar, two pixel columns (at +/-30% of the bar width, clear of the
error bar drawn on the centre line) are scanned upwards from the axis zero while the pixel colour stays within a
tolerance of the bar's fill; the top of that run is the bar top. Bar heights therefore come from the raster, not from
the vector rectangles. Shared with the first method: the tick-label calibration of the value scale and the bar
positions and colours (stated in the report of agreement).

Usage:  python -I src/04_pixel_extract_infaar.py      writes outputs/infaar_pixel_check.csv
"""
import importlib.util
import os
import sys

import fitz
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
PDF = os.path.join(REPO, "data_raw", "infaar_2019_22.pdf")
OUT = os.path.join(REPO, "outputs")
ZOOM = 4.0                       # 288 dpi
COLOUR_TOL = 60                  # Euclidean RGB distance accepted as "bar colour" (anti-aliased edges excluded)
GAP = int(2 * ZOOM)              # bridge up to 2 pt of non-bar pixels (gridlines drawn over bars, axis line at base)

spec = importlib.util.spec_from_file_location("ext", os.path.join(HERE, "01_extract_infaar.py"))
ext = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ext)


def hex_rgb(h):
    return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)], float)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    bars = pd.read_csv(os.path.join(OUT, "infaar_extracted.csv"))
    d = fitz.open(PDF)
    rows = []
    for page, b in bars.groupby("page"):
        pg = d[page - 1]
        pix = pg.get_pixmap(matrix=fitz.Matrix(ZOOM, ZOOM), alpha=False)
        img = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, 3).astype(float)
        axes = ext.axes_on(pg.get_text("words"))
        for _, r in b.iterrows():
            a = min(axes, key=lambda a: abs(a["x1"] - r.axis_x) + abs(a["y_zero"] - r.axis_y0))
            y = int(round(r.axis_y0 * ZOOM)) - 2
            target = hex_rgb(r.fill)
            # An error bar is drawn down each bar's centre line, so two columns at +/-30% of the bar width are
            # scanned. Upwards from just above the axis line, runs of non-bar pixels up to GAP are bridged (gridlines);
            # the bar top is the highest bar-coloured pixel before a longer gap.
            tops = []
            for off in (-0.3, 0.3):
                col = int(round((r.bar_x + off * r.bar_w) * ZOOM))
                last, miss, yy = None, 0, y
                while yy > 0 and miss <= GAP:
                    if np.linalg.norm(img[yy, col] - target) < COLOUR_TOL:
                        last, miss = yy, 0
                    else:
                        miss += 1
                    yy -= 1
                if last is not None:
                    tops.append(last / ZOOM)
            if not tops:
                rows.append({**r.to_dict(), "pixel_value": np.nan, "pixel_note": "bar colour not found above base"})
                continue
            top = float(np.mean(tops))
            rows.append({**r.to_dict(), "pixel_value": round(a["slope"] * top + a["icpt"], 1), "pixel_note": ""})
    res = pd.DataFrame(rows)
    res["abs_diff"] = (res.pixel_value - res.value_pct).abs()
    res.to_csv(os.path.join(OUT, "infaar_pixel_check.csv"), index=False)
    ok = res.dropna(subset=["pixel_value"])
    print(f"bars: {len(res)}; read by pixel method: {len(ok)}; not read: {len(res) - len(ok)}")
    print(f"agreement with vector extraction: mean |diff| {ok.abs_diff.mean():.2f} pp, median {ok.abs_diff.median():.2f}, "
          f"95th pct {ok.abs_diff.quantile(.95):.2f}, max {ok.abs_diff.max():.2f}; within 1 pp {100 * (ok.abs_diff <= 1).mean():.1f}%, "
          f"within 2 pp {100 * (ok.abs_diff <= 2).mean():.1f}%")
    print(res.sort_values("abs_diff", ascending=False).head(8)[["page", "axis_x", "category", "fill", "value_pct", "pixel_value", "pixel_note"]].to_string(index=False))
    if len(res) > len(ok):
        print(res[res.pixel_value.isna()][["page", "axis_x", "bar_x", "fill", "value_pct"]].to_string(index=False))


if __name__ == "__main__":
    main()
