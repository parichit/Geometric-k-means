"""E1/E2 supporting figure: the five geokmeans algorithms on sensor.

Panels:
  (a) table  -- SSE at convergence per algorithm and k
  (b) line   -- mean distance computations vs k, log scale
  (c) line   -- mean runtime vs k

This is NOT Figure 1 of the paper. Figure 1 is the twitter three-way
comparison built by scripts/plot_results.py; this script keeps the older
sensor panels reproducible on their own. See docs/E1_E2_SENSOR.md.

Reads the tidy summary written by report.py and plots means with std error
bars. Nothing is recomputed here.

Usage:
    python scripts/plot_e1e2.py [--summary results/reports/..._exactness.csv]
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ALGO_ORDER = ["lloyd", "hamerly", "annulus", "exponion", "geo"]
LABEL = {"lloyd": "Lloyd", "hamerly": "Hamerly", "annulus": "Annulus",
         "exponion": "Exponion", "geo": "Geometric"}

# Colour-blind-safe qualitative palette; geo last so it reads as the subject.
COLOR = {"lloyd": "#4C72B0", "hamerly": "#DD8452", "annulus": "#55A868",
         "exponion": "#C44E52", "geo": "#8172B3"}
MARKER = {"lloyd": "o", "hamerly": "s", "annulus": "^",
          "exponion": "D", "geo": "*"}


def _algos(df):
    return [a for a in ALGO_ORDER if a in set(df["algorithm"])]


def _line_panel(ax, df, mean_col, std_col, ylabel, title, scale=1.0, log=False):
    for algo in _algos(df):
        sub = df[df.algorithm == algo].sort_values("k")
        ax.errorbar(sub["k"], sub[mean_col] / scale,
                    yerr=sub[std_col] / scale,
                    label=LABEL[algo], color=COLOR[algo], marker=MARKER[algo],
                    markersize=7 if algo != "geo" else 11,
                    linewidth=1.8, capsize=3, capthick=1.2)
    ax.set_xlabel("number of clusters $k$")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, fontweight="bold", loc="left")
    if log:
        ax.set_yscale("log")
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def _table_panel(ax, df, ks):
    ax.axis("off")
    algos = _algos(df)
    cols = ["Algorithm"] + [f"$k$={int(k)}" for k in ks]
    rows = []
    for algo in algos:
        cells = [LABEL[algo]]
        for k in ks:
            r = df[(df.algorithm == algo) & (df.k == k)]
            if r.empty:
                cells.append("---")
            else:
                m = r["sse_mean"].iloc[0]
                sd = r["sse_std"].iloc[0]
                cells.append(f"{m/1e6:.4f}\n± {sd/1e6:.4f}")
        rows.append(cells)

    tbl = ax.table(cellText=rows, colLabels=cols, cellLoc="center",
                   loc="center", bbox=[0, 0, 1, 0.88])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.5)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.6)
        cell.set_edgecolor("#CCCCCC")
        if r == 0:
            cell.set_facecolor("#F0F0F0")
            cell.set_text_props(fontweight="bold")
        elif c == 0:
            cell.set_text_props(color=COLOR[algos[r - 1]], fontweight="bold")
        cell.set_height(0.16)
    ax.set_title("(a)  SSE at convergence  ($\\times 10^{6}$, mean ± std)",
                 fontsize=11, fontweight="bold", loc="left")


def main():
    base = Path(__file__).parent.parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", type=Path,
                    default=base / "results/reports/sensor_highk_exactness.csv")
    ap.add_argument("--out", type=Path,
                    default=base / "results/figures/e1e2_sensor.png")
    ap.add_argument("--dataset", default=None)
    args = ap.parse_args()

    if not args.summary.exists():
        raise SystemExit(f"summary not found: {args.summary}\n"
                         f"Run scripts/report.py first.")
    df = pd.read_csv(args.summary)
    if args.dataset:
        df = df[df.dataset == args.dataset]
    if df.empty:
        raise SystemExit("no rows to plot")

    ks = sorted(df["k"].unique())
    name = ", ".join(sorted(df["dataset"].unique()))
    n_trials = int(df["n_trials"].max())

    fig = plt.figure(figsize=(14, 4.6))
    gs = fig.add_gridspec(1, 6, wspace=1.5)
    ax_t = fig.add_subplot(gs[0, 0:2])
    ax_d = fig.add_subplot(gs[0, 2:4])
    ax_r = fig.add_subplot(gs[0, 4:6])

    _table_panel(ax_t, df, ks)
    _line_panel(ax_d, df, "n_distance_calculations_mean",
                "n_distance_calculations_std",
                "distance computations ($\\times 10^{6}$)",
                "(b)  Distance computations vs $k$", scale=1e6, log=True)
    _line_panel(ax_r, df, "wall_clock_seconds_mean", "wall_clock_seconds_std",
                "runtime (s)", "(c)  Runtime vs $k$")

    ax_d.set_xticks(ks)
    ax_r.set_xticks(ks)
    handles, labels = ax_r.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels),
               frameon=False, bbox_to_anchor=(0.5, -0.04), fontsize=10)
    fig.suptitle(f"geokmeans algorithms on {name}  \u2014  mean \u00b1 std over "
                 f"{n_trials} trials, shared k-means++ initialisation",
                 fontsize=12, fontweight="bold", y=1.02)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(args.out.with_suffix(".pdf"), bbox_inches="tight",
                facecolor="white")
    print(f"\u2713 Generated: {args.out}")
    print(f"\u2713 Generated: {args.out.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
