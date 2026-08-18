"""Single summary figure from the benchmark CSVs.

Panels:
  (a) table  -- SSE at convergence per algorithm and k
  (b) line   -- mean distance computations vs k
  (c) line   -- mean runtime vs k

Reads the tidy summaries written by report.py and plots means with std error
bars. Nothing is recomputed here.

Usage:
    python scripts/plot_results.py [--summary results/reports/..._exactness.csv]
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

# Row 2 (twitter) compares two implementations rather than five algorithms.
IMPL_LABEL = {"geokmeans_geo": "Geometric-k-means",
              "sklearn_lloyd": "scikit-learn (Lloyd)"}
IMPL_COLOR = {"geokmeans_geo": "#8172B3", "sklearn_lloyd": "#937860"}
IMPL_MARKER = {"geokmeans_geo": "*", "sklearn_lloyd": "o"}
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


def _impl_panel(ax, df, mean_col, std_col, ylabel, title, scale=1.0,
                log=False, note=None):
    """Row-2 panel: one line per implementation, mean with std error bars."""
    for impl in ["sklearn_lloyd", "geokmeans_geo"]:
        sub = df[df.implementation == impl].sort_values("k")
        if sub.empty:
            continue
        ax.errorbar(sub["k"], sub[mean_col] / scale, yerr=sub[std_col] / scale,
                    label=IMPL_LABEL[impl], color=IMPL_COLOR[impl],
                    marker=IMPL_MARKER[impl],
                    markersize=11 if impl == "geokmeans_geo" else 7,
                    linewidth=1.8, capsize=3, capthick=1.2)
    ax.set_xlabel("number of clusters $k$")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, fontweight="bold", loc="left")
    if log:
        ax.set_yscale("log")
    if note:
        ax.text(0.5, -0.30, note, transform=ax.transAxes, fontsize=7.5,
                style="italic", color="#555555", ha="center", va="top")
    ax.grid(alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def main():
    base = Path(__file__).parent.parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", type=Path,
                    default=base / "results/reports/sensor_highk_exactness.csv")
    ap.add_argument("--out", type=Path,
                    default=base / "results/figures/benchmark_summary.png")
    ap.add_argument("--dataset", default=None)
    ap.add_argument("--twitter", type=Path,
                    default=base / "results/reports/e3_twitter_summary.csv",
                    help="E3 summary CSV; row 2 is omitted if absent")
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

    tw = pd.read_csv(args.twitter) if args.twitter and args.twitter.exists() else None
    if tw is not None and tw.empty:
        tw = None

    # Six columns so that two panels on row 2 centre cleanly under three on
    # row 1: row 1 takes 0-2, 2-4, 4-6; row 2 takes 1-3 and 3-5.
    if tw is None:
        fig = plt.figure(figsize=(14, 4.6))
        gs = fig.add_gridspec(1, 6, wspace=1.5)
    else:
        fig = plt.figure(figsize=(14, 10.2))
        # generous hspace: row 1 needs room for its x-labels, its own legend,
        # and the row-2 section heading between the two rows
        gs = fig.add_gridspec(2, 6, wspace=1.5, hspace=0.62,
                              top=0.93, bottom=0.07)
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

    if tw is None:
        fig.legend(handles, labels, loc="lower center", ncol=len(labels),
                   frameon=False, bbox_to_anchor=(0.5, -0.04), fontsize=10)
        fig.suptitle(f"geokmeans on {name}  —  mean ± std over {n_trials} "
                     f"trials, shared k-means++ initialisation",
                     fontsize=12, fontweight="bold", y=1.02)
    else:
        fig.legend(handles, labels, loc="upper center", ncol=len(labels),
                   frameon=False, bbox_to_anchor=(0.5, 0.505), fontsize=10)

        ax_w = fig.add_subplot(gs[1, 1:3])
        ax_i = fig.add_subplot(gs[1, 3:5])
        tw_ks = sorted(tw["k"].unique())
        tw_trials = int(tw["n_trials"].max())
        tw_name = ", ".join(sorted(tw["dataset"].unique()))

        _impl_panel(ax_w, tw, "wall_clock_seconds_mean",
                    "wall_clock_seconds_std", "runtime (s)",
                    "(d)  Runtime vs $k$")
        dmax = tw["n_distance_calculations_mean"].max()
        dscale, dexp = (1e9, 9) if dmax >= 1e9 else (1e6, 6)
        _impl_panel(ax_i, tw, "n_distance_calculations_mean",
                    "n_distance_calculations_std",
                    f"distance computations ($\\times 10^{{{dexp}}}$)",
                    "(e)  Distance computations vs $k$",
                    scale=dscale, log=True,
                    note="scikit-learn count is analytic ($n\\,k$ per pass); "
                         "geokmeans is measured")
        for ax in (ax_w, ax_i):
            ax.set_xticks(tw_ks)

        h2, l2 = ax_w.get_legend_handles_labels()
        fig.legend(h2, l2, loc="lower center", ncol=len(l2), frameon=False,
                   bbox_to_anchor=(0.5, -0.02), fontsize=10)

        fig.text(0.5, 0.975,
                 f"(a-c)  geokmeans algorithms on {name}  —  mean ± std over "
                 f"{n_trials} trials, shared k-means++ initialisation",
                 ha="center", fontsize=12.5, fontweight="bold")
        fig.text(0.5, 0.452,
                 f"(d-e)  Geometric-k-means vs scikit-learn Lloyd on {tw_name}"
                 f"  —  mean ± std over {tw_trials} trials, "
                 f"harmonised tolerance",
                 ha="center", fontsize=12.5, fontweight="bold")
        fig.suptitle("")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(args.out.with_suffix(".pdf"), bbox_inches="tight",
                facecolor="white")
    print(f"✓ Generated: {args.out}")
    print(f"✓ Generated: {args.out.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
