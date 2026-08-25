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

The line panels are drawn with seaborn (`sns.lineplot`) on a seaborn theme;
the error-bar caps use matplotlib because seaborn has no API for pre-computed
mean/std columns, and panel (a) is a matplotlib table, which seaborn has no
equivalent of.

Usage:
    python scripts/plot_e1e2.py [--summary results/reports/..._exactness.csv]

To restyle the figure -- titles, labels, font sizes, figure size, colours,
spacing -- edit the CONFIGURATION block below and nothing else.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

# ============================================================================
# CONFIGURATION -- everything a user would want to tweak lives here.
# ============================================================================

# --- which algorithms, in which order ------------------------------------
ALGO_ORDER = ["lloyd", "hamerly", "annulus", "exponion", "geo"]

# Legend / table row names. Edit the right-hand side to rename a series.
LABEL = {"lloyd": "Lloyd", "hamerly": "Hamerly", "annulus": "Annulus",
         "exponion": "Exponion", "geo": "Geometric"}

# Colour-blind-safe qualitative palette; geo last so it reads as the subject.
# Colour follows the algorithm, never its rank, so dropping one from the sweep
# does not repaint the others.
COLOR = {"lloyd": "#4C72B0", "hamerly": "#DD8452", "annulus": "#55A868",
         "exponion": "#C44E52", "geo": "#8172B3"}
MARKER = {"lloyd": "o", "hamerly": "s", "annulus": "^",
          "exponion": "D", "geo": "*"}
# Per-algorithm marker size; geo is larger so the subject stands out.
MARKER_SIZE = {"lloyd": 7, "hamerly": 7, "annulus": 7, "exponion": 7,
               "geo": 11}

# --- seaborn theme -------------------------------------------------------
# context: "paper" | "notebook" | "talk" | "poster" -- scales every element.
# style:   "whitegrid" | "ticks" | "white" | "darkgrid".
# font_scale: multiplies *all* font sizes on top of FONT_SIZE below.
THEME = dict(context="notebook", style="whitegrid", font_scale=1.0)

# Font family. Set to e.g. ["Times New Roman"] for a serif paper figure; a
# missing family silently falls back to the default.
FONT_FAMILY = None          # None = leave seaborn's default alone

FONT_SIZE = {
    "panel_title": 11,      # "(b) Distance computations vs k"
    "axis_label": 10,       # x/y axis labels
    "tick_label": 9,        # numbers along the axes
    "legend": 10,           # the shared legend under the figure
    "suptitle": 12,         # the figure title on top
    "table": 8.5,           # the SSE table in panel (a)
}

# Raw matplotlib rcParams layered on top of the seaborn style.
RC_OVERRIDES = {
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#CCCCCC",   # left/bottom spine colour
    "axes.linewidth": 0.8,
    "grid.color": "#8A8A8A",
    "grid.linewidth": 0.6,
    "grid.alpha": 0.25,
    "savefig.facecolor": "white",
}

GRID_AXIS = "both"          # "y" | "x" | "both" | "none" -- which gridlines show
HIDE_SPINES = ("top", "right")   # spines removed from the line panels

# --- figure geometry -----------------------------------------------------
LAYOUT = dict(
    figsize=(14, 4.6),      # (width, height) in inches for the whole figure
    dpi=200,                # raster resolution of the .png (the .pdf is vector)
    wspace=1.5,             # horizontal gap between the three panels
    suptitle_y=1.02,        # vertical position of the figure title
    legend_y=-0.04,         # vertical position of the shared legend
    legend_ncol=None,       # None = one row; set an int to wrap the legend
)

# --- mark styling --------------------------------------------------------
LINE = dict(
    linewidth=1.8,          # thickness of the connecting line
    marker_edgecolor="none",  # "none" = flat marks; "white" = ringed marks
    err_capsize=3,          # width of the error-bar caps
    err_capthick=1.2,
    err_linewidth=1.2,
)

# --- table styling (panel a) --------------------------------------------
TABLE = dict(
    scale=1e6,              # SSE is divided by this before printing
    decimals=4,             # digits after the point for mean and std
    bbox=[0, 0, 1, 0.88],   # (x, y, width, height) of the table inside the axes
    row_height=0.16,
    header_facecolor="#F0F0F0",
    edgecolor="#CCCCCC",
    edge_linewidth=0.6,
    missing="---",          # printed when a cell has no data
)

# --- all user-visible strings -------------------------------------------
X_LABEL = "number of clusters $k$"

PANEL = {
    "table": dict(
        title="(a)  SSE at convergence  ($\\times 10^{6}$, mean ± std)",
        row_header="Algorithm",
        col_header="$k$={k}",    # {k} is filled in per column
    ),
    "distance": dict(
        title="(b)  Distance computations vs $k$",
        ylabel="distance computations ($\\times 10^{6}$)",
        scale=1e6,               # values are divided by this before plotting
        logy=True,
    ),
    "runtime": dict(
        title="(c)  Runtime vs $k$",
        ylabel="runtime (s)",
        scale=1.0,
        logy=False,
    ),
}

# {dataset} and {n_trials} are filled in from the summary CSV.
SUPTITLE = ("geokmeans algorithms on {dataset}  —  mean ± std over "
            "{n_trials} trials, shared k-means++ initialisation")

# ============================================================================
# END OF CONFIGURATION -- plotting logic below.
# ============================================================================


def _apply_theme():
    """Install the seaborn theme plus our rcParam overrides.

    Called from main() rather than at import time so that importing this
    module does not mutate a caller's global matplotlib state.
    """
    rc = dict(RC_OVERRIDES)
    rc.update({
        "axes.titlesize": FONT_SIZE["panel_title"],
        "axes.labelsize": FONT_SIZE["axis_label"],
        "xtick.labelsize": FONT_SIZE["tick_label"],
        "ytick.labelsize": FONT_SIZE["tick_label"],
        "legend.fontsize": FONT_SIZE["legend"],
        "figure.titlesize": FONT_SIZE["suptitle"],
    })
    if FONT_FAMILY:
        rc["font.family"] = FONT_FAMILY
    sns.set_theme(rc=rc, **THEME)


def _algos(df):
    return [a for a in ALGO_ORDER if a in set(df["algorithm"])]


def _line_panel(ax, df, mean_col, std_col, spec, ks):
    """One line panel: seaborn draws the series, matplotlib the error caps.

    seaborn aggregates raw observations itself and has no parameter for a
    pre-computed std column, so the marks come from `sns.lineplot` with
    `errorbar=None` and the caps are overlaid with `ax.errorbar(fmt="none")`
    in the matching series colour.
    """
    algos = _algos(df)
    if not algos:
        return False
    labels = [LABEL[a] for a in algos]
    palette = {LABEL[a]: COLOR[a] for a in algos}
    markers = {LABEL[a]: MARKER[a] for a in algos}
    scale = spec.get("scale", 1.0)

    data = df[df.algorithm.isin(algos)].copy()
    data["series"] = data["algorithm"].map(LABEL)
    data["_mean"] = data[mean_col] / scale
    data["_std"] = data[std_col] / scale
    data = data.sort_values(["series", "k"])

    sns.lineplot(
        data=data, x="k", y="_mean", hue="series", style="series",
        hue_order=labels, style_order=labels,
        palette=palette, markers=markers, dashes=False,   # all lines solid
        linewidth=LINE["linewidth"],
        markeredgecolor=LINE["marker_edgecolor"],
        errorbar=None, legend=True, zorder=3, ax=ax,
    )

    # lineplot draws one Line2D per series and takes a single markersize, so
    # the per-algorithm sizes are applied afterwards, before the error bars
    # add lines of their own.
    size_by_label = {LABEL[a]: MARKER_SIZE.get(a, 7) for a in algos}
    for line in ax.lines:
        if line.get_label() in size_by_label:
            line.set_markersize(size_by_label[line.get_label()])

    for algo in algos:
        sub = data[data.series == LABEL[algo]]
        ax.errorbar(sub["k"], sub["_mean"], yerr=sub["_std"].fillna(0.0),
                    fmt="none", ecolor=COLOR[algo],
                    elinewidth=LINE["err_linewidth"],
                    capsize=LINE["err_capsize"], capthick=LINE["err_capthick"],
                    zorder=2)

    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(spec["ylabel"])
    ax.set_title(spec["title"], fontsize=FONT_SIZE["panel_title"],
                 fontweight="bold", loc="left")
    if spec.get("logy"):
        ax.set_yscale("log")
    ax.grid(False)
    if GRID_AXIS in ("y", "both"):
        ax.yaxis.grid(True)
    if GRID_AXIS in ("x", "both"):
        ax.xaxis.grid(True)
    ax.set_axisbelow(True)
    ax.set_xticks(ks)
    sns.despine(ax=ax, top="top" in HIDE_SPINES, right="right" in HIDE_SPINES,
                left="left" in HIDE_SPINES, bottom="bottom" in HIDE_SPINES)
    return True


def _table_panel(ax, df, ks):
    """Panel (a): SSE per algorithm and k, as a table.

    Stays plain matplotlib -- seaborn has no table artist. Row labels are
    inked in the algorithm's own colour so the table and the line panels
    share one identity scheme.
    """
    ax.axis("off")
    algos = _algos(df)
    spec = PANEL["table"]
    cols = [spec["row_header"]] + [spec["col_header"].format(k=int(k))
                                   for k in ks]
    rows = []
    for algo in algos:
        cells = [LABEL[algo]]
        for k in ks:
            r = df[(df.algorithm == algo) & (df.k == k)]
            if r.empty:
                cells.append(TABLE["missing"])
            else:
                m = r["sse_mean"].iloc[0] / TABLE["scale"]
                sd = r["sse_std"].iloc[0] / TABLE["scale"]
                d = TABLE["decimals"]
                cells.append(f"{m:.{d}f}\n± {sd:.{d}f}")
        rows.append(cells)

    tbl = ax.table(cellText=rows, colLabels=cols, cellLoc="center",
                   loc="center", bbox=TABLE["bbox"])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(FONT_SIZE["table"])
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(TABLE["edge_linewidth"])
        cell.set_edgecolor(TABLE["edgecolor"])
        if r == 0:                                   # header row
            cell.set_facecolor(TABLE["header_facecolor"])
            cell.set_text_props(fontweight="bold")
        elif c == 0:                                 # algorithm name column
            cell.set_text_props(color=COLOR[algos[r - 1]], fontweight="bold")
        cell.set_height(TABLE["row_height"])
    ax.set_title(spec["title"], fontsize=FONT_SIZE["panel_title"],
                 fontweight="bold", loc="left")


def _shared_legend(fig, ax):
    """Lift one panel's seaborn legend up to a single figure-level legend."""
    handles, labels = ax.get_legend_handles_labels()
    for axis in fig.axes:                     # drop the per-panel legends
        if axis.get_legend() is not None:
            axis.get_legend().remove()
    if not handles:
        return
    fig.legend(handles, labels, loc="lower center",
               ncol=LAYOUT["legend_ncol"] or len(labels), frameon=False,
               bbox_to_anchor=(0.5, LAYOUT["legend_y"]),
               fontsize=FONT_SIZE["legend"])


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

    _apply_theme()
    fig = plt.figure(figsize=LAYOUT["figsize"])
    gs = fig.add_gridspec(1, 6, wspace=LAYOUT["wspace"])
    ax_t = fig.add_subplot(gs[0, 0:2])
    ax_d = fig.add_subplot(gs[0, 2:4])
    ax_r = fig.add_subplot(gs[0, 4:6])

    _table_panel(ax_t, df, ks)
    _line_panel(ax_d, df, "n_distance_calculations_mean",
                "n_distance_calculations_std", PANEL["distance"], ks)
    _line_panel(ax_r, df, "wall_clock_seconds_mean", "wall_clock_seconds_std",
                PANEL["runtime"], ks)

    _shared_legend(fig, ax_r)
    fig.suptitle(SUPTITLE.format(dataset=name, n_trials=n_trials),
                 fontsize=FONT_SIZE["suptitle"], fontweight="bold",
                 y=LAYOUT["suptitle_y"])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=LAYOUT["dpi"], bbox_inches="tight",
                facecolor="white")
    fig.savefig(args.out.with_suffix(".pdf"), bbox_inches="tight",
                facecolor="white")
    print(f"✓ Generated: {args.out}")
    print(f"✓ Generated: {args.out.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
