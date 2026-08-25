"""Figure 1: Geometric-k-means vs scikit-learn Lloyd vs R k-means on twitter.

Three panels, all with k on the x-axis:

  (a) line  -- mean distance computations vs k, log scale
  (b) line  -- mean runtime vs k
  (c) bars  -- mean total energy (RAPL package + DRAM) vs k

Reads only the tidy summary written by bench_twitter.py and plots means with
std error bars. Nothing is recomputed here.

Drawing is done with seaborn (`sns.lineplot` / `sns.barplot`) on top of a
seaborn theme; only the error-bar caps and the annotations fall back to raw
matplotlib, because seaborn has no API for *pre-computed* mean/std columns.

The five-panel sensor figure that used to live here is now
scripts/plot_e1e2.py; the paper's Figure 1 reports twitter only.

Usage:
    python scripts/plot_results.py
    python scripts/plot_results.py --summary results/reports/e3_twitter_summary_dryrun.csv

To restyle the figure -- titles, labels, font sizes, figure size, colours,
spacing -- edit the CONFIGURATION block below and nothing else.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

# ============================================================================
# CONFIGURATION -- everything a user would want to tweak lives here.
# ============================================================================

# --- which implementations, in which order -------------------------------
# Fixed display order. Colour follows the implementation, never its rank, so a
# sweep that omits one of the three does not repaint the survivors.
IMPL_ORDER = ["geokmeans_geo", "sklearn_lloyd", "r_hartigan_wong",
              "r_lloyd", "r_macqueen", "r_forgy"]

# Legend / series names. Edit the right-hand side to rename a series.
# Anything between $...$ is LaTeX-style mathtext.
IMPL_LABEL = {
    "geokmeans_geo": "Geometric-$k$-means",
    "sklearn_lloyd": "scikit-learn (Lloyd)",
    "r_hartigan_wong": "R stats::kmeans (Hartigan-Wong)",
    "r_lloyd": "R stats::kmeans (Lloyd)",
    "r_macqueen": "R stats::kmeans (MacQueen)",
    "r_forgy": "R stats::kmeans (Forgy)",
}

# Okabe-Ito vermillion / blue / bluish-green. Verified with the dataviz
# validator against a white surface: lightness band, chroma floor, all-pairs
# CVD separation (worst 11.0 deutan), normal-vision separation and contrast
# all pass. Markers and hatches carry the same identity, so nothing is
# distinguished by colour alone -- if you swap a colour, re-validate, and keep
# the marker/hatch distinct too.
IMPL_COLOR = {"geokmeans_geo": "#D55E00", "sklearn_lloyd": "#0072B2",
              "r_hartigan_wong": "#009E73", "r_lloyd": "#009E73",
              "r_macqueen": "#009E73", "r_forgy": "#009E73"}
IMPL_MARKER = {"geokmeans_geo": "o", "sklearn_lloyd": "s",     # line panels
               "r_hartigan_wong": "^", "r_lloyd": "^",
               "r_macqueen": "^", "r_forgy": "^"}
IMPL_HATCH = {"geokmeans_geo": "", "sklearn_lloyd": "//",      # bar panel
              "r_hartigan_wong": "xx", "r_lloyd": "xx",
              "r_macqueen": "xx", "r_forgy": "xx"}

# Text and axis-line inks. Keep both dark enough to read on white.
INK, MUTED = "#1A1A1A", "#5A5A5A"

# --- seaborn theme -------------------------------------------------------
# context: "paper" | "notebook" | "talk" | "poster" -- scales every element.
# style:   "whitegrid" | "ticks" | "white" | "darkgrid".
# font_scale: multiplies *all* font sizes on top of FONT_SIZE below; the
#             fastest single knob for "make the text bigger".
THEME = dict(context="notebook", style="whitegrid", font_scale=1.0)

# Font family. Set to e.g. ["Times New Roman"] or ["DejaVu Serif"] for a
# serif paper figure; a missing family silently falls back to the default.
FONT_FAMILY = None          # None = leave seaborn's default alone

# Point sizes for each kind of text. font_scale above multiplies these.
FONT_SIZE = {
    "panel_title": 11,      # "(a) Distance computations vs k"
    "axis_label": 10,       # x/y axis labels
    "tick_label": 9,        # numbers along the axes
    "legend": 10,           # the shared legend under the figure
    "suptitle": 12.5,       # the two-line figure title on top
    "footnote": 8,          # the methodology notes under the figure
    "placeholder": 9,       # "no energy data in this run" message
}

# Raw matplotlib rcParams layered on top of the seaborn style. Grid, spines
# and tick colours are set here.
RC_OVERRIDES = {
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#CCCCCC",   # left/bottom spine colour
    "axes.linewidth": 0.8,
    "grid.color": "#8A8A8A",
    "grid.linewidth": 0.6,
    "grid.alpha": 0.25,
    "axes.labelcolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "savefig.facecolor": "white",
}

GRID_AXIS = "y"             # "y" | "x" | "both" | "none" -- which gridlines show
HIDE_SPINES = ("top", "right")   # spines removed from every panel

# --- figure geometry -----------------------------------------------------
LAYOUT = dict(
    figsize=(15, 4.4),      # (width, height) in inches for the whole figure
    dpi=200,                # raster resolution of the .png (the .pdf is vector)
    wspace=0.32,            # horizontal gap between the three panels
    suptitle_y=1.10,        # vertical position of the figure title (axes-free coords)
    legend_y=-0.10,         # vertical position of the shared legend
    legend_ncol=None,       # None = one row; set an int to wrap the legend
    footnote_xy=(0.075, -0.185),   # (x, y) of the methodology notes, figure coords
)

# --- mark styling --------------------------------------------------------
# Line panels (a) and (b).
LINE = dict(
    linewidth=2.0,          # thickness of the connecting line
    markersize=8,           # data-point marker size (>= 8 keeps it legible)
    marker_edgecolor="none",  # "none" = flat marks; "white" = ringed marks
    err_capsize=3,          # width of the error-bar caps
    err_capthick=1.2,
    err_linewidth=1.2,
)

# Bar panel (c).
BAR = dict(
    group_width=0.78,       # fraction of the k-slot occupied by a whole group
    gap=0.06,               # 6% of each bar's slot becomes a surface gap
    edgecolor="white",      # the 2 px ring that separates adjacent bars
    edge_linewidth=1.2,
    err_capsize=3,
    err_capthick=1.0,
    err_linewidth=1.0,
)

# --- all user-visible strings -------------------------------------------
X_LABEL = "number of clusters $k$"

PANEL = {
    "distance": dict(
        title="(a)  Distance computations vs $k$",
        ylabel="distance computations (log scale)",
        logy=True,          # counts span orders of magnitude
    ),
    "runtime": dict(
        title="(b)  Runtime vs $k$",
        ylabel="runtime (s)",
        logy=False,
    ),
    "energy": dict(
        title="(c)  Total energy vs $k$",
        ylabel="total energy (J)",
        logy=False,
    ),
}

# {dataset} and {n_trials} are filled in from the summary CSV.
SUPTITLE = ("Geometric-$k$-means vs scikit-learn Lloyd vs R stats::kmeans on "
            "{dataset}\nmean ± std over {n_trials} trials, shared k-means++ "
            "initialisation, single-threaded")

# Methodology notes printed under the figure. Drop an entry to remove it; the
# energy note is only shown when the run actually carries energy data.
FOOTNOTE = {
    "distance": ("(a)  distance computations are measured for geokmeans; "
                 "scikit-learn and R expose no counter, so theirs is the "
                 "analytic $(\\mathrm{iter}+1)\\,n\\,k$ — exact for "
                 "Lloyd, an upper bound for Hartigan-Wong."),
    "energy": ("(c)  Intel RAPL, package + DRAM domains, sampled around the "
               "fit only. R is measured inside its own process, so its "
               "startup and data load are excluded as well."),
}

PLACEHOLDER_TEXT = ("no energy data in this run\n\n"
                    "RAPL was unavailable or --no-energy was set;\n"
                    "check with  python scripts/energy.py --check")

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


def _impls(df, mean_col=None):
    """Implementations present in the summary, in IMPL_ORDER.

    With `mean_col`, drops the ones that have no data for that metric (e.g. R
    in a run where energy was not measured). Unrecognised implementations are
    appended rather than silently dropped.
    """
    present = set(df["implementation"])
    known = [i for i in IMPL_ORDER if i in present]
    impls = known + sorted(present - set(known))
    if mean_col is not None:
        impls = [i for i in impls
                 if not df.loc[df.implementation == i, mean_col].isna().all()]
    return impls


def _label(impl):
    return IMPL_LABEL.get(impl, impl)


def _series_style(impls):
    """Palette / marker / hatch dicts keyed by *display label*.

    seaborn keys its `palette`, `markers` and legend by the values in the hue
    column, so the frame carries pretty labels and these follow suit.
    """
    labels = [_label(i) for i in impls]
    palette = {_label(i): IMPL_COLOR.get(i, "#666666") for i in impls}
    markers = {_label(i): IMPL_MARKER.get(i, "o") for i in impls}
    hatches = {_label(i): IMPL_HATCH.get(i, "") for i in impls}
    return labels, palette, markers, hatches


def _style(ax, title, ylabel, xlabel=X_LABEL):
    """Shared per-panel cosmetics: labels, title, grid, spines."""
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=FONT_SIZE["panel_title"], fontweight="bold",
                 loc="left", color=INK)
    ax.grid(False)
    if GRID_AXIS in ("y", "both"):
        ax.yaxis.grid(True)
    if GRID_AXIS in ("x", "both"):
        ax.xaxis.grid(True)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=FONT_SIZE["tick_label"])
    sns.despine(ax=ax, top="top" in HIDE_SPINES, right="right" in HIDE_SPINES,
                left="left" in HIDE_SPINES, bottom="bottom" in HIDE_SPINES)


def _line_panel(ax, df, mean_col, std_col, spec, ks):
    """One line panel: seaborn draws the series, matplotlib the error caps.

    seaborn aggregates raw observations itself and has no parameter for a
    pre-computed std column, so the marks come from `sns.lineplot` with
    `errorbar=None` and the caps are overlaid with `ax.errorbar(fmt="none")`
    in the matching series colour.
    """
    impls = _impls(df, mean_col)
    if not impls:
        return False
    labels, palette, markers, _ = _series_style(impls)

    data = df[df.implementation.isin(impls)].copy()
    data["series"] = data["implementation"].map(_label)
    data = data.sort_values(["series", "k"])

    sns.lineplot(
        data=data, x="k", y=mean_col, hue="series", style="series",
        hue_order=labels, style_order=labels,
        palette=palette, markers=markers, dashes=False,   # all lines solid
        linewidth=LINE["linewidth"], markersize=LINE["markersize"],
        markeredgecolor=LINE["marker_edgecolor"],
        errorbar=None, legend=True, zorder=3, ax=ax,
    )

    for label in labels:
        sub = data[data.series == label]
        ax.errorbar(sub["k"], sub[mean_col], yerr=sub[std_col].fillna(0.0),
                    fmt="none", ecolor=palette[label],
                    elinewidth=LINE["err_linewidth"],
                    capsize=LINE["err_capsize"], capthick=LINE["err_capthick"],
                    zorder=2)

    _style(ax, spec["title"], spec["ylabel"])
    if spec.get("logy"):
        ax.set_yscale("log")
    ax.set_xticks(ks)
    return True


def _bar_panel(ax, df, mean_col, std_col, spec, ks):
    """Grouped bars: one group per k, one bar per implementation.

    `gap` leaves a surface gap between adjacent bars so the group reads as a
    set of marks rather than one striped block. Bar centres are read back off
    the patches seaborn drew, so the error caps land correctly whatever
    seaborn's dodge arithmetic does.
    """
    impls = _impls(df, mean_col)
    if not impls:
        return False
    labels, palette, _, hatches = _series_style(impls)

    data = df[df.implementation.isin(impls)].copy()
    data["series"] = data["implementation"].map(_label)

    sns.barplot(
        data=data, x="k", y=mean_col, hue="series",
        hue_order=labels, order=ks, palette=palette,
        estimator="mean", errorbar=None,          # one row per cell already
        width=BAR["group_width"], gap=BAR["gap"], dodge=True,
        edgecolor=BAR["edgecolor"], linewidth=BAR["edge_linewidth"],
        zorder=3, ax=ax,
    )

    # seaborn emits one BarContainer per hue level, in hue_order.
    for label, container in zip(labels, ax.containers):
        sub = data[data.series == label].set_index("k")
        errs, xs, ys = [], [], []
        for k, patch in zip(ks, container):
            patch.set_hatch(hatches[label])
            mean = sub[mean_col].get(k, np.nan)
            std = sub[std_col].get(k, np.nan)
            if pd.isna(mean):
                continue
            xs.append(patch.get_x() + patch.get_width() / 2)
            ys.append(mean)
            errs.append(0.0 if pd.isna(std) else std)
        if xs:
            ax.errorbar(xs, ys, yerr=errs, fmt="none", ecolor=MUTED,
                        elinewidth=BAR["err_linewidth"],
                        capsize=BAR["err_capsize"],
                        capthick=BAR["err_capthick"], zorder=4)

    _style(ax, spec["title"], spec["ylabel"])
    if spec.get("logy"):
        ax.set_yscale("log")
    # x is categorical here (one slot per k), so the ticks are 0..len(ks)-1
    # and only their labels carry the k values.
    ax.set_xticks(range(len(ks)), labels=[str(int(k)) for k in ks])
    return True


def _energy_unavailable(ax, spec):
    """Stand-in for panel (c) when the run measured no energy."""
    _style(ax, spec["title"], "", xlabel="")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    sns.despine(ax=ax, left=True, bottom=True, top=True, right=True)
    ax.text(0.5, 0.5, PLACEHOLDER_TEXT, transform=ax.transAxes,
            ha="center", va="center", fontsize=FONT_SIZE["placeholder"],
            color=MUTED, style="italic")


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
    ap.add_argument("--summary", "--twitter", dest="summary", type=Path,
                    default=base / "results/reports/e3_twitter_summary.csv",
                    help="tidy mean/std CSV from bench_twitter.py")
    ap.add_argument("--out", type=Path,
                    default=base / "results/figures/benchmark_summary.png")
    args = ap.parse_args()

    if not args.summary.exists():
        raise SystemExit(f"summary not found: {args.summary}\n"
                         f"Run scripts/bench_twitter.py first.")
    df = pd.read_csv(args.summary)
    if df.empty:
        raise SystemExit("no rows to plot")

    ks = sorted(df["k"].unique())
    name = ", ".join(sorted(df["dataset"].unique()))
    n_trials = int(df["n_trials"].max())
    has_energy = ("energy_total_joules_mean" in df.columns
                  and df["energy_total_joules_mean"].notna().any())

    _apply_theme()
    fig, axes = plt.subplots(1, 3, figsize=LAYOUT["figsize"])
    ax_d, ax_r, ax_e = axes

    # Counts span orders of magnitude, so the axis is log and the ticks stay
    # raw powers of ten. Dividing by a scale factor as well would print
    # "2 x 10^0" underneath a "(x 10^6)" axis label.
    _line_panel(ax_d, df, "n_distance_calculations_mean",
                "n_distance_calculations_std", PANEL["distance"], ks)

    _line_panel(ax_r, df, "wall_clock_seconds_mean", "wall_clock_seconds_std",
                PANEL["runtime"], ks)

    if has_energy:
        _bar_panel(ax_e, df, "energy_total_joules_mean",
                   "energy_total_joules_std", PANEL["energy"], ks)
    else:
        _energy_unavailable(ax_e, PANEL["energy"])

    # One legend for the whole figure: identity is never colour-alone, since
    # every series also carries a distinct marker and bar hatch.
    _shared_legend(fig, ax_r)

    notes = [FOOTNOTE["distance"]]
    if has_energy:
        notes.append(FOOTNOTE["energy"])
    fig.text(*LAYOUT["footnote_xy"], "\n".join(notes), ha="left", va="top",
             fontsize=FONT_SIZE["footnote"], style="italic", color=MUTED)

    fig.suptitle(SUPTITLE.format(dataset=name, n_trials=n_trials),
                 fontsize=FONT_SIZE["suptitle"], fontweight="bold",
                 y=LAYOUT["suptitle_y"], color=INK)
    fig.subplots_adjust(wspace=LAYOUT["wspace"])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=LAYOUT["dpi"], bbox_inches="tight",
                facecolor="white")
    fig.savefig(args.out.with_suffix(".pdf"), bbox_inches="tight",
                facecolor="white")
    print(f"✓ Generated: {args.out}")
    print(f"✓ Generated: {args.out.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
