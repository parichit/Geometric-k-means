"""Figure 1: Geometric-k-means vs scikit-learn Lloyd vs R k-means on twitter.

Three panels, all with k on the x-axis:

  (a) line  -- mean distance computations vs k, log scale
  (b) line  -- mean runtime vs k
  (c) bars  -- mean total energy (RAPL package + DRAM) vs k

Reads only the tidy summary written by bench_twitter.py and plots means with
std error bars. Nothing is recomputed here.

The five-panel sensor figure that used to live here is now
scripts/plot_e1e2.py; the paper's Figure 1 reports twitter only.

Usage:
    python scripts/plot_results.py
    python scripts/plot_results.py --summary results/reports/e3_twitter_summary_dryrun.csv
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Fixed display order. Colour follows the implementation, never its rank, so a
# sweep that omits one of the three does not repaint the survivors.
IMPL_ORDER = ["geokmeans_geo", "sklearn_lloyd", "r_hartigan_wong",
              "r_lloyd", "r_macqueen", "r_forgy"]

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
# distinguished by colour alone.
IMPL_COLOR = {"geokmeans_geo": "#D55E00", "sklearn_lloyd": "#0072B2",
              "r_hartigan_wong": "#009E73", "r_lloyd": "#009E73",
              "r_macqueen": "#009E73", "r_forgy": "#009E73"}
IMPL_MARKER = {"geokmeans_geo": "o", "sklearn_lloyd": "s",
               "r_hartigan_wong": "^", "r_lloyd": "^",
               "r_macqueen": "^", "r_forgy": "^"}
IMPL_HATCH = {"geokmeans_geo": "", "sklearn_lloyd": "//",
              "r_hartigan_wong": "xx", "r_lloyd": "xx",
              "r_macqueen": "xx", "r_forgy": "xx"}

INK, MUTED = "#1A1A1A", "#5A5A5A"


def _impls(df):
    present = set(df["implementation"])
    known = [i for i in IMPL_ORDER if i in present]
    # anything unrecognised still gets plotted rather than silently dropped
    return known + sorted(present - set(known))


def _label(impl):
    return IMPL_LABEL.get(impl, impl)


def _style(ax, title, ylabel, ks):
    ax.set_xlabel("number of clusters $k$")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, fontweight="bold", loc="left", color=INK)
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=9)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#CCCCCC")


def _line_panel(ax, df, mean_col, std_col, ylabel, title, ks,
                scale=1.0, log=False, note=None):
    for impl in _impls(df):
        sub = df[df.implementation == impl].sort_values("k")
        if sub.empty or sub[mean_col].isna().all():
            continue
        ax.errorbar(sub["k"], sub[mean_col] / scale,
                    yerr=sub[std_col].fillna(0.0) / scale,
                    label=_label(impl), color=IMPL_COLOR.get(impl, "#666666"),
                    marker=IMPL_MARKER.get(impl, "o"), markersize=8,
                    linewidth=2.0, capsize=3, capthick=1.2, zorder=3)
    _style(ax, title, ylabel, ks)
    if log:
        ax.set_yscale("log")
    ax.set_xticks(ks)
    if note:
        ax.text(0.0, -0.30, note, transform=ax.transAxes, fontsize=7.5,
                style="italic", color=MUTED, ha="left", va="top")


def _bar_panel(ax, df, mean_col, std_col, ylabel, title, ks, note=None):
    """Grouped bars: one group per k, one bar per implementation.

    A 2 px surface gap separates adjacent bars so the group reads as a set of
    marks rather than one striped block.
    """
    impls = [i for i in _impls(df)
             if not df[df.implementation == i][mean_col].isna().all()]
    if not impls:
        return False
    n = len(impls)
    idx = np.arange(len(ks), dtype=float)
    group_w = 0.78
    bar_w = group_w / n
    for j, impl in enumerate(impls):
        sub = df[df.implementation == impl].set_index("k")
        means = [sub[mean_col].get(k, np.nan) for k in ks]
        errs = [sub[std_col].get(k, 0.0) if pd.notna(sub[std_col].get(k, 0.0))
                else 0.0 for k in ks]
        offset = (j - (n - 1) / 2) * bar_w
        ax.bar(idx + offset, means, bar_w * 0.94,  # 6% of the slot = the gap
               yerr=errs, label=_label(impl),
               color=IMPL_COLOR.get(impl, "#666666"),
               hatch=IMPL_HATCH.get(impl, ""), edgecolor="white",
               linewidth=1.2, capsize=3,
               error_kw=dict(ecolor=MUTED, capthick=1.0, lw=1.0), zorder=3)
    _style(ax, title, ylabel, ks)
    ax.set_xticks(idx)
    ax.set_xticklabels([str(int(k)) for k in ks])
    if note:
        ax.text(0.0, -0.30, note, transform=ax.transAxes, fontsize=7.5,
                style="italic", color=MUTED, ha="left", va="top")
    return True


def _energy_unavailable(ax, title):
    _style(ax, title, "", [])
    ax.set_xlabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(False)
    ax.text(0.5, 0.5,
            "no energy data in this run\n\n"
            "RAPL was unavailable or --no-energy was set;\n"
            "check with  python scripts/energy.py --check",
            transform=ax.transAxes, ha="center", va="center",
            fontsize=9, color=MUTED, style="italic")


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

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    ax_d, ax_r, ax_e = axes

    # Counts span orders of magnitude, so the axis is log and the ticks stay
    # raw powers of ten. Dividing by a scale factor as well would print
    # "2 x 10^0" underneath a "(x 10^6)" axis label.
    _line_panel(ax_d, df, "n_distance_calculations_mean",
                "n_distance_calculations_std",
                "distance computations (log scale)",
                "(a)  Distance computations vs $k$", ks, log=True)

    _line_panel(ax_r, df, "wall_clock_seconds_mean", "wall_clock_seconds_std",
                "runtime (s)", "(b)  Runtime vs $k$", ks)

    if has_energy:
        _bar_panel(ax_e, df, "energy_total_joules_mean",
                   "energy_total_joules_std", "total energy (J)",
                   "(c)  Total energy vs $k$", ks)
    else:
        _energy_unavailable(ax_e, "(c)  Total energy vs $k$")

    # One legend for the whole figure: identity is never colour-alone, since
    # every series also carries a distinct marker and bar hatch.
    handles, labels = ax_r.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels),
               frameon=False, bbox_to_anchor=(0.5, -0.10), fontsize=10)

    notes = ["(a)  distance computations are measured for geokmeans; "
             "scikit-learn and R expose no counter, so theirs is the analytic "
             "$(\\mathrm{iter}+1)\\,n\\,k$ \u2014 exact for Lloyd, an upper "
             "bound for Hartigan-Wong."]
    if has_energy:
        notes.append("(c)  Intel RAPL, package + DRAM domains, sampled around "
                     "the fit only. R is measured inside its own process, so "
                     "its startup and data load are excluded as well.")
    fig.text(0.075, -0.185, "\n".join(notes), ha="left", va="top",
             fontsize=8, style="italic", color=MUTED)

    fig.suptitle(
        f"Geometric-$k$-means vs scikit-learn Lloyd vs R stats::kmeans on "
        f"{name}\nmean ± std over {n_trials} trials, shared k-means++ "
        f"initialisation, single-threaded",
        fontsize=12.5, fontweight="bold", y=1.10, color=INK)
    fig.subplots_adjust(wspace=0.32)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(args.out.with_suffix(".pdf"), bbox_inches="tight",
                facecolor="white")
    print(f"✓ Generated: {args.out}")
    print(f"✓ Generated: {args.out.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
