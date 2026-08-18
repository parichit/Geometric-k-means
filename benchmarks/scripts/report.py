"""Generate tables and reports from benchmark results.

Reads raw CSV results from results/raw/ and generates the outputs specified in
config.yaml (LaTeX tables, markdown reports).

Usage:
    python report.py
    python report.py --dry-run   # read the *_dryrun.csv results
"""
import argparse
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Any, Dict, List

# Metrics that are aggregated across trials. Anything else is metadata.
AGGREGATED = [
    "wall_clock_seconds",
    "n_distance_calculations",
    "n_iterations",
    "sse",
    "sse_ratio_vs_lloyd",
]

ALGO_ORDER = ["lloyd", "hamerly", "annulus", "exponion", "geo"]


def _ordered(df: pd.DataFrame, col: str = "algorithm") -> pd.DataFrame:
    """Sort by the canonical algorithm order, unknown names last."""
    rank = {a: i for i, a in enumerate(ALGO_ORDER)}
    return df.assign(_r=df[col].map(lambda a: rank.get(a, len(rank)))) \
             .sort_values(["_r"] + [c for c in ("k",) if c in df.columns]) \
             .drop(columns="_r")


def aggregate_trials(df: pd.DataFrame, group_cols: List[str]) -> pd.DataFrame:
    """Collapse trials to a per-cell median (used for the LaTeX tables)."""
    metrics = [c for c in AGGREGATED if c in df.columns]
    agg = df.groupby(group_cols, as_index=False)[metrics].median()
    counts = df.groupby(group_cols, as_index=False).size().rename(
        columns={"size": "n_trials"}
    )
    return agg.merge(counts, on=group_cols, how="left")


def summarise_mean_std(df: pd.DataFrame, group_cols: List[str]) -> pd.DataFrame:
    """Per-cell mean and sample standard deviation of every metric.

    Emitted as a tidy CSV (one row per cell, `<metric>_mean` / `<metric>_std`
    columns) so downstream plotting can read it without re-deriving anything.
    """
    metrics = [c for c in AGGREGATED + ["ari_vs_lloyd",
                                        "max_relative_sse_deviation"]
               if c in df.columns]
    g = df.groupby(group_cols, as_index=False)
    out = g[metrics].mean().rename(columns={m: f"{m}_mean" for m in metrics})
    std = g[metrics].std(ddof=1).rename(columns={m: f"{m}_std" for m in metrics})
    out = out.merge(std, on=group_cols, how="left")
    for m in metrics:                       # a single trial has no spread
        out[f"{m}_std"] = out[f"{m}_std"].fillna(0.0)
    n = g.size().rename(columns={"size": "n_trials"})
    return out.merge(n, on=group_cols, how="left")


def fmt_mean_std(mean: float, std: float, decimals: int = 3) -> str:
    return f"{mean:.{decimals}f} ± {std:.{decimals}f}"


def format_metric(col: str, val: float) -> str:
    """Render one metric value for a LaTeX cell."""
    if pd.isna(val):
        return "---"
    if col == "wall_clock_seconds":
        return f"{val:.3f}"
    if col == "n_distance_calculations":
        return f"{val / 1e6:.2f}"
    if col == "sse_ratio_vs_lloyd":
        return f"{val:.4f}"
    if col == "n_iterations":
        return f"{val:.0f}"
    return f"{val:.2f}"


def format_column_name(col: str, short: bool = False) -> str:
    """Format a metric name for a LaTeX table header."""
    mapping = {
        "wall_clock_seconds": ("Time (s)", "Time"),
        "n_distance_calculations": ("Distances ($\\times 10^6$)", "Dist."),
        "sse_ratio_vs_lloyd": ("SSE ratio", "SSE"),
        "n_iterations": ("Iters", "Iters"),
        "sse": ("SSE", "SSE"),
    }
    if col not in mapping:
        return col.replace("_", " ")
    return mapping[col][1 if short else 0]


def format_implementation_name(impl: str, short: bool = False) -> str:
    """Format an implementation name for display."""
    if impl.startswith("sklearn"):
        return "sklearn" if short else "scikit-learn (Lloyd)"
    if impl.startswith("geokmeans"):
        algo = impl.replace("geokmeans_", "")
        if short:
            return algo
        display = "Geometric" if algo == "geo" else algo.capitalize()
        return f"geokmeans ({display})"
    return impl


def _estimated_note(used: bool) -> List[str]:
    """Footnote for estimated cells.

    Written with plain LaTeX rather than `threeparttable` so the generated
    tables need nothing beyond `booktabs`.
    """
    if not used:
        return []
    return [
        "\\par\\smallskip",
        "{\\footnotesize $^{\\dagger}$ scikit-learn exposes no distance counter; "
        "this is the analytic count for Lloyd's, not a measurement.}",
    ]


def generate_latex_table_main(df: pd.DataFrame, config: Dict[str, Any],
                              output_config: Dict[str, Any],
                              output_path: Path) -> bool:
    """Generate the main comparison table (Table 1)."""
    dataset = output_config["dataset"]
    k = output_config["k"]
    columns = output_config["columns"]
    rows = output_config["rows"]

    subset = df[(df["dataset"] == dataset) & (df["k"] == k)]
    if subset.empty:
        available = sorted(
            {(d, int(kk)) for d, kk in zip(df["dataset"], df["k"])}
        )
        print(f"⚠️  No data for dataset={dataset!r} k={k}, skipping "
              f"{output_config['id']}. Available (dataset, k): {available}")
        return False

    n_trials = int(subset["n_trials"].max())
    trial_word = "trial" if n_trials == 1 else "trials"

    lines = [
        "\\begin{table}[ht]",
        "\\centering",
        f"\\caption{{Performance comparison on \\texttt{{{dataset}}} "
        f"($k={k}$, median of {n_trials} {trial_word})}}",
        "\\label{tab:comparison}",
        f"\\begin{{tabular}}{{{'l' + 'r' * len(columns)}}}",
        "\\toprule",
        " & ".join(["Implementation"] + [format_column_name(c) for c in columns])
        + " \\\\",
        "\\midrule",
    ]

    used_estimate = False
    for impl in rows:
        row = subset[subset["implementation"] == impl]
        if row.empty:
            print(f"⚠️  Missing data for {impl} on {dataset} k={k}")
            continue
        row = row.iloc[0]

        values = [format_implementation_name(impl)]
        for col in columns:
            cell = format_metric(col, row[col]) if col in row else "---"
            if col == "n_distance_calculations" and row.get("distances_estimated"):
                cell += "$^{\\dagger}$"
                used_estimate = True
            values.append(cell)

        lines.append(" & ".join(values) + " \\\\")
        if impl.startswith("sklearn"):
            lines.append("\\midrule")

    lines += ["\\bottomrule", "\\end{tabular}"]
    lines += _estimated_note(used_estimate)
    lines += ["\\end{table}"]

    write_lines(output_path, lines)
    return True


def generate_latex_table_appendix(df: pd.DataFrame, config: Dict[str, Any],
                                  output_config: Dict[str, Any],
                                  output_path: Path) -> bool:
    """Generate the full appendix table: all datasets x all k values.

    One tabular, with dataset and k as leading columns. The previous version
    emitted several `tabular` environments and bare `\\\\` line breaks directly
    inside one `table` float, which is not valid LaTeX and will not compile.
    """
    columns = output_config["columns"]

    lines = [
        "\\begin{table}[ht]",
        "\\centering",
        "\\caption{Full benchmark results: all datasets and $k$ values "
        "(median across trials)}",
        "\\label{tab:appendix_full}",
        "\\scriptsize",
        f"\\begin{{tabular}}{{ll{'r' * (len(columns) + 1)}}}",
        "\\toprule",
        " & ".join(
            ["Dataset", "$k$", "Impl."]
            + [format_column_name(c, short=True) for c in columns]
        )
        + " \\\\",
        "\\midrule",
    ]

    used_estimate = False
    for dataset in sorted(df["dataset"].unique()):
        for k in sorted(df[df["dataset"] == dataset]["k"].unique()):
            block = df[(df["dataset"] == dataset) & (df["k"] == k)]
            for i, (_, row) in enumerate(block.iterrows()):
                values = [
                    dataset.replace("_", "\\_") if i == 0 else "",
                    f"{int(k)}" if i == 0 else "",
                    format_implementation_name(row["implementation"], short=True),
                ]
                for col in columns:
                    cell = format_metric(col, row[col]) if col in row else "---"
                    if col == "n_distance_calculations" and row.get("distances_estimated"):
                        cell += "$^{\\dagger}$"
                        used_estimate = True
                    values.append(cell)
                lines.append(" & ".join(values) + " \\\\")
            lines.append("\\midrule")

    if lines[-1] == "\\midrule":
        lines.pop()

    lines += ["\\bottomrule", "\\end{tabular}"]
    lines += _estimated_note(used_estimate)
    lines += ["\\end{table}"]

    write_lines(output_path, lines)
    return True


def generate_exactness_report(df: pd.DataFrame, config: Dict[str, Any],
                              output_config: Dict[str, Any],
                              output_path: Path) -> bool:
    """E1 report: mean +/- std per algorithm and k, plus a tidy CSV.

    The CSV alongside the markdown is what the plotting script consumes.
    """
    reference = config["exactness_check"].get("reference", "geokmeans_lloyd")
    summary = summarise_mean_std(df, ["dataset", "algorithm", "k"])

    csv_path = output_path.with_suffix(".csv")
    _ordered(summary).to_csv(csv_path, index=False)
    print(f"✓ Generated: {csv_path}")

    n_trials = int(summary["n_trials"].max())
    lines = [
        "# E1: Exactness Verification Report\n",
        f"**{len(df)} runs across {df['algorithm'].nunique()} algorithms, "
        f"{n_trials} trials per (algorithm, k).**  ",
        f"All algorithms start from identical initial centroids; ratios are "
        f"against `{reference}` **at the same seed**.\n",
        f"Datasets: {', '.join(sorted(df['dataset'].unique()))}  ",
        f"k values: {', '.join(str(int(k)) for k in sorted(df['k'].unique()))}\n",
    ]

    # --- exactness verdict ------------------------------------------------
    cmp = df[df["implementation"] != reference]
    ari_fail = int((cmp["ari_vs_lloyd"] < 0.9999).sum())
    sse_fail = int((cmp["max_relative_sse_deviation"] > 1e-6).sum())
    verdict = "PASS" if (ari_fail == 0 and sse_fail == 0) else "FAIL"
    lines += [
        "## Exactness\n",
        f"- Comparisons against `{reference}`: **{len(cmp)}**",
        f"- ARI < 0.9999: **{ari_fail}**",
        f"- Relative SSE deviation > 1e-6: **{sse_fail}**",
        f"- Worst ARI: {cmp['ari_vs_lloyd'].min():.6f}" if len(cmp) else "",
        f"- Worst SSE deviation: {cmp['max_relative_sse_deviation'].max():.2e}"
        if len(cmp) else "",
        f"\n**Verdict: {verdict}**\n",
    ]

    # --- mean +/- std tables ---------------------------------------------
    for k in sorted(df["k"].unique()):
        sub = _ordered(summary[summary["k"] == k])
        lines += [
            f"\n## k = {int(k)}  (mean ± std over {n_trials} trials)\n",
            "| Algorithm | Iterations | Runtime (s) | SSE ratio | Distances (M) |",
            "|---|---|---|---|---|",
        ]
        for _, r in sub.iterrows():
            lines.append(
                f"| {r['algorithm']} "
                f"| {fmt_mean_std(r['n_iterations_mean'], r['n_iterations_std'], 1)} "
                f"| {fmt_mean_std(r['wall_clock_seconds_mean'], r['wall_clock_seconds_std'], 3)} "
                f"| {fmt_mean_std(r['sse_ratio_vs_lloyd_mean'], r['sse_ratio_vs_lloyd_std'], 6)} "
                f"| {fmt_mean_std(r['n_distance_calculations_mean'] / 1e6, r['n_distance_calculations_std'] / 1e6, 2)} |"
            )

    # --- SSE at convergence ----------------------------------------------
    lines += ["\n## SSE at convergence (mean ± std)\n",
              "| Algorithm | " + " | ".join(f"k={int(k)}"
                                            for k in sorted(df['k'].unique())) + " |",
              "|---" * (1 + df['k'].nunique()) + "|"]
    for algo in [a for a in ALGO_ORDER if a in set(summary["algorithm"])]:
        cells = []
        for k in sorted(df["k"].unique()):
            r = summary[(summary.algorithm == algo) & (summary.k == k)]
            cells.append(f"{r['sse_mean'].iloc[0]:.6e} ± {r['sse_std'].iloc[0]:.1e}"
                         if len(r) else "---")
        lines.append(f"| {algo} | " + " | ".join(cells) + " |")

    failures = cmp[(cmp["ari_vs_lloyd"] < 0.9999)
                   | (cmp["max_relative_sse_deviation"] > 1e-6)]
    if not failures.empty:
        lines.append("\n## Failures\n")
        cols = ["dataset", "algorithm", "k", "seed",
                "ari_vs_lloyd", "max_relative_sse_deviation"]
        lines.append(failures[cols].to_markdown(index=False))

    write_lines(output_path, [ln for ln in lines if ln != ""])
    return True


def write_lines(output_path: Path, lines: List[str]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")
    print(f"✓ Generated: {output_path}")


def read_results(path: Path, label: str) -> pd.DataFrame:
    """Read a results CSV, returning an empty frame if it is missing or blank."""
    if not path.exists():
        print(f"⚠️  {label} results not found at {path}")
        return pd.DataFrame()
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        print(f"⚠️  {label} results at {path} are empty")
        return pd.DataFrame()
    if df.empty:
        print(f"⚠️  {label} results at {path} contain no rows")
    return df


def main():
    parser = argparse.ArgumentParser(description="Generate benchmark reports")
    parser.add_argument("--config", type=Path,
                        default=Path(__file__).parent.parent / "config.yaml")
    parser.add_argument("--dry-run", action="store_true",
                        help="Read the *_dryrun.csv results")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)

    base_dir = args.config.parent.resolve()
    results_dir = base_dir / "results" / "raw"
    suffix = "_dryrun" if args.dry_run else ""

    print("=" * 80)
    print("GENERATING BENCHMARK REPORTS")
    print("=" * 80)

    e1 = read_results(results_dir / f"e1_exactness{suffix}.csv", "E1")
    e2 = read_results(results_dir / f"e2_comparison{suffix}.csv", "E2")

    e2_agg = (
        aggregate_trials(e2, ["dataset", "implementation", "k"])
        if not e2.empty else pd.DataFrame()
    )
    if not e2.empty:
        e2_summary = summarise_mean_std(e2, ["dataset", "algorithm", "k"])
        e2_csv = base_dir / "results" / "raw" / f"e2_summary{suffix}.csv"
        _ordered(e2_summary).to_csv(e2_csv, index=False)
        print(f"✓ Generated: {e2_csv}")

    generated = 0
    for output_spec in config["outputs"]:
        output_path = base_dir / output_spec["output_path"]
        if args.dry_run:
            output_path = output_path.with_name(
                f"{output_path.stem}_dryrun{output_path.suffix}"
            )
        output_type = output_spec["type"]

        if output_type == "latex_table":
            if e2_agg.empty:
                print(f"⚠️  Skipping {output_spec['id']}: no E2 results")
                continue
            if output_spec["id"] == "table1_main":
                ok = generate_latex_table_main(e2_agg, config, output_spec, output_path)
            elif output_spec["id"] == "appendix_full_matrix":
                ok = generate_latex_table_appendix(e2_agg, config, output_spec, output_path)
            else:
                print(f"⚠️  Unknown table id: {output_spec['id']}")
                continue
            generated += bool(ok)

        elif output_type == "markdown_report":
            if e1.empty:
                print(f"⚠️  Skipping {output_spec['id']}: no E1 results")
                continue
            generated += bool(
                generate_exactness_report(e1, config, output_spec, output_path)
            )

        else:
            print(f"⚠️  Unknown output type: {output_type}")

    print("\n" + "=" * 80)
    print(f"REPORT GENERATION COMPLETE ({generated} output(s) written)")
    print("=" * 80)


if __name__ == "__main__":
    main()
