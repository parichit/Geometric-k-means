"""Benchmark runner: executes E1 (exactness) and E2 (comparison) experiments.

Usage:
    python run.py --experiment E1   # Exactness verification
    python run.py --experiment E2   # Comparison against baselines
    python run.py --experiment all  # Both E1 and E2
    python run.py --dry-run         # Tiny smoke run over the real pipeline
"""
import sys
from pathlib import Path

# Must precede numpy/sklearn/geokmeans: see _threads.py for why.
sys.path.insert(0, str(Path(__file__).parent))
from _threads import pin_single_threaded, report_thread_limits  # noqa: E402

pin_single_threaded()

import argparse  # noqa: E402
import copy  # noqa: E402
import yaml  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from typing import Any, Dict, List  # noqa: E402
from tqdm import tqdm  # noqa: E402

from init_generator import generate_and_save_init  # noqa: E402
from utils import (  # noqa: E402
    Timer, compute_sse, adjusted_rand_index, load_dataset,
)

# Metrics recorded for every trial. peak_rss_mb is deliberately absent.
METRIC_COLUMNS = [
    "wall_clock_seconds",
    "n_distance_calculations",
    "n_iterations",
    "sse",
]


def import_geokmeans():
    """Import GeoKMeans with a diagnosis for the two common failure modes."""
    try:
        import geokmeans
    except ImportError as exc:
        raise ImportError(
            "geokmeans package not found. Install with:\n"
            "  cd geokmeans-py && pip install -e ."
        ) from exc

    if getattr(geokmeans, "__file__", None) is None:
        # A directory named `geokmeans` on sys.path (the repo's C++ source tree
        # is one) is picked up as an empty namespace package and shadows the
        # real install. Running from the repo root triggers this.
        raise ImportError(
            "`import geokmeans` resolved to a namespace package at "
            f"{list(getattr(geokmeans, '__path__', []))}, not the installed "
            "library. A directory named 'geokmeans' on sys.path is shadowing "
            "it -- run the benchmarks from the benchmarks/ directory, or "
            "install the package into the active environment:\n"
            "  cd geokmeans-py && pip install -e ."
        )

    from geokmeans import GeoKMeans

    return GeoKMeans


def run_geokmeans(X: np.ndarray, k: int, init_path: Path, algorithm: str,
                  tol: float, max_iter: int) -> Dict[str, Any]:
    """Run a geokmeans algorithm from explicit initial centroids.

    ``init_path`` is passed as a path rather than an array on purpose: the C++
    core reads centroid files directly, so nothing is serialised inside the
    timed region. Passing an array instead would make every geokmeans timing
    include a CSV write that the timed region should not contain.

    Returns
    -------
    result : dict
        Contains: labels, centroids, n_iterations, n_distance_calculations, sse
    """
    GeoKMeans = import_geokmeans()

    kmeans = GeoKMeans(
        n_clusters=k,
        algorithm=algorithm,
        init=str(init_path),
        tol=tol,
        max_iter=max_iter,
        random_state=0,  # unused when init is explicit
        verbose=False,
    )
    with Timer() as timer:
        kmeans.fit(X)

    return {
        "labels": kmeans.labels_,
        "centroids": kmeans.cluster_centers_,
        "n_iterations": int(kmeans.n_iter_),
        "n_distance_calculations": int(kmeans.n_distance_calculations_),
        "sse": compute_sse(X, kmeans.labels_, kmeans.cluster_centers_),
        "wall_clock_seconds": timer.wall_clock_seconds,
    }


def ensure_init(X: np.ndarray, dataset_name: str, k: int, seed: int,
                init_dir: Path) -> Path:
    """Return the shared init file for this (dataset, k, seed), creating it once."""
    init_file = init_dir / f"{dataset_name}_k{k}_seed{seed}.csv"
    if not init_file.exists():
        generate_and_save_init(X, dataset_name, k, seed, init_dir)
    return init_file


def run_single_trial(X: np.ndarray, dataset_name: str, implementation: str,
                     algorithm: str, k: int, seed: int, init_dir: Path,
                     tol: float, max_iter: int) -> Dict[str, Any]:
    """Run one benchmark trial and return its metrics plus metadata."""
    init_file = ensure_init(X, dataset_name, k, seed, init_dir)
    result = run_geokmeans(X, k, init_file, algorithm, tol, max_iter)
    result.update(
        dataset=dataset_name,
        implementation=implementation,
        algorithm=algorithm,
        k=k,
        seed=seed,
    )
    return result


def load_all_datasets(config: Dict[str, Any], data_dir: Path) -> Dict[str, np.ndarray]:
    """Load every configured dataset once, reporting any that are missing.

    Raises if none could be loaded: an empty result set silently produces empty
    CSVs and unreadable tables downstream, which is worse than stopping here.
    """
    loaded, missing = {}, []
    for dataset_config in config["datasets"]:
        name = dataset_config["name"]
        try:
            loaded[name] = load_dataset(dataset_config, data_dir)
            n, d = loaded[name].shape
            print(f"  loaded {name}: {n} rows x {d} features")
        except FileNotFoundError as exc:
            missing.append(name)
            print(f"  ⚠️  skipping {name}: {exc}".splitlines()[0])

    if not loaded:
        raise SystemExit(
            f"\nNo datasets could be loaded (missing: {', '.join(missing)}).\n"
            f"Place the CSV files in {data_dir} and re-run."
        )
    if missing:
        print(f"\n  ⚠️  {len(missing)} dataset(s) missing and excluded from "
              f"results: {', '.join(missing)}")
    return loaded


def run_e1_exactness(config: Dict[str, Any], base_dir: Path,
                     datasets: Dict[str, np.ndarray]) -> pd.DataFrame:
    """E1: verify every algorithm reproduces Lloyd's clustering.

    For each dataset x k x seed, run all geokmeans algorithms from identical
    initial centroids and compare each against the reference implementation.
    """
    print("=" * 80)
    print("E1: EXACTNESS VERIFICATION")
    print("=" * 80)

    controls = config["controls"]
    tol, max_iter = float(controls["tol"]), int(controls["max_iter"])
    k_values = config["exactness_check"]["k_values"]
    seeds = config["exactness_check"]["seeds"]
    reference = config["exactness_check"].get("reference", "geokmeans_lloyd")
    ref_algo = reference.replace("geokmeans_", "")

    algos = list(config["implementations"]["geokmeans"])
    if ref_algo not in algos:
        raise SystemExit(
            f"E1 reference '{reference}' is not in implementations.geokmeans "
            f"({algos}); add it to config.yaml."
        )
    compared = [a for a in algos if a != ref_algo]

    init_dir = base_dir / "inits"
    results: List[Dict[str, Any]] = []
    failures = 0

    # One reference run plus one run per compared algorithm, per (k, seed).
    per_cell = 1 + len(compared)
    total = len(datasets) * len(k_values) * len(seeds) * per_cell

    with tqdm(total=total, desc="E1 Progress") as pbar:
        for dataset_name, X in datasets.items():
            print(f"\nDataset: {dataset_name}")

            for k in k_values:
                for seed in seeds:
                    ref = run_single_trial(
                        X, dataset_name, reference, ref_algo,
                        k, seed, init_dir, tol, max_iter,
                    )
                    ref_sse, ref_labels = ref["sse"], ref["labels"]

                    # The reference is reported alongside the others so that
                    # every algorithm appears in the summary, trivially exact
                    # against itself.
                    ref["sse_ratio_vs_lloyd"] = 1.0
                    ref["max_relative_sse_deviation"] = 0.0
                    ref["ari_vs_lloyd"] = 1.0
                    results.append(ref)
                    pbar.update(1)

                    for algo in compared:
                        result = run_single_trial(
                            X, dataset_name, f"geokmeans_{algo}",
                            algo, k, seed, init_dir, tol, max_iter,
                        )

                        deviation = abs(result["sse"] - ref_sse) / (ref_sse + 1e-12)
                        ari = adjusted_rand_index(ref_labels, result["labels"])

                        result["sse_ratio_vs_lloyd"] = result["sse"] / (ref_sse + 1e-12)
                        result["max_relative_sse_deviation"] = deviation
                        result["ari_vs_lloyd"] = ari

                        results.append(result)
                        pbar.update(1)

                        if deviation > 1e-6 or ari < 0.9999:
                            failures += 1
                            pbar.write(
                                f"  ⚠️  EXACTNESS FAILURE: geokmeans_{algo} "
                                f"{dataset_name} k={k} seed={seed} "
                                f"SSE dev={deviation:.2e} ARI={ari:.6f}"
                            )

    print(f"\nE1 complete: {len(results)} comparisons, {failures} failure(s)")
    return pd.DataFrame(results).drop(columns=["labels", "centroids"])


def run_e2_comparison(config: Dict[str, Any], base_dir: Path,
                      datasets: Dict[str, np.ndarray]) -> pd.DataFrame:
    """E2: time every algorithm across datasets, k values and trials."""
    print("=" * 80)
    print("E2: COMPARISON AGAINST BASELINES")
    print("=" * 80)

    controls = config["controls"]
    tol, max_iter = float(controls["tol"]), int(controls["max_iter"])
    trials = int(controls["trials_per_config"])
    warmups = int(controls.get("warmup_trials", 0))
    k_values = config["k_values"]

    all_impls = [{"name": f"geokmeans_{a}", "algorithm": a}
                 for a in config["implementations"]["geokmeans"]]

    init_dir = base_dir / "inits"
    results: List[Dict[str, Any]] = []

    total = len(datasets) * len(k_values) * trials * len(all_impls)

    with tqdm(total=total, desc="E2 Progress") as pbar:
        for dataset_name, X in datasets.items():
            print(f"\nDataset: {dataset_name}")

            for k in k_values:
                for impl in all_impls:
                    # Untimed fits first, so the first timed trial is not the
                    # one paying for cold caches and lazy library loading.
                    for _ in range(warmups):
                        run_single_trial(
                            X, dataset_name, impl["name"], impl["algorithm"],
                            k, 0, init_dir, tol, max_iter,
                        )

                    for trial in range(trials):
                        result = run_single_trial(
                            X, dataset_name, impl["name"], impl["algorithm"],
                            k, trial, init_dir, tol, max_iter,
                        )
                        result["trial"] = trial
                        results.append(result)
                        pbar.update(1)

    df = pd.DataFrame(results).drop(columns=["labels", "centroids"])
    return add_sse_ratio(df, config)


def sse_reference(config: Dict[str, Any]) -> str:
    """Implementation that sse_ratio_vs_lloyd is measured against."""
    return config.get("controls", {}).get("sse_reference", "geokmeans_lloyd")


def add_sse_ratio(df: pd.DataFrame, config: Dict[str, Any]) -> pd.DataFrame:
    """Add sse_ratio_vs_lloyd, relative to the baseline run at the SAME seed.

    The denominator must be matched on seed, not aggregated across seeds. Each
    seed is a different initialisation and converges to a different local
    optimum, so dividing by a cross-seed median measures seed-to-seed spread
    rather than any difference between implementations -- every implementation
    in a given seed then shows the same misleading offset, and a genuine
    divergence hides inside it.

    Done as a single vectorised join rather than the previous rescan of the
    whole accumulated result list per (dataset, k), which was quadratic.
    """
    baseline = sse_reference(config)
    if baseline not in set(df["implementation"]):
        raise SystemExit(
            f"SSE reference '{baseline}' produced no results, so "
            f"sse_ratio_vs_lloyd cannot be computed. Set controls.sse_reference "
            f"to an implementation that this run actually includes."
        )

    reference = (
        df[df["implementation"] == baseline]
        .groupby(["dataset", "k", "seed"])["sse"]
        .median()  # one row per seed already; median is a no-op guard
        .rename("_ref_sse")
    )
    df = df.merge(reference, on=["dataset", "k", "seed"], how="left")
    df["sse_ratio_vs_lloyd"] = df["sse"] / df["_ref_sse"]
    return df.drop(columns=["_ref_sse"])


def apply_dry_run(config: Dict[str, Any]) -> Dict[str, Any]:
    """Shrink the config to the smallest run that still exercises every path."""
    config = copy.deepcopy(config)

    # Keep the k that table1_main asks for, so a dry run actually exercises
    # main-table generation instead of skipping it for want of that cell.
    table_ks = [
        o["k"] for o in config.get("outputs", [])
        if o.get("type") == "latex_table" and o.get("k") in config["k_values"]
    ]
    config["k_values"] = [table_ks[0]] if table_ks else config["k_values"][:1]
    config["controls"]["trials_per_config"] = 2
    config["controls"]["warmup_trials"] = 1
    config["exactness_check"]["k_values"] = config["exactness_check"]["k_values"][:1]
    config["exactness_check"]["seeds"] = config["exactness_check"]["seeds"][:2]
    return config


def main():
    parser = argparse.ArgumentParser(description="Run geokmeans benchmarks")
    parser.add_argument("--experiment", choices=["E1", "E2", "all"], default="all",
                        help="Which experiment to run")
    parser.add_argument("--config", type=Path,
                        default=Path(__file__).parent.parent / "config.yaml",
                        help="Path to config file")
    parser.add_argument("--dry-run", action="store_true",
                        help="Smallest run that exercises the whole pipeline")
    parser.add_argument("--subsample", type=int, default=None,
                        help="Use only the first N rows of each dataset")
    args = parser.parse_args()

    with open(args.config) as f:
        config = yaml.safe_load(f)
    if args.dry_run:
        config = apply_dry_run(config)
        print("DRY RUN: reduced k values, seeds and trials\n")

    base_dir = args.config.parent.resolve()
    results_dir = base_dir / "results" / "raw"
    results_dir.mkdir(parents=True, exist_ok=True)

    print(f"Thread limits: {report_thread_limits()}")
    print("\nLoading datasets...")
    datasets = load_all_datasets(config, base_dir / "data")
    if args.subsample:
        datasets = {n: X[: args.subsample] for n, X in datasets.items()}
        print(f"  subsampled to first {args.subsample} rows per dataset")

    suffix = "_dryrun" if args.dry_run else ""

    if args.experiment in ("E1", "all"):
        e1 = run_e1_exactness(config, base_dir, datasets)
        out = results_dir / f"e1_exactness{suffix}.csv"
        e1.to_csv(out, index=False)
        print(f"\n✓ E1 results saved to {out}")

    if args.experiment in ("E2", "all"):
        e2 = run_e2_comparison(config, base_dir, datasets)
        out = results_dir / f"e2_comparison{suffix}.csv"
        e2.to_csv(out, index=False)
        print(f"\n✓ E2 results saved to {out}")

    print("\n" + "=" * 80)
    print("BENCHMARK COMPLETE")
    print("=" * 80)
    print("\nNext steps:")
    print("  1. Review raw results in results/raw/")
    print(f"  2. Generate tables and reports: python report.py{' --dry-run' if args.dry_run else ''}")


if __name__ == "__main__":
    main()
