"""E3: Geometric-k-means vs scikit-learn Lloyd on twitter.

Only two implementations, and the comparison is deliberately like-for-like:

  * scikit-learn runs at its native tolerance.
  * geo runs at the HARMONISED threshold sqrt(mean_var * tol / k), which makes
    the geokmeans stopping test (`sqrt(S/k) <= thr`) algebraically identical to
    scikit-learn's (`S <= mean_var * tol`). Without this geo runs a far
    stricter test, does ~40% more iterations, and its wall clock is charged for
    work scikit-learn never does.

Recorded per trial: total runtime, total iterations, time per iteration
(distance computations are kept as an extra column -- free to record, and the
plot can be switched to them without re-running).

Usage:
    python scripts/bench_twitter.py --k 100 300 500 --trials 5
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _threads import pin_single_threaded, report_thread_limits  # noqa: E402

pin_single_threaded()

import argparse  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from init_generator import generate_and_save_init, load_init  # noqa: E402
from utils import Timer, compute_sse, load_dataset  # noqa: E402

TOL, MAX_ITER = 1e-4, 300

METRICS = ["wall_clock_seconds", "n_iterations", "ms_per_iteration",
           "n_distance_calculations"]


def prepare_inits(X, dataset, ks, trials, init_dir, suffix=""):
    """Generate every (k, trial) initialisation once, up front, and cache it.

    k-means++ is O(k n d) per seed, which at k=500 on twitter is substantial
    work in its own right. Doing it here rather than lazily inside the trial
    loop keeps it off the benchmark's critical path, makes its cost visible
    separately, and guarantees that a re-run reuses byte-identical centroids
    instead of regenerating them.

    Returns {(k, seed): path}.
    """
    paths, made = {}, 0
    init_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nPreparing initial centroids for {len(ks)} k values "
          f"x {trials} trials ...")
    with Timer() as t:
        for k in ks:
            for seed in range(trials):
                ip = init_dir / f"{dataset}{suffix}_k{k}_seed{seed}.csv"
                if not ip.exists():
                    generate_and_save_init(X, f"{dataset}{suffix}", k, seed,
                                           init_dir)
                    made += 1
                paths[(k, seed)] = ip
    print(f"  {made} generated, {len(paths) - made} reused from cache "
          f"({t.wall_clock_seconds:.1f}s)")
    return paths


def harmonised_tol(X, k, tol=TOL):
    """geokmeans threshold reproducing scikit-learn's stopping test."""
    mean_var = float(np.mean(np.var(X, axis=0)))
    return float(np.sqrt(mean_var * tol / k)), mean_var


def run_geo(X, k, init_path, thr):
    from geokmeans import GeoKMeans
    km = GeoKMeans(n_clusters=k, algorithm="geo", init=str(init_path),
                   tol=thr, max_iter=MAX_ITER)
    with Timer() as t:
        km.fit(X)
    return dict(implementation="geokmeans_geo", n_iterations=int(km.n_iter_),
                n_distance_calculations=int(km.n_distance_calculations_),
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds)


def run_sklearn(X, k, C0):
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=k, init=C0, n_init=1, algorithm="lloyd",
                tol=TOL, max_iter=MAX_ITER, random_state=0)
    with Timer() as t:
        km.fit(X)
    n_it = int(km.n_iter_)
    return dict(implementation="sklearn_lloyd", n_iterations=n_it,
                # analytic: Lloyd computes every n*k distance per pass, plus
                # the final assignment pass. scikit-learn has no counter.
                n_distance_calculations=(n_it + 1) * X.shape[0] * k,
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="twitter")
    ap.add_argument("--path", default="data/Twitter.csv")
    ap.add_argument("--k", type=int, nargs="+", default=[100, 300, 500])
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--warmup", type=int, default=0,
                    help="untimed fits per (k, impl); 0 by default because at "
                         "this scale each fit is many seconds, so cold-start "
                         "effects are negligible and warm-up would roughly "
                         "double a multi-hour sweep")
    ap.add_argument("--subsample", type=int, default=None,
                    help="use only the first N rows (smoke testing)")
    ap.add_argument("--dry-run", action="store_true",
                    help="tiny end-to-end run: 20k rows, k=[10,20], 2 trials. "
                         "Exercises every code path and the plotting chain in "
                         "under a minute. Writes to *_dryrun.csv.")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    if args.dry_run:
        args.subsample = args.subsample or 20000
        args.k = [10, 20]
        args.trials = 2
        print("DRY RUN: 20k rows, k=[10, 20], 2 trials\n")

    base = Path(__file__).parent.parent
    print(f"Thread limits: {report_thread_limits()}")
    X = load_dataset({"name": args.dataset, "path": args.path}, base / "data")
    if args.subsample:
        X = np.ascontiguousarray(X[:args.subsample])
        print(f"  subsampled to first {args.subsample} rows")
    print(f"{args.dataset}: {X.shape[0]} rows x {X.shape[1]} features")

    init_dir = base / "inits"
    suffix = f"_n{args.subsample}" if args.subsample else ""
    init_paths = prepare_inits(X, args.dataset, args.k, args.trials,
                               init_dir, suffix)
    rows = []

    for k in args.k:
        thr, mean_var = harmonised_tol(X, k)
        print(f"\n{'=' * 70}\nk = {k}   harmonised tol = {thr:.6e} "
              f"(mean_var = {mean_var:.6g})\n{'=' * 70}")

        for seed in range(args.trials):
            ip = init_paths[(k, seed)]
            C0 = load_init(ip)

            for _ in range(args.warmup):
                run_geo(X, k, ip, thr)
                run_sklearn(X, k, C0)

            for r in (run_geo(X, k, ip, thr), run_sklearn(X, k, C0)):
                r["ms_per_iteration"] = (1000 * r["wall_clock_seconds"]
                                         / max(r["n_iterations"], 1))
                r.update(dataset=args.dataset, k=k, seed=seed, trial=seed)
                rows.append(r)
                print(f"  seed {seed}  {r['implementation']:<15} "
                      f"{r['n_iterations']:>4} iters  "
                      f"{r['wall_clock_seconds']:>8.3f}s  "
                      f"{r['ms_per_iteration']:>8.2f} ms/iter")

    df = pd.DataFrame(rows)
    tag = "_dryrun" if args.dry_run else ""
    out = args.out or base / f"results/raw/e3_twitter{tag}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\n✓ saved per-trial results: {out}")

    g = df.groupby(["dataset", "implementation", "k"], as_index=False)
    summary = g[METRICS].mean().rename(columns={m: f"{m}_mean" for m in METRICS})
    std = g[METRICS].std(ddof=1).rename(columns={m: f"{m}_std" for m in METRICS})
    summary = summary.merge(std, on=["dataset", "implementation", "k"])
    for m in METRICS:
        summary[f"{m}_std"] = summary[f"{m}_std"].fillna(0.0)
    summary = summary.merge(g.size().rename(columns={"size": "n_trials"}),
                            on=["dataset", "implementation", "k"])

    sum_out = base / f"results/reports/e3_twitter_summary{tag}.csv"
    sum_out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(sum_out, index=False)
    print(f"✓ saved summary (mean/std): {sum_out}\n")

    print(summary[["implementation", "k", "wall_clock_seconds_mean",
                   "wall_clock_seconds_std", "n_iterations_mean",
                   "n_iterations_std", "ms_per_iteration_mean",
                   "ms_per_iteration_std"]].to_string(index=False))


if __name__ == "__main__":
    main()
