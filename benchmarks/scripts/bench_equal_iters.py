"""Compare scikit-learn Lloyd against C++ geo at an EQUAL iteration count.

The problem
-----------
The two libraries stop on different rules. scikit-learn tests
``S <= mean_var * tol``; the geokmeans core tests ``sqrt(S/k) <= tol``, where S
is the summed squared centroid shift. On real data the geokmeans test is
thousands of times stricter, so it runs far more iterations and its wall-clock
total is charged for work scikit-learn never does. Comparing the totals is
therefore not a comparison of the algorithms.

The fix
-------
Run geo first at its normal tolerance and record N, the iterations it needed.
Then run scikit-learn with ``tol=0.0, max_iter=N``. ``tol=0`` removes the
tolerance exit, so the loop is pinned to exactly N iterations (scikit-learn's
`strict_convergence` label-stability exit can still fire earlier; that is a
genuine fixed point and is reported when it happens).

With the iteration count equalised, **time per iteration** is the comparison
that means something: one Lloyd iteration computes all n*k distances, one geo
iteration computes only what its pruning admits.

Usage:
    python scripts/bench_equal_iters.py --dataset twitter --k 50
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
from utils import Timer, compute_sse, adjusted_rand_index, load_dataset  # noqa: E402

TOL, MAX_ITER = 1e-4, 300


def run_geo(X, k, init_path):
    from geokmeans import GeoKMeans
    km = GeoKMeans(n_clusters=k, algorithm="geo", init=str(init_path),
                   tol=TOL, max_iter=MAX_ITER)
    with Timer() as t:
        km.fit(X)
    return dict(implementation="geokmeans_geo", n_iterations=int(km.n_iter_),
                n_distance_calculations=int(km.n_distance_calculations_),
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds, labels=km.labels_,
                stopping_rule=f"geokmeans tol={TOL}")


def harmonised_tol(X, k, tol=TOL):
    """geokmeans threshold reproducing scikit-learn's stopping test.

    scikit-learn stops when S <= mean_var * tol; the geokmeans core stops when
    sqrt(S/k) <= thr. Setting thr = sqrt(mean_var * tol / k) makes the two the
    same test on the same quantity.
    """
    mean_var = float(np.mean(np.var(X, axis=0)))
    return float(np.sqrt(mean_var * tol / k)), mean_var


def run_geo_tol(X, k, init_path, thr, tag):
    from geokmeans import GeoKMeans
    km = GeoKMeans(n_clusters=k, algorithm="geo", init=str(init_path),
                   tol=thr, max_iter=MAX_ITER)
    with Timer() as t:
        km.fit(X)
    return dict(implementation=tag, n_iterations=int(km.n_iter_),
                n_distance_calculations=int(km.n_distance_calculations_),
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds, labels=km.labels_,
                stopping_rule=f"geokmeans tol={thr:.6e} (harmonised)")


def run_sklearn(X, k, C0, tol, max_iter, tag):
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=k, init=C0, n_init=1, algorithm="lloyd",
                tol=tol, max_iter=max_iter, random_state=0)
    with Timer() as t:
        km.fit(X)
    n_it = int(km.n_iter_)
    return dict(implementation=tag, n_iterations=n_it,
                # Lloyd computes every n*k distance each pass, plus the final
                # assignment pass. Analytic, not measured -- sklearn has no counter.
                n_distance_calculations=(n_it + 1) * X.shape[0] * k,
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds, labels=km.labels_,
                stopping_rule=f"tol={tol}, max_iter={max_iter}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="twitter")
    ap.add_argument("--path", default=None, help="CSV path override")
    ap.add_argument("--k", type=int, nargs="+", default=[50])
    ap.add_argument("--no-warmup", action="store_true",
                    help="skip untimed warm-up fits (cold-start effects are "
                         "negligible for multi-second fits and warm-up doubles "
                         "the cost of a large sweep)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    base = Path(__file__).parent.parent
    print(f"Thread limits: {report_thread_limits()}")

    path = args.path or f"data/{args.dataset}.csv"
    X = load_dataset({"name": args.dataset, "path": path}, base / "data")
    print(f"{args.dataset}: {X.shape[0]} rows x {X.shape[1]} features")

    init_dir = base / "inits"
    rows = []

    for k in args.k:
        print(f"\n{'=' * 70}\nk = {k}\n{'=' * 70}")
        ip = init_dir / f"{args.dataset}_k{k}_seed{args.seed}.csv"
        if not ip.exists():
            print("generating shared k-means++ initialisation ...")
            generate_and_save_init(X, args.dataset, k, args.seed, init_dir)
        C0 = load_init(ip)
        thr, mean_var = harmonised_tol(X, k)

        if not args.no_warmup:
            print("warm-up (untimed) ...")
            run_geo(X, k, ip)
            run_sklearn(X, k, C0, TOL, MAX_ITER, "warmup")

        print("1/4  geo at its own tolerance ...")
        geo = run_geo(X, k, ip)
        n = geo["n_iterations"]
        print(f"     geo used {n} iterations")

        print(f"2/4  sklearn Lloyd pinned to {n} iterations (tol=0) ...")
        sk_eq = run_sklearn(X, k, C0, 0.0, n, "sklearn_lloyd_equal_iters")
        if sk_eq["n_iterations"] < n:
            print(f"     NOTE: sklearn stopped at {sk_eq['n_iterations']} of {n} "
                  f"via strict_convergence (labels stable)")

        print("3/4  sklearn Lloyd at its own tolerance ...")
        sk_nat = run_sklearn(X, k, C0, TOL, MAX_ITER, "sklearn_lloyd_native")

        print(f"4/4  geo at harmonised tolerance {thr:.4e} "
              f"(mean_var={mean_var:.6g}) ...")
        geo_h = run_geo_tol(X, k, ip, thr, "geokmeans_geo_harmonised")

        ref = geo["labels"]
        for r in (geo, sk_eq, sk_nat, geo_h):
            r["ari_vs_geo"] = adjusted_rand_index(ref, r["labels"])
            r["sse_ratio_vs_geo"] = r["sse"] / geo["sse"]
            r["ms_per_iteration"] = (1000 * r["wall_clock_seconds"]
                                     / max(r["n_iterations"], 1))
            r["ns_per_distance"] = (1e9 * r["wall_clock_seconds"]
                                    / r["n_distance_calculations"])
            r.pop("labels")
            r.update(dataset=args.dataset, k=k, seed=args.seed)
            rows.append(r)

        df_k = pd.DataFrame(rows[-4:])
        print()
        print(df_k[["implementation", "n_iterations", "wall_clock_seconds",
                    "ms_per_iteration"]].to_string(index=False))

    df = pd.DataFrame(rows)
    tag = "_".join(str(k) for k in args.k)
    out = args.out or base / f"results/raw/equal_iters_{args.dataset}_k{tag}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\n✓ saved {out}\n")

    cols = ["k", "implementation", "n_iterations", "wall_clock_seconds",
            "ms_per_iteration", "n_distance_calculations", "ns_per_distance",
            "sse_ratio_vs_geo", "ari_vs_geo"]
    print(df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
