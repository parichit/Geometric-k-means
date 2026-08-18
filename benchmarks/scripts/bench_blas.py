"""Compare scikit-learn Lloyd against BLAS-vectorised geo on shared inits.

Also runs the C++ `geo` as a reference row, because without it you cannot tell
whether vectorising actually bought anything.

Usage:
    python scripts/bench_blas.py [--k 50 80 100] [--trials 5]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _threads import pin_single_threaded, report_thread_limits  # noqa: E402

pin_single_threaded()

import argparse  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from geo_blas import GeoKMeansBLAS  # noqa: E402
from init_generator import load_init  # noqa: E402
from utils import Timer, compute_sse, adjusted_rand_index  # noqa: E402

TOL, MAX_ITER = 1e-4, 300


def run_sklearn(X, k, C0):
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=k, init=C0, n_init=1, algorithm="lloyd",
                tol=TOL, max_iter=MAX_ITER, random_state=0)
    with Timer() as t:
        km.fit(X)
    return dict(implementation="sklearn_lloyd", n_iterations=int(km.n_iter_),
                n_distance_calculations=(int(km.n_iter_) + 1) * X.shape[0] * k,
                n_distances_materialised=(int(km.n_iter_) + 1) * X.shape[0] * k,
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds, labels=km.labels_)


def run_geo_blas(X, k, C0):
    km = GeoKMeansBLAS(k, init=C0, tol=TOL, max_iter=MAX_ITER)
    with Timer() as t:
        km.fit(X)
    return dict(implementation="geo_blas", n_iterations=km.n_iter_,
                n_distance_calculations=km.n_distance_calculations_,
                n_distances_materialised=km.n_distances_materialised_,
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds, labels=km.labels_)


def run_geo_cpp(X, k, init_path):
    from geokmeans import GeoKMeans
    km = GeoKMeans(n_clusters=k, algorithm="geo", init=str(init_path),
                   tol=TOL, max_iter=MAX_ITER)
    with Timer() as t:
        km.fit(X)
    return dict(implementation="geo_cpp", n_iterations=int(km.n_iter_),
                n_distance_calculations=int(km.n_distance_calculations_),
                n_distances_materialised=int(km.n_distance_calculations_),
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds, labels=km.labels_)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, nargs="+", default=[50, 80, 100])
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).parent.parent / "results/raw/blas_comparison.csv")
    args = ap.parse_args()

    base = Path(__file__).parent.parent
    print(f"Thread limits: {report_thread_limits()}")
    X = pd.read_csv(base / "data/sensor.csv", header=None).to_numpy(np.float64)
    X = np.ascontiguousarray(X)
    print(f"sensor: {X.shape[0]} rows x {X.shape[1]} features\n")

    rows = []
    for k in args.k:
        for seed in range(args.trials):
            ip = base / f"inits/sensor_k{k}_seed{seed}.csv"
            C0 = load_init(ip)

            # untimed warm-up per implementation
            run_sklearn(X, k, C0); run_geo_blas(X, k, C0); run_geo_cpp(X, k, ip)

            ref_labels = None
            for fn in (lambda: run_sklearn(X, k, C0),
                       lambda: run_geo_blas(X, k, C0),
                       lambda: run_geo_cpp(X, k, ip)):
                r = fn()
                if ref_labels is None:
                    ref_labels = r["labels"]
                    r["ari_vs_sklearn"] = 1.0
                else:
                    r["ari_vs_sklearn"] = adjusted_rand_index(ref_labels, r["labels"])
                r.pop("labels")
                r.update(k=k, seed=seed)
                rows.append(r)
            print(f"  k={k} seed={seed} done")

    df = pd.DataFrame(rows)
    ref = (df[df.implementation == "sklearn_lloyd"]
           .set_index(["k", "seed"])["sse"].rename("_ref"))
    df = df.merge(ref, left_on=["k", "seed"], right_index=True, how="left")
    df["sse_ratio_vs_sklearn"] = df["sse"] / df["_ref"]
    df = df.drop(columns=["_ref"])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\n✓ saved {args.out}")


if __name__ == "__main__":
    main()
