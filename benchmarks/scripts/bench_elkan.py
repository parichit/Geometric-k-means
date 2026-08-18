"""scikit-learn Elkan vs the C++ geo, on shared initial centroids.

scikit-learn exposes no distance counter, and unlike Lloyd's there is no
analytic count for Elkan (it prunes), so no distance column is reported here --
there would be no honest number for the Elkan side.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _threads import pin_single_threaded, report_thread_limits  # noqa: E402

pin_single_threaded()

import argparse  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from init_generator import load_init  # noqa: E402
from utils import Timer, compute_sse, adjusted_rand_index  # noqa: E402

TOL, MAX_ITER = 1e-4, 300


def run_sklearn_elkan(X, k, C0):
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=k, init=C0, n_init=1, algorithm="elkan",
                tol=TOL, max_iter=MAX_ITER, random_state=0)
    with Timer() as t:
        km.fit(X)
    return dict(implementation="sklearn_elkan", n_iterations=int(km.n_iter_),
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds, labels=km.labels_)


def run_geo_cpp(X, k, init_path):
    from geokmeans import GeoKMeans
    km = GeoKMeans(n_clusters=k, algorithm="geo", init=str(init_path),
                   tol=TOL, max_iter=MAX_ITER)
    with Timer() as t:
        km.fit(X)
    return dict(implementation="geokmeans_geo", n_iterations=int(km.n_iter_),
                n_distance_calculations=int(km.n_distance_calculations_),
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds, labels=km.labels_)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, nargs="+", default=[50, 80, 100])
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).parent.parent / "results/raw/elkan_comparison.csv")
    args = ap.parse_args()

    base = Path(__file__).parent.parent
    print(f"Thread limits: {report_thread_limits()}")
    X = np.ascontiguousarray(
        pd.read_csv(base / "data/sensor.csv", header=None).to_numpy(np.float64))
    print(f"sensor: {X.shape[0]} rows x {X.shape[1]} features\n")

    rows = []
    for k in args.k:
        for seed in range(args.trials):
            ip = base / f"inits/sensor_k{k}_seed{seed}.csv"
            C0 = load_init(ip)

            run_sklearn_elkan(X, k, C0)          # untimed warm-up
            run_geo_cpp(X, k, ip)

            a = run_sklearn_elkan(X, k, C0)
            b = run_geo_cpp(X, k, ip)
            b["ari_vs_elkan"] = adjusted_rand_index(a["labels"], b["labels"])
            a["ari_vs_elkan"] = 1.0
            for r in (a, b):
                r.pop("labels")
                r.update(k=k, seed=seed)
                rows.append(r)
            print(f"  k={k} seed={seed} done")

    df = pd.DataFrame(rows)
    ref = (df[df.implementation == "sklearn_elkan"]
           .set_index(["k", "seed"])["sse"].rename("_ref"))
    df = df.merge(ref, left_on=["k", "seed"], right_index=True, how="left")
    df["sse_ratio_vs_elkan"] = df["sse"] / df["_ref"]
    df = df.drop(columns=["_ref"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"\n✓ saved {args.out}")


if __name__ == "__main__":
    main()
