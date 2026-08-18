"""Pre-generate and cache shared k-means++ initial centroids.

One initialisation per (k, trial), written to
``inits/<dataset>_k<K>_seed<S>.csv`` at full float64 precision. Every
implementation in every experiment then reads those same files, so all
algorithms start from identical points and a re-run reuses byte-identical
centroids rather than regenerating them.

Worth running separately before a long benchmark: k-means++ is O(k n d) per
seed, which on twitter at k=500 is substantial work that would otherwise sit
on the benchmark's critical path. Already-present files are left untouched.

Usage:
    python scripts/make_inits.py --dataset twitter --path data/Twitter.csv \
        --k 100 300 500 --trials 5
    python scripts/make_inits.py --dataset sensor --k 50 80 100 --trials 5
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _threads import pin_single_threaded  # noqa: E402

pin_single_threaded()

import argparse  # noqa: E402
import numpy as np  # noqa: E402

from init_generator import generate_and_save_init, load_init  # noqa: E402
from utils import Timer, load_dataset  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--path", default=None,
                    help="CSV path; defaults to data/<dataset>.csv")
    ap.add_argument("--k", type=int, nargs="+", required=True)
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--subsample", type=int, default=None)
    ap.add_argument("--force", action="store_true",
                    help="regenerate even if the file already exists")
    args = ap.parse_args()

    base = Path(__file__).parent.parent
    path = args.path or f"data/{args.dataset}.csv"
    X = load_dataset({"name": args.dataset, "path": path}, base / "data")
    suffix = ""
    if args.subsample:
        X = np.ascontiguousarray(X[:args.subsample])
        suffix = f"_n{args.subsample}"
    print(f"{args.dataset}: {X.shape[0]} rows x {X.shape[1]} features")

    init_dir = base / "inits"
    init_dir.mkdir(parents=True, exist_ok=True)
    name = f"{args.dataset}{suffix}"

    made = reused = 0
    with Timer() as total:
        for k in args.k:
            for seed in range(args.trials):
                ip = init_dir / f"{name}_k{k}_seed{seed}.csv"
                if ip.exists() and not args.force:
                    reused += 1
                    continue
                with Timer() as t:
                    generate_and_save_init(X, name, k, seed, init_dir)
                made += 1
                print(f"  k={k:<4} seed={seed}  {t.wall_clock_seconds:7.2f}s"
                      f"  -> {ip.name}")

    print(f"\n{made} generated, {reused} already present "
          f"({total.wall_clock_seconds:.1f}s total)")

    # sanity: every requested file exists and has the right shape
    bad = []
    for k in args.k:
        for seed in range(args.trials):
            ip = init_dir / f"{name}_k{k}_seed{seed}.csv"
            c = load_init(ip)
            if c.shape != (k, X.shape[1]):
                bad.append((ip.name, c.shape, (k, X.shape[1])))
    if bad:
        for n_, got, want in bad:
            print(f"  ✗ {n_}: shape {got}, expected {want}")
        raise SystemExit("some initialisations are malformed")
    print(f"✓ verified {len(args.k) * args.trials} files in {init_dir}")


if __name__ == "__main__":
    main()
