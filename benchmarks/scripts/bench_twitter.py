"""E3/E4/E5 on twitter: Geometric-k-means vs scikit-learn Lloyd vs R k-means.

One sweep produces all three panels of Figure 1:

  E3  distance computations and runtime, geo vs scikit-learn  (the original)
  E4  the same, with R's default `stats::kmeans` (Hartigan-Wong) added
  E5  total RAPL energy for each of the three, per fit

The three implementations are run interleaved within a (k, seed) cell so none
of them is systematically favoured by drift in machine state, and all three
start from the SAME k-means++ centroids read off disk.

Tolerances are not the same object in the three libraries:

  * scikit-learn stops when  S <= mean_var * tol.
  * geokmeans stops when  sqrt(S/k) <= tol, so it runs at the HARMONISED
    threshold sqrt(mean_var * tol / k), which makes the two tests
    algebraically identical. Without this geo runs a far stricter test, does
    ~40% more iterations, and its wall clock is charged for work scikit-learn
    never does.
  * R's Hartigan-Wong has no tolerance at all -- it runs to a local optimum
    where no single point transfer reduces SSE. That is a strictly stronger
    stopping condition than either of the above and cannot be relaxed, which
    is worth stating alongside the figure rather than trying to hide.

Usage:
    python scripts/bench_twitter.py --k 100 300 500 --trials 5
    python scripts/bench_twitter.py --dry-run
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _threads import pin_single_threaded, report_thread_limits  # noqa: E402

pin_single_threaded()

import argparse  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import energy as energy_mod  # noqa: E402
from init_generator import generate_and_save_init, load_init  # noqa: E402
from utils import Timer, compute_sse, load_dataset  # noqa: E402

TOL, MAX_ITER = 1e-4, 300

IMPL_CHOICES = ("geo", "sklearn", "r")

METRICS = ["wall_clock_seconds", "n_iterations", "ms_per_iteration",
           "n_distance_calculations", "sse",
           *energy_mod.ENERGY_COLUMNS]


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


# --------------------------------------------------------------- runners ---

def run_geo(X, k, init_path, thr, meter):
    from geokmeans import GeoKMeans
    km = GeoKMeans(n_clusters=k, algorithm="geo", init=str(init_path),
                   tol=thr, max_iter=MAX_ITER)
    with meter, Timer() as t:
        km.fit(X)
    return dict(implementation="geokmeans_geo", n_iterations=int(km.n_iter_),
                n_distance_calculations=int(km.n_distance_calculations_),
                distance_count_method="measured",
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds,
                converged=int(km.n_iter_) < MAX_ITER,
                **meter.result())


def run_sklearn(X, k, C0, meter):
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=k, init=C0, n_init=1, algorithm="lloyd",
                tol=TOL, max_iter=MAX_ITER, random_state=0)
    with meter, Timer() as t:
        km.fit(X)
    n_it = int(km.n_iter_)
    return dict(implementation="sklearn_lloyd", n_iterations=n_it,
                # analytic: Lloyd computes every n*k distance per pass, plus
                # the final assignment pass. scikit-learn has no counter.
                n_distance_calculations=(n_it + 1) * X.shape[0] * k,
                distance_count_method="analytic",
                sse=compute_sse(X, km.labels_, km.cluster_centers_),
                wall_clock_seconds=t.wall_clock_seconds,
                converged=n_it < MAX_ITER,
                **meter.result())


class RRunner:
    """Drives `stats::kmeans` in a fresh Rscript process per (k, seed).

    The dataset is handed over once as a column-major float64 blob rather than
    re-parsed from CSV by R: it guarantees R clusters bit-identical data to the
    Python implementations, and turns a multi-hundred-megabyte `read.csv` into
    a `readBin`. A fresh process per cell keeps R's fits interleaved with the
    Python ones instead of batched into one long tail.
    """

    SCRIPT = Path(__file__).parent / "run_r_kmeans.R"

    def __init__(self, X, algorithm, workdir, max_iter=MAX_ITER):
        self.algorithm = algorithm
        self.max_iter = max_iter
        self.n, self.d = X.shape
        self.dir = Path(workdir)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.blob = self.dir / "X.f64"
        # tofile always writes C order, so transposing first puts the bytes on
        # the wire in the column-major layout R's matrix() expects.
        np.ascontiguousarray(X.T).tofile(self.blob)
        self.env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                        OPENBLAS_NUM_THREADS="1", R_ENABLE_JIT="3")

    @staticmethod
    def available():
        return shutil.which("Rscript") is not None

    def run(self, k, seed, init_path):
        out = self.dir / f"r_k{k}_seed{seed}.kv"
        cmd = ["Rscript", "--vanilla", str(self.SCRIPT),
               "--data", str(self.blob), "--n", str(self.n), "--d", str(self.d),
               "--init", str(init_path), "--k", str(k),
               "--max-iter", str(self.max_iter),
               "--algorithm", self.algorithm, "--out", str(out)]
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              env=self.env)
        if proc.returncode != 0 or not out.exists():
            raise RuntimeError(
                f"Rscript failed (exit {proc.returncode}) for k={k} "
                f"seed={seed}\n--- stdout ---\n{proc.stdout}"
                f"\n--- stderr ---\n{proc.stderr}")
        return self._parse(out)

    @staticmethod
    def _parse(path):
        row = {}
        for line in path.read_text().splitlines():
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            row[key] = val
        for key in ("wall_clock_seconds", "sse", "n_distance_calculations",
                    *energy_mod.ENERGY_COLUMNS):
            # R writes NA for an unmeasured energy domain; everything else is
            # a plain number. Anything unparseable becomes NaN rather than
            # taking down a multi-hour sweep.
            raw = row.get(key, "").strip()
            try:
                row[key] = float(raw)
            except ValueError:
                row[key] = float("nan")
        row["n_iterations"] = int(row["n_iterations"])
        row["converged"] = row.get("converged") == "True"
        return row


# ------------------------------------------------------------------ main ---

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="twitter")
    ap.add_argument("--path", default="data/Twitter.csv")
    ap.add_argument("--k", type=int, nargs="+", default=[100, 300, 500])
    ap.add_argument("--trials", type=int, default=5)
    ap.add_argument("--impls", nargs="+", default=list(IMPL_CHOICES),
                    choices=IMPL_CHOICES,
                    help="which implementations to run (default: all three)")
    ap.add_argument("--r-algorithm", default="Hartigan-Wong",
                    choices=["Hartigan-Wong", "Lloyd", "MacQueen", "Forgy"],
                    help="R stats::kmeans algorithm; Hartigan-Wong is R's "
                         "default and the one the paper compares against")
    ap.add_argument("--energy", dest="energy", action="store_true",
                    default=None, help="require RAPL; fail if unreadable")
    ap.add_argument("--no-energy", dest="energy", action="store_false",
                    help="skip energy measurement entirely")
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

    # --- energy availability, decided once and reported --------------------
    rapl_ok, rapl_why = energy_mod.availability()
    if args.energy and not rapl_ok:
        raise SystemExit(f"--energy requested but RAPL is unusable:\n  "
                         f"{rapl_why}\nRun: python scripts/energy.py --check")
    measure_energy = rapl_ok if args.energy is None else args.energy
    print(f"Energy (RAPL): {'ON' if measure_energy else 'OFF'} -- {rapl_why}")
    if not measure_energy:
        print("  energy columns will be NaN; panel (c) will be omitted")

    def new_meter():
        return energy_mod.EnergyMeter() if measure_energy \
            else energy_mod.EnergyMeter(domains=[])

    # --- data --------------------------------------------------------------
    X = load_dataset({"name": args.dataset, "path": args.path}, base / "data")
    if args.subsample:
        X = np.ascontiguousarray(X[:args.subsample])
        print(f"  subsampled to first {args.subsample} rows")
    print(f"{args.dataset}: {X.shape[0]} rows x {X.shape[1]} features")

    init_dir = base / "inits"
    suffix = f"_n{args.subsample}" if args.subsample else ""
    init_paths = prepare_inits(X, args.dataset, args.k, args.trials,
                               init_dir, suffix)

    # --- R --------------------------------------------------------------
    r_runner, r_tmp = None, None
    if "r" in args.impls:
        if not RRunner.available():
            raise SystemExit("Rscript not found on PATH. Install R, or drop "
                             "R from the sweep with --impls geo sklearn")
        r_tmp = tempfile.mkdtemp(prefix="rbench_",
                                 dir=str(base / "results"))
        with Timer() as t:
            r_runner = RRunner(X, args.r_algorithm, r_tmp)
        print(f"R handover blob: {X.nbytes / 1e6:.0f} MB written in "
              f"{t.wall_clock_seconds:.1f}s -> {r_tmp}")

    rows = []
    try:
        for k in args.k:
            thr, mean_var = harmonised_tol(X, k)
            print(f"\n{'=' * 70}\nk = {k}   harmonised tol = {thr:.6e} "
                  f"(mean_var = {mean_var:.6g})\n{'=' * 70}")

            for seed in range(args.trials):
                ip = init_paths[(k, seed)]
                C0 = load_init(ip)

                for _ in range(args.warmup):
                    if "geo" in args.impls:
                        run_geo(X, k, ip, thr, new_meter())
                    if "sklearn" in args.impls:
                        run_sklearn(X, k, C0, new_meter())

                trial_rows = []
                if "geo" in args.impls:
                    trial_rows.append(run_geo(X, k, ip, thr, new_meter()))
                if "sklearn" in args.impls:
                    trial_rows.append(run_sklearn(X, k, C0, new_meter()))
                if "r" in args.impls:
                    trial_rows.append(r_runner.run(k, seed, ip))

                for r in trial_rows:
                    r["ms_per_iteration"] = (1000 * r["wall_clock_seconds"]
                                             / max(r["n_iterations"], 1))
                    r.update(dataset=args.dataset, k=k, seed=seed, trial=seed)
                    rows.append(r)
                    e = r.get("energy_total_joules", float("nan"))
                    e_txt = f"{e:>9.1f} J" if e == e else "        - "
                    print(f"  seed {seed}  {r['implementation']:<22} "
                          f"{r['n_iterations']:>4} iters  "
                          f"{r['wall_clock_seconds']:>8.3f}s  "
                          f"{r['ms_per_iteration']:>8.2f} ms/iter  {e_txt}")
    finally:
        if r_tmp:
            shutil.rmtree(r_tmp, ignore_errors=True)

    df = pd.DataFrame(rows)
    tag = "_dryrun" if args.dry_run else ""
    out = args.out or base / f"results/raw/e3_twitter{tag}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\n✓ saved per-trial results: {out}")

    metrics = [m for m in METRICS if m in df.columns]
    g = df.groupby(["dataset", "implementation", "k"], as_index=False)
    summary = g[metrics].mean().rename(columns={m: f"{m}_mean" for m in metrics})
    std = g[metrics].std(ddof=1).rename(columns={m: f"{m}_std" for m in metrics})
    summary = summary.merge(std, on=["dataset", "implementation", "k"])
    for m in metrics:
        # a single trial has no std; an unmeasured metric stays NaN on purpose
        summary[f"{m}_std"] = summary[f"{m}_std"].where(
            summary[f"{m}_mean"].isna(), summary[f"{m}_std"].fillna(0.0))
    summary = summary.merge(g.size().rename(columns={"size": "n_trials"}),
                            on=["dataset", "implementation", "k"])

    sum_out = base / f"results/reports/e3_twitter_summary{tag}.csv"
    sum_out.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(sum_out, index=False)
    print(f"✓ saved summary (mean/std): {sum_out}\n")

    cols = ["implementation", "k", "n_distance_calculations_mean",
            "wall_clock_seconds_mean", "wall_clock_seconds_std",
            "n_iterations_mean", "energy_total_joules_mean",
            "energy_total_joules_std"]
    print(summary[[c for c in cols if c in summary.columns]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
