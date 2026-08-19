# Geokmeans Benchmark Pipeline

Benchmark suite for the JMLR MLOSS paper on `geokmeans`. Produces empirical
comparisons between geokmeans algorithms and scikit-learn's Lloyd implementation.

## Overview

Five algorithms are benchmarked: **Lloyd's, Hamerly, Annulus, Exponion and
Geometric-k-means**. There is no external baseline -- see *No external baseline*
below.

- **E1 (Exactness Verification):** confirms every algorithm produces the same
  clustering as Lloyd's, given identical initial centroids. Reports
  **mean ± std** of iterations, runtime and SSE ratio.
- **E2 (Comparison):** measures distance computations, iterations, wall-clock
  time and solution quality across the five algorithms.
- **E3/E4/E5 (twitter, Figure 1):** Geometric-k-means vs scikit-learn Lloyd vs
  R's default `stats::kmeans` (Hartigan-Wong) — distance computations, runtime
  and **total RAPL energy** per fit. This is the figure the paper reports; see
  [BENCHMARKING_GUIDE.md](BENCHMARKING_GUIDE.md).

## Quick Start

### 1. Setup

```bash
cd benchmarks

python3 -m venv .venv
.venv/bin/pip install -r requirements.txt

# Build and install geokmeans into the same environment
.venv/bin/pip install ../geokmeans-py
```

The pipeline must run in an environment where `geokmeans` is installed. Note
that the repository root contains a directory literally named `geokmeans/` (the
C++ sources), which Python picks up as an empty namespace package and which
shadows the real install if you run from the repo root. Run from `benchmarks/`.
`run.py` detects this case and says so rather than failing obscurely.

### 2. Prepare Datasets

Place datasets in `data/` as CSV files (rows = samples, columns = features):

```
benchmarks/data/
├── sensor.csv
└── twitter.csv
```

Headers are detected and skipped automatically, and non-numeric columns (labels,
ids) are dropped with a note. A dataset listed in `config.yaml` but missing from
`data/` is skipped with a warning; if *none* load, the run stops rather than
producing empty results.

### 3. Verify initialization works

```bash
.venv/bin/python scripts/test_init.py
```

This gates the benchmark: it checks that explicit initialization is honoured and
that the array and file routes agree. Exits non-zero on failure.

### 4. Run Benchmarks

```bash
# Smallest run that exercises the whole pipeline end to end
.venv/bin/python scripts/run.py --dry-run

# Full run
.venv/bin/python scripts/run.py --experiment all

# Or individually
.venv/bin/python scripts/run.py --experiment E1
.venv/bin/python scripts/run.py --experiment E2
```

Useful flags: `--dry-run` (one k, two seeds, two trials, writes to
`*_dryrun.csv`) and `--subsample N` (first N rows of each dataset).

**Note:** experiments are single-threaded and may take hours on large datasets.

### 5. Generate Reports

```bash
.venv/bin/python scripts/report.py --config config_sensor_highk.yaml
```

Reads `results/raw/` and writes LaTeX tables to `results/tables/`, markdown
reports to `results/reports/`, and **tidy mean/std CSVs** next to them
(`results/reports/*_exactness.csv`, `results/raw/e2_summary.csv`). Those CSVs
are the plotting input -- one row per (dataset, algorithm, k) with
`<metric>_mean` / `<metric>_std` columns.

### 6. Plot

```bash
.venv/bin/python scripts/plot_results.py
```

`plot_results.py` writes Figure 1 to
`results/figures/benchmark_summary.{png,pdf}` -- three panels on twitter, all
with `k` on the x-axis, comparing Geometric-k-means, scikit-learn Lloyd and R
`stats::kmeans`:

- **(a)** mean distance computations vs `k`, log scale
- **(b)** mean runtime vs `k`
- **(c)** mean total energy vs `k`, grouped bars

Lines and bars carry std error bars. Panel (c) falls back to a placeholder if
the run carries no energy data. The E1/E2 sensor panels now live in
`scripts/plot_e1e2.py`. Nothing is recomputed -- the scripts only read the
summary CSVs.

## Configuration

Edit `config.yaml` to change datasets, `k` values, algorithms, controls, and
outputs. `outputs[].columns` is the single source of truth for which metrics
reach the tables.

## Methodology

### Shared Initialization

All implementations in a given (dataset, k, seed) trial start from **identical**
initial centroids:

1. k-means++ centroids are generated once per (dataset, k, seed)
2. saved to `inits/<dataset>_k<K>_seed<S>.csv` at full float64 precision
3. every implementation reads that same file

scikit-learn receives the loaded array; geokmeans receives the **path**, which
its C++ core reads directly. That asymmetry is deliberate: passing an array to
geokmeans makes it serialise the centroids to a temporary CSV inside `fit()`,
which would put a file write inside the timed region that scikit-learn never
pays.

### Controls

- **Single-threaded.** `scripts/_threads.py` sets the OpenMP/BLAS thread
  variables and re-execs the interpreter if numpy is already loaded. This is
  required: those libraries read the variables at load time, so setting them
  after `import numpy` — as an earlier version did — had no effect at all. The
  run log prints the thread limits actually in force.
- **Warm-up.** `controls.warmup_trials` untimed fits precede the timed trials,
  so the first trial is not the one paying for cold caches.
- **Trials.** `controls.trials_per_config` timed trials per cell. The E1 report
  and the plotting CSVs use **mean ± std**; the LaTeX tables use the median.

### Metrics

| Metric | Source |
|---|---|
| `wall_clock_seconds` | `time.perf_counter()` around `fit()` |
| `n_distance_calculations` | measured for geokmeans; **estimated** for scikit-learn |
| `n_iterations` | reported by each library |
| `sse` | recomputed identically for every implementation from returned labels and centroids |
| `sse_ratio_vs_lloyd` | `sse` / `sse` of `controls.sse_reference` **at the same seed** |
| `ari_vs_lloyd` | E1 only, against `exactness_check.reference` |

The SSE-ratio denominator is matched on **seed**, not aggregated across seeds.
Each seed is a different initialisation converging to a different local optimum,
so a cross-seed denominator measures seed spread rather than any difference
between implementations, and hides genuine divergence inside a common offset.

Two further notes:

- **No external baseline.** scikit-learn has been removed from this pipeline
  entirely. It computes in `float64` while the geokmeans core computes in
  `float`, so the two agree only to ~2.4e-08 relative, and its stopping rule
  differs (see *Comparing against scikit-learn*), which makes it unusable as an
  exactness reference. Standalone comparison scripts remain available:
  `scripts/bench_elkan.py` (scikit-learn Elkan vs Geometric) and
  `scripts/bench_blas.py` (scikit-learn Lloyd vs the BLAS prototype).
- **Memory is not measured.** RSS sampled after `fit()` returns is neither a
  peak nor a delta, so rather than print a meaningless column, `peak_rss_mb` was
  removed from the pipeline entirely.

### Comparing against scikit-learn

If you re-add a scikit-learn baseline, two things make the comparison inexact:

1. **Different stopping rules.** scikit-learn scales `tol` by the mean feature
   variance and tests the squared centre shift; the geokmeans core compares
   `sqrt(S/k)` against `tol` directly, where `S` is the summed squared centroid
   shift. The two coincide when the geokmeans threshold is set to
   `sqrt(mean_var * tol / k)`; without that conversion geokmeans runs a
   substantially stricter test and therefore more iterations.
2. **Different precision.** The geokmeans core computes in `float`,
   scikit-learn in `float64`, so the two agree only to ~2.4e-08 relative
   however the tolerances are set.

## Output Files

```
benchmarks/results/
├── raw/
│   ├── e1_exactness.csv
│   └── e2_comparison.csv
├── tables/
│   ├── table1_main.tex        # needs only booktabs
│   └── appendix_matrix.tex
└── reports/
    └── exactness_report.md
```

## Directory Structure

```
benchmarks/
├── config.yaml              # Master configuration
├── requirements.txt
├── scripts/
│   ├── _threads.py          # Thread pinning; must be imported before numpy
│   ├── init_generator.py    # k-means++ initialization
│   ├── utils.py             # Timing, SSE, dataset loading
│   ├── run.py               # Experiment runner (E1, E2)
│   ├── report.py            # Tables, reports, tidy mean/std CSVs
│   ├── plot_results.py      # Figure 1: twitter, 3 panels
│   ├── plot_e1e2.py         # E1/E2 sensor figure (3 panels)
│   ├── make_inits.py        # Pre-stage shared k-means++ centroids
│   ├── bench_twitter.py     # E3/E4/E5: geo vs scikit-learn vs R, + energy
│   ├── run_r_kmeans.R       # One timed stats::kmeans fit, driven by the above
│   ├── energy.py            # RAPL meter (Linux powercap); --check to verify
│   ├── geo_blas.py          # BLAS-vectorised geo (validated reference impl)
│   ├── bench_elkan.py       # Standalone: scikit-learn Elkan vs Geometric
│   ├── bench_blas.py        # Standalone: scikit-learn Lloyd vs geo_blas
│   └── test_init.py         # Pre-flight check on explicit initialization
├── data/                    # Place datasets here (CSV)
├── inits/                   # Generated initial centroids (auto-created)
└── results/                 # Benchmark outputs (auto-created)
    ├── raw/                 # per-trial CSVs + e2_summary.csv
    ├── tables/              # LaTeX
    ├── reports/             # markdown + tidy mean/std CSVs
    └── figures/             # benchmark_summary.{png,pdf}, e1e2_sensor.{png,pdf}
```

Documentation: [BENCHMARKING_GUIDE.md](BENCHMARKING_GUIDE.md) covers Figure 1
end to end; [docs/E1_E2_SENSOR.md](docs/E1_E2_SENSOR.md) covers the older
sensor experiments.

## Troubleshooting

**"resolved to a namespace package"** — a directory named `geokmeans` is
shadowing the install. Run from `benchmarks/`, not the repo root.

**"Dataset file not found"** — put the CSV in `benchmarks/data/` and check the
path in `config.yaml`.

**"No data for dataset=... k=..."** — `outputs.table1_main` names a
dataset/`k` that the run did not cover. The message lists what is available.

**Experiments are slow** — reduce `k_values`, use `--subsample N`, or lower
`controls.trials_per_config` (at a cost to statistical validity).

## Extending the Pipeline

**Add a dataset:** drop the CSV in `data/`, add `{name, path}` to
`config.yaml`.

**Add an output:** add an entry under `outputs:` and re-run `report.py` — no
need to re-run the experiments.

**Add a metric:** compute it in the runner in `scripts/run.py`, add it to
`AGGREGATED` in `scripts/report.py`, and list it in an output's `columns`.

## Citation

```bibtex
@article{sharma2025geometrickmeans,
  title={Geometric-k-means: A Bound Free Approach to Fast and Eco-Friendly k-means},
  author={Sharma, Parichit and Stanislaw, Marcin and Kurban, Hasan and Kulekci, Oguzhan and Dalkilic, Mehmet},
  year={2025},
  eprint={2508.06353},
  archivePrefix={arXiv},
  primaryClass={cs.LG}
}
```

## Contact

- Open an issue on the GitHub repository
- Contact: Parichit Sharma (parishar@iu.edu)
