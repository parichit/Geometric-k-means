# Geokmeans Benchmarking Guide

Running the benchmark experiments for the JMLR MLOSS paper, including on a
remote Linux workstation as a background job.

---

## 0. What the library patches changed

The benchmarks depend on fixes made to the **canonical C++ sources in `src/`**,
which `geokmeans-py/sync_sources.py` vendors into the Python package. `src/` is
shared with the CRAN package, so these change R behaviour too.

| File | Fix |
|---|---|
| `algo_utils.hpp` | `check_convergence` summed the squared difference of only the **last feature dimension** (missing braces on the inner loop). Lloyd's, Geo and Ball converged far too early. |
| `ball_kmeans++_xf.hpp` | Same missing-braces bug in Ball's own `check_convergence`. |
| `lloyd_kmeans.hpp` | Returned the **previous** iteration's centroids whenever it converged (off by one). |
| `exponion.hpp`, `geokmeans.hpp` | Returned the **zeroed accumulator** whenever the loop hit `max_iter` — all-zero centroids, silently. |
| `_core.cpp` | `init` accepts an ndarray *or a file path*; array route writes a real temp file at 17 significant digits (was 6, silently truncating callers' centroids) and cleans up on every exit path (was leaking on exceptions). Previously derived its temp directory from `__FILE__`, the build machine's source path. |
| `estimator.py` | `init` accepts `'random'`, an ndarray, or a path. |

The last three fixes are what make E1 pass. Rebuild after any change to `src/`:

```bash
cd geokmeans-py && python sync_sources.py && pip install .
```

---

## 1. Setup (Linux workstation)

```bash
git clone <repo> && cd Geometric-k-means

python3 -m venv benchmarks/.venv
benchmarks/.venv/bin/pip install -r benchmarks/requirements.txt

# build and install the patched geokmeans into the same environment
cd geokmeans-py && python sync_sources.py
../benchmarks/.venv/bin/pip install .
cd ..
```

Requires a C++ compiler and CMake (scikit-build-core drives the build).

> **Run from `benchmarks/`, not the repo root.** The repo root contains a
> directory literally named `geokmeans/` (the C++ sources) which Python picks up
> as an empty namespace package and which shadows the real install. `run.py`
> detects this and says so rather than failing obscurely.

### Datasets

Place CSVs in `benchmarks/data/`. Headers are detected and skipped; non-numeric
columns are dropped with a note.

```
benchmarks/data/
├── sensor.csv      #  58,509 x 48
└── Twitter.csv     # 583,250 x 78
```

---

## 2. Verify before committing compute

Two fast checks. Run both; each takes seconds.

```bash
cd benchmarks

# 1. explicit initialisation is honoured, and all algorithms agree
.venv/bin/python scripts/test_init.py          # exits non-zero on failure

# 2. whole E1/E2 pipeline end to end, reduced k/seeds/trials
.venv/bin/python scripts/run.py --dry-run --config config_sensor_highk.yaml
.venv/bin/python scripts/report.py --dry-run --config config_sensor_highk.yaml

# 3. whole E3 pipeline: 20k rows, k=[10,20], 2 trials, ~4 seconds
.venv/bin/python scripts/bench_twitter.py --dry-run
.venv/bin/python scripts/plot_results.py \
    --twitter results/reports/e3_twitter_summary_dryrun.csv \
    --out results/figures/benchmark_summary_dryrun.png
```

Dry-run artefacts are suffixed `_dryrun` and never collide with real results.

---

## 3. The experiments

| | What | Algorithms | Dataset |
|---|---|---|---|
| **E1** | Exactness: does every algorithm reproduce Lloyd's clustering from identical initial centroids? Reports **mean ± std** of iterations, runtime, SSE ratio. | Lloyd, Hamerly, Annulus, Exponion, Geo | sensor |
| **E2** | Distance computations, iterations, runtime, solution quality across the same five. | same | sensor |
| **E3** | Geometric-k-means vs **scikit-learn Lloyd**, runtime and iterations vs `k`. | Geo, sklearn Lloyd | twitter |

`k = 50, 80, 100` for E1/E2; `k = 100, 300, 500` for E3. Five trials each.

**Ball k-means and Elkan are excluded** from E1/E2 by config. scikit-learn is
absent from E1/E2 entirely — see *Why no scikit-learn baseline in E1/E2* below.

---

## 4. Running the full benchmarks in the background

### Pre-stage the initial centroids

k-means++ is O(*knd*) per seed — on twitter at k=500 that is minutes of work per
trial. Generate all of them once, up front, so the cost is off the benchmark's
critical path and visible on its own:

```bash
cd benchmarks
mkdir -p logs

.venv/bin/python scripts/make_inits.py \
    --dataset twitter --path data/Twitter.csv --k 100 300 500 --trials 5
.venv/bin/python scripts/make_inits.py \
    --dataset sensor --k 50 80 100 --trials 5
```

One file per (k, trial) in `inits/`, at full float64 precision. Existing files
are reused untouched unless `--force`, so re-runs use byte-identical centroids.
The script verifies every file's shape before exiting. `bench_twitter.py` also
does this pass internally at startup, so this step is optional — but running it
separately keeps the timed section clean and lets you see the cost.

### Launch

```bash
cd benchmarks

# E1 + E2 on sensor
nohup .venv/bin/python scripts/run.py \
      --config config_sensor_highk.yaml --experiment all \
      > logs/e1e2.log 2>&1 &

# E3 on twitter
nohup .venv/bin/python scripts/bench_twitter.py \
      --k 100 300 500 --trials 5 \
      > logs/e3.log 2>&1 &
```

`mkdir -p logs` first. Or use `tmux new -s bench` and run in the foreground.

> **Do not pipe the output through `grep`/`tail`.** A run that was piped
> measured **33% CPU** and produced timings inflated ~40% uniformly, against
> 99% CPU for the same workload writing straight to a file — `tqdm` redraws
> constantly and the pipeline blocks the writer. Redirect to a file instead.

### Validate the timings afterwards

Timing results are only meaningful if the process actually had a core:

```bash
/usr/bin/time -v .venv/bin/python scripts/run.py ...   # Linux
# check "Percent of CPU this job got" -- want >95%
```

or while running:

```bash
ps -o pid=,etime=,time=,%cpu= -p $(pgrep -f scripts/run.py | tail -1)
```

If `%cpu` is well below 95, something else is competing and the runtime column
is not trustworthy. Distance counts, iteration counts and SSE are deterministic
and unaffected.

### Expected cost

E3 dominates: scikit-learn's per-iteration work is O(*nkd*), so at k=500 on
twitter each Lloyd pass is ~22.7 GFLOP. Budget hours, not minutes. `--warmup`
defaults to `0` for E3 for this reason — at this scale each fit is many seconds,
so cold-start effects are negligible and warm-up would roughly double the sweep.

---

## 5. Reports and the figure

```bash
.venv/bin/python scripts/report.py --config config_sensor_highk.yaml
.venv/bin/python scripts/plot_results.py
```

`report.py` writes LaTeX tables, a markdown exactness report, and **tidy
mean/std CSVs** which are the plotting input. `plot_results.py` reads only
those CSVs and recomputes nothing.

```
benchmarks/results/
├── raw/
│   ├── e1_exactness.csv          # per-trial
│   ├── e2_comparison.csv         # per-trial
│   ├── e2_summary.csv            # mean/std
│   └── e3_twitter.csv            # per-trial
├── reports/
│   ├── sensor_highk_exactness.md
│   ├── sensor_highk_exactness.csv    # mean/std  <- plot input (a-c)
│   └── e3_twitter_summary.csv        # mean/std  <- plot input (d-e)
├── tables/                        # LaTeX, needs only booktabs
└── figures/
    └── benchmark_summary.{png,pdf}
```

### The figure

Five panels — three on the top row, two centred beneath:

- **(a)** table of SSE at convergence per algorithm and k (sensor)
- **(b)** mean distance computations vs k, log scale (sensor)
- **(c)** mean runtime vs k (sensor)
- **(d)** mean runtime vs k — Geo vs scikit-learn Lloyd (twitter)
- **(e)** mean distance computations vs k, log scale — Geo vs scikit-learn
  Lloyd (twitter). scikit-learn's count is analytic (`n·k` per pass, it has no
  counter); geokmeans' is measured. The panel is annotated to say so.

All lines carry std error bars. Panels (d-e) are omitted automatically if the
E3 summary CSV is absent, so the figure still builds from E1/E2 alone.

---

## 6. Methodology

### Shared initialisation

k-means++ centroids are generated once per (dataset, k, seed), saved to
`inits/<dataset>_k<K>_seed<S>.csv` at full float64 precision, and every
implementation reads that same file. geokmeans receives the **path**, which its
C++ core reads directly — passing an array instead makes it serialise the
centroids to a temp CSV *inside* `fit()`, putting a file write in the timed
region.

### Harmonised tolerance (E3)

The two libraries stop on different rules:

| | stops when |
|---|---|
| scikit-learn | `S <= mean_var * tol` |
| geokmeans | `sqrt(S/k) <= tol` |

where `S` is the summed squared centroid shift. On real data the geokmeans test
is **thousands of times stricter** (6005x at k=50 on sensor), so it runs far more
iterations and its wall clock is charged for work scikit-learn never does.

E3 therefore runs geo at

```
thr = sqrt(mean_var * tol / k)
```

which makes the two tests algebraically identical. This is why E3 is a fair
runtime comparison and a naive one is not. Measured on twitter at k=50, this
took geo from 264 iterations / 20.7 s to 200 / 16.2 s against scikit-learn's
192 / 19.8 s — i.e. from a 5% loss to a 22% win.

### Controls

- **Single-threaded.** `scripts/_threads.py` sets the OpenMP/BLAS thread
  variables and re-execs the interpreter if numpy is already imported — those
  libraries read the variables at load time, so setting them after
  `import numpy` has no effect. On Linux `threadpoolctl` can verify the limit
  took; on macOS with Accelerate there is no controllable pool to inspect, and
  the run log says so.
- **Trials.** Five per cell. E1 reports and plot CSVs use mean ± std; the LaTeX
  tables use the median.
- **SSE ratio** is measured against `controls.sse_reference` **at the same
  seed**. Aggregating the denominator across seeds measures seed-to-seed spread
  rather than any difference between implementations.
- **Memory is not measured.** RSS sampled after `fit()` returns is neither a
  peak nor a delta.

### Why no scikit-learn baseline in E1/E2

The geokmeans core computes in `float`, scikit-learn in `float64`, so the two
agree only to ~2.4e-08 relative regardless of tolerance. That makes scikit-learn
unusable as an *exactness* reference. It remains a valid *performance*
comparator, which is what E3 is for.

---

## 7. Known issues

- **Tie-breaking.** Hamerly/Annulus/Exponion/Geo occasionally differ from
  Lloyd's at ~1e-5 relative SSE, always *together* and with identical ARI —
  points equidistant from two centroids resolve differently under pruning than
  under an exact argmin. 56/60 comparisons are exact on sensor; the misses fall
  in a single (k, seed) cell.
- **Ball k-means** was excluded after scoring 5/45 exact with ARI as low as
  0.983 and its own divergent iteration counts. Its separate `check_convergence`
  is the suspect. Not yet investigated.
- **`geo_blas.py`** is a numpy reimplementation of Geo, validated bit-exact
  against the C++ (identical iterations, ARI 1.000000, SSE ratio 1.000000000).
  Its BLAS kernels hit 2.9 ns/distance versus ~113 ns for the scalar C++, but
  Python per-cluster loop overhead is ~65% of its runtime, so it is *slower*
  overall. It is useful as an executable specification of Geo's semantics; a
  C++/Eigen version would be needed to test the vectorisation hypothesis
  properly.

---

## 8. Standalone comparison scripts

Not part of E1/E2/E3; each writes its own CSV.

```bash
.venv/bin/python scripts/bench_elkan.py --k 50 80 100    # sklearn Elkan vs Geo
.venv/bin/python scripts/bench_equal_iters.py \
    --dataset twitter --path data/Twitter.csv --k 50     # equal-iteration study
.venv/bin/python scripts/bench_blas.py                   # sklearn Lloyd vs geo_blas
```

`bench_equal_iters.py` pins scikit-learn to geo's iteration count via
`tol=0, max_iter=N` and reports ms/iteration and ns/distance for both, plus geo
at the harmonised tolerance. It flags when scikit-learn exits early via
`strict_convergence` rather than reaching N.
