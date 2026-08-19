# E1 / E2 — exactness and algorithm comparison on `sensor`

These are the original two experiments. They are **not** part of the paper's
Figure 1, which reports twitter only (see `../BENCHMARKING_GUIDE.md`). This
file keeps them reproducible and records the methodology decisions behind
them.

| | What | Algorithms | Dataset |
|---|---|---|---|
| **E1** | Exactness: does every algorithm reproduce Lloyd's clustering from identical initial centroids? Reports **mean ± std** of iterations, runtime, SSE ratio. | Lloyd, Hamerly, Annulus, Exponion, Geo | sensor |
| **E2** | Distance computations, iterations, runtime, solution quality across the same five. | same | sensor |

`k = 50, 80, 100`, five trials each. **Ball k-means and Elkan are excluded** by
config. scikit-learn is absent entirely — see *Why no scikit-learn baseline*
below.

---

## Running

```bash
cd benchmarks

# verification first: explicit initialisation honoured, all algorithms agree
.venv/bin/python scripts/test_init.py            # exits non-zero on failure
.venv/bin/python scripts/run.py    --dry-run --config config_sensor_highk.yaml
.venv/bin/python scripts/report.py --dry-run --config config_sensor_highk.yaml

# pre-stage the shared k-means++ centroids
.venv/bin/python scripts/make_inits.py --dataset sensor --k 50 80 100 --trials 5

# the sweep
mkdir -p logs
nohup .venv/bin/python scripts/run.py \
      --config config_sensor_highk.yaml --experiment all \
      > logs/e1e2.log 2>&1 &

# reports and the supporting figure
.venv/bin/python scripts/report.py    --config config_sensor_highk.yaml
.venv/bin/python scripts/plot_e1e2.py
```

`report.py` writes LaTeX tables, a markdown exactness report, and tidy
mean/std CSVs. `plot_e1e2.py` reads only those CSVs and recomputes nothing.
It writes `results/figures/e1e2_sensor.{png,pdf}` with three panels:

- **(a)** table of SSE at convergence per algorithm and `k`
- **(b)** mean distance computations vs `k`, log scale
- **(c)** mean runtime vs `k`

```
results/
├── raw/
│   ├── e1_exactness.csv          # per-trial
│   ├── e2_comparison.csv         # per-trial
│   └── e2_summary.csv            # mean/std
├── reports/
│   ├── sensor_highk_exactness.md
│   └── sensor_highk_exactness.csv    # mean/std  <- plot input
├── tables/                        # LaTeX, needs only booktabs
└── figures/e1e2_sensor.{png,pdf}
```

---

## Library fixes these experiments depend on

E1 only passes because of fixes to the canonical C++ sources in `src/`, which
`geokmeans-py/sync_sources.py` vendors into the Python package. `src/` is
shared with the CRAN package, so these changed R behaviour too.

| File | Fix |
|---|---|
| `algo_utils.hpp` | `check_convergence` summed the squared difference of only the **last feature dimension** (missing braces on the inner loop). Lloyd's, Geo and Ball converged far too early. |
| `ball_kmeans++_xf.hpp` | Same missing-braces bug in Ball's own `check_convergence`. |
| `lloyd_kmeans.hpp` | Returned the **previous** iteration's centroids whenever it converged (off by one). |
| `exponion.hpp`, `geokmeans.hpp` | Returned the **zeroed accumulator** whenever the loop hit `max_iter` — all-zero centroids, silently. |
| `_core.cpp` | `init` accepts an ndarray *or a file path*; the array route writes a real temp file at 17 significant digits (was 6, silently truncating callers' centroids) and cleans up on every exit path (was leaking on exceptions). Previously derived its temp directory from `__FILE__`, the build machine's source path. |
| `estimator.py` | `init` accepts `'random'`, an ndarray, or a path. |

Rebuild after any change to `src/`:

```bash
cd geokmeans-py && python sync_sources.py && pip install .
```

---

## Methodology

### Shared initialisation

k-means++ centroids are generated once per (dataset, k, seed), saved to
`inits/<dataset>_k<K>_seed<S>.csv` at full float64 precision, and every
implementation reads that same file. geokmeans receives the **path**, which
its C++ core reads directly — passing an array instead makes it serialise the
centroids to a temp CSV *inside* `fit()`, putting a file write in the timed
region.

### Controls

- **Single-threaded.** `scripts/_threads.py` sets the OpenMP/BLAS thread
  variables and re-execs the interpreter if numpy is already imported — those
  libraries read the variables at load time, so setting them after
  `import numpy` has no effect.
- **Trials.** Five per cell. E1 and the plot CSVs use mean ± std; the LaTeX
  tables use the median.
- **SSE ratio** is measured against `controls.sse_reference` **at the same
  seed**. Aggregating the denominator across seeds measures seed-to-seed
  spread rather than any difference between implementations.
- **Memory is not measured.** RSS sampled after `fit()` returns is neither a
  peak nor a delta.

### Why no scikit-learn baseline

The geokmeans core computes in `float`, scikit-learn in `float64`, so the two
agree only to ~2.4e-08 relative regardless of tolerance. That makes
scikit-learn unusable as an *exactness* reference. It remains a valid
*performance* comparator, which is what the twitter experiments are for.

---

## Known issues

- **Tie-breaking.** Hamerly/Annulus/Exponion/Geo occasionally differ from
  Lloyd's at ~1e-5 relative SSE, always *together* and with identical ARI —
  points equidistant from two centroids resolve differently under pruning than
  under an exact argmin. 56/60 comparisons are exact on sensor; the misses
  fall in a single (k, seed) cell.
- **Ball k-means** was excluded after scoring 5/45 exact with ARI as low as
  0.983 and its own divergent iteration counts. Its separate
  `check_convergence` is the suspect. Not yet investigated.
- **`geo_blas.py`** is a numpy reimplementation of Geo, validated bit-exact
  against the C++ (identical iterations, ARI 1.000000, SSE ratio
  1.000000000). Its BLAS kernels hit 2.9 ns/distance versus ~113 ns for the
  scalar C++, but Python per-cluster loop overhead is ~65% of its runtime, so
  it is *slower* overall. It is useful as an executable specification of Geo's
  semantics; a C++/Eigen version would be needed to test the vectorisation
  hypothesis properly.

---

## Standalone comparison scripts

Not part of any experiment; each writes its own CSV.

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
