# Geokmeans Benchmarking Guide

Everything needed to reproduce **Figure 1** of the JMLR MLOSS paper on a Linux
workstation. Figure 1 reports the **twitter** dataset only.

Older experiments — E1 (exactness) and E2 (algorithm comparison), both on
`sensor` — are documented separately in [`docs/E1_E2_SENSOR.md`](docs/E1_E2_SENSOR.md).
They still run, but nothing in Figure 1 comes from them.

---

## The experiments

One sweep, `scripts/bench_twitter.py`, produces all three panels:

| | What | Panel |
|---|---|---|
| **E3** | Distance computations and runtime, Geometric-*k*-means vs scikit-learn Lloyd | (a), (b) |
| **E4** | The same two metrics for **R's default `stats::kmeans`** (Hartigan-Wong) | (a), (b) |
| **E5** | **Total energy** per fit for all three, measured with Intel RAPL | (c) |

All three implementations run interleaved inside each (k, seed) cell, from the
same k-means++ centroids read off disk. `k = 100, 300, 500`, five trials each.

Figure 1, all panels with `k` on the x-axis:

- **(a)** distance computations vs `k`, log scale — line plot
- **(b)** runtime vs `k` — line plot
- **(c)** total energy vs `k` — grouped bars

---

## 1. Setup

```bash
git clone <repo> && cd Geometric-k-means

python3 -m venv benchmarks/.venv
benchmarks/.venv/bin/pip install -r benchmarks/requirements.txt

# build and install the patched geokmeans into the same environment
cd geokmeans-py && python sync_sources.py
../benchmarks/.venv/bin/pip install .
cd ..
```

Requires a C++ compiler, CMake (scikit-build-core drives the build), and **R**
on `PATH` — only base R is used, no packages beyond `stats`.

> **Run from `benchmarks/`, not the repo root.** The repo root contains a
> directory literally named `geokmeans/` (the C++ sources) which Python picks
> up as an empty namespace package and which shadows the real install.

### Dataset

```
benchmarks/data/Twitter.csv      # 583,250 x 78
```

Headers are detected and skipped; non-numeric columns are dropped with a note.

### RAPL access

Since CVE-2020-8694 the energy counters are root-readable only. Grant read
access once per boot:

```bash
sudo chmod -R a+r /sys/class/powercap/intel-rapl*
benchmarks/.venv/bin/python benchmarks/scripts/energy.py --check
```

`--check` lists the domains, says which are counted, and prints a one-second
idle sample. It exits non-zero if the counters are unusable. If no domains
appear at all, try `sudo modprobe intel_rapl_common`.

Energy is opt-out: the sweep measures it whenever RAPL is readable. Pass
`--energy` to make an unreadable counter a hard failure instead of a silent
downgrade, or `--no-energy` to skip it.

---

## 2. Verify before committing compute

```bash
cd benchmarks

# ~1 minute: 20k rows, k=[10,20], 2 trials, all three implementations
.venv/bin/python scripts/bench_twitter.py --dry-run
.venv/bin/python scripts/plot_results.py \
    --summary results/reports/e3_twitter_summary_dryrun.csv \
    --out results/figures/benchmark_summary_dryrun.png
```

Dry-run artefacts are suffixed `_dryrun` and never collide with real results.

Sanity check on the dry-run CSV: R's SSE should match scikit-learn's to a few
parts in 10⁴ or better. A larger gap means R did not receive the same data or
the same initial centroids.

---

## 3. Run

### Pre-stage the initial centroids

k-means++ is O(*knd*) per seed — on twitter at k=500 that is minutes per trial.
Generate them once, up front, so the cost is off the benchmark's critical path:

```bash
cd benchmarks && mkdir -p logs
.venv/bin/python scripts/make_inits.py \
    --dataset twitter --path data/Twitter.csv --k 100 300 500 --trials 5
```

One file per (k, trial) in `inits/`, at full float64 precision. Existing files
are reused untouched unless `--force`. `bench_twitter.py` does this pass
internally at startup too, so the step is optional.

### Launch

```bash
nohup .venv/bin/python scripts/bench_twitter.py \
      --k 100 300 500 --trials 5 \
      > logs/figure1.log 2>&1 &
```

> **Do not pipe the output through `grep`/`tail`.** A piped run measured **33%
> CPU** and produced timings inflated ~40% uniformly, against 99% CPU for the
> same workload writing straight to a file — `tqdm` redraws constantly and the
> pipeline blocks the writer. Redirect to a file instead.

**Leave the machine otherwise idle.** RAPL measures the whole package, not the
process, so anything else running is charged to whichever fit is in flight.

### Expected cost

Budget hours, not minutes. scikit-learn's per-iteration work is O(*nkd*), so
at k=500 each Lloyd pass is ~22.7 GFLOP; R's Hartigan-Wong runs to a full local
optimum and is slower still. `--warmup` defaults to `0` because at this scale
each fit is many seconds, so cold-start effects are negligible.

### Validate the timings afterwards

```bash
/usr/bin/time -v .venv/bin/python scripts/bench_twitter.py ...
# check "Percent of CPU this job got" -- want >95%
```

or while running:

```bash
ps -o pid=,etime=,time=,%cpu= -p $(pgrep -f bench_twitter.py | tail -1)
```

Below ~95% something is competing and both runtime and energy are untrustworthy.
Distance counts, iteration counts and SSE are deterministic and unaffected.

---

## 4. The figure

```bash
.venv/bin/python scripts/plot_results.py
```

Reads only the tidy summary and recomputes nothing.

```
benchmarks/results/
├── raw/e3_twitter.csv                    # per-trial, every column
├── reports/e3_twitter_summary.csv        # mean/std  <- the plot's only input
└── figures/benchmark_summary.{png,pdf}
```

Panel (c) is replaced by a placeholder if the run carries no energy data, so
the figure still builds.

---

## 5. Methodology

### Shared initialisation

k-means++ centroids are generated once per (k, seed), saved at full float64
precision, and read by all three implementations. geokmeans receives the
**path**, which its C++ core reads directly — passing an array instead makes it
serialise the centroids to a temp CSV *inside* `fit()`, putting a file write in
the timed region.

### Getting the same data into R

R is handed the dataset as a raw column-major float64 blob written by the
Python side *after* it has dropped non-numeric columns and applied
`--subsample`. So R clusters bit-identical data, and a `readBin` replaces a
multi-hundred-megabyte `read.csv`. A fresh `Rscript` runs per (k, seed) rather
than one batched pass, which keeps R's fits interleaved with the Python ones.

### The three stopping rules are not the same object

| | stops when |
|---|---|
| scikit-learn | `S <= mean_var * tol` |
| geokmeans | `sqrt(S/k) <= tol` |
| R (Hartigan-Wong) | no point transfer reduces SSE — no tolerance at all |

`S` is the summed squared centroid shift. On real data the geokmeans test is
thousands of times stricter (6005× at k=50 on sensor), so geo is run at

```
thr = sqrt(mean_var * tol / k)
```

which makes its test algebraically identical to scikit-learn's. Measured on
twitter at k=50 this took geo from 264 iterations / 20.7 s to 200 / 16.2 s
against scikit-learn's 192 / 19.8 s — from a 5% loss to a 22% win.

**R cannot be harmonised.** Hartigan-Wong runs to a local optimum, a strictly
stronger condition than either tolerance test, and `stats::kmeans` exposes no
knob for it. Its runtime therefore includes work the other two never do; state
this alongside the figure rather than papering over it. `iter.max` is set to
300 to match (R's own default is 10, which would stop the fit early).

### Distance computations

Measured for geokmeans, which counts them. Neither scikit-learn nor R exposes a
counter, so both use the analytic `(iter + 1) * n * k` — one full assignment
pass per iteration plus the initial one. That is **exact** for Lloyd and an
**upper bound** for Hartigan-Wong, whose optimal-transfer stage scans only live
clusters and whose quick-transfer stage touches two centres per point. The
column `distance_count_method` records which applies to each row, and the
figure is annotated.

### Energy

`scripts/energy.py` reads the powercap sysfs tree around the fit and nothing
else. Three things a naive read gets wrong, all handled:

- **Wraparound** — `energy_uj` is free-running and rolls over roughly every 60 s
  under load, so deltas are taken modulo `max_energy_range_uj + 1`.
- **Double counting** — `core`/`uncore` are subdomains *of* `package`, and
  `psys` covers the whole platform including it. Only `package-*` and `dram`
  are summed into `energy_total_joules`; the rest are reported for information.
- **Scope** — R samples RAPL *inside its own process*, around the `kmeans()`
  call only, so R's startup and data load are excluded exactly as they are for
  the Python implementations.

Recent AMD parts publish the same interface under the same `intel-rapl` names,
so nothing here is Intel-specific in practice.

### Controls

- **Single-threaded.** `scripts/_threads.py` sets the OpenMP/BLAS thread
  variables and re-execs the interpreter if numpy is already imported — those
  libraries read the variables at load time, so setting them after
  `import numpy` has no effect. The same variables are passed to `Rscript`.
- **Trials.** Five per cell; the summary and the figure use mean ± std.
- **Memory is not measured.** RSS sampled after `fit()` returns is neither a
  peak nor a delta.

---

## 6. Troubleshooting

**"resolved to a namespace package"** — a directory named `geokmeans` is
shadowing the install. Run from `benchmarks/`, not the repo root.

**"Rscript not found on PATH"** — install R, or drop it from the sweep with
`--impls geo sklearn`.

**R fails to start: `error while loading shared libraries: lib*.so`** — the R
on `PATH` was built against libraries this environment no longer provides.
Common on a cluster after a module or conda change. `ldd $(command -v R)`-style
diagnosis aside, the durable fix without root is to build your own:

```bash
scripts/install_r_local.sh --prefix /nobackup/$USER/r-local
source /nobackup/$USER/r-local/env.sh
```

It probes for R's dependencies, builds only the missing ones, and links
everything with an RPATH into the prefix — so the resulting R resolves its own
libraries and keeps working inside `nohup`/batch shells that never read your
profile. Re-running is safe; each component is stamped and skipped. Only base R
is needed, so the defaults (no X11, no Tcl/Tk, no readline, no ICU) are fine.
It deliberately does **not** use `--enable-R-shlib` or an external threaded
BLAS: the first costs a few percent on every call and the second would
reintroduce parallelism the benchmark's single-thread control assumes away.

**"--energy requested but RAPL is unusable"** — run
`python scripts/energy.py --check` and follow the fix it prints.

**"initial centers must be distinct" from R** — two k-means++ centroids
coincide, which means the dataset has duplicate rows and `k` is large relative
to the number of distinct ones. Regenerate with `make_inits.py --force`.

**Energy numbers look implausibly large** — something else was running. RAPL is
package-wide, not per-process.
