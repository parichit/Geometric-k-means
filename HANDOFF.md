# Handoff: geokmeans Benchmark Pipeline

## 1. Context

`geokmeans` is a library implementing seven exact k-means algorithms
(Lloyd's, Elkan, Hamerly, Annulus, Exponion, Ball k-means,
Geometric-k-means) behind one C++ core, with R (CRAN) and Python (PyPI)
front ends. Every algorithm returns the same clustering solution as
Lloyd's given identical initialization; they differ only in how many
distance computations they perform to get there.

This benchmark pipeline produces the empirical results for the "Empirical
Comparison" section of a JMLR MLOSS software paper. It needs to generate:
one compact main table (a representative dataset/k, all implementations),
one full appendix table (all datasets × all k values), an exactness
verification report, and an energy summary. **The exact target table
format is in Section 5 below — match it precisely**, since the paper's
`\label{tab:comparison}` and column structure are already written and
waiting for these numbers.

## 2. Goal

Three experiments, in priority order:

- **E1 — Exactness verification.** Confirm all seven `geokmeans`
  algorithms produce the identical clustering solution as Lloyd's, given
  identical initial centroids. This is the paper's central claim and
  must be verified, not assumed.
- **E2 — Comparison against incumbent tools.** Benchmark `geokmeans`
  against scikit-learn's `KMeans` (Lloyd and Elkan) and R's
  `stats::kmeans` — the tools practitioners already use — on wall-clock
  time, distance computations, memory, and solution quality.
- **E3 — Energy.** Wrap the E2 runs with RAPL energy instrumentation
  (package + DRAM). No additional runs needed — same executions,
  additional measurement.

## 3. Methodology (must follow precisely)

### 3.1 Shared initialization (the hard part)

All implementations in a given (dataset, $k$, seed) trial **must start
from the identical set of initial centroids** — not the same random
seed, since different libraries consume RNG state differently and will
diverge even with "the same" seed.

Implement this as:
1. For each (dataset, $k$, seed), run k-means++ **once** to generate an
   initial centroid array.
2. Save it to disk (e.g. `inits/{dataset}_k{k}_seed{seed}.csv`).
3. Every implementation in the trial reads this same file and is
   initialized from it directly, not from a seed or a fresh k-means++
   call.

Both scikit-learn (`KMeans(init=<ndarray>, n_init=1)`) and R
(`stats::kmeans(x, centers=<matrix>)`) accept an explicit centroid array
directly — confirm this works as expected with a small smoke test before
running the full sweep, rather than assuming it.

### 3.2 Tie-breaking

Ties (a point exactly equidistant from two or more centroids) must be
broken the same way across all seven `geokmeans` algorithms: lowest
centroid index. **Verify this is actually true in the current codebase**
before relying on it — if any algorithm breaks ties differently, either
fix it or scope the exactness claim accordingly and flag it back to me.

### 3.3 Convergence tolerance

scikit-learn's `tol` is relative inertia change; `geokmeans`'s stopping
criterion is centroid movement. Don't assume matching the numeric
tolerance value produces matching stopping behavior — validate by
comparing iteration counts across implementations on a few trial runs,
and adjust tolerances until they agree, or document the discrepancy.

### 3.4 Controls

- Single-threaded: `OMP_NUM_THREADS=1`, and equivalent for any other
  implementation that parallelizes by default (scikit-learn does).
- `n_init=1` everywhere (we're supplying explicit init, per 3.1).
- CPU governor pinned to `performance`; record full machine spec
  (CPU model, RAM, OS, compiler/library versions) to
  `results/machine_spec.md`.
- 5 trials per (dataset, algorithm, $k$) configuration; report
  median and IQR.
- E1 uses its own denser sweep (Section 3.5) — separate from E2/E3's
  5-trial default.

### 3.5 E1 specifics

For each dataset × $k \in \{5, 25, 100\}$ × 10 seeds: run all seven
`geokmeans` algorithms from the identical initial centroids (Section
3.1). Record final SSE, iteration count, and the full label vector for
each. Compute, per algorithm relative to Lloyd's:
- max relative SSE deviation (expect ~1e-12, floating-point noise)
- Adjusted Rand Index between label vectors (expect 1.0)

If either deviates meaningfully from the expected value on any run, stop
and flag it — don't average it away or silently exclude it.

## 4. Configuration

Everything below — datasets, metrics, algorithms, $k$ values, and
**outputs** — must be driven by a single editable config file, not
hardcoded. Use the schema below as a starting point; adjust as needed,
but preserve the principle: I edit `config.yaml`, not pipeline code, to
change what's measured or what's produced.

```yaml
# benchmarks/config.yaml

datasets:
  - name: breast_cancer
    path: data/breast_cancer.csv
  - name: credit_risk
    path: data/credit_risk.csv
  - name: sensor
    path: data/sensor.csv
  - name: mnist784
    path: data/mnist784.csv
  - name: covertype
    path: data/covertype.csv
  - name: gaussian_1m_d50
    generator: gaussian_mixture
    n: 1000000
    d: 50
    k_true: 50
    seed: 0
  # EDIT ME: add, remove, or replace datasets freely.
  # `path` datasets load from CSV; `generator` datasets are synthesized.

k_values: [10, 50, 200]

implementations:
  baselines:
    - {name: sklearn_lloyd,   library: scikit-learn, algorithm: lloyd}
    - {name: sklearn_elkan,   library: scikit-learn, algorithm: elkan}
    - {name: r_stats_kmeans,  library: stats,        algorithm: hartigan-wong}
  geokmeans:
    [lloyd, elkan, hamerly, annulus, exponion, ball, geometric]

metrics:
  # Default set. Add/remove freely; each must have a computable
  # definition in the pipeline (see Section 6).
  - wall_clock_seconds
  - n_distance_calculations
  - peak_rss_mb
  - sse
  - sse_ratio_vs_lloyd
  - n_iterations
  - energy_joules          # RAPL, Linux only
  - ari_vs_lloyd           # E1 only
  - max_relative_sse_deviation   # E1 only

controls:
  threads: 1
  cpu_governor: performance
  trials_per_config: 5
  aggregate: median_iqr

exactness_check:
  k_values: [5, 25, 100]
  seeds: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
  tie_break_rule: lowest_centroid_index   # CONFIRM against actual code

outputs:
  # This list is the extension point. Each entry has a `type`; the
  # pipeline should dispatch on `type` via a small registry, so adding
  # a new output doesn't require touching experiment-running code.

  - id: table1_main
    type: latex_table
    dataset: mnist784        # EDIT: pick the representative dataset
    k: 50                    # EDIT: pick the representative k
    columns: [wall_clock_seconds, n_distance_calculations, peak_rss_mb, sse_ratio_vs_lloyd, energy_joules]
    rows:
      [sklearn_lloyd, sklearn_elkan, r_stats_kmeans,
       geokmeans_lloyd, geokmeans_elkan, geokmeans_hamerly,
       geokmeans_annulus, geokmeans_exponion, geokmeans_ball,
       geokmeans_geometric]
    output_path: results/tables/table1_main.tex

  - id: appendix_full_matrix
    type: latex_table
    datasets: all
    k_values: all
    columns: [wall_clock_seconds, n_distance_calculations, peak_rss_mb, sse_ratio_vs_lloyd, energy_joules]
    output_path: results/tables/appendix_matrix.tex

  - id: exactness_report
    type: markdown_report
    contents: [ari_vs_lloyd, max_relative_sse_deviation]
    output_path: results/reports/exactness_report.md

  - id: energy_summary_figure
    type: figure
    kind: bar_chart
    metric: energy_joules
    group_by: implementation
    dataset: mnist784
    output_path: results/figures/energy_summary.pdf
```

## 5. Required output: Table 1 exact format

This must match what's already drafted in the paper. Columns, in order:

| Column | Meaning |
|---|---|
| Implementation | e.g. "scikit-learn (Lloyd)", "geokmeans (Hamerly)" |
| Wall-clock (s) | median of trials |
| Distances (×10⁶) | median distance computations, in millions |
| Peak RSS (MB) | median peak resident memory |
| SSE ratio | final SSE ÷ that implementation's own Lloyd's-equivalent SSE (pinned to 1.000 for the Lloyd baselines) |
| Energy (J) | median RAPL package+DRAM energy |

Row order: the three baselines first (scikit-learn Lloyd, scikit-learn
Elkan, R `stats::kmeans`), then a rule, then all seven `geokmeans`
algorithms in the order listed in Section 4. Output as a LaTeX
`tabular` (booktabs style: `\toprule`/`\midrule`/`\bottomrule`), ready to
paste into the paper — match the structure of the placeholder already in
`geokmeans_mloss.tex` (`\label{tab:comparison}`).

## 6. Architecture requirement: separate "run" from "report"

Split the pipeline into two independent stages:

1. **`run.py`** (or equivalent) — executes all configured trials,
   writes raw results to `results/raw/*.csv` (one row per trial: dataset,
   implementation, $k$, seed, and every metric value). This is the
   expensive step.
2. **`report.py`** (or equivalent) — reads `results/raw/*.csv` and the
   `outputs:` section of the config, and generates every requested
   table/figure/report. This is cheap and should be safely re-runnable
   any time I edit the `outputs:` config, without re-running experiments.

This is the mechanism that satisfies "I need to be able to add/edit
outputs": I edit `config.yaml`'s `outputs:` list and re-run `report.py`
only.

## 7. Out of scope (don't build these)

- GPU comparisons (faiss, cuML) — different hardware class.
- Accuracy/quality benchmarks beyond exactness (NMI vs. ground truth,
  etc.) — redundant, since exactness guarantees identical partitions.
- Comparison against the individual algorithms' original research code
  (Hamerly's `fast-kmeans`, `eakmeans`, etc.) — already done in the
  method paper; this pipeline compares against scikit-learn and R only.
- Multithreaded scaling curves.
- Sweeps over tolerance/initialization schemes beyond what's needed for
  Section 3.3's validation.

## 8. Open items to confirm during implementation

These are unresolved facts from the paper-writing side that affect this
pipeline directly — confirm each rather than assuming:

- **Tie-breaking rule**: confirmed as lowest-centroid-index across all
  seven algorithms? (Section 3.2)
- **R dispatcher signature**: is it `geokmeans(X, k, algorithm)` or does
  it use `centers` (matching `stats::kmeans`'s convention) or another
  name? Needed to call the R package correctly from the harness.
- **R baseline choice**: `stats::kmeans` defaults to Hartigan-Wong: confirm
  this (not Lloyd or MacQueen) is the intended baseline variant.
- **Representative dataset/$k$ for Table 1**: currently placeholder
  `mnist784` / `k=50` in the config — replace with the actual choice.
- **Machine**: which physical machine this runs on, and whether RAPL is
  accessible (requires Linux, appropriate permissions, and Intel/AMD
  hardware support).

Flag anything else that seems ambiguous or under-specified rather than
guessing.
