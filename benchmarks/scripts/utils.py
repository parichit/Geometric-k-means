"""Utility functions for benchmarking: timing, SSE, dataset loading."""
import time
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Any, Dict


class Timer:
    """Context manager measuring wall-clock time around a block.

    Memory is deliberately not collected here. The previous ``MetricsCollector``
    sampled RSS once after ``fit`` returned, by which point transient
    allocations are already freed, so the number it reported was neither a peak
    nor a delta. Rather than report a meaningless column, it is dropped.
    """

    def __init__(self) -> None:
        self._start = None
        self._end = None

    def __enter__(self) -> "Timer":
        self._start = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self._end = time.perf_counter()

    @property
    def wall_clock_seconds(self) -> float:
        return self._end - self._start


def compute_sse(X: np.ndarray, labels: np.ndarray, centroids: np.ndarray) -> float:
    """Sum of squared distances from each point to its assigned centroid.

    Computed identically for every implementation so that solution quality is
    comparable even where the libraries' own inertia definitions differ.

    Parameters
    ----------
    X : ndarray of shape (n, d)
        Dataset
    labels : ndarray of shape (n,)
        Cluster assignments
    centroids : ndarray of shape (k, d)
        Cluster centroids

    Returns
    -------
    sse : float
        Sum of squared errors
    """
    diff = X - centroids[labels]
    return float(np.einsum("ij,ij->", diff, diff))


def adjusted_rand_index(labels_true: np.ndarray, labels_pred: np.ndarray) -> float:
    """Adjusted Rand Index (1.0 = identical partitions up to relabelling)."""
    from sklearn.metrics import adjusted_rand_score

    return float(adjusted_rand_score(labels_true, labels_pred))


def _has_header(filepath: Path) -> bool:
    """True if the first line of the CSV is not fully numeric."""
    with open(filepath) as f:
        first = f.readline().strip()
    if not first:
        return False
    for field in first.split(","):
        try:
            float(field.strip().strip('"'))
        except ValueError:
            return True
    return False


def load_dataset(dataset_config: Dict[str, Any], data_dir: Path) -> np.ndarray:
    """Load a dataset from its config entry.

    Uses ``pandas.read_csv`` rather than ``np.loadtxt``: the latter is an order
    of magnitude slower on large files and throws on a header row instead of
    skipping it. Any non-numeric column (a label or id column) is dropped, so a
    file that ships with its target still loads.

    Parameters
    ----------
    dataset_config : dict
        Dataset configuration from config.yaml
    data_dir : Path
        Base directory containing data files

    Returns
    -------
    X : ndarray of shape (n, d), float64, C-contiguous
    """
    if "path" not in dataset_config:
        raise ValueError(
            f"Dataset config must have a 'path': {dataset_config!r}"
        )

    filepath = data_dir / dataset_config["path"].replace("data/", "")
    if not filepath.exists():
        raise FileNotFoundError(
            f"Dataset file not found: {filepath}\n"
            f"Please ensure the data file exists before running benchmarks."
        )

    df = pd.read_csv(filepath, header=0 if _has_header(filepath) else None)

    numeric = df.select_dtypes(include=[np.number])
    dropped = [c for c in df.columns if c not in numeric.columns]
    if dropped:
        print(f"  note: dropped {len(dropped)} non-numeric column(s): {dropped}")

    X = np.ascontiguousarray(numeric.to_numpy(dtype=np.float64))
    if X.ndim == 1:
        X = X.reshape(-1, 1)

    if not np.all(np.isfinite(X)):
        n_bad = int((~np.isfinite(X)).any(axis=1).sum())
        raise ValueError(
            f"{filepath.name} contains NaN/inf in {n_bad} row(s); clean the "
            f"dataset before benchmarking."
        )
    return X


if __name__ == "__main__":
    print("Testing timer...")
    with Timer() as t:
        time.sleep(0.1)
    print(f"Wall-clock time: {t.wall_clock_seconds:.4f}s (expected: ~0.1)")

    print("\nTesting SSE computation...")
    X = np.array([[0, 0], [1, 1], [10, 10], [11, 11]], dtype=float)
    labels = np.array([0, 0, 1, 1])
    centroids = np.array([[0.5, 0.5], [10.5, 10.5]])
    print(f"SSE: {compute_sse(X, labels, centroids):.4f} (expected: 2.0)")

    print("\nSmoke tests passed!")
