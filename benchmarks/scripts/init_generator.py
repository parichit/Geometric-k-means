"""Generate shared initial centroids using k-means++.

For each (dataset, k, seed), generate initial centroids ONCE and save to disk.
All implementations in E1 and E2 read these same centroids, so every algorithm
starts from an identical point and the comparison is fair.
"""
import numpy as np
from pathlib import Path


def kmeans_plus_plus_init(X: np.ndarray, k: int, seed: int) -> np.ndarray:
    """Generate k initial centroids using k-means++.

    Keeps a running minimum squared distance instead of recomputing distances
    to every previously chosen centroid on each step, which takes the cost from
    O(k^2 n d) down to O(k n d).

    Parameters
    ----------
    X : ndarray of shape (n, d)
        Dataset
    k : int
        Number of clusters
    seed : int
        Random seed

    Returns
    -------
    centroids : ndarray of shape (k, d)
        Initial centroids
    """
    if not 1 <= k <= X.shape[0]:
        raise ValueError(f"k={k} must be in [1, n_samples={X.shape[0]}]")

    rng = np.random.RandomState(seed)
    n_samples, n_features = X.shape

    centroids = np.empty((k, n_features), dtype=np.float64)
    centroids[0] = X[rng.randint(n_samples)]

    # Squared distance from every point to the nearest centroid chosen so far.
    closest_sq = np.sum((X - centroids[0]) ** 2, axis=1)

    for i in range(1, k):
        total = closest_sq.sum()
        if total <= 0:
            # Fewer than i distinct points: fall back to uniform sampling so we
            # still return k rows rather than dividing by zero.
            centroids[i] = X[rng.randint(n_samples)]
        else:
            cumprobs = np.cumsum(closest_sq / total)
            idx = int(np.searchsorted(cumprobs, rng.rand()))
            idx = min(idx, n_samples - 1)  # guard against fp overshoot
            centroids[i] = X[idx]

        np.minimum(closest_sq, np.sum((X - centroids[i]) ** 2, axis=1),
                   out=closest_sq)

    return centroids


def generate_and_save_init(
    X: np.ndarray,
    dataset_name: str,
    k: int,
    seed: int,
    output_dir: Path,
) -> Path:
    """Generate initial centroids and save to CSV at full float64 precision.

    Returns
    -------
    Path
        Path to saved initialization file
    """
    centroids = kmeans_plus_plus_init(X, k, seed)

    output_dir.mkdir(parents=True, exist_ok=True)
    filepath = output_dir / f"{dataset_name}_k{k}_seed{seed}.csv"
    np.savetxt(filepath, centroids, delimiter=",", fmt="%.18e")

    return filepath


def load_init(filepath: Path) -> np.ndarray:
    """Load initial centroids from CSV, always as a 2-D (k x d) array."""
    centroids = np.loadtxt(filepath, delimiter=",", ndmin=2)
    return np.ascontiguousarray(centroids, dtype=np.float64)


if __name__ == "__main__":
    # Smoke test
    from sklearn.datasets import make_blobs

    X, _ = make_blobs(n_samples=100, n_features=5, centers=3, random_state=42)

    print("Testing k-means++ initialization...")
    centroids = kmeans_plus_plus_init(X, k=3, seed=0)
    print(f"Generated centroids shape: {centroids.shape}")

    print("Checking determinism...")
    assert np.array_equal(centroids, kmeans_plus_plus_init(X, k=3, seed=0))
    print("  same seed -> same centroids ✓")

    output_dir = Path(__file__).parent.parent / "inits"
    filepath = generate_and_save_init(X, "test", k=3, seed=0, output_dir=output_dir)
    print(f"Saved to: {filepath}")

    loaded = load_init(filepath)
    print(f"Loaded centroids shape: {loaded.shape}")
    print(f"Round-trip exact: {np.array_equal(centroids, loaded)}")

    filepath.unlink()
    print("Smoke test passed!")
