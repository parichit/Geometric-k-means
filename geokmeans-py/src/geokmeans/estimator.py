"""scikit-learn-style estimator for the Geometric-k-means family."""
from __future__ import annotations

from typing import Optional, Union

import numpy as np

from . import _core

__all__ = ["GeoKMeans", "ALGORITHMS"]

ALGORITHMS = _core.ALGORITHMS  # ("geo", "lloyd", "elkan", "hamerly", ...)

# Accept a few friendly aliases in addition to the canonical names above.
_ALIASES = {
    "geometric": "geo",
    "geo-kmeans": "geo",
    "geokmeans": "geo",
    "kmeans": "lloyd",
    "lloyd's": "lloyd",
    "ball-kmeans": "ball",
    "ballkmeans": "ball",
}

_SEED_MAX = 2_147_483_647  # fits in C++ int; avoids the Ball -1 (time-seeded) path


class GeoKMeans:
    """k-means clustering with a choice of acceleration algorithm.

    All seven algorithms compute the *same* Lloyd's k-means solution; they
    differ only in how aggressively they prune distance computations. ``"geo"``
    (Geometric-k-means) is the bound-free method from Sharma et al. (2026).

    Parameters
    ----------
    n_clusters : int, default=8
        Number of clusters ``k``.
    algorithm : str, default="geo"
        One of ``geokmeans.ALGORITHMS``: ``"geo"``, ``"lloyd"``, ``"elkan"``,
        ``"hamerly"``, ``"annulus"``, ``"exponion"``, ``"ball"``.
    max_iter : int, default=300
        Maximum number of iterations.
    tol : float, default=1e-4
        Convergence threshold on centroid movement between iterations.
    random_state : int or None, default=None
        Seed for centroid initialisation. ``None`` draws a fresh random seed.
    verbose : bool, default=False
        If True, stream the C++ routine's progress messages to stdout.

    Attributes
    ----------
    cluster_centers_ : ndarray of shape (n_clusters, n_features)
    labels_ : ndarray of shape (n_samples,)
        0-based cluster index for each training point.
    n_iter_ : int
        Iterations run before convergence.
    n_distance_calculations_ : int
        Point-to-centroid distance computations performed (the efficiency
        metric the algorithms optimise).

    Examples
    --------
    >>> import numpy as np
    >>> from geokmeans import GeoKMeans
    >>> X = np.random.default_rng(0).random((100, 5))
    >>> km = GeoKMeans(n_clusters=3, algorithm="geo", random_state=0).fit(X)
    >>> km.labels_.shape
    (100,)
    """

    def __init__(
        self,
        n_clusters: int = 8,
        *,
        algorithm: str = "geo",
        max_iter: int = 300,
        tol: float = 1e-4,
        random_state: Optional[int] = None,
        verbose: bool = False,
    ) -> None:
        self.n_clusters = n_clusters
        self.algorithm = algorithm
        self.max_iter = max_iter
        self.tol = tol
        self.random_state = random_state
        self.verbose = verbose

    # ------------------------------------------------------------------ utils
    def _resolve_algorithm(self) -> str:
        a = str(self.algorithm).lower().strip()
        a = _ALIASES.get(a, a)
        if a not in ALGORITHMS:
            raise ValueError(
                f"algorithm must be one of {ALGORITHMS} (or an alias); got "
                f"{self.algorithm!r}"
            )
        return a

    def _resolve_seed(self) -> int:
        if self.random_state is None:
            return int(np.random.default_rng().integers(1, _SEED_MAX))
        seed = int(self.random_state)
        if not (0 < seed <= _SEED_MAX):
            # Map any int deterministically into the valid positive range.
            seed = (abs(seed) % (_SEED_MAX - 1)) + 1
        return seed

    @staticmethod
    def _check_array(X) -> np.ndarray:
        X = np.ascontiguousarray(X, dtype=np.float64)
        if X.ndim != 2:
            raise ValueError(f"Expected a 2-D array, got shape {X.shape}")
        if not np.all(np.isfinite(X)):
            raise ValueError("Input contains NaN or infinite values")
        return X

    # -------------------------------------------------------------------- API
    def fit(self, X, y=None) -> "GeoKMeans":
        """Compute clustering and store results on ``self``."""
        X = self._check_array(X)
        n_samples = X.shape[0]
        if not (1 <= self.n_clusters <= n_samples):
            raise ValueError(
                f"n_clusters={self.n_clusters} must be in [1, n_samples="
                f"{n_samples}]"
            )

        result = _core.run(
            algorithm=self._resolve_algorithm(),
            data=X,
            num_clusters=int(self.n_clusters),
            threshold=float(self.tol),
            num_iterations=int(self.max_iter),
            seed=self._resolve_seed(),
            verbose=bool(self.verbose),
        )

        self.cluster_centers_ = np.asarray(result["centroids"], dtype=np.float64)
        self.labels_ = np.asarray(result["labels"], dtype=np.int64)
        self.n_iter_ = int(result["n_iter"])
        self.n_distance_calculations_ = int(result["distance_calculations"])
        return self

    def predict(self, X) -> np.ndarray:
        """Assign each row of ``X`` to the nearest fitted centroid."""
        if not hasattr(self, "cluster_centers_"):
            raise RuntimeError("This GeoKMeans instance is not fitted yet.")
        X = self._check_array(X)
        # (n, 1, d) - (k, d) -> (n, k) squared distances
        d = ((X[:, None, :] - self.cluster_centers_[None, :, :]) ** 2).sum(axis=2)
        return d.argmin(axis=1).astype(np.int64)

    def fit_predict(self, X, y=None) -> np.ndarray:
        """Fit, then return the training labels."""
        return self.fit(X).labels_

    def __repr__(self) -> str:
        return (
            f"GeoKMeans(n_clusters={self.n_clusters}, algorithm={self.algorithm!r}, "
            f"max_iter={self.max_iter}, tol={self.tol}, "
            f"random_state={self.random_state})"
        )
