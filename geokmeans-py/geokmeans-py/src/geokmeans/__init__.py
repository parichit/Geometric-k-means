"""geokmeans — fast, eco-friendly k-means clustering.

A Python interface to the C++ Geometric-k-means library: seven k-means
acceleration algorithms (Geometric, Lloyd, Elkan, Hamerly, Annulus, Exponion,
Ball k-means++) behind one scikit-learn-style estimator.

Based on Sharma et al. (2026), "Geometric-k-means: A Bound Free Approach to
Fast and Eco-Friendly k-means", Machine Learning 115(2):30.
"""
from __future__ import annotations

from .estimator import ALGORITHMS, GeoKMeans

try:  # populated from package metadata when installed
    from importlib.metadata import version as _version

    __version__ = _version("geokmeans")
except Exception:  # pragma: no cover - source tree / editable edge cases
    __version__ = "0.0.0"

__all__ = ["GeoKMeans", "ALGORITHMS", "__version__"]
