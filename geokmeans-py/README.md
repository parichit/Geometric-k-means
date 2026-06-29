# geokmeans

Fast, eco-friendly **k-means clustering** for Python. A thin, scikit-learn-style
wrapper over a C++ library implementing seven k-means algorithms — including
**Geometric-k-means**, the bound-free method from Sharma et al. (2026).

All seven algorithms return the *same* Lloyd's k-means solution; they differ in
how aggressively they prune distance computations, so you get identical results
with far less work.

| Algorithm | `algorithm=` |
|-----------|--------------|
| Geometric-k-means | `"geo"` (default) |
| Lloyd's | `"lloyd"` |
| Elkan | `"elkan"` |
| Hamerly | `"hamerly"` |
| Annulus | `"annulus"` |
| Exponion | `"exponion"` |
| Ball k-means++ | `"ball"` |

## Installation

```bash
pip install geokmeans
```

Prebuilt wheels are published for Linux, macOS, and Windows (CPython 3.9–3.13),
so no compiler is required on common platforms.

## Quickstart

```python
import numpy as np
from geokmeans import GeoKMeans

X = np.random.default_rng(0).random((1000, 20))

km = GeoKMeans(n_clusters=8, algorithm="geo", random_state=0).fit(X)

km.labels_                     # cluster index per point  (n_samples,)
km.cluster_centers_            # centroids                (n_clusters, n_features)
km.n_iter_                     # iterations to convergence
km.n_distance_calculations_    # distance computations performed

km.predict(np.random.random((5, 20)))   # assign new points
```

The API mirrors scikit-learn (`fit`, `predict`, `fit_predict`,
`cluster_centers_`, `labels_`, `n_iter_`), so it drops into existing pipelines.

## Long-running jobs

The compute runs in C++ with the GIL released, so a long `fit` won't freeze the
interpreter and you can run it off the main thread:

```python
from concurrent.futures import ProcessPoolExecutor

def cluster(X):
    return GeoKMeans(n_clusters=50, algorithm="geo", random_state=0).fit(X)

with ProcessPoolExecutor() as pool:
    result = pool.submit(cluster, X).result()
```

Pass `verbose=True` to stream convergence progress.

## Citation

If you use this package, please cite:

> Sharma, P., Stanislaw, M., Kurban, H., Kulekci, O., & Dalkilic, M. (2026).
> Geometric-k-means: A Bound Free Approach to Fast and Eco-Friendly k-means.
> *Machine Learning*, 115(2), 30. https://doi.org/10.1007/s10994-025-06891-1

## Documentation

Full documentation: https://parichit.github.io/Geometric-k-means/

## License

Apache-2.0. See [LICENSE](LICENSE).
