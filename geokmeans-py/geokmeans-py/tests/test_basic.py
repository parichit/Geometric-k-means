"""Tests for the geokmeans package."""
import numpy as np
import pytest

from geokmeans import ALGORITHMS, GeoKMeans


@pytest.fixture
def blobs():
    """Three well-separated 2-D Gaussian blobs, 200 points each."""
    rng = np.random.default_rng(0)
    X = np.vstack([rng.normal(m, 0.25, (200, 2))
                   for m in ([0, 0], [6, 6], [0, 6])])
    return np.ascontiguousarray(X, dtype=float)


def test_algorithms_tuple():
    assert set(ALGORITHMS) == {
        "geo", "lloyd", "elkan", "hamerly", "annulus", "exponion", "ball"
    }


@pytest.mark.parametrize("algorithm", ALGORITHMS)
def test_fit_shapes_and_attributes(blobs, algorithm):
    km = GeoKMeans(n_clusters=3, algorithm=algorithm, random_state=1).fit(blobs)
    assert km.cluster_centers_.shape == (3, 2)
    assert km.labels_.shape == (600,)
    assert set(np.unique(km.labels_)).issubset({0, 1, 2})
    assert km.n_iter_ >= 1
    assert km.n_distance_calculations_ > 0


@pytest.mark.parametrize("algorithm", ALGORITHMS)
def test_recovers_separated_clusters(blobs, algorithm):
    # With a good seed every algorithm should recover the 200/200/200 split.
    km = GeoKMeans(n_clusters=3, algorithm=algorithm, random_state=1).fit(blobs)
    sizes = sorted(np.bincount(km.labels_, minlength=3))
    assert sizes == [200, 200, 200]


def test_all_algorithms_agree(blobs):
    # Same seed -> identical clustering across algorithms (up to label perm).
    def canonical(labels):
        # relabel by first-appearance order so permutations compare equal
        mapping, out, nxt = {}, [], 0
        for v in labels:
            if v not in mapping:
                mapping[v] = nxt
                nxt += 1
            out.append(mapping[v])
        return tuple(out)

    ref = canonical(GeoKMeans(n_clusters=3, algorithm="lloyd",
                              random_state=1).fit(blobs).labels_)
    for algo in ALGORITHMS:
        got = canonical(GeoKMeans(n_clusters=3, algorithm=algo,
                                  random_state=1).fit(blobs).labels_)
        assert got == ref, f"{algo} disagrees with lloyd"


def test_geo_does_fewer_distances_than_lloyd(blobs):
    lloyd = GeoKMeans(n_clusters=3, algorithm="lloyd", random_state=1).fit(blobs)
    geo = GeoKMeans(n_clusters=3, algorithm="geo", random_state=1).fit(blobs)
    assert geo.n_distance_calculations_ < lloyd.n_distance_calculations_


def test_reproducible_with_seed(blobs):
    a = GeoKMeans(n_clusters=3, random_state=42).fit(blobs)
    b = GeoKMeans(n_clusters=3, random_state=42).fit(blobs)
    np.testing.assert_array_equal(a.labels_, b.labels_)
    np.testing.assert_allclose(a.cluster_centers_, b.cluster_centers_)


def test_predict_matches_fit_labels(blobs):
    km = GeoKMeans(n_clusters=3, random_state=1).fit(blobs)
    np.testing.assert_array_equal(km.predict(blobs), km.labels_)


def test_fit_predict(blobs):
    km = GeoKMeans(n_clusters=3, random_state=1)
    labels = km.fit_predict(blobs)
    np.testing.assert_array_equal(labels, km.labels_)


def test_aliases(blobs):
    for alias in ("geometric", "geo-kmeans", "GEO"):
        km = GeoKMeans(n_clusters=3, algorithm=alias, random_state=1).fit(blobs)
        assert km.cluster_centers_.shape == (3, 2)


def test_invalid_algorithm(blobs):
    with pytest.raises(ValueError, match="algorithm must be one of"):
        GeoKMeans(n_clusters=3, algorithm="nope").fit(blobs)


def test_invalid_n_clusters(blobs):
    with pytest.raises(ValueError):
        GeoKMeans(n_clusters=0).fit(blobs)
    with pytest.raises(ValueError):
        GeoKMeans(n_clusters=10_000).fit(blobs)


def test_rejects_nan():
    X = np.ones((10, 3))
    X[0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        GeoKMeans(n_clusters=2).fit(X)


def test_predict_before_fit():
    with pytest.raises(RuntimeError, match="not fitted"):
        GeoKMeans(n_clusters=2).predict(np.ones((3, 2)))
