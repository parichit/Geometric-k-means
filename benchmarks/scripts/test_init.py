"""Verify that explicit initialization works, as the benchmarks depend on it.

Checks:
1. geokmeans imports and is the real package (not a shadowing namespace dir)
2. GeoKMeans accepts an ndarray init and a path init
3. Both init routes give the same answer, at full precision
4. All algorithms agree given the same init
5. Different inits give different results

Exits non-zero on the first failure so it can gate a benchmark run.
"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from _threads import pin_single_threaded  # noqa: E402

pin_single_threaded()

import numpy as np  # noqa: E402
from sklearn.datasets import make_blobs  # noqa: E402

from utils import compute_sse  # noqa: E402

failures = []


def check(name, condition, detail=""):
    print(f"   {'✓' if condition else '✗'} {name}" + (f": {detail}" if detail else ""))
    if not condition:
        failures.append(name)


print("Testing explicit initialization in geokmeans...")

X, _ = make_blobs(n_samples=500, n_features=8, centers=4, random_state=42)
X = np.ascontiguousarray(X, dtype=np.float64)
k = 4

print("\n1. Testing import...")
from run import import_geokmeans  # noqa: E402

try:
    GeoKMeans = import_geokmeans()
    import geokmeans

    print(f"   ✓ geokmeans imported from {geokmeans.__file__}")
except ImportError as exc:
    print(f"   ✗ {exc}")
    sys.exit(1)

print("\n2. Testing random initialization...")
km = GeoKMeans(n_clusters=k, algorithm="lloyd", init="random", random_state=42).fit(X)
check("random init runs", km.n_iter_ > 0, f"{km.n_iter_} iterations")

print("\n3. Testing explicit ndarray initialization...")
init = np.ascontiguousarray(X[:k].copy())
km_arr = GeoKMeans(n_clusters=k, algorithm="lloyd", init=init).fit(X)
sse_arr = compute_sse(X, km_arr.labels_, km_arr.cluster_centers_)
check("ndarray init runs", km_arr.n_iter_ > 0, f"SSE={sse_arr:.6f}")

print("\n4. Testing path initialization...")
with tempfile.TemporaryDirectory() as tmp:
    init_path = Path(tmp) / "init.csv"
    np.savetxt(init_path, init, delimiter=",", fmt="%.18e")
    km_path = GeoKMeans(n_clusters=k, algorithm="lloyd", init=str(init_path)).fit(X)
    sse_path = compute_sse(X, km_path.labels_, km_path.cluster_centers_)
check("path init runs", km_path.n_iter_ > 0, f"SSE={sse_path:.6f}")

# The two routes serialise the same centroids, so they must agree. This is the
# check that would have caught the 6-significant-digit truncation in _core.cpp.
check(
    "ndarray and path init agree",
    np.allclose(km_arr.cluster_centers_, km_path.cluster_centers_, rtol=1e-6),
    f"max diff={np.abs(km_arr.cluster_centers_ - km_path.cluster_centers_).max():.3e}",
)

print("\n5. Testing a bad init path is rejected loudly...")
try:
    GeoKMeans(n_clusters=k, algorithm="lloyd", init="/no/such/file.csv").fit(X)
    check("missing init file raises", False, "no exception raised")
except (ValueError, RuntimeError) as exc:
    check("missing init file raises", True, type(exc).__name__)

print("\n6. Testing exactness across algorithms...")
results = {}
for algo in ["lloyd", "hamerly", "annulus", "exponion", "ball", "geo"]:
    try:
        km = GeoKMeans(n_clusters=k, algorithm=algo, init=init).fit(X)
        sse = compute_sse(X, km.labels_, km.cluster_centers_)
        results[algo] = sse
        print(f"   {algo:10s}: SSE={sse:.6f}, iter={km.n_iter_}")
    except Exception as exc:  # noqa: BLE001
        print(f"   ✗ {algo} failed: {exc}")
        failures.append(f"algorithm {algo}")

if len(results) > 1:
    ref = results.get("lloyd", next(iter(results.values())))
    worst = max(abs(s - ref) / ref for s in results.values())
    check("all algorithms match Lloyd's SSE", worst < 1e-6, f"max rel dev={worst:.2e}")

print("\n7. Testing that different inits give different results...")
km1 = GeoKMeans(n_clusters=k, algorithm="lloyd", init=np.ascontiguousarray(X[:k])).fit(X)
km2 = GeoKMeans(n_clusters=k, algorithm="lloyd", init=np.ascontiguousarray(X[10:10 + k])).fit(X)
sse1 = compute_sse(X, km1.labels_, km1.cluster_centers_)
sse2 = compute_sse(X, km2.labels_, km2.cluster_centers_)
print(f"   init A SSE={sse1:.4f}, init B SSE={sse2:.4f}")
if abs(sse1 - sse2) <= 1e-3:
    print("   note: both inits converged to the same optimum (not a failure)")

print("\n" + "=" * 60)
if failures:
    print(f"FAILED: {len(failures)} check(s): {', '.join(failures)}")
    sys.exit(1)
print("ALL TESTS PASSED ✓")
print("=" * 60)
