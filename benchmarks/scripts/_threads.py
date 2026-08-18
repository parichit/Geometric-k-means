"""Single-thread pinning. Import this before numpy/sklearn/geokmeans.

OpenMP and the BLAS backends read their thread-count environment variables
when the shared library is loaded, which happens on the first ``import numpy``.
Setting the variables after that point -- as the previous ``ensure_single_threaded``
did -- has no effect, so the benchmarks were not in fact single-threaded.

The only reliable fix in-process is to set the variables and re-exec the
interpreter before any of those libraries are imported.
"""
import os
import sys

_THREAD_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)
_SENTINEL = "GEOKMEANS_BENCH_SINGLE_THREADED"


def pin_single_threaded() -> None:
    """Force one thread everywhere, re-execing this process if necessary."""
    if os.environ.get(_SENTINEL) == "1":
        return

    for var in _THREAD_VARS:
        os.environ[var] = "1"
    os.environ[_SENTINEL] = "1"

    # numpy already loaded means a BLAS is already loaded with its own thread
    # pool, so the variables above came too late. Restart with them in place.
    if "numpy" in sys.modules:
        os.execv(sys.executable, [sys.executable] + sys.argv)


def report_thread_limits() -> str:
    """Describe the thread pools actually in effect, for the run log."""
    try:
        import threadpoolctl
    except ImportError:
        return "threadpoolctl not installed; thread limits unverified"

    # Belt and braces: cap any pool that ignored the environment.
    threadpoolctl.threadpool_limits(limits=1)
    pools = threadpoolctl.threadpool_info()
    if not pools:
        return "no native thread pools detected"
    return ", ".join(
        f"{p.get('user_api') or p.get('internal_api')}={p.get('num_threads')}"
        for p in pools
    )
