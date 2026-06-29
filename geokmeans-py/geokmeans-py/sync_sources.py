#!/usr/bin/env python3
"""Vendor the canonical C++ sources into the Python package tree.

The single source of truth for the algorithm code lives at the repository
root in ``src/`` (shared with the CRAN package). A PyPI sdist, however, must
be self-contained: ``pip install`` from source on a user's machine cannot
reach back up to ``../../src``. So before building an sdist/wheel we copy the
headers (and a minimal Eigen) into ``src/geokmeans/_cpp/``.

Run this whenever the root C++ sources change, and always before ``build``:

    python sync_sources.py

Root stays canonical; the vendored copy is a build artifact.
"""
from __future__ import annotations

import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent          # geokmeans-py/
REPO_ROOT = HERE.parent                          # repository root
SRC = REPO_ROOT / "src"                          # canonical C++ library
EIGEN = REPO_ROOT / "eigen" / "Eigen"            # bundled Eigen headers
DEST = HERE / "src" / "geokmeans" / "_cpp"       # vendored destination

# The Ball k-means header has '++' in its filename upstream; rename on copy so
# the include path is portable (matches what the CRAN package did).
RENAME = {"ball_kmeans++_xf.hpp": "ball_kmeans_xf.hpp"}


def main() -> None:
    if not SRC.is_dir():
        raise SystemExit(f"Cannot find canonical C++ sources at {SRC}")
    if not EIGEN.is_dir():
        raise SystemExit(f"Cannot find bundled Eigen at {EIGEN}")

    if DEST.exists():
        shutil.rmtree(DEST)
    DEST.mkdir(parents=True)

    # 1. Copy every header from the root C++ library.
    headers = sorted(SRC.glob("*.hpp"))
    if not headers:
        raise SystemExit(f"No .hpp headers found in {SRC}")
    for h in headers:
        target = DEST / RENAME.get(h.name, h.name)
        shutil.copy2(h, target)
        print(f"  header  {h.name} -> {target.name}")

    # 2. Vendor Eigen (Ball k-means hard-requires it). Self-contained sdist,
    #    no network needed at install time.
    shutil.copytree(EIGEN, DEST / "eigen" / "Eigen")
    print(f"  eigen   {EIGEN} -> {DEST / 'eigen' / 'Eigen'}")

    print(f"Synced {len(headers)} headers + Eigen into {DEST}")


if __name__ == "__main__":
    main()
