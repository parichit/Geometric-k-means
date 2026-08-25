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

import os
import shutil
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent          # geokmeans-py/
REPO_ROOT = HERE.parent                          # repository root
SRC = REPO_ROOT / "src"                          # canonical C++ library
EIGEN = REPO_ROOT / "eigen" / "Eigen"            # bundled Eigen headers
DEST = HERE / "src" / "geokmeans" / "_cpp"       # vendored destination

# The Ball k-means header has '++' in its filename upstream; rename on copy so
# the include path is portable (matches what the CRAN package did).
RENAME = {"ball_kmeans++_xf.hpp": "ball_kmeans_xf.hpp"}


def clear_dest(dest: Path) -> None:
    """Remove the vendored tree, tolerating NFS silly-rename files.

    On a network filesystem -- this repo is often built on /nobackup or a
    similar cluster scratch mount -- deleting a file another process still
    holds open does not unlink it. NFS renames it to .nfsXXXXXXXX and keeps it
    alive, so rmtree empties a directory, calls rmdir, and finds that file
    sitting in it: OSError 39, "Directory not empty".

    Retrying wins once the holder exits. If it does not, renaming the tree
    aside always works -- rename has no such restriction -- so the build can
    proceed regardless and the leftovers are reported rather than silently
    left behind.
    """
    for attempt in range(3):
        try:
            shutil.rmtree(dest)
            return
        except FileNotFoundError:
            return
        except OSError as exc:
            if attempt == 2:
                break
            print(f"  {dest.name}: {exc.strerror}; retrying in 1s "
                  f"(NFS may still be holding a deleted file open)")
            time.sleep(1.0)

    stale = dest.with_name(f"{dest.name}.stale.{os.getpid()}")
    dest.rename(stale)
    print(f"  could not delete {dest.name}; moved it aside as {stale.name}")
    shutil.rmtree(stale, ignore_errors=True)
    if stale.exists():
        print(f"  NOTE: {stale} is still on disk because something holds a "
              f"file in it open.\n        Find it with: lsof +D {stale}\n"
              f"        Then: rm -rf {stale}")


def main() -> None:
    if not SRC.is_dir():
        raise SystemExit(f"Cannot find canonical C++ sources at {SRC}")
    if not EIGEN.is_dir():
        raise SystemExit(f"Cannot find bundled Eigen at {EIGEN}")

    if DEST.exists():
        clear_dest(DEST)
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
