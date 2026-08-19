"""RAPL energy measurement (Linux only).

Intel RAPL exposes cumulative energy counters through the powercap sysfs tree:

    /sys/class/powercap/intel-rapl:0/name              -> "package-0"
    /sys/class/powercap/intel-rapl:0/energy_uj         -> microjoules, wraps
    /sys/class/powercap/intel-rapl:0/intel-rapl:0:0/   -> "core"
    /sys/class/powercap/intel-rapl:0/intel-rapl:0:1/   -> "dram"

Recent AMD parts (Zen and later) publish the same interface under the same
``intel-rapl`` names, so nothing here is Intel-specific in practice.

Three things make a naive read wrong, and each is handled below:

* **Wraparound.** ``energy_uj`` is a free-running counter that rolls over at
  ``max_energy_range_uj`` -- roughly every 60 s on a package domain under load.
  Deltas are taken modulo ``max + 1``.
* **Double counting.** ``core``/``uncore`` are *subdomains of* ``package``, and
  ``psys`` (where present) covers the whole platform including the package.
  Only ``package-*`` and ``dram`` are summed into the total; the rest are
  reported for information.
* **Permissions.** Since CVE-2020-8694 the counters are root-readable only
  (mode 0400). See ``--check`` output for the fix.

The meter is a context manager and never raises on an unmeasurable system: it
returns NaN, so a benchmark run on a laptop without RAPL still produces the
same columns.

Usage:
    python scripts/energy.py --check          # domains, permissions, 1 s sample
"""
from __future__ import annotations

import math
import os
import time
from pathlib import Path
from typing import Dict, List, Optional

#: GEOKMEANS_RAPL_ROOT relocates the powercap tree. It exists so the energy
#: path can be exercised against a synthetic tree on a machine without RAPL,
#: and is honoured by scripts/run_r_kmeans.R too so both sides agree.
RAPL_ROOT = Path(os.environ.get("GEOKMEANS_RAPL_ROOT", "/sys/class/powercap"))

#: Columns the meter always produces, so the CSV schema does not depend on the
#: machine the sweep happened to run on.
ENERGY_COLUMNS = ("energy_total_joules", "energy_pkg_joules",
                  "energy_dram_joules", "power_watts")

_NAN = float("nan")


class RaplDomain:
    """One powercap domain: a package, or a subdomain such as dram."""

    def __init__(self, path: Path, parent: Optional[str] = None) -> None:
        self.path = path
        self.parent = parent
        self.name = (path / "name").read_text().strip()
        self.max_uj = int((path / "max_energy_range_uj").read_text().strip())
        self._energy_file = path / "energy_uj"

    @property
    def kind(self) -> str:
        """'package', 'dram', 'psys' or 'other' -- decides what is summed."""
        if self.name.startswith("package"):
            return "package"
        if self.name == "dram":
            return "dram"
        if self.name == "psys":
            return "psys"
        return "other"

    def read_uj(self) -> int:
        return int(self._energy_file.read_text().strip())

    def delta_uj(self, before: int, after: int) -> int:
        """After-minus-before, corrected for counter wraparound."""
        return (after - before) % (self.max_uj + 1)

    def readable(self) -> bool:
        try:
            self.read_uj()
            return True
        except OSError:
            return False


def discover_domains(root: Path = RAPL_ROOT) -> List[RaplDomain]:
    """Every RAPL domain on the machine, packages first then their children.

    ``intel-rapl-mmio:*`` is deliberately not walked: on the parts that expose
    it, it duplicates the MSR-backed ``intel-rapl:*`` package counter.
    """
    domains: List[RaplDomain] = []
    if not root.is_dir():
        return domains
    for pkg_dir in sorted(root.glob("intel-rapl:[0-9]*")):
        try:
            pkg = RaplDomain(pkg_dir)
        except OSError:
            continue
        domains.append(pkg)
        for sub_dir in sorted(pkg_dir.glob("intel-rapl:*:*")):
            try:
                domains.append(RaplDomain(sub_dir, parent=pkg.name))
            except OSError:
                continue
    return domains


def availability(root: Path = RAPL_ROOT) -> tuple[bool, str]:
    """(usable, human-readable reason). Never raises."""
    if not root.is_dir():
        return False, (f"{root} does not exist -- RAPL is a Linux powercap "
                       f"interface; not available on this platform")
    domains = discover_domains(root)
    if not domains:
        return False, (f"no intel-rapl domains under {root} -- the CPU or "
                       f"kernel does not expose RAPL "
                       f"(try: sudo modprobe intel_rapl_common)")
    unreadable = [d.name for d in domains if not d.readable()]
    if len(unreadable) == len(domains):
        return False, ("RAPL counters are present but unreadable (mode 0400 "
                       "since CVE-2020-8694). Fix with:\n"
                       "    sudo chmod -R a+r /sys/class/powercap/intel-rapl*\n"
                       "or run the benchmark as root.")
    if unreadable:
        return True, f"readable, except: {', '.join(unreadable)}"
    return True, f"{len(domains)} domain(s): " + ", ".join(d.name for d in domains)


class EnergyMeter:
    """Context manager accumulating RAPL energy across a block.

    Reading the counters costs a few sysfs reads (~10 us), so wrapping a fit
    that takes seconds adds nothing measurable.

    >>> with EnergyMeter() as m:
    ...     do_work()
    >>> m.result()["energy_total_joules"]

    On a machine without usable counters every value is NaN and ``available``
    is False -- no exception, no missing columns.
    """

    def __init__(self, domains: Optional[List[RaplDomain]] = None) -> None:
        self.domains = discover_domains() if domains is None else domains
        self.domains = [d for d in self.domains if d.readable()]
        self.available = bool(self.domains)
        self._before: Dict[Path, int] = {}
        self._joules: Dict[str, float] = {}
        self._seconds = _NAN

    def __enter__(self) -> "EnergyMeter":
        self._t0 = time.perf_counter()
        self._before = {d.path: d.read_uj() for d in self.domains}
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        after = {d.path: d.read_uj() for d in self.domains}
        self._seconds = time.perf_counter() - self._t0
        self._joules = {}
        for d in self.domains:
            self._joules[self._key(d)] = (
                d.delta_uj(self._before[d.path], after[d.path]) / 1e6)

    @staticmethod
    def _key(d: "RaplDomain") -> str:
        return d.name if d.parent is None else f"{d.parent}/{d.name}"

    def per_domain(self) -> Dict[str, float]:
        """Joules per domain, keyed 'package-0', 'package-0/dram', ..."""
        return dict(self._joules)

    def result(self) -> Dict[str, float]:
        """The four benchmark columns. NaN throughout if RAPL is unusable."""
        if not self.available or not self._joules:
            return {c: _NAN for c in ENERGY_COLUMNS}
        pkg = sum(self._joules[self._key(d)]
                  for d in self.domains if d.kind == "package")
        dram = sum(self._joules[self._key(d)]
                   for d in self.domains if d.kind == "dram")
        total = pkg + dram
        return {
            "energy_total_joules": total,
            "energy_pkg_joules": pkg,
            "energy_dram_joules": dram,
            "power_watts": total / self._seconds if self._seconds > 0 else _NAN,
        }


def nan_result() -> Dict[str, float]:
    """Placeholder row for an implementation that was not energy-measured."""
    return {c: _NAN for c in ENERGY_COLUMNS}


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check", action="store_true",
                    help="report domains, permissions and a 1 s idle sample")
    ap.add_argument("--seconds", type=float, default=1.0)
    args = ap.parse_args()

    ok, reason = availability()
    print(f"RAPL available: {ok}\n  {reason}\n")
    if not ok:
        raise SystemExit(1)

    domains = discover_domains()
    print(f"{'domain':<24}{'kind':<10}{'counted':<10}max_energy_range_uj")
    for d in domains:
        key = d.name if d.parent is None else f"{d.parent}/{d.name}"
        counted = "yes" if d.kind in ("package", "dram") else "no"
        print(f"{key:<24}{d.kind:<10}{counted:<10}{d.max_uj}")

    if args.check:
        print(f"\nSampling {args.seconds:g}s idle ...")
        m = EnergyMeter()
        with m:
            time.sleep(args.seconds)
        for key, j in m.per_domain().items():
            print(f"  {key:<24}{j:10.3f} J")
        r = m.result()
        print(f"\n  total {r['energy_total_joules']:.3f} J "
              f"({r['power_watts']:.2f} W idle baseline)")
        if not math.isfinite(r["energy_total_joules"]):
            raise SystemExit(1)


if __name__ == "__main__":
    main()
