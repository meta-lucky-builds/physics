"""Small helper shared by this section's verify_*.py scripts.

Puts the repo root on sys.path, exposes the shared constants as K, and keeps a
tally of checks so that a script exits non-zero when any check fails.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from common import constants as K  # noqa: E402

_failed = []
_count = [0]


def section(title):
    print(f"\n--- {title}")


def check(name, ok, detail=""):
    """Record one pass/fail test and print it."""
    _count[0] += 1
    if not ok:
        _failed.append(name)
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


def finish():
    print(f"\n{_count[0] - len(_failed)}/{_count[0]} checks passed")
    if _failed:
        print("FAILED:", "; ".join(_failed))
        sys.exit(1)
