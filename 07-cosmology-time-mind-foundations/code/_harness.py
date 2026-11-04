"""Small shared test harness for this section's verify_*.py scripts."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


class Checks:
    def __init__(self, title):
        self.failed = 0
        self.count = 0
        print(title)
        print("=" * len(title))

    def section(self, text):
        print(f"\n-- {text}")

    def check(self, label, ok):
        self.count += 1
        if not ok:
            self.failed += 1
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")

    def finish(self):
        print(f"\n{self.count - self.failed}/{self.count} checks passed")
        sys.exit(1 if self.failed else 0)
