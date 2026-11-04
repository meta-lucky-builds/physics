"""Check that every question is answered and every verification script passes.

Run from anywhere: python tools/check.py
"""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LABELS = {"Derived", "Code", "Established", "Estimate", "Open", "Conceptual", "Speculative"}
SECTIONS = {
    "F": "01-forces-and-rigid-bodies",
    "E": "02-elasticity-and-collisions",
    "T": "03-chemistry-thermo-radiation",
    "L": "04-light-relativity-quantum",
    "A": "05-antimatter-engine",
    "B": "06-black-holes-and-spacetime",
    "C": "07-cosmology-time-mind-foundations",
}


def question_ids():
    text = (ROOT / "QUESTIONS.md").read_text(encoding="utf-8")
    return re.findall(r"^- \*\*([A-Z]\d{2})\*\*", text, flags=re.M)


def check_answers(ids):
    problems = []
    for letter, folder in SECTIONS.items():
        readme = ROOT / folder / "README.md"
        wanted = [i for i in ids if i[0] == letter]
        if not readme.exists():
            problems.append(f"{folder}: README.md is missing")
            continue
        text = readme.read_text(encoding="utf-8")
        blocks = re.split(r"^### (?=[A-Z]\d{2}\b)", text, flags=re.M)[1:]
        found = {block[:3]: block for block in blocks}
        for qid in wanted:
            block = found.get(qid)
            if block is None:
                problems.append(f"{qid}: no answer heading")
                continue
            status = re.search(r"\*\*Status:\*\*\s*(.+)", block)
            if status is None:
                problems.append(f"{qid}: no status")
            else:
                used = {s.strip() for s in status.group(1).split("+")}
                unknown = used - LABELS
                if unknown:
                    problems.append(f"{qid}: unknown status {sorted(unknown)}")
                elif "Code" in used and "code/" not in block:
                    problems.append(f"{qid}: status says Code but links no script")
            if "**Short answer:**" not in block:
                problems.append(f"{qid}: no short answer")
        extra = sorted(set(found) - set(wanted))
        if extra:
            problems.append(f"{folder}: answers for unknown IDs {extra}")
        for link in re.findall(r"\]\(([^)#\s]+)\)", text):
            if not link.startswith("http") and not (readme.parent / link).exists():
                problems.append(f"{folder}: broken link {link}")
    return problems


def run_scripts():
    problems = []
    for script in sorted(ROOT.glob("[0-9]*/code/verify_*.py")):
        name = script.relative_to(ROOT).as_posix()
        try:
            done = subprocess.run(
                [sys.executable, str(script)], cwd=ROOT, capture_output=True,
                text=True, encoding="utf-8", errors="replace", timeout=180,
            )
        except subprocess.TimeoutExpired:
            problems.append(f"{name}: timed out")
            print(f"TIMEOUT {name}")
            continue
        passed = done.returncode == 0
        print(f"{'PASS' if passed else 'FAIL'} {name}")
        if not passed:
            problems.append(f"{name}: exit {done.returncode}\n{done.stderr[-800:]}")
    return problems


def main():
    ids = question_ids()
    only = sys.argv[1].upper() if len(sys.argv) > 1 else None
    problems = check_answers(ids)
    if "--no-run" not in sys.argv:
        problems += run_scripts()
    if only in SECTIONS:
        problems = [p for p in problems if p.startswith(only) or SECTIONS[only] in p]
    print(f"\n{len(ids)} questions, {len(problems)} problems")
    for problem in problems:
        print(" -", problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
