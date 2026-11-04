"""F02: which way does friction act when you press?

A block on a surface of any orientation is simulated with a penalty normal contact and an
elastic-slip tangential contact (a spring that slips once its force reaches mu N). Nothing
in the contact law is told the direction friction 'should' take.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
fails = []
g = 9.81


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        fails.append(name)


def simulate(F_ext, normal, mu=0.5, m=1.0, T=1.5, dt=4e-5, gravity=True):
    """Return time-average of last 10% of contact forces and the tangential drift.

    normal: unit vector out of the surface, pointing towards the block.
    F_ext : applied force vector (N) in world axes (x horizontal, y up).
    """
    n = np.array(normal, float) / np.linalg.norm(normal)
    t = np.array([n[1], -n[0]])
    grav = np.array([0.0, -m * g]) if gravity else np.zeros(2)
    F = np.array(F_ext, float) + grav
    kn, cn = 1e6, 2 * np.sqrt(1e6 * m)
    kt, ct = 2e5, 2 * np.sqrt(2e5 * m)
    x = np.zeros(2)     # displacement; x.n is height above the surface
    v = np.zeros(2)
    anchor = 0.0        # tangential position of the sticking spring's anchor
    N_hist, f_hist = [], []
    steps = int(T / dt)
    for i in range(steps):
        xn, xt = x @ n, x @ t
        vn, vt = v @ n, v @ t
        N = max(0.0, -kn * xn - cn * vn) if xn < 0 else 0.0
        f = -kt * (xt - anchor) - ct * vt
        limit = mu * N
        if abs(f) > limit:  # slip: move the anchor so the spring force sits at the limit
            f = np.sign(f) * limit
            anchor = xt + (f + ct * vt) / kt
        a = (F + N * n + f * t) / m
        v = v + dt * a
        x = x + dt * v
        if i > 0.9 * steps:
            N_hist.append(N)
            f_hist.append(f * t)
    drift = (x @ t)
    return np.mean(N_hist), np.mean(f_hist, axis=0), drift, v @ t


floor = (0, 1)
wall = (1, 0)  # block to the right of a vertical wall

print("Case 1: press straight down on a floor with 20 N")
N, f, d, v = simulate((0, -20), floor)
print(f"    N = {N:.2f} N up, friction = {f} N")
check("no friction at all; only the normal force goes up, N = m g + F",
      np.linalg.norm(f) < 1e-6 and abs(N - (g + 20)) < 0.01)

print("Case 2: press down 20 N and push right 5 N (limit mu N = 14.7 N)")
N, f, d, v = simulate((5, -20), floor)
print(f"    N = {N:.2f} N, friction = {np.round(f, 3)} N, drift = {d * 1e3:.3f} mm")
check("friction balances only the sideways part and points left", abs(f[0] + 5) < 0.01 and abs(f[1]) < 1e-6)

print("Case 3: pressing harder raises the limit: push right 8 N with 0 N press, then with 20 N press")
_, _, d0, v0 = simulate((8, 0), floor)
_, _, d1, v1 = simulate((8, -20), floor)
print(f"    no press: speed {abs(v0):.2f} m/s; with press: speed {abs(v1):.4f} m/s")
check("the press stops the slide by raising mu N", abs(v0) > 1 and abs(v1) < 1e-3)

print("Case 4: slide-free block pressed against a wall with 30 N, gravity down")
N, f, d, v = simulate((-30, 0), wall)
print(f"    N = {N:.2f} N, friction = {np.round(f, 3)} N")
check("friction points up on a wall, equal to m g", abs(f[1] - g) < 0.01 and abs(f[0]) < 1e-6)

print("Case 5: on that wall add a 5 N downward force")
N, f, d, v = simulate((-30, -5), wall)
print(f"    N = {N:.2f} N, friction = {np.round(f, 3)} N")
check("downward force on a vertical surface: friction does act upward", abs(f[1] - (g + 5)) < 0.01)

print("Case 6: friction cone. Weightless block, floor, mu = 0.5, cone half-angle 26.57 deg")
for scale in (10, 1000):
    N, f, d, v = simulate((scale * np.sin(np.radians(25)), -scale * np.cos(np.radians(25))), floor,
                          gravity=False, T=0.5)
    print(f"    {scale:5d} N at 25 deg from the normal: slip speed {abs(v):.4f} m/s")
    check(f"inside the cone it never slides ({scale} N)", abs(v) < 1e-3)
    N, f, d, v = simulate((scale * np.sin(np.radians(28)), -scale * np.cos(np.radians(28))), floor,
                          gravity=False, T=0.5)
    print(f"    {scale:5d} N at 28 deg from the normal: slip speed {abs(v):.3f} m/s")
    check(f"just outside the cone it slides ({scale} N)", abs(v) > 0.05)

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All F02 checks passed.")
