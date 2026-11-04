"""F07: tipping two 60 kg people with a horizontal forehead push.

The person is a rigid block that can pivot about the edge of the stance. Thresholds are found
by bisection on a simulated rocking motion, not from the formula they are compared with.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
fails = []
g = 9.81


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        fails.append(name)


def falls(F, duration, m, h, b, H, depth=None, t_max=30.0):
    """Simulate rotation about the rear stance edge. True if the COM passes over the edge.

    theta is the tilt about the pivot. Position vectors from the pivot at theta = 0:
    COM (-b, h), forehead (-b, H). The ground holds the block at theta = 0 until the net
    moment lifts it.
    """
    I = m * (4 / 3) * (h ** 2 + b ** 2)  # uniform block of height 2h, depth 2b, about an edge
    theta_c = np.arctan2(b, h)

    def moment(t, th):
        push = F if t <= duration else 0.0
        return push * (H * np.cos(th) + b * np.sin(th)) - m * g * (b * np.cos(th) - h * np.sin(th))

    def rhs(t, y):
        th, w = y
        acc = moment(t, th) / I
        if th <= 0 and acc < 0:  # still sitting on the ground
            return [0.0, 0.0] if w <= 0 else [w, acc]
        return [w, acc]

    def over(t, y):
        return y[0] - theta_c
    over.terminal = True
    over.direction = 1

    def landed(t, y):
        return y[0] + 1e-9
    landed.terminal = True
    landed.direction = -1

    t_end = min(t_max, duration + 8.0)
    sol = solve_ivp(rhs, (0, duration), [0.0, 0.0], events=over, rtol=1e-8, atol=1e-10, max_step=0.02)
    if sol.t_events[0].size:
        return True
    sol2 = solve_ivp(rhs, (duration, t_end), sol.y[:, -1], events=[over, landed], rtol=1e-8,
                     atol=1e-10, max_step=0.02)
    return bool(sol2.t_events[0].size)


def threshold(duration, m, h, b, H, hi=2000.0):
    lo = 0.0
    for _ in range(24):
        mid = 0.5 * (lo + hi)
        if falls(mid, duration, m, h, b, H):
            hi = mid
        else:
            lo = mid
    return hi


m = 60.0
people = {"2 m person": 2.0, "1 m person": 1.0}
print("Sustained push held for 5 s; forehead at 0.93 of height; COM at 0.50 of height\n")
print(f"{'case':42s} {'simulated F (N)':>16s} {'m g b / H (N)':>14s}")
results = {}
for label, b_of in (("same stance, b = 0.10 m", lambda s: 0.10), ("stance scales with height, b = 0.05 s", lambda s: 0.05 * s)):
    for name, s in people.items():
        h, H, b = 0.5 * s, 0.93 * s, b_of(s)
        F = threshold(5.0, m, h, b, H)
        results[(label, name)] = F
        pred = m * g * b / H
        print(f"{name + ', ' + label:42s} {F:16.2f} {pred:14.2f}")
        check(f"static threshold equals m g b / H ({name}, {label})", abs(F / pred - 1) < 0.01)

a = results[("same stance, b = 0.10 m", "1 m person")] / results[("same stance, b = 0.10 m", "2 m person")]
c = results[("stance scales with height, b = 0.05 s", "1 m person")] / results[("stance scales with height, b = 0.05 s", "2 m person")]
check("same stance: the short person needs twice the force", abs(a - 2) < 0.02, f"ratio {a:.3f}")
check("stance scaling with height: same force for both", abs(c - 1) < 0.01, f"ratio {c:.3f}")

# COM height has no effect on the static threshold.
print("\nCOM height varied at fixed forehead height 0.93 m, b = 0.10 m:")
vals = []
for h in (0.35, 0.50, 0.55, 0.70):
    F = threshold(5.0, m, h, 0.10, 0.93)
    vals.append(F)
    print(f"    h = {h:.2f} m -> {F:.2f} N")
check("static threshold does not depend on COM height", max(vals) / min(vals) - 1 < 0.01)

# Short push: 0.1 s shove. Now inertia and COM height matter.
print("\nShort shove lasting 0.10 s, threshold force and impulse:")
imp = {}
for label, b_of in (("b = 0.10 m", lambda s: 0.10), ("b = 0.05 s", lambda s: 0.05 * s)):
    for name, s in people.items():
        h, H, b = 0.5 * s, 0.93 * s, b_of(s)
        F = threshold(0.10, m, h, b, H, hi=20000.0)
        imp[(label, name)] = F * 0.10
        print(f"    {name}, {label}: F = {F:7.1f} N, impulse = {F * 0.10:6.1f} N s")
r1 = imp[("b = 0.10 m", "1 m person")] / imp[("b = 0.10 m", "2 m person")]
r2 = imp[("b = 0.05 s", "1 m person")] / imp[("b = 0.05 s", "2 m person")]
print(f"    impulse ratio short/tall: same stance {r1:.2f}; scaled stance {r2:.2f}")
check("short shove, same stance: short person is harder", r1 > 1.2)
check("short shove, scaled stance: tall person is harder", r2 < 0.9)

# Energy barrier and geometry that do depend on COM height.
print("\nGeometry with b = 0.10 m:")
for name, s in people.items():
    h = 0.5 * s
    theta = np.degrees(np.arctan2(0.10, h))
    barrier = m * g * (np.hypot(h, 0.10) - h)
    tau = np.sqrt((4 / 3) * (h ** 2 + 0.01) / (g * np.hypot(h, 0.10)))
    print(f"    {name}: tips past {theta:.1f} deg, energy barrier {barrier:.1f} J, fall time scale {tau:.2f} s")

# Sliding versus tipping, with an assumed shoe friction of 0.5.
mu = 0.5
print(f"\nSliding needs mu m g = {mu * m * g:.0f} N with mu = {mu}, far above the tipping force. "
      "Tipping comes first when b/H < mu.")
check("tipping precedes sliding for both people", 0.10 / 0.93 < mu and 0.10 / 1.86 < mu)

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All F07 checks passed.")
