"""F08: can mu exceed 1, and does atomically flat mean more friction?

Checks: (1) Tabor's junction growth by solving for the sliding force at fixed load, (2) adhesion
makes F/N diverge at small load, (3) commensurate vs incommensurate atomic contact, rigid and
with an elastic chain (Frenkel-Kontorova) relaxed under a rising force.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        fails.append(name)


# 1. Junction growth. Yield: p^2 + alpha s^2 = p0^2, p = W/A, s = F/A. Interface shear limit
# s_i = c s0 with s0 = p0 / sqrt(alpha) the bulk shear strength. Solve for F at fixed load.
p0, W = 1.0, 1.0


def area_for(F, alpha):
    """Area that satisfies the yield condition for tangential force F at load W."""
    return np.sqrt(W ** 2 + alpha * F ** 2) / p0


def mu_numeric(c, alpha):
    s0 = p0 / np.sqrt(alpha)
    # Interface slides when F/A(F) reaches c s0. F/A rises monotonically towards s0.
    target = c * s0
    return brentq(lambda F: F / area_for(F, alpha) - target, 1e-9, 1e6, xtol=1e-12) / W


print("Junction growth, mu at gross sliding (solved numerically):")
print(f"{'c = s_i/s0':>10} " + " ".join(f"alpha={a:<4}" for a in (3, 9, 12)))
for c in (0.2, 0.5, 0.8, 0.95, 0.99, 0.999):
    print(f"{c:10.3f} " + " ".join(f"{mu_numeric(c, a):10.3f}" for a in (3, 9, 12)))
worst = max(abs(mu_numeric(c, a) * np.sqrt(a * (1 - c ** 2)) / c - 1)
            for c in (0.2, 0.5, 0.95, 0.99) for a in (3, 9, 12))
check("numerical mu matches c / sqrt(alpha (1 - c^2))", worst < 1e-6, f"worst fractional error {worst:.1e}")
m95 = [mu_numeric(0.95, a) for a in (3, 9, 12)]
check("5% weaker interface gives mu of order unity, as Tabor (1959) reports",
      all(0.5 < x < 2 for x in m95), f"{np.round(m95, 2)}")
check("a nearly clean interface gives mu well above 1", mu_numeric(0.999, 9) > 5, f"{mu_numeric(0.999, 9):.1f}")

# 2. Adhesion: real area A = A_0 + W/p0 under a constant adhesive contact A_0. F = s_i A.
s_i, A0, H = 200e6, 1e-8, 1e9   # Pa, m^2, Pa; soft-metal values for illustration
print("\nAdhesive contact: F = tau (A0 + N/H), tau = 200 MPa, A0 = 0.01 mm^2, H = 1 GPa")
vals = []
for N in (1e-3, 1e-1, 1.0, 10.0, 100.0, 1000.0):
    F = s_i * (A0 + N / H)
    vals.append(F / N)
    print(f"    N = {N:8.3f} N: F = {F:9.3f} N, F/N = {F / N:9.3f}")
check("F/N exceeds 1 at low load and settles to tau/H = 0.2 at high load",
      vals[0] > 1 and abs(vals[-1] - s_i / H) < 0.01)

# van der Waals pressure between flat surfaces, P = A_H / (6 pi D^3).
A_H, D = 1e-19, 0.3e-9
P = A_H / (6 * np.pi * D ** 3)
print(f"\nvan der Waals pressure at D = 0.3 nm, A_H = 1e-19 J: {P / 1e6:.0f} MPa; over 1 cm^2: {P * 1e-4 / 1e3:.1f} kN")

# 3a. Rigid chain of N atoms on a sinusoidal substrate of period 1. Force needed to start sliding.


def rigid_force_per_atom(N, spacing):
    x0 = np.linspace(0, 1, 2001)
    pos = np.arange(N) * spacing
    total = np.abs(np.sin(2 * np.pi * (x0[:, None] + pos)).sum(axis=1))
    return total.max() / N


golden = (np.sqrt(5) - 1) / 2
print("\nRigid chain, peak lateral force per atom (units of the single-atom force):")
print(f"{'N':>7} {'aligned':>10} {'golden-ratio spacing':>22}")
vals = []
for N in (10, 100, 1000, 10000):
    a, b = rigid_force_per_atom(N, 1.0), rigid_force_per_atom(N, golden)
    vals.append(b)
    print(f"{N:7d} {a:10.4f} {b:22.6f}")
check("aligned contact: force per atom does not fall with size", rigid_force_per_atom(10000, 1.0) > 0.99)
check("misaligned contact: force per atom falls roughly as 1/N", vals[-1] < vals[0] / 100,
      f"{vals[0]:.4f} to {vals[-1]:.6f}")

# 3b. Elastic chain (Frenkel-Kontorova). Overdamped dynamics under a rising force.


def depinning_force(K, spacing_num, spacing_den, V=1.0):
    """Periodic chain of spacing_den atoms over spacing_num substrate periods.

    Returns force per atom, in units of the single-atom force amplitude V/(2 pi), at which the
    chain slides steadily. Equation of motion: dx/dt = K (x[i+1] + x[i-1] - 2 x[i])
    - amp sin(2 pi x[i]) + F, overdamped, with unit mobility.
    """
    n = spacing_den
    a = spacing_num / spacing_den
    amp = V / (2 * np.pi)
    x = np.arange(n) * a + 0.1

    def f(x, F):
        right = np.roll(x, -1)
        right[-1] += spacing_num
        left = np.roll(x, 1)
        left[0] -= spacing_num
        return K * (right + left - 2 * x) - amp * np.sin(2 * np.pi * x) + F

    dt = min(0.2 / (4 * K + amp * 2 * np.pi), 0.2)
    steps = int(300.0 / dt)     # 300 time units at each force level
    mark = 2 * steps // 3
    # Sliding test: mean speed over the last third compared with the free-chain speed F.
    for F in amp * np.logspace(np.log10(0.005), np.log10(1.25), 25):
        for i in range(steps):
            x = x + dt * f(x, F)
            if i == mark:
                x_mark = x.mean()
        speed = (x.mean() - x_mark) / (dt * (steps - mark - 1))
        if speed > 0.1 * F:
            return F / amp
    return 1.25


print("\nElastic chain, force per atom to slide / single-atom force, N = 55 on 34 periods (a near the golden ratio):")
results = {}
for K in (0.02, 0.05, 0.3, 1.0, 3.0):
    results[K] = depinning_force(K, 34, 55)
    print(f"    spring K = {K:5.2f}: {results[K]:.3f}")
aligned = depinning_force(0.3, 55, 55)
print(f"    aligned (a = 1), K = 0.30: {aligned:.3f}")
check("aligned chain needs the full single-atom force", aligned > 0.9)
check("stiff misaligned chain (K = 3) slides almost freely", results[3.0] < 0.05)
check("soft misaligned chain pins: the Aubry transition", results[0.02] > 10 * results[3.0])

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All F08 checks passed.")
