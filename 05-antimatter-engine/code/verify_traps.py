"""A13, A17, A21(E): what really confines and sorts electrons and positrons.

Orbits are integrated with the Boris scheme using only the Lorentz force. Nothing about
"which sign is trapped" is assumed; each case is classified from the trajectory.
Dimensionless units: particle mass 1, charge +1 or -1, reference cyclotron frequency 1.
"""
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from common.constants import c, k_B, m_e  # noqa: E402
from scipy.constants import alpha, epsilon_0, physical_constants  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


def boris(q, x, v, field, dt, steps, escape=200.0, record=False):
    """Non-relativistic Boris pusher in plain floats. field(x, y, z) gives (Ex,Ey,Ez,Bx,By,Bz)."""
    x0, x1, x2 = map(float, x)
    v0, v1, v2 = map(float, v)
    r_max, z_max, escaped = 0.0, 0.0, False
    track = []
    half = 0.5 * dt * q
    for _ in range(steps):
        Ex, Ey, Ez, Bx, By, Bz = field(x0, x1, x2)
        um0, um1, um2 = v0 + half * Ex, v1 + half * Ey, v2 + half * Ez
        t0, t1, t2 = half * Bx, half * By, half * Bz
        f = 2 / (1 + t0 * t0 + t1 * t1 + t2 * t2)
        s0, s1, s2 = f * t0, f * t1, f * t2
        up0 = um0 + um1 * t2 - um2 * t1
        up1 = um1 + um2 * t0 - um0 * t2
        up2 = um2 + um0 * t1 - um1 * t0
        uq0 = um0 + up1 * s2 - up2 * s1
        uq1 = um1 + up2 * s0 - up0 * s2
        uq2 = um2 + up0 * s1 - up1 * s0
        v0, v1, v2 = uq0 + half * Ex, uq1 + half * Ey, uq2 + half * Ez
        x0, x1, x2 = x0 + dt * v0, x1 + dt * v1, x2 + dt * v2
        r_max = max(r_max, math.hypot(x0, x1))
        z_max = max(z_max, abs(x2))
        if record:
            track.append((x0, x1, x2))
        if r_max > escape or z_max > escape:
            escaped = True
            break
    if record:
        return r_max, z_max, escaped, np.array(track)
    return r_max, z_max, escaped


print("1. Penning trap: uniform B plus quadrupole E. Try every sign combination.")
print("   Potential phi = s (w_z^2/2)(z^2 - r^2/2), so s=+1 is a well for positive charge.")
print("     q    B     s     r_max     z_max     result")
w_z = 0.1


def penning(B0, s):
    k = s * w_z**2

    def f(x, y, z):
        return k * x / 2, k * y / 2, -k * z, 0.0, 0.0, B0
    return f


outcome = {}
for q in (+1, -1):
    for B0 in (+1.0, -1.0):
        for s in (+1, -1):
            r_m, z_m, esc = boris(q, [1.0, 0.0, 1.0], [0.0, 0.1, 0.0], penning(B0, s), 0.05, 40000)
            outcome[(q, B0, s)] = not esc
            print(f"   {q:+d} {B0:+5.1f} {s:+3d} {r_m:10.2f} {z_m:9.2f}   {'confined' if not esc else 'ESCAPES'}")
check("each sign of charge is confined only when s matches its sign",
      all(outcome[(q, B0, s)] == (q == s) for (q, B0, s) in outcome))
check("flipping B never changes which charge is confined",
      all(outcome[(q, 1.0, s)] == outcome[(q, -1.0, s)] for q in (1, -1) for s in (1, -1)))

r_m, _, esc = boris(+1, [1.0, 0.0, 1.0], [0.0, 0.1, 0.0], penning(0.1, +1), 0.05, 40000)
crit = math.sqrt(2) * w_z
print(f"   weak field test: B = 0.1 is below sqrt(2) w_z = {crit:.3f}, escapes = {esc}")
check("below the Penning stability limit w_c^2 > 2 w_z^2 even the right charge escapes", esc)
_, _, esc2 = boris(+1, [1.0, 0.0, 1.0], [0.0, 0.1, 0.0], penning(0.2, +1), 0.05, 40000)
print(f"   B = 0.2 is above it, escapes = {esc2}")
check("just above the limit the correct charge is held", not esc2)

*_, track = boris(+1, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], penning(1.0, +1), 0.05, 20000, record=True)
phase = np.unwrap(np.arctan2(track[:, 1], track[:, 0]))
*_, track2 = boris(+1, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0], penning(-1.0, +1), 0.05, 20000, record=True)
phase2 = np.unwrap(np.arctan2(track2[:, 1], track2[:, 0]))
print(f"   B=+1 net rotation {phase[-1] - phase[0]:+.2f} rad, B=-1 {phase2[-1] - phase2[0]:+.2f} rad:"
      " polarity only reverses the sense of rotation")
check("reversing B reverses the magnetron rotation", phase[-1] * phase2[-1] < 0)

print("\n2. Magnetic bottle: B alone, mirror ratio 5")
L = 10.0
B0 = 1.0


def bottle(sign):
    def f(x, y, z):
        Bz = sign * B0 * (1 + 4 * (z / L) ** 2)
        Br_over_r = -sign * B0 * 4 * z / L**2
        return 0.0, 0.0, 0.0, Br_over_r * x, Br_over_r * y, Bz
    return f


loss_cone = math.degrees(math.asin(math.sqrt(1 / 5)))
print(f"   mirror ratio 5, predicted loss-cone half-angle {loss_cone:.1f} deg from the axis")
print("   pitch[deg]  +,B+    +,B-    -,B+    -,B-")
agreement = True
for pitch in (20, 35, 60, 80):
    row = []
    v = 0.5
    vperp, vz = v * math.sin(math.radians(pitch)), v * math.cos(math.radians(pitch))
    for q, sign in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        start_v = [0.0, vperp, vz]
        # Start the gyro-orbit near the axis at the field minimum.
        r_m, z_m, esc = boris(q, [0.5, 0.0, 0.0], start_v, bottle(sign), 0.02, 20000, escape=L)
        row.append(not esc)
    print(f"   {pitch:8d}   " + "  ".join("held " if t else "lost " for t in row))
    agreement &= len(set(row)) == 1
    expected = pitch > loss_cone
    check(f"pitch {pitch} deg behaves as the loss cone predicts", row[0] == expected)
check("a magnetic bottle treats both charges and both polarities alike", agreement)

print("\n3. Nested well: one axial potential holding opposite charges in different places")
print("   phi(z) = A[(z/a)^4 - 2 (z/a)^2]. Positive charge sits at z = +-a, negative near z = 0.")
A, a = 1.0, 2.0


def nested_force(z, q):
    return -q * A * (4 * z**3 / a**4 - 4 * z / a**2)


def run_1d(q, z0, v0, dt=0.002, steps=100000):
    z, v, lo, hi = z0, v0, z0, z0
    for _ in range(steps):
        v += 0.5 * dt * nested_force(z, q)
        z += dt * v
        v += 0.5 * dt * nested_force(z, q)
        lo, hi = min(lo, z), max(hi, z)
        if abs(z) > 50:
            return lo, hi, True
    return lo, hi, False


lo_p, hi_p, esc_p = run_1d(+1, a + 0.2, 0.0)
lo_n, hi_n, esc_n = run_1d(-1, 0.2, 0.0)
print(f"   cold positive stays in z = [{lo_p:.2f}, {hi_p:.2f}], cold negative in z = [{lo_n:.2f}, {hi_n:.2f}]")
check("both species are axially confined at once", not esc_p and not esc_n)
check("cold species do not overlap", lo_p > hi_n)
lo_h, hi_h, esc_h = run_1d(+1, a, math.sqrt(2 * (A + 0.3)))
print(f"   positive given energy above the central hill passes z = 0: range [{lo_h:.2f}, {hi_h:.2f}]")
check("an energised positive crosses through the negative species' region", lo_h < 0 < hi_h and not esc_h)

print("\n4. Crossed E and B fields (A17). Does a 'mass spectrometer style' filter sort charge?")
E_y, B_z = 0.1, 1.0
v_d = E_y / B_z


def crossed(x, y, z):
    return 0.0, E_y, 0.0, 0.0, 0.0, B_z


def e_only(x, y, z):
    return 0.0, E_y, 0.0, 0.0, 0.0, 0.0


def b_only(x, y, z):
    return 0.0, 0.0, 0.0, 0.0, 0.0, B_z


T = 50.0
dt = 0.005
n = int(T / dt)
*_, tp = boris(+1, [0, 0, 0], [v_d, 0, 0], crossed, dt, n, escape=1e9, record=True)
*_, tn = boris(-1, [0, 0, 0], [v_d, 0, 0], crossed, dt, n, escape=1e9, record=True)
sep = np.linalg.norm(tp[-1] - tn[-1])
print(f"   at the balance speed v = E/B: both go straight; separation after t={T:.0f}: {sep:.2e}")
check("at the balance speed a Wien filter does not separate e+ from e-", sep < 1e-8)

T_long = 2000.0
n_long = int(T_long / 0.02)
*_, tp2 = boris(+1, [0, 0, 0], [0.5, 0, 0], crossed, 0.02, n_long, escape=1e9, record=True)
*_, tn2 = boris(-1, [0, 0, 0], [0.5, 0, 0], crossed, 0.02, n_long, escape=1e9, record=True)
xp, xn = tp2, tn2
sep_all = np.linalg.norm(xp - xn, axis=1)
drift_p, drift_n = xp[-1, 0] / T_long, xn[-1, 0] / T_long
print(f"   off balance, v = 0.5: drift speed e+ {drift_p:.4f}, e- {drift_n:.4f}, E/B = {v_d}")
print(f"   separation stays bounded: max {sep_all.max():.2f} over t = {T_long:.0f}; gyroradius-scale"
      f" 2(v-E/B)/w_c x 2 = {4 * (0.5 - v_d):.2f}")
check("long crossed fields drift both charges the same way", abs(drift_p - v_d) < 2e-3 and abs(drift_n - v_d) < 2e-3)
check("in long crossed fields the separation never grows", sep_all.max() < 1.05 * 4 * (0.5 - v_d) + 1e-6)

*_, ta = boris(+1, [0, 0, 0], [0.5, 0, 0], e_only, dt, n, escape=1e9, record=True)
*_, tb = boris(-1, [0, 0, 0], [0.5, 0, 0], e_only, dt, n, escape=1e9, record=True)
sep_e = np.linalg.norm(ta[-1] - tb[-1])
print(f"   E field alone: separation after t={T:.0f} is {sep_e:.1f}, prediction E t^2 = {E_y * T**2:.1f}")
check("an E field alone separates the charges and keeps doing so", abs(sep_e - E_y * T**2) / (E_y * T**2) < 1e-2)
half_turn = int(math.pi / dt)
*_, tc = boris(+1, [0, 0, 0], [0.5, 0, 0], b_only, dt, half_turn, escape=1e9, record=True)
*_, td = boris(-1, [0, 0, 0], [0.5, 0, 0], b_only, dt, half_turn, escape=1e9, record=True)
y_turn = (tc[-1][1], td[-1][1])
print(f"   B field alone, half turn: y = {y_turn[0]:+.3f} for e+, {y_turn[1]:+.3f} for e- (opposite bends)")
check("a B-only spectrometer bends the two charges to opposite sides", y_turn[0] * y_turn[1] < 0)

print("\n5. Stored density, energy and annihilation time (SI)")
r_e = physical_constants["classical electron radius"][0]
sv = math.pi * r_e**2 * c
print(f"   rate coefficient sigma v = pi r_e^2 c = {sv:.3e} m^3/s")
print("    n [1/m^3]     annihilation time (no Coulomb focusing)")
for n_den in (1e15, 1e16, 1e20):
    print(f"   {n_den:9.0e} {1 / sv / n_den:18.2e} s")
B_T = 5.0
n_b = epsilon_0 * B_T**2 / (2 * m_e)
E_pair = 2 * m_e * c**2
print(f"   Brillouin limit at {B_T} T: n = {n_b:.2e} 1/m^3")
print(f"   stored rest energy if all paired at that density: {n_b * E_pair:.2e} J/m^3 "
      f"(petrol ~3e10 J/m^3, general knowledge)")
check("even the Brillouin-limit positron plasma stores far less than chemical fuel per volume",
      n_b * E_pair < 1e9)
N_lab = 1e10
print(f"   1e10 positrons (Danielson et al. storage scale): rest-energy of pairs {N_lab * E_pair:.2e} J")
for T_K in (4.0, 300.0):
    v_rel = math.sqrt(8 * k_B * T_K / (math.pi * 0.5 * m_e))
    eta_s = 2 * math.pi * alpha / (v_rel / c)
    S = eta_s / (1 - math.exp(-eta_s))
    print(f"   Coulomb (Sommerfeld) enhancement at {T_K:.0f} K relative speed: {S:.0f}x")

print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILED: ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
