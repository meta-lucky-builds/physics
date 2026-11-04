"""F06: pen in soil, tower overturning, temperature at foundation depth, heavy ball at the base.

The stability of a body on a spherical base is tested by integrating its rolling equations,
generated symbolically from the Lagrangian, then bisecting for the stability boundary.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
fails = []
g = 9.81


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        fails.append(name)


# 1. Short pile rotating about its toe. Passive pressure per unit width rises with depth.
gamma = 16e3   # N/m^3, soil unit weight, assumed
phi = np.radians(30)  # friction angle, assumed
Kp = (1 + np.sin(phi)) / (1 - np.sin(phi))


def capacity(L, d, factor=3.0):
    """Moment about the pile toe with Broms' 3 Kp gamma z d resistance."""
    val, _ = quad(lambda z: factor * Kp * gamma * z * d * (L - z), 0, L)
    return val


d_pen = 0.01
print(f"Kp = {Kp:.2f}. Pen d = {d_pen * 100:.0f} cm in soil, capacity about the toe:")
for L in (0.005, 0.01, 0.02, 0.04):
    print(f"    depth {L * 100:4.1f} cm: {capacity(L, d_pen) * 1e3:7.3f} mN m")
check("capacity scales as L^3 (doubling depth gives 8x)",
      abs(capacity(0.04, d_pen) / capacity(0.02, d_pen) - 8) < 1e-6)
check("capacity equals (1/2) Kp gamma d L^3", abs(capacity(0.04, d_pen) / (0.5 * Kp * gamma * d_pen * 0.04 ** 3) - 1) < 1e-6)

# Pen load: its weight tilted by 10 deg plus a light finger touch at the top.
m_pen, h_pen = 0.010, 0.14
for depth in (0.005, 0.02, 0.05):
    arm = (h_pen - depth) / 2 + depth
    weight_moment = m_pen * g * arm * np.sin(np.radians(10))
    touch = 0.05 * (h_pen - depth + depth)  # 0.05 N at the top, about the toe
    M_load = weight_moment + touch
    print(f"    depth {depth * 100:.1f} cm: load {M_load * 1e3:6.2f} mN m vs capacity {capacity(depth, d_pen) * 1e3:7.2f} mN m "
          f"-> {'stands' if capacity(depth, d_pen) > M_load else 'falls'}")
L_needed = (m_pen * g * 0.07 * np.sin(np.radians(10)) + 0.05 * 0.14) / (0.5 * Kp * gamma * d_pen)
L_needed = L_needed ** (1 / 3)
print(f"    depth needed for that load: about {L_needed * 100:.1f} cm")
check("shallow pen falls, deep pen stands, in the model",
      capacity(0.005, d_pen) < 0.05 * 0.14 and capacity(0.05, d_pen) > 0.05 * 0.14)

# 2. Scale: geometric enlargement by s.
print("\nGeometric scaling by s: wind moment ~ s^3, weight moment ~ s^4, bearing pressure ~ s.")

# 3. Tower overturning. Burj Khalifa order of magnitude.
H = 828.0
W = 500e6 * g  # 500 000 t, official figure via fact sheet
half_width = 35.0  # m, estimated from raft/pile layout (see README)
q = 0.5 * 1.2 * 50.0 ** 2  # 50 m/s wind
Cd = 1.3
width = lambda z: 60.0 - 50.0 * z / H  # assumed taper from 60 m to 10 m
F_wind, _ = quad(lambda z: q * Cd * width(z), 0, H)
M_wind, _ = quad(lambda z: q * Cd * width(z) * z, 0, H)
M_rest = W * half_width
print(f"\nTower: W = {W:.2e} N. 50 m/s wind: F = {F_wind:.2e} N, M = {M_wind:.2e} N m. "
      f"Restoring W b = {M_rest:.2e} N m. Factor {M_rest / M_wind:.1f} (Poulos & Bunce quote about 5).")
check("gravity moment beats wind overturning by the same order as the published factor",
      2 < M_rest / M_wind < 25)
for area in (4000, 8000):
    print(f"    mean bearing pressure on {area} m^2: {W / area / 1e6:.2f} MPa")

# 4. Temperature at foundation depth.
for label, depth in (("Burj Khalifa piles", 50), ("Petronas foundations", 120), ("Kola borehole", 12262)):
    print(f"    {label}: {depth} m -> {depth * 25 / 1000:.1f} to {depth * 30 / 1000:.1f} K warmer at 25-30 K/km")
check("deepest building foundation warms by under 4 K at 30 K/km", 120 * 30 / 1000 < 4)
print(f"    Petronas depth is {100 * 120 / 12262:.1f}% of the Kola depth")

# 5. Roly-poly. Rolling sphere of radius R, COM offset d below the sphere's centre.
th, thd, R_s, d_s, m_s, Ic, g_s, t = sp.symbols("theta thetadot R d m I_c g t")
theta = sp.Function("theta")(t)
xc = R_s * theta - d_s * sp.sin(theta)
yc = R_s - d_s * sp.cos(theta)
Lag = sp.Rational(1, 2) * m_s * (sp.diff(xc, t) ** 2 + sp.diff(yc, t) ** 2) \
    + sp.Rational(1, 2) * Ic * sp.diff(theta, t) ** 2 - m_s * g_s * yc
eq = sp.diff(sp.diff(Lag, sp.diff(theta, t)), t) - sp.diff(Lag, theta)
acc_sym = sp.solve(eq, sp.diff(theta, t, 2))[0].subs({sp.diff(theta, t): thd, theta: th})
acc_num = sp.lambdify((th, thd, R_s, d_s, m_s, Ic, g_s), acc_sym, "numpy")


def stays_upright(d, R, m=1.0, k=0.5, theta0=0.02, T=None):
    """Integrate from a small tilt. True if |theta| stays under 10 x the start for T."""
    I_c = m * (k * R) ** 2
    T = T or 40 * np.sqrt(R / g)

    def rhs(t_, y):
        return [y[1], acc_num(y[0], y[1], R, d, m, I_c, g)]

    def blown(t_, y):
        return abs(y[0]) - 10 * theta0
    blown.terminal = True
    sol = solve_ivp(rhs, (0, T), [theta0, 0.0], events=blown, rtol=1e-8, atol=1e-10)
    return sol.t_events[0].size == 0


def boundary(R):
    lo, hi = -0.3 * R, 0.3 * R
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if stays_upright(mid, R):
            hi = mid
        else:
            lo = mid
    return hi


print("\nBody on a sphere: stability boundary by bisection (simulated). h_c = R - d is the COM height.")
for R in (1.0, 50.0):
    b = boundary(R)
    print(f"    R = {R:5.1f} m: boundary at d = {b / R:+.4f} R, that is h_c = {(R - b) / R:.4f} R")
    check(f"upright is stable only for h_c < R (R = {R})", abs(b / R) < 0.01)

# Uniform ball with any tower on top: COM above the sphere's centre, so d < 0.
print("\nUniform ball of radius 100 m, ball mass 100 x tower, tower H = 828 m, on the ball:")
R_b, H_t, mt = 100.0, 828.0, 1.0
mb = 100.0
d_uni = -(mt * (R_b + H_t / 2)) / (mb + mt)
print(f"    COM is {-d_uni:.2f} m above the ball's centre")
k_sys = np.sqrt(0.4 * mb * R_b ** 2 / (mb + mt)) / R_b
check("a uniform ball with a tower on top falls over in the simulation",
      not stays_upright(d_uni, R_b, k=k_sys))
check("a bare uniform ball (d = 0) has zero angular acceleration at any tilt",
      all(abs(acc_num(a, 0.0, R_b, 0.0, 1.0, 0.4 * R_b ** 2, g)) < 1e-12 for a in (0.1, 0.5, 1.5)))

# Bottom-heavy ball: only the lower half filled, density rho. m delta = rho pi R^4 / 4 with delta=3R/8.
M_t = 500e6
M_w = M_wind
print("\nBottom-heavy hemisphere under a 500 000 t, 828 m tower, top half empty:")
for rho, mat in ((7850.0, "steel"), (2400.0, "concrete")):
    def excess(R):
        return rho * np.pi * R ** 4 / 4 - M_t * (R + H / 2)
    from scipy.optimize import brentq
    R_min = brentq(excess, 1, 2000)
    # Tilt limit H/500 under the 50 m/s wind moment: S sin(theta) = M_w
    theta_lim = 1 / 500
    def excess2(R):
        return (rho * np.pi * R ** 4 / 4 - M_t * (R + H / 2)) * g * theta_lim - M_w
    R_lim = brentq(excess2, 1, 2000)
    mass_lim = rho * (2 / 3) * np.pi * R_lim ** 3
    wt = (mass_lim + M_t) * g
    print(f"    {mat}: neutral at R = {R_min:.0f} m ({rho * (2 / 3) * np.pi * R_min ** 3 / M_t:.0f} x tower mass); "
          f"tilt H/500 at R = {R_lim:.0f} m, ball {mass_lim:.2e} kg = {mass_lim / M_t:.0f} x tower")
    for q_b in (1e6, 10e6):
        area = wt / q_b
        a_cap = np.sqrt(area / np.pi)
        if a_cap < R_lim:
            sink = R_lim - np.sqrt(R_lim ** 2 - a_cap ** 2)
            print(f"        ground bearing {q_b / 1e6:.0f} MPa: contact radius {a_cap:.0f} m, sinks {sink:.0f} m")
        else:
            print(f"        ground bearing {q_b / 1e6:.0f} MPa: contact radius {a_cap:.0f} m > R, ball buries past its equator")
        check(f"{mat}: ball cannot rest on a point at {q_b / 1e6:.0f} MPa", a_cap > 0.2 * R_lim)

# Flat base has zero tilt until W b; sphere tilts under any moment.
print(f"\nFlat rigid base: no tilt until M = W b = {M_rest:.2e} N m. Sphere with restoring S: tilt = M / S for any M > 0.")

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All F06 checks passed.")
