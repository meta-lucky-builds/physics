"""F10: drilling through the Earth and falling through the hole.

Checks: geometry of a line through a ball, the shell theorem by direct integration,
PREM totals against independently measured Earth values, fall times against Klotz (2015),
the Coriolis deflection in a simulated drop, rod stress and the plumb-line miss.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common.constants import G, M_earth  # noqa: E402
import prem  # noqa: E402

fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        fails.append(name)


R = prem.R_MEAN
omega = 2 * np.pi / (23.9345 * 3600)  # NASA fact sheet sidereal period

# 1. A straight line entering a ball always leaves it, at the antipode if through the centre.
rng = np.random.default_rng(1)
worst = 0.0
for _ in range(20000):
    p = rng.normal(size=3)
    p *= R / np.linalg.norm(p)
    inner = rng.normal(size=3)
    inner *= R * rng.uniform(0, 0.999) / np.linalg.norm(inner)
    d = inner - p
    t2 = -2 * np.dot(p, d) / np.dot(d, d)  # second root of |p + t d| = R (first is t=0)
    exit_pt = p + t2 * d
    worst = max(worst, abs(np.linalg.norm(exit_pt) - R) / R)
    assert t2 > 1  # the exit lies beyond the interior point on the same straight line
check("line through interior exits on the far surface (20000 random lines)", worst < 1e-12,
      f"worst radius error fraction {worst:.1e}")
lat, lon = 20.0, 78.0
print(f"    antipode of central India ({lat}N, {lon}E) is ({lat}S, {180 - lon}W), in the eastern Pacific")

# 2. Shell theorem by direct integration over a thin shell (no use of the theorem).


def shell_force(a, R_s=1.0, m=1.0):
    """Net attraction towards the shell centre on a unit mass at radius a (G=1)."""
    def integrand(u):  # u = cos(theta)
        s2 = a * a + R_s * R_s - 2 * a * R_s * u
        return 0.5 * m * (R_s * u - a) / s2 ** 1.5
    pts = [1.0] if abs(a - R_s) < 0.2 else None
    val, _ = quad(integrand, -1, 1, limit=400, points=pts)
    return val


inside = [shell_force(a) for a in (0.05, 0.3, 0.6, 0.9)]
outside = [shell_force(a) * (-a * a) for a in (1.2, 2.0, 5.0)]
check("shell theorem: zero force inside a shell", max(abs(x) for x in inside) < 1e-8,
      f"max |F| = {max(abs(x) for x in inside):.1e} (G m = 1, R = 1)")
check("shell theorem: outside it acts like a point mass", max(abs(x - 1) for x in outside) < 1e-8,
      f"F r^2 / (G m) = {[round(x, 10) for x in outside]}")

# 3. PREM integrated, against measurements that were not used to build this code.
P = prem.profile()
factor = P["I"] / (P["M"] * R ** 2)
check("PREM mass matches NASA 5.9722e24 kg to 0.1%", abs(P["M"] / M_earth - 1) < 1e-3,
      f"M = {P['M']:.4e} kg")
check("PREM I/(M R^2) matches NASA 0.3308 to 0.1%", abs(factor / 0.3308 - 1) < 1e-3,
      f"factor = {factor:.4f}")
g_fn = prem.gravity_function()
check("PREM surface gravity near 9.8 m/s^2", abs(g_fn(R) - 9.82) < 0.05, f"g(R) = {g_fn(R):.3f} m/s^2")
i_peak = np.argmax(P["g"])
print(f"    gravity peaks at r = {P['r'][i_peak] / 1e3:.0f} km with g = {P['g'][i_peak]:.2f} m/s^2")
# Hydrostatic central pressure from dP/dr = -rho g.
r_all, rho_all, g_all = P["r"], P["rho"], P["g"]
p_centre = np.sum(0.5 * (rho_all[1:] * g_all[1:] + rho_all[:-1] * g_all[:-1]) * np.diff(r_all))
print(f"    hydrostatic central pressure = {p_centre / 1e9:.0f} GPa")

# 4. Fall times by integrating the equation of motion.


def fall_time(g_of_r):
    def rhs(t, y):
        return [y[1], -g_of_r(abs(y[0])) * np.sign(y[0])]

    def at_centre(t, y):
        return y[0]
    at_centre.terminal = True
    sol = solve_ivp(rhs, (0, 1e5), [R, 0.0], events=at_centre, rtol=1e-9, atol=1e-6, max_step=5.0)
    return 2 * sol.t_events[0][0], abs(sol.y_events[0][0][1])


g0 = G * P["M"] / R ** 2
t_uni, v_uni = fall_time(lambda r: g0 * r / R)
t_con, v_con = fall_time(lambda r: g0)
t_prem, v_prem = fall_time(g_fn)


def mmss(t):
    return f"{int(t // 60)} min {t % 60:04.1f} s"


print(f"    uniform density : {mmss(t_uni)}, peak speed {v_uni / 1e3:.2f} km/s")
print(f"    constant gravity: {mmss(t_con)}, peak speed {v_con / 1e3:.2f} km/s")
print(f"    PREM            : {mmss(t_prem)}, peak speed {v_prem / 1e3:.2f} km/s")
check("uniform-density time is pi sqrt(R/g) and near Klotz's 42 min 12 s",
      abs(t_uni - (42 * 60 + 12)) < 5 and abs(t_uni - np.pi * np.sqrt(R / g0)) < 0.1)
check("constant-gravity time near Klotz's 38 min 0 s", abs(t_con - 38 * 60) < 5)
check("PREM time near Klotz's 38 min 11 s", abs(t_prem - (38 * 60 + 11)) < 5)

# 5. Coriolis: drop down an equatorial diameter. Integrate in the inertial frame.


def drop_3d(hole_radius):
    """Return time and depth when the falling mass reaches the wall of a hole that
    corotates with the Earth, and its sideways offset at the centre if no wall existed."""
    def rhs(t, y):
        pos = y[:2]
        r = np.hypot(*pos)
        acc = -g_fn(r) * pos / r if r > 0 else np.zeros(2)
        return [y[2], y[3], acc[0], acc[1]]
    y0 = [R, 0.0, 0.0, omega * R]  # released at rest relative to the rotating surface
    sol = solve_ivp(rhs, (0, 1200), y0, rtol=1e-10, atol=1e-6, max_step=1.0, dense_output=True)
    ts = np.linspace(0, 1200, 24001)
    x, y = sol.sol(ts)[:2]
    # Unit vector along the corotating tunnel and sideways offset from its axis.
    ax, ay = np.cos(omega * ts), np.sin(omega * ts)
    along = x * ax + y * ay
    side = -x * ay + y * ax
    hit = np.argmax(np.abs(side) >= hole_radius)
    return ts, along, side, ts[hit], R - along[hit]


ts, along, side, t_hit, depth_hit = drop_3d(5.0)
print(f"    5 m radius hole at the equator: wall reached after {t_hit:.1f} s at depth {depth_hit / 1e3:.2f} km")
cubic = (3 * 5.0 / (omega * g_fn(R))) ** (1 / 3)
check("wall-hit time matches the Coriolis estimate (1/3) omega g t^3 within 3%",
      abs(t_hit / cubic - 1) < 0.03, f"estimate {cubic:.1f} s")
i_c = np.argmax(along <= 0)
print(f"    sideways offset by the time the depth reaches the Earth's centre: {abs(side[i_c]) / 1e3:.0f} km")
check("offset at the centre is hundreds of km", 100e3 < abs(side[i_c]) < 1000e3)

# 6. Rod stress: a rod reaching the centre carries rho_rod * (potential drop).
dphi = np.sum(0.5 * (g_all[1:] + g_all[:-1]) * np.diff(r_all))
rho_steel, yield_steel = 7850.0, 250e6  # typical structural steel, general knowledge
sigma = rho_steel * dphi
print(f"    gravitational potential drop surface to centre: {dphi:.3e} J/kg")
print(f"    steel rod stress at the end of a surface-to-centre rod: {sigma / 1e9:.0f} GPa "
      f"= {sigma / yield_steel:.0f} x yield; self-supporting length at surface g: "
      f"{yield_steel / (rho_steel * 9.81) / 1e3:.1f} km")
check("peak speed from energy equals the integrated peak speed",
      abs(np.sqrt(2 * dphi) / v_prem - 1) < 1e-3, f"sqrt(2 dphi) = {np.sqrt(2 * dphi) / 1e3:.2f} km/s")

# 7. A plumb line does not point at the centre on an oblate Earth.
a, b = 6378.137e3, 6356.752e3  # NASA fact sheet radii
phi_gd = np.radians(np.linspace(0, 90, 18001))
phi_gc = np.arctan((b / a) ** 2 * np.tan(np.clip(phi_gd, 0, np.radians(89.999))))
diff = np.degrees(phi_gd - phi_gc)
miss = R * np.sin(np.radians(diff.max()))
print(f"    plumb line vs geocentric direction: up to {diff.max():.3f} deg at latitude "
      f"{np.degrees(phi_gd[np.argmax(diff)]):.1f}; straight-line miss at the centre {miss / 1e3:.1f} km")
check("plumb-line angle near 0.19 deg", abs(diff.max() - 0.192) < 0.003)

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All F10 checks passed.")
