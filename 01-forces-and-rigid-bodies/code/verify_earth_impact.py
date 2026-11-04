"""F09: spinning ball hit at its edge, then Earth hit by asteroids.

(1) A disc of particles joined by springs is spun, then hit by a projectile that sticks. Its new
spin is measured from how fast particles turn over time, and compared with angular momentum
bookkeeping. (2) Earth: real numbers for Chicxulub-scale impacts, as a vector calculation.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common.constants import G  # noqa: E402
import prem  # noqa: E402

fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        fails.append(name)


# 1. Particle disc: hex lattice, nearest-neighbour springs with internal dampers.
SPACING = 0.15


def build_disc(R=1.0, spacing=SPACING):
    pts = []
    for j in range(-14, 15):
        for i in range(-14, 15):
            x = spacing * (i + 0.5 * (j % 2))
            y = spacing * j * np.sqrt(3) / 2
            if x * x + y * y <= R * R:
                pts.append((x, y))
    pts = np.array(pts)
    d = np.linalg.norm(pts[:, None] - pts[None, :], axis=2)
    ii, jj = np.where((d > 0) & (d < 1.1 * spacing))
    keep = ii < jj
    return pts, ii[keep], jj[keep]


def simulate(omega0, v_proj, b, m_proj=1.0, k=4e4, c=30.0, T=2.0, dt=4e-5):
    """Projectile moving in +x at impact parameter b (y offset). It sticks to the three nearest particles."""
    pos, ii, jj = build_disc()
    n = len(pos)
    m_each = 10.0 / n
    mass = np.full(n + 1, m_each)
    mass[-1] = m_proj
    rest = np.linalg.norm(pos[ii] - pos[jj], axis=1)
    vel = omega0 * np.column_stack((-pos[:, 1], pos[:, 0]))
    # Start the projectile just left of the disc on its line.
    xc = -np.sqrt(max(1.0 - b ** 2, 0.0)) - 0.02
    pos = np.vstack((pos, [xc, b]))
    vel = np.vstack((vel, [v_proj, 0.0]))
    stuck = -1
    pairs_i, pairs_j = ii, jj

    def accel(p, v):
        d = p[pairs_j] - p[pairs_i]
        Ln = np.linalg.norm(d, axis=1)
        u = d / Ln[:, None]
        vrel = ((v[pairs_j] - v[pairs_i]) * u).sum(1)
        f = (k * (Ln - rest_now) + c * vrel)[:, None] * u
        nn = len(p)
        ax = np.bincount(pairs_i, f[:, 0], nn) - np.bincount(pairs_j, f[:, 0], nn)
        ay = np.bincount(pairs_i, f[:, 1], nn) - np.bincount(pairs_j, f[:, 1], nn)
        return np.column_stack((ax, ay)) / mass[:, None]

    rest_now = rest
    a = accel(pos, vel)
    com0 = (mass[:, None] * pos).sum(0) / mass.sum()
    v0 = (mass[:, None] * vel).sum(0) / mass.sum()
    rel, vr = pos - com0, vel - v0
    L0 = (mass * (rel[:, 0] * vr[:, 1] - rel[:, 1] * vr[:, 0])).sum()
    steps = int(T / dt)
    track, inertias, ref = [], [], None
    for s in range(steps):
        vel += 0.5 * dt * a
        pos += dt * vel
        if stuck < 0:
            dd = np.linalg.norm(pos[:-1] - pos[-1], axis=1)
            if dd.min() < SPACING:
                # Stick: bond to the three nearest particles so the projectile cannot swing freely.
                near = np.argsort(dd)[:3]
                pairs_i = np.concatenate((ii, near))
                pairs_j = np.concatenate((jj, np.full(near.size, n)))
                rest_now = np.concatenate((rest, np.maximum(dd[near], 0.05)))
                stuck = s
        a = accel(pos, vel)
        vel += 0.5 * dt * a
        if s % 250 == 0 and s * dt > T * 0.5:
            # Whole-body orientation: least-squares rotation from a reference snapshot.
            com = (mass[:, None] * pos).sum(0) / mass.sum()
            rel_now = pos - com
            if ref is None:
                ref = rel_now.copy()
            cross = (mass * (ref[:, 0] * rel_now[:, 1] - ref[:, 1] * rel_now[:, 0])).sum()
            dot = (mass * (ref * rel_now).sum(1)).sum()
            track.append((s * dt, np.arctan2(cross, dot)))
            inertias.append((mass * (rel_now ** 2).sum(1)).sum())
    tt, ang = np.array(track).T
    rate = np.polyfit(tt, np.unwrap(ang), 1)[0]
    com = (mass[:, None] * pos).sum(0) / mass.sum()
    v_f = (mass[:, None] * vel).sum(0) / mass.sum()
    rel, vr = pos - com, vel - v_f
    L_f = (mass * (rel[:, 0] * vr[:, 1] - rel[:, 1] * vr[:, 0])).sum()
    I_f = float(np.mean(inertias))
    return rate, L0, L_f, I_f, np.linalg.norm(v_f - v0)


print("Spinning particle disc (10 kg, radius 1 m, 2 rad/s anticlockwise), hit by a 1 kg sticking projectile at 4 m/s")
print("Positive offset b means the projectile passes above the centre (y > 0), moving +x, a clockwise push.")
print(f"{'b (m)':>7} {'tracked spin rad/s':>19} {'L_before/I_after':>17} {'|L_f-L_0|/|L_0|':>17}")
measured = {}
for b in (0.0, 0.7, -0.7):
    rate, L0, Lf, If, dv = simulate(2.0, 4.0, b)
    measured[b] = rate
    print(f"{b:7.2f} {rate:19.4f} {L0 / If:17.4f} {abs(Lf - L0) / abs(L0):17.1e}")
    check(f"spin rate from particle tracking equals L_before/I_after for b={b}", abs(rate - L0 / If) < 0.02)
    check(f"total angular momentum conserved for b={b}", abs(Lf - L0) / abs(L0) < 1e-3)
check("hit against the spin slows it", measured[0.7] < measured[0.0])
check("hit along the spin speeds it up", measured[-0.7] > measured[0.0])

# 2. Earth, real numbers.
P = prem.profile()
M = 5.9722e24
Re = prem.R_MEAN
I_e = 0.3308 * M * Re ** 2     # NASA fact sheet factor
T_sid = 23.9345 * 3600
om = 2 * np.pi / T_sid
L_e = I_e * om
E_rot = 0.5 * I_e * om ** 2
print(f"\nEarth: I = {I_e:.3e} kg m^2, omega = {om:.4e} rad/s, L = {L_e:.3e} kg m^2/s, rotational KE = {E_rot:.2e} J")
check("PREM moment of inertia agrees with the NASA factor within 0.2%", abs(P["I"] / I_e - 1) < 2e-3)
i_cmb = np.argmin(np.abs(P["r"] - 3480e3))
frac_core = P["I_enc"][i_cmb] / P["I"]
print(f"    core holds {100 * frac_core:.1f}% of the moment of inertia (PREM), mantle+crust {100 * (1 - frac_core):.1f}%")


def outcome(m, v, lat_deg, lon_deg, angle_deg, heading_deg):
    """Stick mass m at the surface. heading: compass direction of horizontal motion (90 = east)."""
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    east = np.array([-np.sin(lon), np.cos(lon), 0.0])
    north = np.cross(up, east)
    h = np.radians(heading_deg)
    a = np.radians(angle_deg)
    vel = v * (np.cos(a) * (np.sin(h) * east + np.cos(h) * north) - np.sin(a) * up)
    r = Re * up
    L_imp = m * np.cross(r, vel)
    L_tot = np.array([0.0, 0.0, L_e]) + L_imp
    # Inertia tensor after adding a point mass at r.
    It = np.eye(3) * I_e + m * (Re ** 2 * np.eye(3) - np.outer(r, r))
    om_new = np.linalg.solve(It, L_tot)
    T_new = 2 * np.pi / np.linalg.norm(om_new)
    tilt = np.arccos(om_new[2] / np.linalg.norm(om_new))
    return {"L": np.linalg.norm(L_imp), "dT_ms": (T_new - T_sid) * 1e3, "tilt_mas": np.degrees(tilt) * 3.6e6,
            "E": 0.5 * m * v ** 2, "dv": m * v / M}


def sphere_mass(d_km, rho=2630.0):
    return rho * np.pi / 6 * (d_km * 1e3) ** 3


scen = [
    ("10 km, 20 km/s (Schulte 2010 size)", sphere_mass(10), 20e3),
    ("17 km, 12 km/s (Collins 2020 sim.)", sphere_mass(17), 12e3),
    ("D-M 2014 largest, 4.6e17 kg, 15.9 km/s", 4.6e17, 15.9e3),
]
print("\nBest case for spin-up: 45 deg impact on the equator, heading east. Then the worst case, heading west:")
print(f"{'scenario':>42} {'m (kg)':>9} {'E (J)':>9} {'E (Mt TNT)':>11} {'L/L_earth':>10} {'dT east (ms)':>13} {'dT west (ms)':>13}")
res = {}
for name, m, v in scen:
    east = outcome(m, v, 0.0, 0.0, 45.0, 90.0)
    west = outcome(m, v, 0.0, 0.0, 45.0, 270.0)
    res[name] = (east, west)
    print(f"{name:>42} {m:9.2e} {east['E']:9.2e} {east['E'] / 4.184e15:11.2e} {east['L'] / L_e:10.2e} "
          f"{east['dT_ms']:13.3f} {west['dT_ms']:13.3f}")

name = scen[0][0]
east, west = res[name]
check("east-heading impact shortens the day, west-heading lengthens it", east["dT_ms"] < 0 < west["dT_ms"])
# Rigid-body estimate and vertical limit.
m0, v0 = scen[0][1], scen[0][2]
est = -T_sid * (m0 * v0 * Re * np.cos(np.radians(45))) / L_e * 1e3
print(f"    estimate dT = -T L_imp / L_earth = {est:.3f} ms for the first scenario")
check("vector result matches the simple estimate within 10%", abs(east["dT_ms"] / est - 1) < 0.1)
vert = outcome(m0, v0, 0.0, 0.0, 90.0, 90.0)
print(f"    head-on (vertical) impact: L_imp = {vert['L']:.1e}; dT = {vert['dT_ms']:+.4f} ms, from added mass alone")
check("vertical impact gives no spin torque; only a small lengthening from added inertia",
      vert["L"] < 1e-6 * east["L"] and 0 < vert["dT_ms"] < 0.2)

# Chicxulub geometry per Collins: 21.38 N, 89.52 W, from the northeast (heading SW = 225), 45-60 deg.
print("\nChicxulub site with the Collins et al. trajectory (heading SW, modern latitude):")
for name, m, v in scen[:2]:
    for ang in (45.0, 60.0):
        o = outcome(m, v, 21.38, -89.52, ang, 225.0)
        print(f"    {name}, {ang:.0f} deg: dT = {o['dT_ms']:+.3f} ms, axis shift {o['tilt_mas']:.0f} mas, Earth dv = {o['dv'] * 1e6:.1f} micron/s")
o = outcome(scen[0][1], scen[0][2], 21.38, -89.52, 45.0, 225.0)
check("this trajectory has a westward part, so it slowed the spin", o["dT_ms"] > 0)

# Energy shared. Spin energy change vs impact energy.
for name, m, v in scen[:1]:
    o = outcome(m, v, 0.0, 0.0, 45.0, 90.0)
    dE = om * o["L"]
    print(f"\nEnergy for the first scenario: impact {o['E']:.2e} J; spin energy change omega L = {dE:.2e} J "
          f"= {100 * dE / o['E']:.1f}% of it; ratio equals 2 (omega R) cos(angle) / v")
    check("spin energy share is small and equals 2 omega R cos(angle)/v",
          abs(dE / o["E"] - 2 * om * Re * np.cos(np.radians(45)) / v) < 1e-6)

# Comparisons.
print("\nFor scale: tidal slowing +1.8 ms per century (Stephenson et al. 2016);")
print("2004 Sumatra earthquake about -2.68 microseconds (NASA).")

# What would it take?
target = 1 / 24
m_need = target * L_e / (20e3 * Re * np.cos(np.radians(45)))
d_need = (6 * m_need / (np.pi * 2630)) ** (1 / 3)
print(f"To shorten the day by 1 hour at 20 km/s, 45 deg: m = {m_need:.2e} kg, diameter {d_need / 1e3:.0f} km, "
      f"energy {0.5 * m_need * 20e3 ** 2:.1e} J = {0.5 * m_need * 20e3 ** 2 / east['E']:.0f} x the first scenario")
check("that impactor is over 1000 km across and carries more energy than the Earth's spin",
      d_need > 1000e3 and 0.5 * m_need * 20e3 ** 2 > E_rot)

# Earth-Moon angular momentum.
m_moon, a_moon = 0.07346e24, 0.3844e9
L_moon = m_moon * np.sqrt(G * M * a_moon)
print(f"Moon's orbital L = {L_moon:.2e}, Earth's spin L = {L_e:.2e}, total {L_moon + L_e:.2e} kg m^2/s; "
      f"total / spin = {(L_moon + L_e) / L_e:.1f}")
L_theia = 6.4e23 * 10e3 * Re
print(f"Mars-mass impactor (6.4e23 kg, general knowledge) at 10 km/s, grazing: L = {L_theia:.1e} kg m^2/s, "
      f"{L_theia / (L_moon + L_e):.1f} x the present system")
check("a Mars-sized grazing impact carries Earth-Moon scale angular momentum",
      0.5 < L_theia / (L_moon + L_e) < 2)

# Speed bounds.
v_esc = np.sqrt(2 * G * M / Re)
print(f"\nImpact speed bounds: escape speed {v_esc / 1e3:.1f} km/s; solar-escape 42.1 + orbital 29.8 km/s head-on combined with Earth's well: "
      f"{np.sqrt((42.1e3 + 29.8e3) ** 2 + v_esc ** 2) / 1e3:.0f} km/s")

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All F09 checks passed.")
