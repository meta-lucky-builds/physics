"""C35-C37: replaced theories keep their accuracy where they were tested; laws compress data."""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from _harness import Checks
from common.constants import G, c, year, M_sun, M_earth

C = Checks("Effective theories: Newton inside Einstein, and why patterns matter")

C.section("Newtonian limit by Taylor expansion")
v, m, cc = sp.symbols("v m c", positive=True)
E = m * cc**2 / sp.sqrt(1 - v**2 / cc**2)
ser = sp.series(E, v, 0, 7).removeO()
print(f"  E = {sp.simplify(ser)}")
C.check("E = m c^2 + m v^2/2 + (3/8) m v^4/c^2 + ...",
        sp.simplify(ser - (m * cc**2 + m * v**2 / 2 + sp.Rational(3, 8) * m * v**4 / cc**2
                           + sp.Rational(5, 16) * m * v**6 / cc**4)) == 0)
x, t, u = sp.symbols("x t u")
gam = 1 / sp.sqrt(1 - u**2 / cc**2)
xp = gam * (x - u * t)
tp = gam * (t - u * x / cc**2)
C.check("Lorentz transform becomes the Galilean transform as c -> infinity",
        sp.limit(xp, cc, sp.oo) == x - u * t and sp.limit(tp, cc, sp.oo) == t)

C.section("GPS: relativistic rate offset, from the orbit itself")
GM = 3.986004418e14
a = 26561.75e3  # semi-major axis, m (from Ashby 2003, value recalled; 12 sidereal-hour orbit)
Re = 6378137.0
J2 = 1.0826e-3
w = 7.2921151e-5
sid = 86164.1
print(f"  period check: a = {a/1e3:.1f} km gives T = {2*np.pi*np.sqrt(a**3/GM)/3600:.3f} h")
Phi_geoid = -(GM / Re) * (1 + J2 / 2) - 0.5 * (w * Re)**2  # equatorial surface effective potential
Phi_sat = -GM / a
v2 = GM / a
grav = (Phi_sat - Phi_geoid) / c**2
vel = -v2 / (2 * c**2)
net = grav + vel
day = 86400
print(f"  gravity: {grav:+.4e} = {grav*day*1e6:+.2f} us/day; speed: {vel:+.4e} = {vel*day*1e6:+.2f} us/day")
print(f"  net rate offset {net:+.4e} = {net*day*1e6:+.2f} us/day (Physics Today, Ashby 2002: 4.4645e-10)")
print(f"  ignoring this: {net*day*c/1e3:.1f} km of ranging error per day")
C.check("computed net offset agrees with the published 4.4645e-10 within 0.2%", abs(net - 4.4645e-10) / 4.4645e-10 < 2e-3)
# Independent route: integrate an actual eccentric orbit with Newtonian gravity, accumulate proper
# time from the weak-field metric, and read the average rate and the periodic eccentricity term.
ecc = 0.01
rp = a * (1 - ecc)
vp = np.sqrt(GM * (1 + ecc) / rp)


def rhs(tt, y):
    r3 = (y[0]**2 + y[1]**2)**1.5
    r = np.sqrt(y[0]**2 + y[1]**2)
    vv = y[2]**2 + y[3]**2
    return [y[2], y[3], -GM * y[0] / r3, -GM * y[1] / r3, (-GM / r - vv / 2 - Phi_geoid) / c**2]


T = 2 * np.pi * np.sqrt(a**3 / GM)
sol = solve_ivp(rhs, [0, 10 * T], [rp, 0, 0, vp, 0], method="DOP853", rtol=1e-11, atol=[1e-3] * 4 + [1e-18],
                dense_output=True)
drift = sol.y[4, -1] / (10 * T)
print(f"  simulated orbit (e = 0.01): mean clock rate offset = {drift:+.4e}")
ts = np.linspace(0, T, 2001)
tau = sol.sol(ts)[4] - drift * ts
amp = (tau.max() - tau.min()) / 2
amp_formula = 2 * np.sqrt(GM * a) * ecc / c**2
print(f"  periodic term amplitude: simulated {amp*1e9:.1f} ns; 2 sqrt(GM a) e/c^2 = {amp_formula*1e9:.1f} ns")
C.check("simulated mean rate agrees with the closed form to 1%", abs(drift - net) / net < 0.01)
C.check("simulated eccentricity term matches formula to 2%", abs(amp - amp_formula) / amp_formula < 0.02)

C.section("Mercury's perihelion: integrate the relativistic orbit, not the formula")
a_m = 57.909e9  # Park et al. and NASA fact sheet give the same rounded values
e_m = 0.2056
T_m = 87.969  # days
GMs = G * M_sun


def perihelion_shift(eps_factor, n_orbits=20):
    """u'' + u = GM/L^2 + eps_factor * 3 GM u^2/c^2 with u = 1/r, phi as the clock."""
    Lsq = GMs * a_m * (1 - e_m**2)
    u0 = (1 + e_m) / (a_m * (1 - e_m**2))
    # Scale u by 1/a_m for numerical conditioning.
    k1 = GMs * a_m / Lsq
    k2 = eps_factor * 3 * GMs / (c**2 * a_m)

    def f(phi, y):
        return [y[1], k1 + k2 * y[0]**2 - y[0]]

    def ev(phi, y):
        return y[1]
    ev.direction = -1
    s = solve_ivp(f, [0, 2 * np.pi * (n_orbits + 0.2)], [u0 * a_m, 0.0], method="DOP853", events=ev, rtol=1e-13,
                  atol=1e-15)
    ph = s.t_events[0]
    ph = ph[ph > 1.0]
    return (ph[-1] - 2 * np.pi * len(ph)) / len(ph)


shift_0 = perihelion_shift(0.0)
shift_gr = perihelion_shift(1.0)
form = 6 * np.pi * GMs / (c**2 * a_m * (1 - e_m**2))
arcsec = 180 / np.pi * 3600
per_cy = 36525 / T_m
print(f"  Newton (no 3GMu^2/c^2 term): shift per orbit = {shift_0:.2e} rad")
print(f"  relativistic orbit integrated: {shift_gr:.4e} rad/orbit = {shift_gr*arcsec*per_cy:.2f} arcsec/century")
print(f"  formula 6 pi G M/(c^2 a (1-e^2)) = {form:.4e} rad/orbit = {form*arcsec*per_cy:.2f} arcsec/century (Will 2014: 42.98)")
C.check("Newtonian control gives no shift", abs(shift_0) < 1e-10)
C.check("integrated value reproduces 42.98 arcsec/century within 0.1%", abs(shift_gr * arcsec * per_cy - 42.98) < 0.05)
total = 575.31
print(f"  measured total {total} arcsec/century (Park et al. 2017): Newtonian perturbations keep {total-42.98:.2f}, or {100*(total-42.98)/total:.1f}%")

C.section("What the change from Newton to Einstein actually touched")
rows = [
    ("Car at 30 m/s, kinetic energy", 0.75 * (30 / c)**2),
    ("Airliner at 250 m/s, kinetic energy", 0.75 * (250 / c)**2),
    ("Earth orbit, v^2/c^2", (29.78e3 / c)**2),
    ("Sun's surface potential GM/(R c^2)", G * M_sun / (6.957e8 * c**2)),
    ("Earth's surface potential GM/(R c^2)", G * M_earth / (6.371e6 * c**2)),
    ("Hydrogen electron, (alpha)^2", (1 / 137.036)**2),
    ("GPS clocks, net offset", net),
    ("Proton at the LHC, 6.8 TeV, gamma - 1", 6.8e12 / 938.272e6 - 1),
]
print(f"  {'quantity':42s} {'relativistic correction':>24s}")
small = 0
for name, val in rows:
    print(f"  {name:42s} {val:24.2e}")
    small += val < 1e-3
print(f"  {small}/{len(rows)} of these hand-picked rows change by under 0.1%. This sample shows the mechanism;")
print("  it is not a census of physics and gives no percentage for 'all of physics'.")
C.check("the same replacement is negligible for slow, weak-field systems and huge at the LHC",
        rows[0][1] < 1e-14 and rows[-1][1] > 1e3)

C.section("A law versus a lookup table: predicting data held back")
names = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]
dist = np.array([57.9, 108.2, 149.6, 228.0, 778.5, 1432.0, 2867.0, 4515.0])  # 1e6 km, NASA fact sheet
per = np.array([88.0, 224.7, 365.2, 687.0, 4331, 10747, 30589, 59800])  # days, NASA fact sheet
train = np.arange(4)
test = np.arange(4, 8)
slope, inter = np.polyfit(np.log(dist[train]), np.log(per[train]), 1)
pred = np.exp(inter + slope * np.log(dist[test]))
err = np.abs(pred - per[test]) / per[test]
print(f"  fit on the 4 inner planets only: exponent = {slope:.4f} (Kepler's 3/2)")
for n_, p_, o_, e_ in zip([names[i] for i in test], pred, per[test], err):
    print(f"    {n_:8s}: predicted {p_:9.0f} d, observed {o_:9.0f} d, error {100*e_:.2f}%")
lookup_mean = np.exp(np.mean(np.log(per[train])))
err_lookup = np.abs(lookup_mean - per[test]) / per[test]
print(f"  'every planet has its own law' has no prediction; using the inner mean gives errors {100*err_lookup.min():.0f}-{100*err_lookup.max():.0f}%")
C.check("a law fit on inner planets predicts the outer ones to under 2%", err.max() < 0.02)
C.check("the fitted exponent is 3/2 within 1%", abs(slope - 1.5) < 0.015)
# Shuffled control: a pattern must not appear in data that has none.
rng = np.random.default_rng(0)
errs = []
for _ in range(500):
    shuffled = rng.permutation(per)
    sl, it = np.polyfit(np.log(dist[train]), np.log(shuffled[train]), 1)
    errs.append(np.max(np.abs(np.exp(it + sl * np.log(dist[test])) - shuffled[test]) / shuffled[test]))
errs = np.array(errs)
print(f"  control, periods shuffled among planets: median worst error {np.median(errs)*100:.0f}%, "
      f"fraction beating the real law {np.mean(errs < err.max()):.4f}")
C.check("control: shuffled data almost never predicts this well", np.mean(errs < err.max()) < 0.01)
C.finish()
