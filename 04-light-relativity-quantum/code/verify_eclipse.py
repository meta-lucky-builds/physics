"""L16: eclipse light deflection and the size of claimed gravity anomalies.

- Numerical ray bending by three independent routes: Newtonian orbit, Schwarzschild
  null geodesic, and ray tracing through an effective refractive index.
- 1919 data and reanalysis, compared with each prediction.
- Newtonian tide over the hours around an eclipse against the reported anomaly.
- Moon's shielding of the Sun's pull implied by the anomaly against the laser ranging limit.
"""
import numpy as np
from scipy.integrate import solve_ivp

from checklib import K, check, finish, section

GM = K.G * K.M_sun
c = K.c
R_sun = 6.957e8  # m, IAU nominal
ARCSEC = np.degrees(1) * 3600
k = GM / c**2
print(f"GM/c^2 for the Sun = {k:.1f} m, Sun radius = {R_sun:.4g} m, ratio {k / R_sun:.3e}")


def deflection_orbit(mu_rhs, b=R_sun):
    """Integrate u'' + u = rhs(u) from closest approach to u = 0; deflection = 2 phi_end - pi."""

    def rhs(phi, y):
        return [y[1], mu_rhs(y[0]) - y[0]]

    def hit(phi, y):
        return y[0]

    hit.terminal = True
    hit.direction = -1
    sol = solve_ivp(rhs, [0, np.pi], [1 / b, 0.0], events=hit, method="DOP853", rtol=1e-13, atol=1e-13 / b)
    return 2 * sol.t_events[0][0] - np.pi


section("Three ways of bending a ray that grazes the Sun")
b = R_sun
d_newton = deflection_orbit(lambda u: k / b**2)  # Newton: h = c b, u'' + u = GM/h^2
d_gr = deflection_orbit(lambda u: 3 * k * u**2)
print(f"  Newtonian particle at c: {d_newton * ARCSEC:.4f} arcsec")
print(f"  Schwarzschild null geodesic: {d_gr * ARCSEC:.4f} arcsec")
check("Newton gives 2GM/(c^2 b) = 0.875 arcsec", abs(d_newton / (2 * k / b) - 1) < 1e-4)
check("GR gives 4GM/(c^2 b) = 1.75 arcsec", abs(d_gr / (4 * k / b) - 1) < 1e-4)
check("GR is twice Newton", abs(d_gr / d_newton - 2) < 2e-4, f"ratio {d_gr / d_newton:.5f}")


def deflection_index(gamma_ppn, b=R_sun, span=4000.0):
    """Trace a ray through n(r) = 1 + (1 + gamma) GM/(r c^2); d(n dr/ds)/ds = grad n."""

    def n(r):
        return 1 + (1 + gamma_ppn) * k / r

    def rhs(s, y):
        x, yy, px, py = y
        r = np.hypot(x, yy)
        nn = n(r)
        dn = -(1 + gamma_ppn) * k / r**3
        return [px / nn, py / nn, dn * x, dn * yy]

    L = span * b
    n0 = n(np.hypot(L, b))
    sol = solve_ivp(rhs, [0, 2 * L], [-L, b, n0, 0.0], method="DOP853", rtol=1e-13, atol=[1e-3, 1e-3, 1e-16, 1e-16])
    return -np.arctan2(sol.y[3, -1], sol.y[2, -1])


d_time = deflection_index(0.0)
d_full = deflection_index(1.0)
print(f"  Ray through an index with gamma = 0 (time curvature only): {d_time * ARCSEC:.4f} arcsec")
print(f"  Ray through an index with gamma = 1 (time + space curvature): {d_full * ARCSEC:.4f} arcsec")
check("index route with gamma = 1 matches the geodesic", abs(d_full / d_gr - 1) < 1e-3)
check("gamma = 0 matches the Newtonian half: space curvature provides the other half", abs(d_time / d_newton - 1) < 1e-3)

section("1919 and the 1979 reanalysis (Kennefick 2009 table, Dyson et al. 1920)")
pred = {"Newton": d_newton * ARCSEC, "Einstein": d_gr * ARCSEC}
data = [
    ("Sobral 4-inch, 1919 reduction", 1.98, 0.18),
    ("Sobral 4-inch, 1979 re-reduction", 1.90, 0.11),
    ("Principe, 1919 reduction", 1.61, 0.30),
    ("Sobral astrographic, 1979 re-reduction", 1.55, 0.34),
]
print("  measurement                                value    sigma from Einstein   sigma from Newton")
zs = {}
for name, val, err in data:
    ze = (val - pred["Einstein"]) / err
    zn = (val - pred["Newton"]) / err
    zs[name] = (ze, zn)
    print(f"  {name:42s} {val:.2f}+-{err:.2f}  {ze:+.2f}                {zn:+.2f}")
print("  Sobral astrographic, 1919 reduction: 0.93, no error stated; dropped by the team for a suspected")
print("  scale change. Its reanalysis moves to 1.55, as the table above shows.")
check("Sobral 4-inch excludes Newton's half value, both reductions, over 5 sigma", abs(zs[data[0][0]][1]) > 5 and abs(zs[data[1][0]][1]) > 5)
check("Principe alone sits within 3 sigma of both predictions", abs(zs[data[2][0]][0]) < 3 and abs(zs[data[2][0]][1]) < 3)
g_vlbi = 0.99983
print(f"  Modern radio VLBI: gamma = {g_vlbi} (Shapiro et al. 2004), i.e. {(1 + g_vlbi) / 2 * pred['Einstein']:.5f} arcsec at the limb")
print(f"  Cassini: gamma - 1 = (2.1 +- 2.3)e-5 (Bertotti et al. 2003)")

section("Newtonian tide around an eclipse against the reported anomaly")
M_moon = 7.346e22
D_sun = 1.496e11
D_moon = 3.844e8
Re = K.R_earth
Om_e = 2 * np.pi / 86164.0
Om_m = 2 * np.pi / (27.32 * 86400)
Om_s = 2 * np.pi / K.year


def tide_up(GMb, Dv, rv):
    """Upward tidal acceleration at Earth-surface position rv from a body at Dv."""
    s = Dv - rv
    a = GMb * (s / np.linalg.norm(s) ** 3 - Dv / np.linalg.norm(Dv) ** 3)
    return a @ (rv / np.linalg.norm(rv))


ts = np.linspace(-3 * 3600, 3 * 3600, 721)
vals = []
vals_nomoon = []
for t in ts:
    ph = Om_e * t
    site = Re * np.array([np.cos(ph), np.sin(ph), 0.0])
    sun = D_sun * np.array([np.cos(Om_s * t), np.sin(Om_s * t), 0.0])
    moon = D_moon * np.array([np.cos(Om_m * t), np.sin(Om_m * t), 0.0])
    vals.append(tide_up(GM, sun, site) + tide_up(K.G * M_moon, moon, site))
    vals_nomoon.append(tide_up(GM, sun, site))
vals = np.array(vals)
print(f"  Sun + Moon vertical tide at a site under both at mid-eclipse: {vals.max() * 1e8:.1f} uGal, range over +-3 h {np.ptp(vals) * 1e8:.1f} uGal (1 uGal = 1e-8 m/s^2)")
wang, wang_err = 7.0e-8, 2.7e-8
print(f"  Wang et al. 2000 anomaly: ({wang * 1e8:.1f} +- {wang_err * 1e8:.1f}) uGal, {wang / wang_err:.1f} sigma, {wang / np.ptp(vals) * 100:.1f}% of that tide's 6 h range")
check("reported anomaly is a small fraction of the tide that must be removed", wang / np.ptp(vals) < 0.2)

section("Moon as a gravity shield?")
g_sun = GM / D_sun**2
frac_needed = wang / g_sun
h_limit = 1.0e-21 * 1e-4 / 1e-3  # Eckhardt 1990: 1e-21 cm^2/g converted to m^2/kg
rho_moon = 3344.0
col = rho_moon * 2 * 1.7374e6
frac_limit = h_limit * col
print(f"  Sun's pull at Earth {g_sun:.3e} m/s^2; a {wang * 1e8:.0f} uGal effect needs shielding fraction {frac_needed:.2e}")
print(f"  Lunar laser ranging limit h < {h_limit:.0e} m^2/kg with Moon column {col:.2e} kg/m^2: fraction < {frac_limit:.2e}")
print(f"  So simple Majorana-type shielding would be too weak by a factor {frac_needed / frac_limit:.1e}")
check("shielding of that size is excluded by over 1e5 in that specific model", frac_needed / frac_limit > 1e5)

finish()
