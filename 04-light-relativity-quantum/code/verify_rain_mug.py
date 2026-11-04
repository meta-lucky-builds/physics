"""L13: how much rain does a mug catch on a trip around the Earth at nearly c?

The count of drops caught is computed in two separate ways by Monte Carlo:
- Earth frame: drops at rest density n fall at u, mug length-contracted along its motion.
- Ship frame: drops have density gamma n, velocity from velocity addition, mug at rest.
Neither calculation uses the other's events, and they must agree, because a count of
events does not depend on frame. Then the real numbers.
"""
import numpy as np

from checklib import K, check, finish, section

rng = np.random.default_rng(17)
c = K.c


def gamma(b):
    return 1 / np.sqrt(1 - b**2)


def forward_earth(n, a, Ltrip, b, u):
    """Disc faces the direction of motion (radius a, unaffected transversely)."""
    T = Ltrip / b
    y_ext = u * T + a
    vol = Ltrip * (2 * y_ext + 2 * a) * (2 * a)
    Np = rng.poisson(n * vol)
    x0 = rng.uniform(0, Ltrip, Np)
    y0 = rng.uniform(-2 * a, 2 * y_ext, Np)
    z0 = rng.uniform(-a, a, Np)
    t_hit = x0 / b
    y_hit = y0 - u * t_hit
    return np.count_nonzero(y_hit**2 + z0**2 < a**2)


def forward_ship(n, a, Ltrip, b, u):
    """Same in the ship frame: gamma n, drops move at -b along x and -u/gamma along y."""
    g = gamma(b)
    Lc = Ltrip / g
    n_p = g * n
    uy = u / g
    T = Lc / b
    y_ext = uy * T + a
    vol = Lc * (2 * y_ext + 2 * a) * (2 * a)
    Np = rng.poisson(n_p * vol)
    x0 = rng.uniform(0, Lc, Np)
    y0 = rng.uniform(-2 * a, 2 * y_ext, Np)
    z0 = rng.uniform(-a, a, Np)
    t_hit = x0 / b
    y_hit = y0 - uy * t_hit
    return np.count_nonzero(y_hit**2 + z0**2 < a**2)


def up_earth(n, a, Ltrip, b, u):
    """Disc faces straight up. In the Earth frame it is an ellipse, a/gamma along x."""
    g = gamma(b)
    T = Ltrip / b
    vol = (Ltrip + 2 * a) * (u * T) * (2 * a)
    Np = rng.poisson(n * vol)
    x0 = rng.uniform(-a, Ltrip + a, Np)
    y0 = rng.uniform(0, u * T, Np)
    z0 = rng.uniform(-a, a, Np)
    t = y0 / u
    xrel = x0 - b * t
    return np.count_nonzero((xrel * g) ** 2 + z0**2 < a**2)


def up_ship(n, a, Ltrip, b, u):
    """Disc at rest, radius a. Drops at gamma n, velocity (-b, -u/gamma)."""
    g = gamma(b)
    n_p = g * n
    uy = u / g
    T = Ltrip / g / b
    Lx = Ltrip / g + 2 * a
    vol = Lx * (uy * T) * (2 * a)
    Np = rng.poisson(n_p * vol)
    x0 = rng.uniform(-a, Lx - a, Np)
    y0 = rng.uniform(0, uy * T, Np)
    z0 = rng.uniform(-a, a, Np)
    t = y0 / uy
    xrel = x0 - b * t
    return np.count_nonzero(xrel**2 + z0**2 < a**2)


section("Monte Carlo counts in both frames (test parameters, c = 1)")
a = 1.0
Ltrip = 500.0
for b in [0.6, 0.99, 0.9999]:
    g = gamma(b)
    u = 0.05
    n = 5.0
    Ne, Ns = forward_earth(n, a, Ltrip, b, u), forward_ship(n, a, Ltrip, b, u)
    pred = n * np.pi * a**2 * Ltrip
    tol = 5 * np.sqrt(pred)
    check(f"bucket facing forward, gamma = {g:.1f}: count = n A L in both frames",
          abs(Ne - pred) < tol and abs(Ns - pred) < tol, f"Earth {Ne}, ship {Ns}, n A L = {pred:.0f}")
for b in [0.6, 0.99]:
    g = gamma(b)
    u = 0.3
    n = 20.0
    Ne, Ns = up_earth(n, a, Ltrip, b, u), up_ship(n, a, Ltrip, b, u)
    pred = n * u * np.pi * a**2 * Ltrip / b / g
    tol = 5 * np.sqrt(pred)
    check(f"mug facing up, gamma = {g:.2f}: count = n u A T / gamma in both frames",
          abs(Ne - pred) < tol and abs(Ns - pred) < tol, f"Earth {Ne}, ship {Ns}, formula {pred:.0f}")

section("Real rain, real mug")
rain_rate = 10e-3 / 3600  # m/s, 10 mm/h, moderate-to-heavy rain
u_drop = 6.5  # m/s, 2 mm drop (general knowledge; Gunn & Kinzer 1949 data)
d_drop = 2e-3
V_drop = np.pi / 6 * d_drop**3
m_drop = 1000 * V_drop
frac = rain_rate / u_drop
n_drops = frac / V_drop
D_mug = 0.08
A = np.pi * D_mug**2 / 4
L_trip = 2 * np.pi * K.R_earth
T = L_trip / c
print(f"  Rain: {rain_rate * 3.6e6:.0f} mm/h, volume fraction of water in air {frac:.2e}, {n_drops:.0f} drops/m^3")
print(f"  Mug mouth {A * 1e4:.1f} cm^2, trip length {L_trip / 1e3:.0f} km")
print(f"  Trip time in Earth frame (v -> c): {T:.4f} s = {T * 1e9:.3e} ns")
check("trip time is about 0.13 s", 0.13 < T < 0.14)
print(f"  On board, 'less than a nanosecond' needs gamma > {T / 1e-9:.2e}")
V_sweep = A * L_trip
forward = frac * V_sweep
print(f"  Bucket facing forward: swept volume {V_sweep:.3e} m^3, water {forward * 1e3:.1f} L, {forward / V_drop:.2e} drops")
up0 = rain_rate * A * T
print(f"  Mug facing up, no contraction (gamma = 1 formula): {up0 * 1e9:.2f} uL = {up0 / V_drop:.2f} drops")
print("  gamma     onboard time      facing forward    facing up")
for g in [1.2, 22.4, 1000, 1e6, 1e9]:
    print(f"  {g:<8.3g}  {T / g:.3e} s     {forward * 1e3:6.1f} L          {up0 / g * 1e9:.3e} uL")
print(f"  Standing still in the same rain for 0.134 s catches {rain_rate * A * T * 1e9:.2f} uL, "
      "the gamma = 1 value of the facing-up formula.")
check("forward catch exceeds the facing-up catch by more than 1e7", forward / up0 > 1e7, f"ratio {forward / up0:.1e}")

section("Why this can't actually happen: energy and stopping")
rho_air = 1.2
m_air = rho_air * A * L_trip
m_water = forward * 1000
m_mug = 0.3
print(f"  Column along the path: air {m_air:.2e} kg, water {m_water:.0f} kg, mug {m_mug} kg")
print("  This assumes sea-level air and rain along the entire path; real rain covers far less of it.")
kt = 4.184e12
print("  gamma    one drop KE         air column KE")
for g in [22.4, 1000]:
    KE_drop = (g - 1) * m_drop * c**2
    KE_air = (g - 1) * m_air * c**2
    print(f"  {g:<7.4g}  {KE_drop:.2e} J = {KE_drop / kt:.2g} kt  {KE_air:.2e} J")
KE_22 = (22.4 - 1) * m_drop * c**2
check("one drop at gamma 22 carries about a kiloton of energy", 1 < KE_22 / kt < 5, f"{KE_22 / kt:.1f} kt TNT")
per_m = rho_air * A * (22.4 - 1) * c**2
print(f"  In the mug's frame each metre of air brings {per_m:.1e} J at gamma 22.4; vaporising it needs about 3e6 J")
print("  If the mug simply keeps the matter it hits, v/c ~ gamma m/(gamma m + swept mass):")
for g in [22.4, 1000, 1e6]:
    s_half = g * m_mug / (rho_air * A)
    print(f"  gamma {g:<8.3g}: speed halves after {s_half / 1e3:.3g} km of air; "
          f"{'trip not completed' if s_half < L_trip else 'mug could finish the distance in principle'}")
gam_need = m_air / m_mug
print(f"  Needed gamma to still be above c/2 after the trip: {gam_need:.1e}")
g_c = c**2 / K.R_earth
print(f"  Following the curved path takes Earth-frame acceleration v^2/R = {g_c:.2e} m/s^2 = {g_c / K.g_0:.1e} g,")
print(f"  and proper acceleration gamma^2 times that, e.g. {22.4**2 * g_c / K.g_0:.1e} g at gamma 22.4")

finish()
