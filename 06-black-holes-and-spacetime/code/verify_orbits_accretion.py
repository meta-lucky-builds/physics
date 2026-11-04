"""Orbits and accretion: ISCO, efficiency, capture, viscous disk, Eddington limit, mergers."""
import numpy as np
from scipy import constants as sc
from scipy.integrate import solve_ivp
from scipy.linalg import solve_banded
from scipy.optimize import brentq, minimize_scalar

from bh_tools import C, M_SGRA, Report, r_g, t_g

R = Report("verify_orbits_accretion: ISCO, capture, disks, Eddington, mergers")
sigma_T = sc.physical_constants["Thomson cross section"][0]

# ------------------------------------------------------------------ ISCO from the potential
R.section("1. ISCO from the effective potential (no formula assumed; M=1)")


def V2(r, L):
    return (1 - 2 / r) * (1 + L**2 / r**2)


RG = np.geomspace(2.05, 400, 40001)


def extrema(L):
    d = np.gradient(V2(RG, L), RG)
    idx = np.where(d[:-1] * d[1:] < 0)[0]
    return RG[idx]


lo, hi = 3.0, 4.2   # bisect on "potential has a local minimum at all"
for _ in range(28):
    mid = (lo + hi) / 2
    if len(extrema(mid)) >= 2:
        hi = mid
    else:
        lo = mid
Lcrit = hi
rc = extrema(Lcrit)
risco = rc.mean()
Eisco = np.sqrt(V2(risco, Lcrit))
print(f"  smallest L with both a max and a min: L = {Lcrit:.4f} (sqrt 12 = {np.sqrt(12):.4f}); radius there r = {risco:.3f} M")
print(f"  orbital energy there E = {Eisco:.5f} (sqrt(8/9) = {np.sqrt(8 / 9):.5f}); radiated fraction = {1 - Eisco:.4f}")
R.check("ISCO at 6M", abs(risco - 6) < 0.05)
R.check("Efficiency 1 - sqrt(8/9) = 5.72%", abs((1 - Eisco) - (1 - np.sqrt(8 / 9))) < 1e-4, f"{100 * (1 - Eisco):.3f}%")


def orbit(r0, kick=1e-3, tmax=4000.0):
    """Integrate the exact radial geodesic equation from a slightly perturbed circular orbit."""
    L = np.sqrt(r0 / (1 - 3 / r0))  # circular angular momentum
    r0p = r0 * (1 + kick)

    def rhs(_, y):
        r, rd = y
        return [rd, -1 / r**2 + L**2 / r**3 - 3 * L**2 / r**4]

    ev = lambda t, y: y[0] - 2.0
    ev.terminal = True
    s = solve_ivp(rhs, [0, tmax], [r0p, 0], events=ev, rtol=1e-10, atol=1e-12)
    return s.y[0].min(), s.t[-1]


rmin7, _ = orbit(7.0, kick=-1e-2)
rmin55, t55 = orbit(5.5, kick=-1e-2)
R.check("Perturbed orbit at r=7M stays bound", rmin7 > 6.0, f"r_min = {rmin7:.3f}")
R.check("Perturbed orbit at r=5.5M plunges through the horizon", rmin55 <= 2.0 + 1e-6, f"crossed at tau = {t55:.0f} M")


def E_kerr(r, a):
    return (r**1.5 - 2 * r**0.5 + a) / (r**0.75 * np.sqrt(r**1.5 - 3 * r**0.5 + 2 * a))


print("  Kerr, prograde equatorial circular-orbit energy of Bardeen-Press-Teukolsky 1972, minimised numerically:")
print("  a/M     r_isco/M   efficiency")
eff = {}
for a in (0.0, 0.5, 0.9, 0.998, 0.99999):
    rmin_ = 1 + np.sqrt(1 - a**2) + 0.0005 if a > 0 else 3.01
    res = minimize_scalar(lambda r: E_kerr(r, a), bounds=(max(rmin_, 1.0005), 12), method="bounded",
                          options={"xatol": 1e-7})
    eff[a] = 1 - res.fun
    print(f"  {a:7.5f}  {res.x:7.3f}  {100 * eff[a]:7.2f}%")
R.check("a=0 recovers 5.72%", abs(eff[0.0] - (1 - np.sqrt(8 / 9))) < 1e-5)
R.check("Rapid spin: efficiency approaches 1-sqrt(1/3) = 42.3% from below", 0.30 < eff[0.998] < eff[0.99999] < 1 - np.sqrt(1 / 3))

# ------------------------------------------------------------------ capture
R.section("2. Black holes do not suck: capture needs a small impact parameter")


def captured(L, E=1.0, kappa=1.0, r0=2000.0):
    """Infall with energy E and angular momentum L; kappa=1 massive, 0 photon. Returns True if r<2M reached."""
    def rhs(_, y):
        r, rd = y
        return [rd, -kappa / r**2 + L**2 / r**3 - 3 * L**2 / r**4]
    rd2 = E**2 - (1 - 2 / r0) * (kappa + L**2 / r0**2)
    ev1 = lambda t, y: y[0] - 2.0001
    ev1.terminal = True
    ev2 = lambda t, y: y[1]
    ev2.terminal = True
    ev2.direction = 1
    s = solve_ivp(rhs, [0, 1e7], [r0, -np.sqrt(rd2)], events=[ev1, ev2], rtol=1e-10, atol=1e-11)
    return len(s.t_events[0]) > 0


Lc = brentq(lambda L: (1.0 if captured(L) else -1.0), 3.5, 4.5, xtol=1e-3)
R.check("Massive particle with E=1 is captured only for L < 4M", abs(Lc - 4) < 0.01, f"threshold L = {Lc:.3f}")
bc = brentq(lambda b: (1.0 if captured(1.0 * b, kappa=0.0) else -1.0), 4.5, 6.0, xtol=1e-3)
R.check("Light is captured for b < 3 sqrt(3) M", abs(bc - 3 * np.sqrt(3)) < 0.01, f"b_c = {bc:.3f} (3 sqrt 3 = {3 * np.sqrt(3):.3f})")
R.note(f"Photon capture area = 27 pi M^2 = {27 / 16:.4f} x horizon area; used as a greybody estimate.")

AU = 1.495978707e11
a_sun = C.G * C.M_sun / AU**2
print(f"  At 1 AU: acceleration toward a 1 M_sun hole = {a_sun:.4e} m/s^2, the same as toward the Sun (shell theorem / Birkhoff).")
print(f"  r_s of 1 M_sun = {2 * r_g(C.M_sun) / 1e3:.3f} km; GR correction at 1 AU 3GM/(rc^2) = {3 * r_g(C.M_sun) / AU:.2e}")
v = 30e3
b_kill = 4 * r_g(C.M_sun) * C.c / v
print(f"  Capture impact parameter at v=30 km/s: b = 4GM/(vc) = {b_kill / 1e3:.0f} km = {b_kill / AU:.2e} AU")

# ------------------------------------------------------------------ angular momentum
R.section("3. Why gas must lose angular momentum, and what carries it away")
ell_pc = np.sqrt(C.G * M_SGRA * 3.0857e16)
ell_isco = np.sqrt(12) * C.G * M_SGRA / C.c
print(f"  Sgr A*: Keplerian specific ang. mom. at 1 pc = {ell_pc:.3e} m^2/s; at ISCO = {ell_isco:.3e} m^2/s; ratio {ell_pc / ell_isco:.0f}")
R.check("Gas from 1 pc must shed > 99% of its angular momentum to reach the ISCO", ell_isco / ell_pc < 0.01)
# two rings exchanging angular momentum at fixed total L: dE = (Om1 - Om2) dL < 0 if inner gives L to outer


def ring_E(r):
    return -0.5 / r


def ring_L(r):
    return np.sqrt(r)


r1, r2 = 10.0, 20.0
dL = 1e-4
E0 = ring_E(r1) + ring_E(r2)
r1n = (ring_L(r1) - dL) ** 2
r2n = (ring_L(r2) + dL) ** 2
R.check("Inner ring giving L to outer ring lowers total energy: inner falls, outer rises",
        ring_E(r1n) + ring_E(r2n) < E0 and r1n < r1 and r2n > r2,
        f"dE = {ring_E(r1n) + ring_E(r2n) - E0:.3e}, (Om1-Om2)dL = {(r1**-1.5 - r2**-1.5) * dL * -1:.3e}")

# viscous ring diffusion, Lynden-Bell & Pringle form, nu = const, x = 2 sqrt r, S = Sigma x
R.note("Viscous disk: dS/dt = (12 nu / x^2) d2S/dx2 with zero torque at the inner edge.")
x_in, x_out, N = 2 * np.sqrt(6.0), 2 * np.sqrt(4e5), 1500
x = np.linspace(x_in, x_out, N)
dx = x[1] - x[0]
nu = 1.0
x0 = 2 * np.sqrt(100.0)
S = np.exp(-0.5 * ((x - x0) / 0.6) ** 2)
S[0] = 0
S[-1] = 0


def mass(S_):
    return (np.pi / 4) * np.trapezoid(x**2 * S_, x)


def angmom(S_):
    return (np.pi / 8) * np.trapezoid(x**3 * S_, x)


M0, L0 = mass(S), angmom(S)
t_visc = 100.0**2 / nu  # r^2/nu
coef = 12 * nu / x**2
# Implicit Crank-Nicolson so that the run stays well inside the time limit
dtI = t_visc / 300.0
steps = 900
r_ = coef[1:-1] * dtI / dx**2 / 2
ab = np.zeros((3, N - 2))
ab[0, 1:] = -r_[:-1]
ab[1, :] = 1 + 2 * r_
ab[2, :-1] = -r_[1:]
Sint = S[1:-1].copy()
report = {}
for n in range(steps):
    rhs_ = Sint + r_ * (np.r_[Sint[1:], 0] - 2 * Sint + np.r_[0, Sint[:-1]])
    Snew = solve_banded((1, 1), ab, rhs_)
    Sint = Snew
    if n + 1 in (100, 300, 900):
        Sfull = np.r_[0, Sint, 0]
        report[n + 1] = (mass(Sfull), angmom(Sfull))
Sfull = np.r_[0, Sint, 0]
Mf, Lf = report[900]
print("  ring starts at r=100 (M=1); time in units of (100)^2/nu:")
for k, (m_, l_) in report.items():
    print(f"   t = {k / 300:.2f} r^2/nu : mass left {m_ / M0:6.3f}, angular momentum left {l_ / L0:6.3f}")
gone = 1 - Mf / M0
ell0 = L0 / M0
ell_in = x_in / 2
lost_L = (M0 - Mf) * ell_in
R.check("Most mass drains inward while most angular momentum stays in the disk", gone > 0.5 and Lf / L0 > 0.8,
        f"mass gone {100 * gone:.0f}%, L kept {100 * Lf / L0:.0f}%")
R.check("Angular momentum budget closes: L_final + L swallowed = L_0 within 2%",
        abs((Lf + lost_L) / L0 - 1) < 0.02, f"ratio {(Lf + lost_L) / L0:.4f}")
Sm = Sfull / (x + 1e-30)
r_arr = x**2 / 4
imed = np.cumsum(x**3 * Sfull)
r_L_half = r_arr[np.searchsorted(imed, imed[-1] / 2)]
print(f"  half of the remaining angular momentum now sits beyond r = {r_L_half:.0f} (started at r = 100): a little gas moved far out.")
R.check("A small mass carries the angular momentum outward", r_L_half > 300)

# ------------------------------------------------------------------ Eddington
R.section("4. Eddington limit from force balance on ionised hydrogen")
L_Edd = 4 * np.pi * C.G * C.M_sun * C.m_p * C.c / sigma_T
print(f"  L_Edd = 4 pi G M m_p c / sigma_T = {L_Edd:.3e} W per M_sun")
r_test = 1e12
frad = L_Edd * sigma_T / (4 * np.pi * r_test**2 * C.c)
fgrav = C.G * C.M_sun * C.m_p / r_test**2
R.check("At L_Edd radiation force on an electron equals gravity on a proton", abs(frad / fgrav - 1) < 1e-12)
eps = 0.1
t_s = (eps / (1 - eps)) * sigma_T * C.c / (4 * np.pi * C.G * C.m_p)
print(f"  e-folding (Salpeter) time for efficiency {eps}: {t_s / C.year / 1e6:.1f} Myr")

def reach_1e9(t, m):
    return m[0] - 1e9


reach_1e9.terminal = True
sol = solve_ivp(lambda t, m: m / t_s, [0, 3e16], [10.0], rtol=1e-10, events=reach_1e9)
tg9 = sol.t[-1] / C.year / 1e6
print(f"  10 M_sun to 1e9 M_sun at the Eddington rate takes {tg9:.0f} Myr")
R.check("Numerical growth time = t_S ln(1e8)", abs(tg9 * 1e6 * C.year / (t_s * np.log(1e8)) - 1) < 1e-3)
Mrho = np.sqrt(3 * C.c**6 / (32 * np.pi * C.G**3 * 1000.0)) / C.M_sun
print(f"  Mean density inside r_s equals water at M = {Mrho:.2e} M_sun (rho ~ M^-2)")

# ------------------------------------------------------------------ mergers
R.section("5. Mergers: orbital decay and the area theorem")
m1, m2 = 36 * C.M_sun, 29 * C.M_sun
k = (64 / 5) * C.G**3 * m1 * m2 * (m1 + m2) / C.c**5


def t_merge_num(a0):
    def merged(t, y):
        return y[0] - 0.001 * a0
    merged.terminal = True
    s = solve_ivp(lambda t, y: [-k / y[0] ** 3], [0, 2 * a0**4 / (4 * k)], [a0], method="LSODA", rtol=1e-9,
                  events=merged)
    return s.t[-1]


print("  Peters 1964 circular-orbit decay da/dt = -(64/5) G^3 m1 m2 (m1+m2) / (c^5 a^3), m = 36+29 M_sun:")
for a0, nm in ((1e9, "1e6 km"), (AU, "1 AU")):
    tn = t_merge_num(a0)
    tc = a0**4 / (4 * k)
    print(f"   from {nm:7s}: numeric {tn / C.year:.3e} yr, a^4/(4k) = {tc / C.year:.3e} yr")
    R.check(f"ODE matches closed form from {nm}", abs(tn / tc - 1) < 1e-3)
Ai = 36**2 + 29**2  # in units of 16 pi M_sun^2
print(f"  Hawking area theorem, non-spinning start: max radiated fraction = 1 - sqrt(m1^2+m2^2)/(m1+m2) = {1 - np.sqrt(Ai) / 65:.3f}")
print(f"  GW150914 measured: 3.0/65 = {3.0 / 65:.3f}")
Mf_, chi = 62, 0.67
Af = 0.5 * Mf_**2 * (1 + np.sqrt(1 - chi**2))
print(f"  Area, central values, ignoring initial spins: A_final/A_initial = {Af / Ai:.2f}")
R.check("GW150914 central values respect the area bound", 3.0 / 65 < 1 - np.sqrt(Ai) / 65 and Af > Ai)

R.section("6. Perseus 'note'")
f_mid = 261.63
for nm, f0 in (("B-flat below middle C", 233.08), ("B-flat above middle C", 466.16)):
    fq = f0 / 2**57
    print(f"  {nm}, 57 octaves down: f = {fq:.2e} Hz, period {1 / fq / C.year / 1e6:.1f} Myr")

R.close()
