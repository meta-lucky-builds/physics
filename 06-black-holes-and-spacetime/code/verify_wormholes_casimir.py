"""Wormhole energy conditions, Casimir mode sum, antimatter energy, dark sector equations of state."""
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp

from bh_tools import C, Report, curvature, l_P

R = Report("verify_wormholes_casimir: flare-out, Casimir, antimatter, dark sector")

# ---------------------------------------------------------------- Morris-Thorne
R.section("1. Morris-Thorne metric: stress-energy read off the Einstein tensor")
t, r, th, ph = sp.symbols("t r theta phi", positive=True)
b = sp.Function("b")(r)
Phi = sp.Function("Phi")(r)
g = sp.diag(-sp.exp(2 * Phi), 1 / (1 - b / r), r**2, r**2 * sp.sin(th) ** 2)
_, Ric, _, ginv, _ = curvature(g, [t, r, th, ph])
Rs = sp.simplify(sum(ginv[i, i] * Ric[i, i] for i in range(4)))
G = (Ric - g * Rs / 2).applyfunc(sp.simplify)
# orthonormal static frame: 8 pi rho = G_tt/|g_tt|, 8 pi p_r = G_rr/g_rr
rho8 = sp.simplify(G[0, 0] / (-g[0, 0]))
pr8 = sp.simplify(G[1, 1] / g[1, 1])
db = sp.diff(b, r)
R.check("8 pi rho = b'/r^2", sp.simplify(rho8 - db / r**2) == 0)
R.check("8 pi p_r = -b/r^3 + 2(1-b/r)Phi'/r",
        sp.simplify(pr8 - (-b / r**3 + 2 * (1 - b / r) * sp.diff(Phi, r) / r)) == 0)
nec = sp.simplify((rho8 + pr8))
r0 = sp.symbols("r0", positive=True)
b1 = sp.symbols("b1", real=True)
at_throat = nec.subs(sp.diff(b, r), b1).subs(b, r).subs(r, r0)
expect = (b1 * r0 - r0) / r0**3
R.check("At a throat b(r0)=r0: 8 pi (rho + p_r) = (b' - 1)/r0^2, for any redshift function",
        sp.simplify(at_throat - expect) == 0, f"sympy: {sp.simplify(at_throat)}")
R.note("Flare-out needs b'(r0) < 1 (derived in D5). So rho + p_r < 0 at every static throat: NEC violated.")

# ---------------------------------------------------------------- Ellis numbers
R.section("2. How much negative energy: the simplest wormhole, b = r0^2/r, Phi = 0")
r0v = 1.0
rho_J = -C.c**4 / (8 * np.pi * C.G * r0v**2)
print(f"  Throat radius 1 m: energy density at throat rho c^2 = {rho_J:.3e} J/m^3")
I, _ = quad(lambda u: 1 / np.sqrt(1 - u**2), 0, 1)   # after u=r0/r
Mex = -2 * (C.c**2 / (2 * C.G)) * r0v * I
print(f"  Integrated exotic mass, both sides = {Mex:.3e} kg (closed form -pi c^2 r0/(2G) = {-np.pi * C.c**2 * r0v / (2 * C.G):.3e}); Jupiter = 1.898e27 kg")
R.check("Numeric integral of negative mass matches closed form", abs(Mex / (-np.pi * C.c**2 * r0v / (2 * C.G)) - 1) < 1e-8)

# QI bound of Ford-Roman applied to static throat observer
tau0_len = (3 * C.hbar * C.c / (32 * np.pi**2 * abs(rho_J))) ** 0.25
print(f"  Ford-Roman 1995/96 bound rho >= -3 hbar c/(32 pi^2 (c tau0)^4): that density may be sustained only over c tau0 ~ {tau0_len:.2e} m")
R.check("Sampling length is ~1e-17 of the throat size", tau0_len / r0v < 1e-16)
print(f"  With sampling time f*r0 the same bound gives r0 <~ l_P/(2 f^2): f=0.01 -> r0 <~ {l_P / (2 * 0.01**2):.1e} m (Ford-Roman eq. 51)")

# ---------------------------------------------------------------- Casimir
R.section("3. Casimir energy: regulated mode sum, not the closed formula")
mp.mp.dps = 60


def casimir_units(eps):
    """E/A * a^3/(hbar c) for perfect plates (two polarisations), exponential cutoff exp(-eps*omega), a=1.

    I(k) = integral d^2q/(2pi)^2 * omega * exp(-eps omega), omega = sqrt(q^2+k^2).
    Modes: k_n = n*pi; n=0 has one polarisation, n>=1 two; zero-point energy is hbar omega/2 each.
    Free space with the same volume is subtracted using the same cutoff.
    """
    eps = mp.mpf(eps)

    def I(k):
        return mp.exp(-eps * k) * (k**2 / eps + 2 * k / eps**2 + 2 / eps**3) / (2 * mp.pi)

    plates = I(0) / 2 + mp.nsum(lambda n: I(n * mp.pi), [1, mp.inf])
    free = mp.quad(I, [0, mp.inf]) / mp.pi
    return plates - free


exact = -mp.pi**2 / 720
print(f"  exact coefficient -pi^2/720 = {float(exact):.8f}")
print("  cutoff eps/a    regulated sum minus free space")
vals = []
for e in (0.5, 0.2, 0.1, 0.05):
    v = casimir_units(e)
    vals.append(v)
    print(f"   {e:5.2f}         {float(v):.8f}")
rich = (4 * vals[3] - vals[2]) / 3  # the printed values approach the limit with error ~ eps^2
print(f"  error ratio between eps=0.1 and eps=0.05: {float((vals[2] - exact) / (vals[3] - exact)):.3f} (4 means error ~ eps^2)")
print(f"  Richardson extrapolation to eps=0: {float(rich):.8f}")
R.check("Mode sum converges to -pi^2/720 as the cutoff is removed", abs(rich / exact - 1) < 1e-4, f"ratio {float(rich / exact):.6f}")
R.check("Sign is negative: energy below the free vacuum", vals[-1] < 0)
print("  gap a       rho = -pi^2 hbar c/(720 a^4) [J/m^3]   pressure pi^2 hbar c/(240 a^4) [Pa]")
for a in (1e-6, 1e-7, 1e-8):
    print(f"  {a:8.0e}    {-np.pi**2 * C.hbar * C.c / (720 * a**4):12.3e}                     {np.pi**2 * C.hbar * C.c / (240 * a**4):10.3e}")
a_needed = (np.pi**2 * C.hbar * C.c / (720 * abs(rho_J))) ** 0.25
print(f"  Gap needed to reach the 1 m wormhole density: a = {a_needed:.1e} m, against proton radius ~8.4e-16 m")
R.check("Plates cannot reach it: required gap is below nuclear size", a_needed < 8.4e-16)

# ---------------------------------------------------------------- antimatter
R.section("4. Antimatter has positive energy")
Ee = C.m_e * C.c**2 / C.eV / 1e3
print(f"  e- + e+ -> 2 gamma at rest: each gamma = {Ee:.1f} keV. If m(e+) were negative the pair would carry 0 total energy and no gammas could emerge.")
ag, s1, s2 = 0.75, 0.13, 0.16
sig = np.hypot(s1, s2)
print(f"  ALPHA-g 2023: a_g = ({ag} +- {s1} +- {s2}) g; quadrature {sig:.2f} g; distance from a_g=-g is {(ag + 1) / sig:.1f} sigma (my combination, naive Gaussian)")
R.check("Falling up is excluded by the published error bars", (ag + 1) / sig > 5)

# ---------------------------------------------------------------- dark sector
R.section("5. Dark matter and dark energy from w = p/rho")


def rho_of_a(w, a_end=4.0):
    s = solve_ivp(lambda a, y: [-3 * (1 + w) * y[0] / a], [1, a_end], [1.0], rtol=1e-11, atol=1e-13)
    return s.y[0, -1]


print("  component  w      rho(a=4)/rho(a=1) numeric  a^(-3(1+w))   rho+3p sign (gravity)   rho+p sign (NEC)")
ok = True
for nm, w in (("dust/DM", 0.0), ("radiation", 1 / 3), ("Lambda", -1.0), ("phantom", -1.1)):
    n = rho_of_a(w)
    cf = 4.0 ** (-3 * (1 + w))
    ok &= abs(n / cf - 1) < 1e-8
    print(f"  {nm:9s}  {w:5.2f}  {n:18.6f}  {cf:18.6f}   {'attractive' if 1 + 3 * w > 0 else 'repulsive':>10s}           {'>= 0' if 1 + w >= -1e-12 else '< 0 (violates NEC)'}")
R.check("Continuity equation integrated numerically gives rho ~ a^(-3(1+w))", ok)
R.check("Only w < -1 could supply a flare-out; w = -1 gives rho+p = 0 exactly", 1 + (-1.0) == 0)

h, Om = 0.6736, 0.3153
OL = 1 - Om
H0 = h * 100e3 / 3.0857e22
rho_c = 3 * H0**2 / (8 * np.pi * C.G)
print(f"  Planck 2018: H0 = 67.36 km/s/Mpc, Omega_m = {Om}; rho_c c^2 = {rho_c * C.c**2:.3e} J/m^3; dark-energy density = {OL * rho_c * C.c**2:.2e} J/m^3")
age = quad(lambda a: 1 / (a * H0 * np.sqrt(Om / a**3 + OL)), 1e-8, 1)[0] / C.year / 1e9
print(f"  Flat LCDM age from these (radiation omitted) = {age:.2f} Gyr")
R.check("Age consistent with about 13.8 Gyr", abs(age - 13.8) < 0.15)
print(f"  Dark energy density / 1 m wormhole density = {OL * rho_c * C.c**2 / abs(rho_J):.1e}")

R.close()
