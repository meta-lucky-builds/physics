"""C09-C17: Goedel metric, Kerr, Gott strings, Morris-Thorne wormhole, time-shifted mouths.

These tests check mathematical properties of solutions. They cannot show any of the
systems can be built or exists.
"""
import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize_scalar
from _harness import Checks
from common.constants import G, c, year, ly

C = Checks("Solutions of GR with closed timelike curves")
Gly = 1e9 * ly

# ---------------------------------------------------------------------------- Goedel
C.section("Goedel metric: symbolic field-equation check (signature -+++, units c = G = 1)")
t, r, ph, z = sp.symbols("t r phi z", real=True)
a = sp.symbols("a", positive=True)
s = sp.sinh(r)
coords = [t, r, ph, z]
g = 4 * a**2 * sp.Matrix([
    [-1, 0, sp.sqrt(2) * s**2, 0],
    [0, 1, 0, 0],
    [sp.sqrt(2) * s**2, 0, s**2 - s**4, 0],
    [0, 0, 0, 1],
])
ginv = sp.simplify(g.inv())
n = 4
Gam = [[[sp.simplify(sum(ginv[i, l] * (sp.diff(g[l, j], coords[k]) + sp.diff(g[l, k], coords[j])
                                      - sp.diff(g[j, k], coords[l])) for l in range(n)) / 2)
         for k in range(n)] for j in range(n)] for i in range(n)]


def ricci(j, k):
    e = 0
    for i in range(n):
        e += sp.diff(Gam[i][j][k], coords[i]) - sp.diff(Gam[i][j][i], coords[k])
        for l in range(n):
            e += Gam[i][i][l] * Gam[l][j][k] - Gam[i][k][l] * Gam[l][j][i]
    return sp.simplify(e)


Ric = sp.Matrix(n, n, lambda j, k: ricci(j, k))
Rs = sp.simplify(sum(ginv[j, k] * Ric[j, k] for j in range(n) for k in range(n)))
u_up = sp.Matrix([1 / (2 * a), 0, 0, 0])
u_low = g * u_up
print(f"  u^mu u_mu = {sp.simplify((u_up.T * u_low)[0])}, Ricci scalar R = {Rs}")
Lam, rho = sp.symbols("Lambda rho")
Ein = Ric - Rs * g / 2 + Lam * g - 8 * sp.pi * rho * u_low * u_low.T
Ein = Ein.applyfunc(sp.simplify)
sols = sp.solve([Ein[1, 1], Ein[0, 0]], [Lam, rho], dict=True)[0]
print(f"  dust plus Lambda solves the field equations with Lambda = {sols[Lam]}, 8 pi G rho = {sp.simplify(8*sp.pi*sols[rho])}")
resid = Ein.subs(sols).applyfunc(sp.simplify)
C.check("all 16 field-equation components vanish", resid == sp.zeros(n, n))
C.check("Lambda is negative: this is not the positive Lambda we measure", sols[Lam].subs(a, 1) < 0)

# Vorticity of the dust
du = sp.Matrix(n, n, lambda i, j: (sp.diff(u_low[j], coords[i]) - sp.diff(u_low[i], coords[j])) / 2)
w2 = sp.simplify(sum(du[i, j] * du[k, l] * ginv[i, k] * ginv[j, l]
                     for i in range(n) for j in range(n) for k in range(n) for l in range(n)) / 2)
print(f"  vorticity^2 = {w2}; 4 pi G rho = {sp.simplify(4*sp.pi*sols[rho])}")
C.check("Goedel's relation omega^2 = 4 pi G rho = -Lambda holds",
        sp.simplify(w2 - 4 * sp.pi * sols[rho]) == 0 and sp.simplify(w2 + sols[Lam]) == 0)

C.section("Goedel: sign change of g_phiphi and the tilting light cone")
gpp = g[2, 2]
rc = np.arcsinh(1.0)
print(f"  g_phiphi = 4 a^2 sinh^2 r (1 - sinh^2 r): zero at sinh r = 1, r_c = ln(1+sqrt 2) = {rc:.6f}")
gfun = sp.lambdify((r, a), gpp)
root = brentq(lambda x: gfun(x, 1.0), 0.5, 1.5)
C.check("numerical root of g_phiphi equals ln(1+sqrt 2)", abs(root - np.log(1 + np.sqrt(2))) < 1e-10)
gmat = sp.lambdify((r, a), g, "numpy")
print(f"  {'r/r_c':>7} {'sinh r':>7} {'g_phiphi/4a^2':>14} {'dt/dphi cone edges':>26}  phi circle")
tilt_ok = True
for rr in [0.3, 0.6, 0.9, 1.0, 1.1, 1.5, 2.0]:
    x = rr * rc
    sh = np.sinh(x)
    gtt, gtp, gppv = -1.0, np.sqrt(2) * sh**2, sh**2 - sh**4
    # null condition: -dt^2 + 2 gtp dt dphi + gpp dphi^2 = 0 gives dt/dphi = gtp +/- sqrt(gtp^2 + gpp)
    m1, m2 = gtp - np.sqrt(gtp**2 + gppv), gtp + np.sqrt(gtp**2 + gppv)
    kind = "spacelike" if gppv > 1e-12 else ("null" if abs(gppv) <= 1e-12 else "TIMELIKE")
    print(f"  {rr:7.2f} {sh:7.3f} {gppv:14.4f} {m1:12.4f}, {m2:12.4f}  {kind}")
    tilt_ok &= (m1 < 1e-9) if rr <= 1 else (m1 > 0)
C.check("lower cone edge crosses dt/dphi = 0 exactly at r_c: the cone has tipped past the phi circle", bool(tilt_ok))

C.section("Goedel: an explicit closed curve, tested for timelike and future-directed")
a_val = 1.0
r0 = 1.5 * rc
G0 = np.array(gmat(r0, a_val), dtype=float)
N_pts = 720
phis = np.linspace(0, -2 * np.pi, N_pts + 1)  # walking in the -phi direction
tang = np.array([0.0, 0.0, -1.0, 0.0])
u = np.array([1 / (2 * a_val), 0, 0, 0])
norm = tang @ G0 @ tang
futur = tang @ G0 @ u
pos_start = np.array([0.0, r0, phis[0], 0.0])
pos_end = np.array([0.0, r0, phis[-1] % (2 * np.pi), 0.0])
print(f"  r = 1.5 r_c, curve t = const, phi from 0 to -2 pi: g(T,T) = {norm:.4f} (negative = timelike), "
      f"g(T,u) = {futur:.4f} (negative = future-directed)")
C.check("the closed curve is timelike everywhere", norm < 0)
C.check("it is future-directed: the traveller always moves forward in their own light cone", futur < 0)
C.check("it closes: same t, r, z and phi mod 2 pi", np.allclose(pos_start, pos_end))
# At r < r_c the same curve fails.
G1 = np.array(gmat(0.8 * rc, a_val), dtype=float)
C.check("the same loop at r = 0.8 r_c is spacelike, so no such loop there", tang @ G1 @ tang > 0)

# Acceleration budget for constant-r loops. T = d/dphi at fixed t is a Killing orbit.
C.section("Goedel: numbers if our universe's density were used")
rho_crit = 3 * (67.4e3 / 3.0856775814913673e22)**2 / (8 * np.pi * G)
omega = np.sqrt(4 * np.pi * G * rho_crit)
a_si = 1 / (np.sqrt(float(sp.simplify(w2.subs(a, 1)))) ) / omega  # w^2 = k/a^2 in units c=1
a_m = a_si * c
r_proper = 2 * a_m * rc
print(f"  density = critical = {rho_crit:.2e} kg/m^3 -> omega = {omega:.2e} rad/s = {omega*year:.2e} rad/yr")
print(f"  rotation period 2 pi/omega = {2*np.pi/omega/year/1e9:.0f} Gyr")
print(f"  Goedel scale a = {a_m/Gly:.1f} Gly; proper CTC radius 2 a r_c = {r_proper/Gly:.1f} Gly")
ratio_h0 = omega / (67.4e3 / 3.0856775814913673e22)
print(f"  omega/H0 for this hypothetical = {ratio_h0:.2f}; Planck 2015 limit (omega/H)_0 < 7.6e-10; "
      f"observed ratio is below hypothetical by {ratio_h0/7.6e-10:.1e}")
H0_si = 67.4e3 / 3.0856775814913673e22
print(f"  that CMB limit in today's rate: omega < {7.6e-10*H0_si*year:.1e} rad/yr (Bianchi VII_h fit, model dependent)")
C.check("a Goedel rotation at the real density is excluded by CMB limits by a factor > 1e8",
        ratio_h0 / 7.6e-10 > 1e8)


def ta(u2):
    """Total integrated acceleration of the constant-r t=const loop, from Killing-orbit formula."""
    sh = np.sqrt(u2)
    x = np.arcsinh(sh)
    N = lambda y: np.sqrt(np.sinh(y)**4 - np.sinh(y)**2)
    dN = (N(x + 1e-6) - N(x - 1e-6)) / 2e-6
    return 2 * np.pi * dN  # |a| = (dN/dr)/(N * 2a) and length = 2 pi * 2a N, so a cancels


res = minimize_scalar(ta, bounds=(1.02, 5), method="bounded")
print(f"  cheapest constant-r t=const loop: sinh^2 r = {res.x:.3f}, total integrated acceleration "
      f"= {res.fun:.2f} (units of c)")
print(f"  ideal photon-rocket mass ratio = exp(TA) = {np.exp(res.fun):.2e} "
      "(Natario 2012 gives a lower-cost non-circular candidate)")
sh = np.sqrt(res.x)
L = 2 * np.pi * 2 * a_m * np.sqrt(sh**4 - sh**2)
print(f"  the traveller's own clock runs {L/Gly:.0f} Gyr per loop at that r")
C.check("even an ideal rocket faces a large mass ratio", np.exp(res.fun) > 1e5)

# ---------------------------------------------------------------------------- Kerr
C.section("Kerr (M = 1, spin a): horizons and where g_phiphi < 0 on the equatorial plane")
for spin in [0.5, 0.9, 0.998]:
    rp = 1 + np.sqrt(1 - spin**2)
    rm = 1 - np.sqrt(1 - spin**2)
    gpp_eq = lambda x: x**2 + spin**2 + 2 * spin**2 / x
    r_ctc = brentq(gpp_eq, -5, -1e-6)
    print(f"  a = {spin}: r+ = {rp:.3f}, r- = {rm:.3f}; g_phiphi < 0 for {r_ctc:.3f} < r < 0")
    C.check(f"a = {spin}: CTC region lies at r < 0, inside r- and r+, behind the horizons",
            r_ctc < 0 and gpp_eq(r_ctc / 2) < 0 and gpp_eq(0.5 * rm) > 0)
print("  Over-spun a > M: 1 - a^2 < 0, no real horizon, so the same region would be exposed.")

# ---------------------------------------------------------------------------- Gott
C.section("Gott strings: holonomy test in 2+1 dimensions")
# Each string gives a rotation by deficit angle D = 8 pi G mu in its rest frame. Going round
# both strings gives the product of the two boosted rotations. CTCs in Gott's construction are
# possible only if the product is a boost (hyperbolic): trace of SO(2,1) matrix > 3.


def boost(b):
    gm = 1 / np.sqrt(1 - b**2)
    return np.array([[gm, gm * b, 0], [gm * b, gm, 0], [0, 0, 1]])


def rot(D):
    return np.array([[1, 0, 0], [0, np.cos(D), -np.sin(D)], [0, np.sin(D), np.cos(D)]])


def trace_pair(b, D):
    r1 = boost(b) @ rot(D) @ boost(-b)
    r2 = boost(-b) @ rot(D) @ boost(b)
    return np.trace(r1 @ r2)


for gmu in [1e-1, 1e-2, 1e-3]:
    D = 8 * np.pi * gmu
    thr = brentq(lambda b: trace_pair(b, D) - 3, 0.0, 1 - 1e-12, xtol=1e-14)
    print(f"  G mu = {gmu}: holonomy turns boost-like at beta = {thr:.8f}; cos(4 pi G mu) = {np.cos(D/2):.8f}")
    C.check(f"G mu = {gmu}: numerical threshold reproduces Gott's beta > cos(4 pi G mu)",
            abs(thr - np.cos(D / 2)) < 1e-6)
for gmu in [1e-7]:
    x = 4 * np.pi * gmu
    print(f"  G mu = {gmu}: need gamma > 1/sin(4 pi G mu) = {1/np.sin(x):.2e}, 1-beta < {x**2/2:.1e}")
gm_needed = 1 / np.sin(4 * np.pi * 1e-7)
E_ratio = gm_needed
print(f"  each string's energy must exceed its rest energy by {E_ratio:.1e} (Planck limit G mu < 1.3e-7 used)")

# ---------------------------------------------------------------------------- Wormhole
C.section("Morris-Thorne wormhole: symbolic energy condition at the throat")
tt, l, th, phi = sp.symbols("t l theta phi", real=True)
R = sp.Function("R")(l)
co = [tt, l, th, phi]
gw = sp.diag(-1, 1, R**2, R**2 * sp.sin(th)**2)
gwi = gw.inv()
Gw = [[[sum(gwi[i, m] * (sp.diff(gw[m, j], co[k]) + sp.diff(gw[m, k], co[j]) - sp.diff(gw[j, k], co[m]))
            for m in range(4)) / 2 for k in range(4)] for j in range(4)] for i in range(4)]


def ric_w(j, k):
    e = 0
    for i in range(4):
        e += sp.diff(Gw[i][j][k], co[i]) - sp.diff(Gw[i][j][i], co[k])
        for m in range(4):
            e += Gw[i][i][m] * Gw[m][j][k] - Gw[i][k][m] * Gw[m][j][i]
    return sp.simplify(e)


Rw = sp.Matrix(4, 4, lambda j, k: ric_w(j, k))
Rsw = sp.simplify(sum(gwi[j, j] * Rw[j, j] for j in range(4)))
Gtt = sp.simplify(Rw[0, 0] - Rsw * gw[0, 0] / 2)
Gll = sp.simplify(Rw[1, 1] - Rsw * gw[1, 1] / 2)
nec = sp.simplify(Gtt + Gll)
print(f"  8 pi G (rho + p_radial) = G_tt + G_ll = {nec}")
C.check("rho + p_radial = -(2 R''/R)/(8 pi G): flare-out R'' > 0 forces NEC violation",
        sp.simplify(nec + 2 * sp.diff(R, l, 2) / R) == 0)
b0 = sp.symbols("b0", positive=True)
Rcat = sp.sqrt(b0**2 + l**2)
val = nec.subs(R, Rcat).doit().subs(l, 0)
val = sp.simplify(val)
print(f"  Ellis throat R = sqrt(b0^2 + l^2): at l = 0, 8 pi G (rho + p) = {val}")
C.check("the sign is negative for a concrete smooth throat", val.subs(b0, 1) < 0)
print(f"  for a 1 m throat in SI: rho + p/c^2 = -c^2/(4 pi G b0^2) = {-c**2/(4*np.pi*G):.2e} kg/m^3 "
      "(compare water, 1e3 kg/m^3)")

C.section("Making a time machine from a wormhole: moving one mouth")
beta = 0.9995
gmm = 1 / np.sqrt(1 - beta**2)
for T_ext in [10.0, 100.0]:
    shift = T_ext * (1 - 1 / gmm)
    print(f"  mouth round trip {T_ext:.0f} yr outside at 0.9995 c: mouth clock {T_ext/gmm:.2f} yr; time shift {shift:.2f} yr")
print("  CTCs appear once shift exceeds the external separation D/c: shift 9.68 yr and D < 9.68 ly.")

C.finish()
