"""Schwarzschild interior: curvature, free fall, proper-time bound, tides, depth, volume.

Geometric units G=c=1 with M=1 unless a line says otherwise.
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from scipy.optimize import minimize_scalar

from bh_tools import C, Report, M_SGRA, r_g, t_P, t_g, curvature

R = Report("verify_geometry: Schwarzschild curvature, free fall, tides, depth")

# ---------------------------------------------------------------- symbolic curvature
R.section("1. Curvature computed symbolically from the metric")
t, r, th, ph, M = sp.symbols("t r theta phi M", positive=True)
f = 1 - 2 * M / r
gS = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th) ** 2)
_, Ric, K, _, _ = curvature(gS, [t, r, th, ph])
R.check("Schwarzschild is vacuum (all Ricci components vanish)", Ric == sp.zeros(4, 4))
R.check("Kretschmann K = 48 M^2 / r^6", sp.simplify(K - 48 * M**2 / r**6) == 0, f"sympy: K = {K}")

# Painleve-Gullstrand form, used for the falling-frame tidal tensor
T = sp.symbols("T", positive=True)
s = sp.sqrt(2 * M / r)
gPG = sp.Matrix([[-(1 - s**2), s, 0, 0], [s, 1, 0, 0], [0, 0, r**2, 0],
                 [0, 0, 0, r**2 * sp.sin(th) ** 2]])
Rd, RicPG, KPG, _, _ = curvature(gPG, [T, r, th, ph])
R.check("PG metric is the same geometry (Ricci = 0, same K)",
        RicPG == sp.zeros(4, 4) and sp.simplify(KPG - 48 * M**2 / r**6) == 0)
R.check("PG slices T=const are flat: dr^2 + r^2 dOmega^2", gPG[1, 1] == 1 and gPG[2, 2] == r**2)


def contract(Rm, a, b, c, d):
    return sp.simplify(sum(Rm[i][j][k][l] * a[i] * b[j] * c[k] * d[l] for i in range(4)
                           for j in range(4) for k in range(4) for l in range(4)))


u = [1, -s, 0, 0]          # four-velocity of a faller with E=1 in PG coordinates
er = [0, 1, 0, 0]          # radial unit vector of the falling frame
eth = [0, 0, 1 / r, 0]     # transverse unit vector
Rr = contract(Rd, er, u, er, u)
Rt = contract(Rd, eth, u, eth, u)
R.check("Falling-frame radial tidal term: d2xi/dtau2 = +2M/r^3 xi", sp.simplify(-Rr - 2 * M / r**3) == 0)
R.check("Falling-frame transverse tidal term: d2xi/dtau2 = -M/r^3 xi", sp.simplify(-Rt + M / r**3) == 0)

# ---------------------------------------------------------------- Kasner
R.section("2. Near r=0 the interior is a Kasner cosmology, not a point")
tau, x, y, z = sp.symbols("tau x y z", positive=True)
p = sp.Rational(-1, 3), sp.Rational(2, 3), sp.Rational(2, 3)
gK = sp.diag(-1, tau ** (2 * p[0]), tau ** (2 * p[1]), tau ** (2 * p[2]))
_, RicK, KK, _, _ = curvature(gK, [tau, x, y, z])
R.check("Kasner (-1/3,2/3,2/3) is an exact vacuum solution", RicK == sp.zeros(4, 4))
R.check("Kasner K = 64/(27 tau^4)", sp.simplify(KK - sp.Rational(64, 27) / tau**4) == 0, f"sympy: K = {KK}")


def tau_E(rr, E=0.0):
    """Proper time from radius rr down to r=0 for energy E (rr<=2)."""
    # substitution r = rr*w^2 removes the sqrt endpoint behaviour at r=0
    def g(w):
        q = rr * w**2
        return 2 * rr * w / np.sqrt(E**2 - 1 + 2 / q)
    return quad(g, 0, 1, limit=400, epsabs=1e-13, epsrel=1e-12)[0]


rs = np.array([1e-5, 2e-5])
ts = np.array([tau_E(q) for q in rs])
slope_r = np.log(rs[1] / rs[0]) / np.log(ts[1] / ts[0])
ax = np.sqrt(2 / rs - 1)
slope_x = np.log(ax[1] / ax[0]) / np.log(ts[1] / ts[0])
R.check("Exact interior: sphere radius ~ tau^(2/3)", abs(slope_r - 2 / 3) < 1e-3, f"fit {slope_r:.4f}")
R.check("Exact interior: t-direction scale ~ tau^(-1/3) (stretches)", abs(slope_x + 1 / 3) < 1e-3, f"fit {slope_x:.4f}")
R.note("Sum of exponents = 1 means the 3-volume ~ tau^1 -> 0 while one length -> infinity.")

# K as a function of remaining proper time is independent of M
R.section("3. Curvature against proper time remaining: the mass drops out")
print("  M/M_sun     E     tau_left(s)   K*c^4*tau^4 (limit 64/27 = %.4f)" % (64 / 27))
ok = True
for Ms in (10.0, 4.15e6, 6.5e9):
    tg = t_g(Ms * C.M_sun)
    for E in (1.0, 0.3):
        rr = 1e-5
        tl = tau_E(rr, E)
        val = 48 / rr**6 * tl**4
        print(f"  {Ms:8.3g}  {E:4.1f}  {tl * tg:11.3e}   {val:.4f}")
        ok &= abs(val - 64 / 27) < 1e-3
R.check("K -> 64/(27 c^4 tau^4) for every mass and energy", ok)
tau_pl = (64 / 27) ** 0.25 * t_P
R.note(f"Curvature reaches K = 1/l_P^4 when tau_left = (64/27)^(1/4) t_P = {tau_pl:.2e} s, for ANY mass.")

# ---------------------------------------------------------------- free fall
R.section("4. Free fall: geodesic equation integrated in ingoing EF coordinates")


def ef_rhs(_, y_):
    v, r_, vd, rd = y_
    F = 1 - 2 / r_
    vdd = -(1 / r_**2) * vd**2
    rdd = F * vdd + (2 / r_**2) * rd * vd
    return [vd, rd, vdd, rdd]


def fall_time(r0):
    F0 = 1 - 2 / r0
    y0 = [0, r0, 1 / np.sqrt(F0), 0.0]
    ev = lambda t_, y_: y_[1] - 0.02
    ev.terminal = True
    sol = solve_ivp(ef_rhs, [0, 1e7], y0, events=ev, rtol=1e-11, atol=1e-13, method="DOP853")
    E = np.sqrt(F0)
    return sol.t[-1] + tau_E(0.02, E)  # analytic-quadrature tail to r=0


print("  r0/M    tau(r0 to 0) numeric   cycloid (pi/2) r0 sqrt(r0/2)")
ok = True
for r0 in (2.2, 3.0, 6.0, 20.0, 200.0):
    a, b = fall_time(r0), (np.pi / 2) * r0 * np.sqrt(r0 / 2)
    print(f"  {r0:6.1f}  {a:18.6f}  {b:18.6f}")
    ok &= abs(a / b - 1) < 2e-5
R.check("Second-order geodesic integration matches the cycloid formula", ok)

R.section("5. Horizon to singularity: the proper-time bound pi*M")
print("  E        tau(2M to 0)/M   (E=1 exact: 4/3; E->0 limit: pi)")
tvals = []
for E in (2.0, 1.0, 0.7, 0.3, 0.1, 0.01, 0.0):
    tv = tau_E(2.0, E)
    tvals.append(tv)
    print(f"  {E:5.2f}  {tv:14.6f}")
R.check("E=1 gives 4/3", abs(tvals[1] - 4 / 3) < 1e-8, f"{tvals[1]:.8f}")
R.check("E=0 gives pi", abs(tvals[-1] - np.pi) < 1e-5, f"{tvals[-1]:.8f}")
R.check("Proper time rises monotonically as E falls to 0", all(np.diff(tvals) > 0))

# Random non-geodesic timelike paths: rockets and wiggles.  With r = M(1+cos eta):
# dtau^2 = M^2 (1+cos eta)^2 d eta^2 - g dt^2 - r^2 dphi^2,  g = tan^2(eta/2).
rng = np.random.default_rng(6)
eta = np.linspace(0, np.pi, 4001)
best = 0
n_ok = 0
for _ in range(400):
    A, k, B, q = rng.uniform(0, 1.2), rng.integers(1, 6), rng.uniform(0, 0.5), rng.integers(1, 6)
    rr = 1 + np.cos(eta)
    drd = 1 + np.cos(eta)
    g = np.tan(eta / 2 * (1 - 1e-12)) ** 2
    tp = A * np.sin(k * eta) * np.cos(eta / 2) ** 2
    pp = B * np.sin(q * eta)
    rad = drd**2 - g * tp**2 - (rr * pp) ** 2 * 0.25
    if rad.min() < -1e-12:
        continue
    n_ok += 1
    best = max(best, np.trapezoid(np.sqrt(np.maximum(rad, 0)), eta))
R.check(f"{n_ok} random accelerated timelike paths all stay below pi", n_ok > 50 and best < np.pi,
        f"longest {best:.5f}")
R.check("Grid resolves the bound itself", abs(np.trapezoid(1 + np.cos(eta), eta) - np.pi) < 1e-6)

yr = C.year
M1 = yr / (np.pi * t_g(C.M_sun))
print(f"\n  One year from horizon to singularity: at the bound pi*M, M = {M1:.3e} M_sun; "
      f"falling from far away (4M/3), M = {yr / ((4 / 3) * t_g(C.M_sun)):.3e} M_sun")
R.check("Mass for one year at the bound is about 2e12 M_sun", 1.9e12 < M1 < 2.2e12)

print("\n  M/M_sun      tau_max=pi*M(s)  tau(E=1)=4M/3(s)")
for Ms in (1, 10, 4.15e6, 6.5e9, M1):
    print(f"  {Ms:9.3g}  {np.pi * t_g(Ms * C.M_sun):12.3e}  {(4 / 3) * t_g(Ms * C.M_sun):12.3e}")

# ---------------------------------------------------------------- tides
R.section("6. Tides on a falling body, checked by differencing two trajectories")


def rhs(_, y_):
    return [y_[1], -1 / y_[0] ** 2]


ok = True
for r0 in (10.0, 2.0, 0.5):
    d = 1e-5
    sol1 = solve_ivp(rhs, [0, 0.05], [r0, -np.sqrt(2 / r0)], rtol=1e-13, atol=1e-15, dense_output=True)
    sol2 = solve_ivp(rhs, [0, 0.05], [r0 + d, -np.sqrt(2 / (r0 + d))], rtol=1e-13, atol=1e-15, dense_output=True)
    h = 0.02
    sep = lambda q: sol2.sol(q)[0] - sol1.sol(q)[0]
    acc = (sep(2 * h) - 2 * sep(h) + sep(0)) / h**2 / sep(h)
    pred = 2 / sol1.sol(h)[0] ** 3
    print(f"  r={r0:4.1f}M  measured stretching rate {acc:.5f}  2M/r^3 = {pred:.5f}")
    ok &= abs(acc / pred - 1) < 0.02
R.check("Separation of two infallers grows at 2M/r^3, inside and outside the horizon", ok)

print("\n  Assumption (an estimate, mine): body length L=2 m, lethal differential acceleration 100 g.")
L, lim = 2.0, 100 * C.g_0
print("  M/M_sun    r_s(km)    tide at horizon(m/s^2)   r_lethal/r_s   tau left at r_lethal(s)   at horizon")
for Ms in (1, 10, 1e3, 1e4, 4.15e6, 6.5e9):
    Mk = Ms * C.M_sun
    rs_ = 2 * r_g(Mk)
    tide = 2 * C.G * Mk * L / rs_**3
    rl = (2 * C.G * Mk * L / lim) ** (1 / 3)
    tl = (2 / 3) * rl**1.5 / np.sqrt(2 * C.G * Mk)
    print(f"  {Ms:8.3g}  {rs_ / 1e3:9.3g}  {tide:18.3e}  {rl / rs_:14.3g}  {tl:16.3e}          "
          f"{'survives' if tide < lim else 'already torn apart'}")
M100 = np.sqrt(C.c**6 * L / (4 * C.G**2 * lim)) / C.M_sun
M1g = np.sqrt(C.c**6 * L / (4 * C.G**2 * C.g_0)) / C.M_sun
print(f"  Horizon tide = c^6 L/(4 G^2 M^2): equals 100 g at M = {M100:.3g} M_sun, 1 g at M = {M1g:.3g} M_sun")
R.check("Survival threshold lies between stellar and supermassive holes", 10 < M100 < 4e6)
# Sun-like star: tidal radius R*(M/M*)^(1/3) equals r_s
Mh = (C.c**2 * 6.96e8 / (2 * C.G * C.M_sun ** (1 / 3))) ** 1.5 / C.M_sun
print(f"  Sun-like star: tidal radius R(M/M_sun)^(1/3) = r_s at M = {Mh:.2e} M_sun")

# ---------------------------------------------------------------- depth
R.section("7. 'Depth': what is invariant and what is not")
print("  M/M_sun   area(m^2)    K(horizon)(m^-4)  tide/L at horizon(s^-2)  tau_max(s)  mean density(kg/m^3)")
rows = {}
for Ms in (10, 4.15e6, 6.5e9):
    Mk = Ms * C.M_sun
    rs_ = 2 * r_g(Mk)
    A = 4 * np.pi * rs_**2
    Kh = 48 * r_g(Mk) ** 2 / rs_**6
    rho = Mk / (4 / 3 * np.pi * rs_**3)
    rows[Ms] = (A, Kh)
    print(f"  {Ms:8.3g}  {A:10.3e}  {Kh:14.3e}  {C.c**6 / (4 * C.G**2 * Mk**2):18.3e}  "
          f"{np.pi * t_g(Mk):10.3e}  {rho:12.3e}")
ratio = rows[6.5e9][0] / rows[10][0]
R.check("Area scales as M^2", abs(ratio / (6.5e8) ** 2 - 1) < 1e-9)
R.check("Horizon curvature scales as M^-4", abs(rows[10][1] / rows[6.5e9][1] / (6.5e8) ** 4 - 1) < 1e-9)
zdepth = 2 * np.sqrt(2 * (10 - 2))
lprop = quad(lambda q: 1 / np.sqrt(1 - 2 / q), 2, 10)[0]
print(f"  Between r=2M and r=10M: Flamm funnel depth z = {zdepth:.3f} M and radial length {lprop:.3f} M on a "
      f"t=const slice; on a PG slice depth 0 and length 8.000 M")
R.check("Same spacetime, different slice, different 'depth' (8M funnel vs flat)", abs(zdepth - 8) < 1e-12)

R.section("8. Interior volume (Christodoulou-Rovelli maximal slice)")
res = minimize_scalar(lambda q: -q**2 * np.sqrt(2 / q - 1), bounds=(0.1, 1.99), method="bounded",
                      options={"xatol": 1e-9})
R.check("r^2 sqrt(2M/r-1) peaks at r = 3M/2", abs(res.x - 1.5) < 1e-6, f"{res.x:.6f}")
R.check("Peak value is (3 sqrt 3/4) M^2", abs(-res.fun - 3 * np.sqrt(3) / 4) < 1e-9)
for name, Mk in (("Sgr A*", M_SGRA), ("10 M_sun", 10 * C.M_sun)):
    Vol = 3 * np.sqrt(3) * np.pi * r_g(Mk) ** 2 * C.c * 13.8e9 * C.year
    Vout = 4 / 3 * np.pi * (2 * r_g(Mk)) ** 3
    print(f"  {name}: V(after v=13.8 Gyr) = {Vol:.3e} m^3, while (4/3)pi r_s^3 = {Vout:.3e} m^3, ratio {Vol / Vout:.2e}")

R.close()
