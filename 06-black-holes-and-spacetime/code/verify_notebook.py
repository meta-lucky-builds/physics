"""The notebook paradox: what a falling or hovering notebook can actually see.

Geometric units G=c=1, M=1 unless stated.  Outside clocks are labelled by the ingoing
null coordinate v = t + r*: light emitted inward at v reaches anything at the same v.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.integrate import quad, solve_ivp  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

from bh_tools import C, M_M87, M_SGRA, Report, T_hawking, l_P, r_g, t_evap, t_g  # noqa: E402

R = Report("verify_notebook: static vs falling clocks, signals received, evaporation map")


def f(r):
    return 1 - 2 / r


def rate_fall(r, E=1.0):
    """d(v of far clock seen)/d(tau) for a radial free faller, from the EF derivation."""
    return 1 / (E + np.sqrt(E**2 - f(r)))


# ------------------------------------------------------------------ 1. dilation table
R.section("1. Rate at which a far clock is seen: hovering vs falling (formulas of D8)")
print("  r/M     hover: 1/sqrt(f)     faller E=1: dt/dtau = E/f     faller E=1: seen rate dv/dtau")
for r in (100.0, 10.0, 3.0, 2.1, 2.001, 2.0, 1.0, 0.1, 0.001):
    hov = f"{1 / np.sqrt(f(r)):12.4f}" if r > 2 else ("    infinite" if r == 2 else "  impossible")
    dtd = f"{1 / f(r):14.4f}" if r != 2 else "      infinite"
    print(f"  {r:6.3f}  {hov}        {dtd}                  {rate_fall(r):.4f}")
R.check("Faller crossing the horizon sees far clocks at exactly half rate (E=1)", abs(rate_fall(2.0) - 0.5) < 1e-15)

# ------------------------------------------------------------------ 2. ray tracing in PG
R.section("2. Independent test: light rays traced in Painleve-Gullstrand coordinates")
R0 = 100.0
T_sing = (2 / 3) * R0**1.5 / np.sqrt(2)  # from dr/dT=-sqrt(2/r) with T=tau for the E=1 faller


def r_faller(T):
    return (R0**1.5 - (3 / np.sqrt(2)) * T) ** (2 / 3)


def trace(T_e):
    """Emit an ingoing radial ray from r=R0 at PG time T_e. Return (T_arrive, r_arrive) or None."""
    def rhs(T, y):
        return [-1 - np.sqrt(2 / max(y[0], 1e-12))]

    def catch(T, y):
        return y[0] - r_faller(min(T, T_sing - 1e-12))

    catch.terminal = True
    catch.direction = -1

    def end(T, y):
        return y[0] - 1e-5

    end.terminal = True
    # the faller is already below R0 at T_e>0, so the event fires when the ray catches it
    s = solve_ivp(rhs, [T_e, T_sing - 1e-9], [R0], events=[catch, end], rtol=1e-11, atol=1e-12)
    if len(s.t_events[0]) == 0:
        return None
    return s.t_events[0][0], s.y_events[0][0][0]


print("  T_emit   arrives at r/M   measured dT_emit/dtau   formula 1/(1+sqrt(2/r))")
ok = True
d = 1e-4
for T_e in (10.0, 200.0, 380.0, 390.0, 390.7, 391.0, 391.2, 391.3):
    a1 = trace(T_e)
    a2 = trace(T_e + d)
    ok &= a1 is not None and a2 is not None
    if a1 is None or a2 is None:
        continue
    meas = d / (a2[0] - a1[0])
    form = 1 / (1 + np.sqrt(2 / a1[1]))
    print(f"  {T_e:7.1f}  {a1[1]:12.5f}   {meas:14.6f}        {form:14.6f}")
    ok &= abs(meas / form - 1) < 2e-3
R.check("Ray-traced clock rate agrees with dv/dtau = 1/(E+sqrt(E^2-f)), through the horizon", ok)

lo, hi = 0.0, T_sing + 200
for _ in range(60):
    mid = (lo + hi) / 2
    if trace(mid) is not None:
        lo = mid
    else:
        hi = mid
T_last = lo
dv_form = quad(lambda r: 1 / (np.sqrt(2 / r) * (1 + np.sqrt(2 / r))), 0, R0, limit=400)[0]
print(f"  Last ray that catches the faller before r=0: emitted at T = {T_last:.5f} M after release at r=100M")
print(f"  EF integral of dv along the worldline: {dv_form:.5f} M; proper time of the whole fall: {T_sing:.5f} M")
R.check("Bisection on traced rays reproduces the EF last-signal time", abs(T_last - dv_form) < 2e-3)

print("  Non-radial light: extra v accumulated between r=100M and r=0 (impact parameter b, in M)")
vals = []
for b in (0.0, 1.0, 3.0, 5.0):
    def dvdr(r):
        w = np.sqrt(1 - f(r) * b**2 / r**2)
        return (b**2 / r**2) / (w * (1 + w))
    g = quad(dvdr, 0, R0, limit=400)[0]
    vals.append(g)
    print(f"   b = {b:3.1f}: extra v = {g:.5f}")
R.check("Radial rays are the fastest: any other ray must leave earlier to arrive", vals[0] == 0 and all(np.diff(vals) > 0))
R.note("Hence v at the moment of impact bounds ALL signals the notebook can ever receive.")


# ------------------------------------------------------------------ 3. energy table
def dv_tau_in(E):
    """(delta v, delta tau) from horizon to singularity, using s = sqrt(E^2-f)."""
    k = 1 - E**2

    def dv(s):
        return 4 / ((s**2 + k) ** 2 * (E + s))

    def dt(s):
        return 4 / (s**2 + k) ** 2
    pts = [p for p in (E, 10 * E, 100 * E, 1.0) if E < p < 1e3]
    a = quad(dv, E, 1e3, points=sorted(set(pts)) or None, limit=800)[0] if E < 1e3 else 0
    b = quad(dt, E, 1e3, limit=800)[0]
    return a, b


def dv_tau_out(E):
    """(delta v, delta tau) from rest at the radius where f=E^2 down to the horizon (E<1)."""
    k = 1 - E**2
    a = quad(lambda s: 4 / ((s**2 + k) ** 2 * (E + s)), 0, E, limit=400)[0]
    b = quad(lambda s: 4 / (s**2 + k) ** 2, 0, E, limit=400)[0]
    return a, b


R.section("3. How much outside time can be seen, by energy")
print("  E          released at r/M   d_tau(in)/M  d_v(in)/M    d_tau(out)/M   d_v(out)/M")
for E in (1.0, 0.9, 0.5, 0.1, 1e-2, 1e-4, 1e-10, 1e-30):
    vi, ti = dv_tau_in(E)
    if E < 1:
        vo, to = dv_tau_out(E)
        rr = 2 / (1 - E**2)
        print(f"  {E:9.2e}  {rr:14.6f}  {ti:10.5f}  {vi:10.4f}    {to:12.5g}  {vo:10.5f}")
    else:
        print(f"  {E:9.2e}      from far      {ti:10.5f}  {vi:10.4f}")
vi1, ti1 = dv_tau_in(1.0)
R.check("E=1 closed form 4M(5/6 - ln 2) = 0.5607 M after the horizon", abs(vi1 - 4 * (5 / 6 - np.log(2))) < 1e-8, f"{vi1:.6f}")
vo_small = dv_tau_out(1e-8)[0]
R.check("Release near the horizon gains 4M ln 2 before crossing (Pavlov 2013, eq. 14 text)",
        abs(vo_small - 4 * np.log(2)) < 1e-5, f"{vo_small:.6f}")
Es = np.array([1e-6, 1e-12, 1e-20, 1e-30])
Vs = np.array([dv_tau_in(e)[0] for e in Es])
slope = np.polyfit(np.log(1 / Es), Vs, 1)
print(f"  Fit d_v(in) = {slope[0]:.4f} M ln(1/E) + {slope[1]:.4f} M")
R.check("Growth with smaller E is only logarithmic, slope 4M", abs(slope[0] - 4) < 1e-3)
C_in = slope[1]

# ------------------------------------------------------------------ 4. physical numbers
R.section("4. Physical masses: seen outside time vs evaporation time")
M1 = C.year / (np.pi * t_g(C.M_sun)) * C.M_sun
print(f"  Mass with pi*GM/c^3 = 1 yr: {M1 / C.M_sun:.3e} M_sun; its t_evap (blackbody) = {t_evap(M1) / C.year:.3e} yr")
print("  Object       tau(E=1)     v seen after horizon (E=1)   v seen, released at d=l_P      t_evap        fraction seen")
ok = True
for nm, M in (("1 M_sun", C.M_sun), ("Sgr A*", M_SGRA), ("M87*", M_M87), ("1-yr hole", M1)):
    tg = t_g(M)
    Emin = l_P / (4 * r_g(M))
    vmax = (4 * np.log(1 / Emin) + C_in + 4 * np.log(2)) * tg
    tv = vi1 * tg
    frac = vmax / t_evap(M)
    print(f"  {nm:10s}  {ti1 * tg:9.3e} s  {tv:12.3e} s               {vmax:12.3e} s       {t_evap(M) / C.year:9.2e} yr  {frac:9.1e}")
    ok &= frac < 1e-50
R.check("Even the most favourable release sees < 1e-50 of the hole's lifetime", ok)
tg1 = t_g(M1)
print(f"  1-yr hole: falling from far sees {vi1 * tg1 / 86400:.0f} days of outside time after crossing, in {ti1 * tg1 / 86400:.0f} days of proper time")
print(f"  Mass lost to Hawking radiation during that interval: dM/M = dt/(3 t_evap) = {vi1 * tg1 / (3 * t_evap(M1)):.1e}")
dil = 1e120
print(f"  Hover dilation {dil:.0e} at this mass needs d = 4GM/(c^2 D) = {4 * r_g(M1) / dil:.1e} m = {4 * r_g(M1) / dil / l_P:.1e} l_P")
print(f"  Largest dilation with d >= l_P: {4 * r_g(M1) / l_P:.1e}")
E_to_see = t_evap(M1) / (4 * tg1)
print(f"  To see through all of evaporation by release alone you need ln(1/E) = {E_to_see:.1e}; the Planck-length cut allows ln(1/E) = {np.log(4 * r_g(M1) / l_P):.0f}")
R.check("The student's 1e120 dilation is sub-Planckian", 4 * r_g(M1) / dil < l_P)

# ------------------------------------------------------------------ 5. Vaidya
R.section("5. Evaporating hole: ingoing Vaidya model with exaggerated rate")
M0 = 1.0


def make_m(Tev):
    def m(v):
        return M0 * max(1 - v / Tev, 0.0) ** (1 / 3)

    def dm(v):
        x = max(1 - v / Tev, 1e-300)
        return -(M0 / (3 * Tev)) * x ** (-2 / 3)
    return m, dm


def vaidya_drop(v0, r0, Tev, tmax=1e7):
    m, dm = make_m(Tev)

    def rhs(tau, y):
        v, r, vd, rd = y
        mv = m(v)
        F = 1 - 2 * mv / r
        Fr = 2 * mv / r**2
        Fv = -2 * dm(v) / r
        vdd = -(mv / r**2) * vd**2
        rdd = F * vdd + Fr * rd * vd + 0.5 * Fv * vd**2
        return [vd, rd, vdd, rdd]

    F0 = 1 - 2 * m(v0) / r0
    y0 = [v0, r0, 1 / np.sqrt(F0), 0.0]

    def hit(tau, y):
        return y[1] - 1e-3 * max(m(y[0]), 1e-3)

    hit.terminal = True

    def gone(tau, y):
        return y[0] - Tev * (1 - 1e-9)

    gone.terminal = True
    s = solve_ivp(rhs, [0, tmax], y0, events=[hit, gone], rtol=1e-9, atol=1e-11, method="DOP853")
    outcome = "singularity" if len(s.t_events[0]) else ("hole gone" if len(s.t_events[1]) else "?")
    norm = -(1 - 2 * np.array([m(q) for q in s.y[0]]) / s.y[1]) * s.y[2] ** 2 + 2 * s.y[2] * s.y[3]
    return outcome, s.t[-1], s.y[0], s.y[1], np.abs(norm + 1).max()


oc, tau_s, vv, rr, nrm = vaidya_drop(0.0, 3.0, 1e8)
E3 = np.sqrt(1 - 2 / 3)
vexp = dv_tau_in(E3)[0] + dv_tau_out(E3)[0]
R.check("Geodesic equations keep u.u = -1", nrm < 1e-6, f"max error {nrm:.1e}")
R.check("Negligible evaporation reproduces the static result", oc == "singularity" and abs((vv[-1] - vexp) / vexp) < 1e-3,
        f"v_end {vv[-1]:.4f} vs {vexp:.4f}")
Tev = 2000.0
print(f"  m(v) = M0 (1 - v/T)^(1/3), T = {Tev:.0f} M0 (a real 1 M_sun hole has T = {t_evap(C.M_sun) / t_g(C.M_sun):.1e} M0; this is a cartoon rate)")
print("  v0      m(v0)    r0      outcome       v_end      T - v_end   proper time")
endings = {}
for v0, r0 in ((0, 3.0), (500, None), (1000, None), (1500, None), (1900, None), (1990, None), (1999, None)):
    mm = make_m(Tev)[0](v0)
    rr0 = 3 * mm if r0 is None else r0
    oc, ta, vv, rr, _ = vaidya_drop(float(v0), rr0, Tev)
    endings[v0] = (oc, vv, rr)
    print(f"  {v0:6.0f}  {mm:6.4f}  {rr0:6.3f}  {oc:12s}  {vv[-1]:9.3f}  {Tev - vv[-1]:9.3f}  {ta:9.3f}")
    assert oc == "singularity"
R.check("Every notebook dropped from nearby hits the singularity before the end of evaporation",
        all(e[0] == "singularity" and e[1][-1] < Tev for e in endings.values()))
print("  Now a notebook hovering at fixed r=6 M0 and released late:")
lates = {}
for v0 in (1000.0, 1900.0, 1960.0, 1975.0, 1990.0):
    oc, ta, vv, rr, _ = vaidya_drop(v0, 6.0, Tev, tmax=1e5)
    lates[v0] = (oc, vv, rr)
    print(f"  v0={v0:7.1f}: {oc:12s} at v={vv[-1]:9.3f}, r={rr[-1]:.4f}")
R.check("Released early: singularity. Released late enough: the hole is gone first",
        lates[1000.0][0] == "singularity" and lates[1990.0][0] == "hole gone")


# Event horizon: outgoing null ray that reaches r=0 exactly at v=T, integrated backward
def event_horizon(Tev):
    m, _ = make_m(Tev)
    s0 = 1e-4
    r_start = np.sqrt(2 * (3 / 4) * M0 * Tev ** (-1 / 3) * s0 ** (4 / 3))

    def rhs(v, y):
        return [0.5 * (1 - 2 * m(v) / y[0])]
    sol = solve_ivp(rhs, [Tev - s0, 0.0], [r_start], rtol=1e-9, atol=1e-11, method="LSODA", dense_output=True)
    return sol


eh = event_horizon(Tev)
vtest = np.linspace(0, Tev - 1, 2000)
reh = eh.sol(vtest)[0]
mv = np.array([make_m(Tev)[0](q) for q in vtest])
R.check("Event horizon lies inside the apparent horizon r=2m(v) while the hole shrinks", np.all(reh < 2 * mv))
print(f"  at v=0: r_EH = {eh.sol(0.0)[0]:.5f}, 2m = 2.00000; difference {2 - eh.sol(0.0)[0]:.2e}")


def out_ray(v0, r0):
    m, _ = make_m(Tev)

    def rhs(v, y):
        return [0.5 * (1 - 2 * m(v) / y[0])]

    def hit(v, y):
        return y[0] - 1e-3

    hit.terminal = True

    def far(v, y):
        return y[0] - 100

    far.terminal = True
    s = solve_ivp(rhs, [v0, Tev + 400], [r0], events=[hit, far], rtol=1e-10, atol=1e-12, method="LSODA")
    return "falls" if len(s.t_events[0]) else "escapes"


a = out_ray(100.0, float(eh.sol(100.0)[0]) * (1 + 1e-3))
b = out_ray(100.0, float(eh.sol(100.0)[0]) * (1 - 1e-3))
R.check("Forward check at v=100: ray 0.1% outside r_EH escapes, 0.1% inside falls", a == "escapes" and b == "falls", f"{a}, {b}")

# ------------------------------------------------------------------ 6. inner horizon
R.section("6. The one place the student's intuition is right: inner horizons (Reissner-Nordstrom, M=1)")
Q = 0.8
rp, rm = 1 + np.sqrt(1 - Q**2), 1 - np.sqrt(1 - Q**2)
kap = (rp - rm) / (2 * rm**2)


def fq(r):
    return 1 - 2 / r + Q**2 / r**2


def v_to_inner(E, q, delta, r_start=10.0):
    def eps(r):
        return E - q * Q / r

    def integrand(u):
        r = rm + np.exp(u)
        e = eps(r)
        ff = fq(r)
        s = np.sqrt(e**2 - ff)
        den = e + s if e >= 0 else -ff / (s - e)
        return np.exp(u) / (s * den)

    def tint(u):
        r = rm + np.exp(u)
        return np.exp(u) / np.sqrt(eps(r) ** 2 - fq(r))
    a_ = np.log(delta)
    b_ = np.log(r_start - rm)
    return quad(integrand, a_, b_, limit=800)[0], quad(tint, a_, b_, limit=800)[0]


print(f"  Q={Q}: r+ = {rp:.3f}, r- = {rm:.3f}, kappa_- = {kap:.3f}")
print("  delta=r-r_-   neutral E=1: dv, dtau      charged q=2,E=1.5 (eps(r-)<0): dv, dtau")
neu, cha, chat = [], [], []
deltas = (1e-2, 1e-4, 1e-6, 1e-8, 1e-10)
for de in deltas:
    vn, tn = v_to_inner(1.0, 0.0, de)
    vc, tc = v_to_inner(1.5, 2.0, de)
    neu.append(vn)
    cha.append(vc)
    chat.append(tc)
    print(f"  {de:9.0e}     {vn:10.5f}  {tn:10.5f}         {vc:10.4f}  {tc:10.5f}")
R.check("Neutral positive-energy infaller: v stays finite at r_-", abs(neu[-1] - neu[-2]) < 1e-5)
sl = (cha[-1] - cha[-2]) / np.log(deltas[-2] / deltas[-1])
R.check("Infaller with negative effective energy at r_-: v -> infinity with slope 1/kappa_-", abs(sl * kap - 1) < 0.01,
        f"slope {sl:.5f}, 1/kappa_- = {1 / kap:.5f}")
R.check("...while its proper time stays finite", abs(chat[-1] - chat[-2]) < 1e-5, f"tau = {chat[-1]:.5f}")
print("  Exponential blueshift of outside signals near that branch: dv/dtau ~ 1/(delta), delta ~ exp(-kappa_- v).")
a_k = 0.9
gpp = lambda r: r**2 + a_k**2 + 2 * a_k**2 / r
print(f"  Kerr a=0.9, equator: g_phi_phi at r=+0.01: {gpp(0.01):.3f}, at r=-0.01: {gpp(-0.01):.3f}, at r=-1: {gpp(-1.0):.3f}, at r=-10: {gpp(-10.0):.3f}")
rneg = brentq(gpp, -1.0, -1e-3)
print(f"  closed phi circles are timelike for {rneg:.3f} < r < 0 (behind the ring, r_- = {1 - np.sqrt(1 - a_k**2):.3f}, r+ = {1 + np.sqrt(1 - a_k**2):.3f})")
R.check("Kerr CTCs exist only at r<0: g_phi_phi>0 for all r>0", gpp(0.01) > 0 and gpp(-0.01) < 0)

# ------------------------------------------------------------------ 7. hover
R.section("7. Hovering: exact acceleration, dilation and heat")
print("  r-2M(M)   a*M exact      1/d (near-horizon)   D=1/sqrt(f)     4M/d")
ok = True
for xx in (1e-2, 1e-6, 1e-10):
    rr_ = 2 + xx
    sqf = np.sqrt(xx / (2 + xx))  # sqrt(f) without round-off
    a_ex = 1 / (rr_**2 * sqf)
    dd = quad(lambda w: 2 * np.sqrt(2 + w**2), 0, np.sqrt(xx))[0]  # q = 2 + w^2 removes the endpoint singularity
    D = 1 / sqf
    print(f"  {xx:8.0e}  {a_ex:12.4f}  {1 / dd:14.4f}      {D:12.4f}  {4 / dd:12.4f}")
    ok &= abs(a_ex * dd - 1) < 0.01 and abs(D * dd / 4 - 1) < 0.01
R.check("Near the horizon a = c^2/d and D = 4GM/(c^2 d), with d the proper distance", ok)
print(f"  Unruh/Tolman heat at proper distance d: T = hbar c/(2 pi k d) = {C.hbar * C.c / (2 * np.pi * C.k_B):.3e} m K / d; 300 K at d = {C.hbar * C.c / (2 * np.pi * C.k_B * 300):.2e} m; Planck-distance T = {C.hbar * C.c / (2 * np.pi * C.k_B * l_P):.2e} K")
print("  Hover at a survivable 1 g (exact):")
for nm, M in (("10 M_sun", 10 * C.M_sun), ("M87*", M_M87), ("1-yr hole", M1), ("1e20 M_sun", 1e20 * C.M_sun)):
    Rg = r_g(M)

    def fun(lx):
        rr_ = 2 + np.exp(lx)
        return C.c**2 / Rg / (rr_**2 * np.sqrt(np.exp(lx) / rr_)) - C.g_0

    lx = brentq(fun, -200, 60)
    rr_ = 2 + np.exp(lx)
    print(f"   {nm:10s}: r = {rr_ / 2:.4g} r_s, clock ratio D = {1 / np.sqrt(np.exp(lx) / rr_):.4g}")
Dcent = 100 * C.year
dc = 4 * r_g(M_SGRA) / Dcent
print(f"  'One second = one century' hover (D = {Dcent:.2e}), near Sgr A*: d = {dc:.1f} m, acceleration c^2/d = {C.c**2 / dc:.1e} m/s^2")
print(f"  The same D at only 1 g needs M = D c^4/(4 G g) = {Dcent * C.c**4 / (4 * C.G * C.g_0) / C.M_sun:.1e} M_sun (non-rotating)")
print("  Near the end of evaporation (blackbody photon-only power):")
print("  M(kg)      r_s(m)     T(K)      P(W)       flux at 1 m (W/m^2)   time left")
for M in (6.7e25, 1e20, 1e15, 1.4e14, 1e12):
    P = C.hbar * C.c**6 / (15360 * np.pi * C.G**2 * M**2)
    print(f"  {M:9.2e}  {2 * r_g(M):9.2e}  {T_hawking(M):9.2e}  {P:9.2e}  {P / (4 * np.pi):12.2e}         {t_evap(M) / C.year:.2e} yr")

# ------------------------------------------------------------------ 8. simultaneity
R.section("8. 'At the same moment' has no invariant meaning")
vi0, ti0 = dv_tau_in(1.0)


def r_star(r):
    return r + 2 * np.log(np.abs(r / 2 - 1))


print("  Along the E=1 worldline (v=0 at r=100M):")
print("  r/M        tau       v         t = v - r*      PG time T=tau")
vstart = 0.0


def v_of_r(r):
    return quad(lambda q: 1 / (np.sqrt(2 / q) * (1 + np.sqrt(2 / q))), r, 100, limit=400, points=[2.0] if r < 2 < 100 else None)[0]


t0 = -r_star(100.0)
for r in (100.0, 10.0, 2.1, 2.0001, 1.9999, 1.5, 0.5, 0.001):
    tau_ = (2 / 3) * (100**1.5 - r**1.5) / np.sqrt(2)
    vv_ = v_of_r(r)
    print(f"  {r:7.4f}  {tau_:9.3f}  {vv_:9.3f}  {vv_ - r_star(r):13.3f}   {tau_:9.3f}")
tA = v_of_r(2.0001) - r_star(2.0001)
tB = v_of_r(1.5) - r_star(1.5)
R.check("Schwarzschild t rises to +infinity at the horizon, then DEcreases inside", tA > tB and v_of_r(1.9999) - r_star(1.9999) > tB)
v_end = v_of_r(1e-6)


def spacelike_chain(v_far):
    """Interval along r~0 -> r=1.5M -> keep r, raise v -> out to r=1000M. All segments must be spacelike."""
    segs = [(1e-3, 1e-3, 1.5), (v_far - v_end, 1.5, 1.5), (1.0, 1.5, 1000.0)]  # (dv, r1, r2)
    worst = np.inf
    for dv, r1, r2 in segs:
        for w in np.linspace(0.0005, 0.9995, 2000):
            rr_ = r1 + w * (r2 - r1)
            ds2 = -f(rr_) * dv**2 + 2 * dv * (r2 - r1)
            worst = min(worst, ds2)
    return worst > 0


vfar_vals = (v_end + 1.0, 1e6, 1e100)
oks = [bool(spacelike_chain(q)) for q in vfar_vals]
print(f"  v at impact = {v_end:.3f} M. Spacelike curves exist from the impact event to outside events at v = {vfar_vals[0]:.1f}, 1e6, 1e100 M: {oks}")
R.check("The impact can be joined by a spacelike slice to ANY outside time after v_end", all(oks))

# ------------------------------------------------------------------ 9. far observer
R.section("9. What the outside sees of the falling notebook")
rs_ = np.array([2 + 1e-3, 2 + 1e-5, 2 + 1e-7, 2 + 1e-9])
us = np.array([v_of_r(q) - 2 * r_star(q) for q in rs_])
rates = np.array([(1 + np.sqrt(2 / q)) / f(q) for q in rs_])   # du/dtau
sl = np.polyfit(us, np.log(rates), 1)[0]
print(f"  d ln(du/dtau)/du = {sl:.5f} /M -> e-folding time {1 / sl:.3f} M (expected 4M)")
R.check("Signals from the notebook redshift exponentially with e-fold time 4GM/c^3", abs(1 / sl - 4) < 0.01)
print(f"  e-fold time: 1 M_sun {4 * t_g(C.M_sun):.2e} s; Sgr A* {4 * t_g(M_SGRA):.1f} s; M87* {4 * t_g(M_M87) / 3600:.1f} h")

R.section("10. Causal trapping: outgoing light rays dr/dv = f/2")
for r0 in (1.999, 1.5, 2.001):
    def hit(v, y):
        return y[0] - 1e-3

    hit.terminal = True

    def far(v, y):
        return y[0] - 1000

    far.terminal = True
    s = solve_ivp(lambda v, y: [0.5 * f(y[0])], [0, 1e5], [r0], events=[hit, far], rtol=1e-10, atol=1e-12, method="LSODA")
    print(f"  'outgoing' ray from r={r0}: ends at r={s.y[0, -1]:.3f} after dv={s.t[-1]:.2f} M")
    if r0 < 2:
        R.check(f"Ray from r={r0} moves inward", np.all(np.diff(s.y[0]) < 0))
    else:
        R.check(f"Ray from r={r0} escapes", s.y[0, -1] > 999)

# ------------------------------------------------------------------ figure
R.section("11. Figure")
BLUE, ORANGE, AQUA, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#c9c8c2"
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.8), dpi=100, facecolor="#fcfcfb",
                               gridspec_kw={"width_ratios": [1, 1.25]})
for ax in (ax1, ax2):
    ax.set_facecolor("#fcfcfb")
TP = 120.0                 # plotted evaporation time: faster still, so that the infall paths are visible
mP = make_m(TP)[0]
vg = np.linspace(0, TP * (1 - 1e-7), 4000)
ax1.plot(2 * np.array([mP(q) for q in vg]), vg, color=MUTED, lw=2, ls="--", label="apparent horizon, r = 2m(v)")
ehP = event_horizon(TP)
vgh = np.linspace(0, TP - 1e-3, 2000)
ax1.plot(ehP.sol(vgh)[0], vgh, color=INK, lw=2, label="event horizon")
ax1.plot([0, 0], [0, TP], color=ORANGE, lw=5, solid_capstyle="butt", label="singularity")
first = True
pv = {}
for v0 in (0.0, 25.0, 50.0, 75.0, 100.0, 112.0, 117.0):
    oc, _, vv_, rr_p, _ = vaidya_drop(v0, 3 * mP(v0), TP)
    pv[v0] = vv_[-1]
    ax1.plot(rr_p, vv_, color=BLUE, lw=2, label="dropped from r = 3m(v0)" if first else None)
    first = False
first = True
for v0 in (20.0, 60.0, 90.0, 100.0, 108.0):
    oc, _, vv_, rr_p, _ = vaidya_drop(v0, 6.0, TP, tmax=1e4)
    ax1.plot(rr_p, vv_, color=AQUA, lw=2, label="hovered at r = 6, then released" if first else None)
    first = False
ax1.plot([0, 7], [pv[25.0]] * 2, color=GRID, lw=1.5)
ax1.text(7, pv[25.0] + 1.5, "last outside light the v0 = 25 notebook receives", ha="right", fontsize=9, color=MUTED)
ax1.plot([0, 7], [TP, TP], color=GRID, lw=1.5)
ax1.text(7, TP + 1.5, "hole gone, v = T", ha="right", fontsize=9, color=MUTED)
ax1.set_xlim(-0.1, 7)
ax1.set_ylim(0, TP * 1.08)
ax1.set_xlabel("r / M0", color=INK)
ax1.set_ylabel("v / M0  (ingoing light rays are horizontal)", color=INK)
ax1.set_title("Computed: Vaidya model, T = 120 M0 (cartoon rate)", fontsize=11, color=INK, loc="left")
ax1.legend(loc="center right", frameon=False, fontsize=9, bbox_to_anchor=(1.0, 0.72))
ax1.tick_params(colors=MUTED)
for sp_ in ("top", "right"):
    ax1.spines[sp_].set_visible(False)

# Schematic Penrose diagram of an evaporating hole, drawn by hand from the standard structure
ax2.plot([0, 0], [-2, 1], color=MUTED, lw=1.5)                # regular centre before collapse
ax2.plot([1, 1], [1, 3], color=MUTED, lw=1.5)                 # regular centre after evaporation
ax2.plot([1, 3], [3, 1], color=INK, lw=1.5)                   # future null infinity
ax2.plot([3, 0], [1, -2], color=INK, lw=1.5)                  # past null infinity
ax2.plot([0, 1], [1, 1], color=ORANGE, lw=5, solid_capstyle="butt")
ax2.plot([0, 1], [0, 1], color=INK, lw=2)                     # event horizon
ax2.plot([1, 2], [1, 2], color=GRID, lw=1.5)                  # last Hawking light
ax2.plot([1.5, 0], [-0.5, 1], color=MUTED, lw=2)              # collapsing shell
ta = np.linspace(0, 1, 400)
ox = ta + 1.6 * np.sin(np.pi * ta)                            # outside observer
oy = -2 + 5 * ta
ax2.plot(ox, oy, color=AQUA, lw=2)
ax2.plot([1.35, 0.3], [-0.25, 1.0], color=BLUE, lw=2)         # early notebook (steeper than light)
ax2.plot([1.2, 0.82], [0.45, 1.0], color=BLUE, lw=2)          # later notebook
xe = 0.3
xs_ = np.linspace(xe, 2.4, 800)
ys_ = 1 - (xs_ - xe)
_, xh, yh = min((np.min((ox - xq) ** 2 + (oy - yq) ** 2), xq, yq) for xq, yq in zip(xs_, ys_))
ax2.plot([xe, xh], [1, yh], color=BLUE, lw=1, ls=":")
ax2.plot([xh], [yh], "o", color=BLUE, ms=6)
_, qq = min((np.min((ox - q) ** 2 + (oy - q) ** 2), q) for q in np.linspace(1, 2, 800))
ax2.plot([qq], [qq], "o", color=ORANGE, ms=6)
tags = (("1", 1.42, -0.27), ("2", 1.28, 0.45), ("3", xh + 0.13, yh - 0.1), ("4", qq + 0.13, qq + 0.02),
        ("H", 0.32, 0.17), ("C", 0.84, 0.0), ("O", 1.98, -0.45))
for tx, xt, yt in tags:
    ax2.text(xt, yt, tx, color=INK, fontsize=10, fontweight="bold", ha="center", va="center")
ax2.text(0.5, 1.13, "singularity", color=INK, fontsize=9, ha="center")
ax2.text(1.1, 1.72, "flat space", color=MUTED, fontsize=9, ha="left")
key = ("1  early notebook", "2  later notebook", "3  last outside event notebook 1 can ever see",
       "4  outside observer receives the last Hawking light", "H  event horizon", "C  collapsing star",
       "O  outside observer")
for i, line in enumerate(key):
    ax2.text(3.2, 2.9 - 0.32 * i, line, color=INK, fontsize=9, va="center")
ax2.text(3.2, 0.2, "Light moves at 45 degrees.\nNotebook 1 and notebook 2 end long\nbefore event 4, which neither can see.",
         color=MUTED, fontsize=9, va="top")
ax2.set_xlim(-0.2, 6.7)
ax2.set_ylim(-2.2, 3.25)
ax2.set_aspect("equal")
ax2.axis("off")
ax2.set_title("Schematic Penrose diagram, not to scale", fontsize=11, color=INK, loc="left")
fig.tight_layout()
out = Path(__file__).resolve().parents[1] / "figures" / "notebook_evaporation.png"
out.parent.mkdir(exist_ok=True)
fig.savefig(out, facecolor=fig.get_facecolor())
kb = out.stat().st_size / 1024
print(f"  wrote figures/{out.name}")
R.check("Figure under 300 KB", kb < 300)

R.close()
