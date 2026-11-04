"""E02, E03: why things bounce, and whether contact area decides it.

Simulated, not assumed:
  1. Coefficient of restitution e of a mass on a tension-free spring-damper, against the
     closed form obtained from the separation condition, with an energy audit.
  2. Contact-area test: same material, different contact area and contact law.
  3. Contact time against relaxation time: a standard linear solid, and a putty-like one.
  4. Drop test in air of ball, pen, sponge and pen-shaped sponge: the air-flow damping
     estimate, then a contact simulation including drag.
  5. Rigid-body drops of a pen, a sponge block and a ball at random orientation, with
     sequential impacts at the corners.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import g_0 as g  # noqa: E402

FAILS = []
RHO_AIR, MU_AIR = 1.225, 1.8e-5


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- 1. spring-damper
def e_kelvin_sim(zeta, rtol=1e-10, m=1.0, k=1.0, v0=1.0):
    """Mass hits a spring-damper that cannot pull. Returns e and the energy audit error."""
    w = np.sqrt(k / m)
    c = 2 * zeta * m * w

    def rhs(t, y):
        d, v, q = y  # compression, compression rate, heat in the damper
        return [v, -(k * d + c * v) / m, c * v**2]

    def lift_off(t, y):
        return k * y[0] + c * y[1]

    lift_off.terminal = True
    lift_off.direction = -1
    # The force starts positive (c*v0 > 0); start just after contact for zeta = 0 safety.
    sol = solve_ivp(rhs, (0, 1e4 / w), [0.0, v0, 0.0], events=lift_off, rtol=rtol,
                    atol=1e-13, method="DOP853", first_step=1e-4 / w)
    d, v, q = sol.y[:, -1]
    e = -v / v0
    audit = 0.5 * m * v0**2 - (0.5 * m * v**2 + 0.5 * k * d**2 + q)
    return e, audit / (0.5 * m * v0**2), d


def e_kelvin_exact(z):
    """Closed form from x(t) = exp(-z w t) sin(wd t) with lift-off at x'' = 0."""
    if z >= 1:
        raise ValueError
    s = np.sqrt(1 - z**2)
    # lift-off phase phi = w_d t*: x''=0 <=> tan(phi) = 2 z s / (2 z^2 - 1) on (0, pi)
    phi = np.arctan2(2 * z * s, 2 * z**2 - 1)
    t = phi / s
    # velocity at t*: x' = exp(-z t)(cos(phi) - z/s sin(phi))
    return -np.exp(-z * t) * (np.cos(phi) - z / s * np.sin(phi))


# ---------------------------------------------------------------- 2. contact law
def e_hunt_crossley(alpha_v, K, n, m=1.0, v0=1.0):
    """F = K d^n (1 + 1.5 alpha d'), with alpha_v = alpha * v0."""
    alpha = alpha_v / v0

    def rhs(t, y):
        d = max(y[0], 0.0)
        return [y[1], -K * d**n * (1 + 1.5 * alpha * y[1]) / m]

    def out(t, y):
        return y[0]

    out.terminal = True
    out.direction = -1
    T = 100 * (m / K) ** (1 / (n + 1))
    sol = solve_ivp(rhs, (1e-9, T), [1e-9 * v0, v0], events=out, rtol=1e-9, atol=1e-12)
    return -sol.y[1, -1] / v0, sol.t[-1]


# ---------------------------------------------------------------- 3. standard linear solid
def e_sls(k_inf, k1, tau, m=1.0, v0=1.0):
    """Spring k_inf in parallel with a Maxwell branch (k1, relaxation time tau)."""

    def rhs(t, y):
        d, v, s = y
        return [v, -(k_inf * d + k1 * s) / m, v - s / tau]

    def lift(t, y):
        return k_inf * y[0] + k1 * y[2]

    lift.terminal = True
    lift.direction = -1
    w = np.sqrt((k_inf + k1) / m)
    T = 400 / w if k_inf == 0 else 4e4 / w
    sol = solve_ivp(rhs, (0, T), [0.0, v0, 0.0], events=lift, method="Radau",
                    rtol=1e-8, atol=1e-11, first_step=1e-5 / w)
    if sol.t_events[0].size == 0:  # the force never returned to zero: the mass never left
        return 0.0, np.inf
    return -sol.y[1, -1] / v0, sol.t[-1]


# ---------------------------------------------------------------- 4. drop test in air
def drop_test(m, k, zeta, area, Cd, h0, peak_only=True):
    """Fall, contact, rebound, peak: air drag in flight, spring-damper in contact."""
    c = 2 * zeta * np.sqrt(k * m)
    b = 0.5 * RHO_AIR * Cd * area

    def rhs(t, y):
        z, v = y
        F = -m * g - b * v * abs(v)
        if z < 0:
            F += max(0.0, -k * z - c * v)
        return [v, F / m]

    def down(t, y):
        return y[0]

    down.terminal = True
    down.direction = -1
    sol = solve_ivp(rhs, (0, 20), [h0, 0.0], events=down, rtol=1e-9, atol=1e-11, method="DOP853",
                    max_step=0.01)
    v_in = sol.y[1, -1]

    def up(t, y):
        return y[0]

    up.terminal = True
    up.direction = 1
    w = np.sqrt(k / m)
    sol2 = solve_ivp(rhs, (0, 5), [-1e-9, v_in], events=up, rtol=1e-9, atol=1e-12,
                     method="LSODA", max_step=0.05 / w)
    v_out = sol2.y[1, -1]
    tc = sol2.t[-1]

    def top(t, y):
        return y[1]

    top.terminal = True
    top.direction = -1
    sol3 = solve_ivp(rhs, (0, 5), [1e-9, v_out], events=top, rtol=1e-9, atol=1e-12, method="DOP853",
                     max_step=0.01)
    peak = sol3.y[0, -1]
    return dict(v_in=-v_in, v_out=v_out, e=v_out / -v_in, tc=tc, peak=peak)


# ---------------------------------------------------------------- 5. rigid polygon drops
def polygon_drop(corners, I_over_m, e, h, theta0, omega0=0.0, T=1.0, max_events=200, n_grid=1200):
    """Frictionless planar rigid body, instantaneous Newton impacts at corners, exact flights.

    corners: body-frame corner offsets (x, z) from the centre of mass.
    h: height of the lowest corner above the floor at release.
    Returns (highest COM rise above its height at first contact, reached in any later flight;
    first-impact COM velocity ratio; energy in rotation after the first impact as a share of
    the incoming energy; number of impacts).
    """
    cx, cz = corners[:, 0], corners[:, 1]

    def heights(z, th):
        return z + cx * np.sin(th) + cz * np.cos(th)

    def levers(th):
        return cx * np.cos(th) - cz * np.sin(th)  # d(height)/d(theta)

    th = theta0
    z = h - np.min(cx * np.sin(th) + cz * np.cos(th))
    vz, om, t = 0.0, omega0, 0.0
    z_first_contact = None
    peak = -np.inf
    first_ratio, rot_share = None, None
    events = 0
    grid = np.linspace(0, 1, n_grid + 1)[1:]
    while t < T and events < max_events:
        span = T - t
        tau = grid * span
        zt = z + vz * tau - 0.5 * g * tau**2
        tht = th + om * tau
        hts = zt[:, None] + np.outer(np.sin(tht), cx) + np.outer(np.cos(tht), cz)
        below = hts < 0
        if not below.any():
            if events and vz > 0:
                peak = max(peak, z + vz**2 / (2 * g))
            break
        idx = np.argmax(below.any(axis=1))
        corner = np.argmax(below[idx])
        lo = tau[idx - 1] if idx > 0 else 0.0
        hi = tau[idx]

        def h_c(s):
            return (z + vz * s - 0.5 * g * s**2 + cx[corner] * np.sin(th + om * s)
                    + cz[corner] * np.cos(th + om * s))

        s = brentq(h_c, lo, hi, xtol=1e-12) if h_c(lo) > 0 else lo
        if events and vz > 0 and vz / g < s:  # apex inside this flight
            peak = max(peak, z + vz**2 / (2 * g))
        z = z + vz * s - 0.5 * g * s**2
        vz = vz - g * s
        th += om * s
        t += s
        lev = levers(th)[corner]
        vc = vz + om * lev
        if vc >= 0:
            break
        if z_first_contact is None:
            z_first_contact = z
            v_before = vz
            ke_before = 0.5 * vz**2 + 0.5 * I_over_m * om**2
        J = -(1 + e) * vc / (1 + lev**2 / I_over_m)
        vz += J
        om += J * lev / I_over_m
        if first_ratio is None:
            first_ratio = vz / -v_before
            rot_share = (0.5 * I_over_m * om**2) / ke_before
        events += 1
    return max(peak - z_first_contact, 0.0), first_ratio, rot_share, events


def rectangle(a, b):
    """a: long side, b: short side."""
    return np.array([[a / 2, b / 2], [a / 2, -b / 2], [-a / 2, b / 2], [-a / 2, -b / 2]]), (a**2 + b**2) / 12


def main():
    print("1. Spring-damper that cannot pull: e versus damping ratio")
    print("     zeta    e(sim)   e(closed form)  naive exp(-pi z/s)  rebound height e^2   energy audit   body still compressed at lift-off (units v0/w)")
    ok_form, ok_audit = True, True
    for z in [0.001, 0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]:
        e, audit, d = e_kelvin_sim(z)
        ex = e_kelvin_exact(z)
        naive = np.exp(-np.pi * z / np.sqrt(1 - z**2))
        print(f"     {z:5.3f}  {e:7.4f}   {ex:7.4f}        {naive:7.4f}            {e**2:7.4f}         {audit:+.1e}     {d:.3f}")
        ok_form &= abs(e - ex) < 1e-6
        ok_audit &= abs(audit) < 1e-6
    for z in [1.0, 2.0, 5.0]:
        e, audit, d = e_kelvin_sim(z)
        print(f"     {z:5.1f}  {e:7.4f}   (overdamped)                        {e**2:7.4f}         {audit:+.1e}     {d:.3f}")
        ok_audit &= abs(audit) < 1e-6
    check("simulated e matches the closed form for zeta < 1", ok_form)
    check("energy audit closes to 1e-6", ok_audit)
    e1, _, _ = e_kelvin_sim(0.2, rtol=1e-7)
    e2, _, _ = e_kelvin_sim(0.2, rtol=1e-11)
    check("result converged against integrator tolerance", abs(e1 - e2) < 1e-5, f"|de| = {abs(e1 - e2):.1e}")
    es = [e_kelvin_sim(z)[0] for z in (0.01, 0.1, 0.3, 0.7, 1.0, 2.0)]
    check("e falls monotonically with damping", all(a > b for a, b in zip(es, es[1:])))

    print("\n2. Contact-area test")
    print("   (a) Fixed material and mass; contact area scales both k and c (Kelvin solid), zeta ~ sqrt(A):")
    base_z = 0.03
    for A in (0.01, 0.1, 1, 10, 100):
        z = base_z * np.sqrt(A)
        e, _, _ = e_kelvin_sim(z)
        print(f"       area x{A:<6g} zeta {z:6.4f} e {e:.3f}, contact time x{1 / np.sqrt(A):.2f}")
    print("   (b) Fixed material, loss proportional to elastic stress (Hunt-Crossley); area sets K and n:")
    results = []
    for K, n in ((1.0, 1.0), (100.0, 1.0), (1e4, 1.0), (1.0, 1.5), (100.0, 1.5), (1e4, 1.5)):
        e, tc = e_hunt_crossley(0.1, K, n)
        results.append(e)
        print(f"       K {K:8g}  n {n:.1f} (1 = flat pad, 1.5 = Hertz ball): e {e:.4f}, contact time {tc:.3f}")
    spread = max(results) - min(results)
    check("(b) e does not depend on area when loss follows elastic stress", spread < 1e-3,
          f"spread {spread:.1e}")
    e_a_small = e_kelvin_sim(base_z * 0.1)[0]
    e_a_big = e_kelvin_sim(base_z * 10)[0]
    check("(a) smaller contact bounces better in a Kelvin solid (the student's direction)",
          e_a_small > e_a_big, f"e {e_a_small:.3f} vs {e_a_big:.3f}")
    print("   (c) Same area, different material: see part 4, where only the material changes")

    print("\n3. Contact time against relaxation time (standard linear solid, k_inf = k1 = 1)")
    print("     tau/t_contact   e        | putty-like solid (k_inf = 0): tau*w   e")
    taus = np.logspace(-3, 3, 13)
    e_tau, de_tau = [], []
    e_p, de_p = [], []
    for tau in taus:
        e, tc = e_sls(1.0, 1.0, tau)
        ep, _ = e_sls(0.0, 1.0, tau)
        e_tau.append(e); de_tau.append(tau / tc)
        e_p.append(ep); de_p.append(tau)
        print(f"     {tau / tc:12.4f}  {e:6.3f}   | {tau:29.4f}  {ep:6.3f}")
    e_tau = np.array(e_tau)
    i_min = np.argmin(e_tau)
    print(f"     minimum e = {e_tau[i_min]:.3f} at tau/t_contact = {de_tau[i_min]:.2f}")
    check("elastic in both limits, lossy when contact time ~ relaxation time",
          e_tau[0] > 0.9 and e_tau[-1] > 0.9 and e_tau[i_min] < 0.75 and 0.05 < de_tau[i_min] < 5)
    check("putty-like solid: no lift-off when relaxation is quick, bounces when contact is fast",
          e_p[0] == 0.0 and e_p[-1] > 0.9, f"e {e_p[0]:.2f} -> {e_p[-1]:.2f}")
    # Derived: the force obeys m F'' + (m/tau) F' + k F = 0, so it reaches zero only if tau*w > 1/2.
    below, above = e_sls(0.0, 1.0, 0.45)[0], e_sls(0.0, 1.0, 0.55)[0]
    print(f"     putty threshold, derived tau*w = 1/2: e(0.45) = {below:.3f}, e(0.55) = {above:.3f}")
    check("putty lift-off threshold at tau*w = 1/2", below == 0.0 and above > 0.0)

    print("\n4. Drop from 1 m in air (inputs are assumed typical values, see README)")
    # Air-flow damping of an open-cell foam under uniform strain rate, air leaving over path L:
    # Darcy: p_avg = mu v L / (3 kappa); zeta_air = mu L / (6 kappa sqrt(E rho)).
    E_f, rho_f = 30e3, 30.0
    print("   foam air-flow damping estimate, E = 30 kPa, rho = 30 kg/m^3, L = air escape path:")
    for kappa in (1e-9, 3e-10, 1e-10):
        for L in (0.005, 0.03):
            za = MU_AIR * L / (6 * kappa * np.sqrt(E_f * rho_f))
            tp = MU_AIR * L**2 / (kappa * E_f)
            print(f"     permeability {kappa:.0e} m^2, L {L * 100:.1f} cm: zeta_air {za:.3f},"
                  f" poroelastic relaxation {1e3 * tp:.3f} ms")
    zeta_h = 0.1  # polymer skeleton hysteresis, assumed
    kappa = 3e-10
    z_block = zeta_h + MU_AIR * 0.03 / (6 * kappa * np.sqrt(E_f * rho_f))
    z_penspg = zeta_h + MU_AIR * 0.005 / (6 * kappa * np.sqrt(E_f * rho_f))

    objects = [
        # name, mass, k, zeta, frontal area, Cd
        ("rubber ball, r 2.5 cm", 1100 * 4 / 3 * np.pi * 0.025**3, 2e5, 0.03, np.pi * 0.025**2, 0.5),
        ("pen on its tip (rigid plastic)", 0.010, 1e5, 0.06, np.pi * 0.005**2, 1.0),
        ("sponge block flat, 10x6x3 cm", rho_f * 10e-2 * 6e-2 * 3e-2, E_f * 60e-4 / 0.03, z_block, 60e-4, 1.1),
        ("pen-shaped sponge on its tip", rho_f * np.pi * 0.005**2 * 0.14, E_f * np.pi * 0.005**2 / 0.14,
         z_penspg, np.pi * 0.005**2, 1.0),
        ("pen-shaped sponge falling flat", rho_f * np.pi * 0.005**2 * 0.14, E_f * 1.4e-3 / 0.01,
         z_penspg, 1.4e-3, 1.1),
    ]
    print("     object                            m(g)     k(N/m)   zeta   v_in   e_contact  contact(ms)  rebound (cm)  no-drag e^2 h (cm)")
    res = {}
    for name, m, k, z, A, Cd in objects:
        r = drop_test(m, k, z, A, Cd, 1.0)
        res[name] = r
        print(f"     {name:32s} {1e3 * m:6.2f}  {k:9.0f}  {z:5.3f}  {r['v_in']:5.2f}   {r['e']:6.3f}   {1e3 * r['tc']:8.2f}    {100 * r['peak']:7.1f}     {100 * e_kelvin_exact(z) ** 2:7.1f}")
    check("rubber ball rebounds far higher than the flat sponge",
          res["rubber ball, r 2.5 cm"]["peak"] > 2 * res["sponge block flat, 10x6x3 cm"]["peak"])
    check("same pen geometry: rigid pen rebounds higher than pen-shaped sponge on its tip",
          res["pen on its tip (rigid plastic)"]["peak"] > 1.5 * res["pen-shaped sponge on its tip"]["peak"])
    # Can a foam pen stay straight when it lands on its tip? Euler load, taken generously
    # as pin-pin (L_eff = L).
    r_tip, L_pen = 0.005, 0.14
    I_sec = np.pi * r_tip**4 / 4
    F_cr = np.pi**2 * E_f * I_sec / L_pen**2
    r_sp = res["pen-shaped sponge on its tip"]
    m_sp, k_sp = objects[3][1], objects[3][2]
    F_peak = k_sp * r_sp["v_in"] * np.sqrt(m_sp / k_sp)
    print(f"     foam pen on its tip: Euler load {1e3 * F_cr:.1f} mN, own weight {1e3 * m_sp * g:.1f} mN,"
          f" elastic peak force {1e3 * F_peak:.0f} mN = {F_peak / F_cr:.0f} x Euler load")
    check("foam pen must buckle on landing, so the straight-column rebound above overstates it",
          F_peak > 5 * F_cr)
    E_plastic = 3e9  # assumed; hollow barrel ignored, solid rod of the same diameter
    print(f"     solid plastic pen for comparison: Euler load {np.pi**2 * E_plastic * I_sec / L_pen**2:.0f} N")
    # pure effect of drag
    v_flat = res["pen-shaped sponge falling flat"]["v_in"]
    print(f"     pen-shaped sponge falling flat hits at only {v_flat:.2f} m/s after 1 m;"
          f" free fall would give {np.sqrt(2 * g):.2f} m/s")

    print("\n5. Rigid-body drop from 0.5 m at random orientation, frictionless floor")
    # formula check: rod at angle, COM rebound ratio after one impact
    corners, Iom = rectangle(0.14, 1e-6)
    ok = True
    for deg in (0, 10, 30, 60, 90):
        # the formula from impulse theory is compared against the simulator's first impact
        th = np.radians(deg)
        # corner at the +x end is lowest for negative theta in this convention
        _, ratio, rs, _ = polygon_drop(corners, Iom, 0.8, 0.5, -th if deg else -1e-9)
        pred = -1 + 1.8 / (1 + 3 * np.cos(th) ** 2)
        print(f"     thin rod tilted {deg:2d} deg from flat, e = 0.8: COM velocity after first hit"
              f" {ratio:+.3f} v, derived {pred:+.3f} v, rotation takes {100 * rs:.0f} % of the energy")
        ok &= abs(ratio - pred) < 1e-3
    check("simulated first impact matches the derived rod formula", ok)

    rng = np.random.default_rng(7)
    N = 300
    cases = [
        ("pen, 14 x 1 cm, e_tip = 0.7", *rectangle(0.14, 0.01), 0.7),
        ("pen, 14 x 1 cm, e_tip = 0.5", *rectangle(0.14, 0.01), 0.5),
        ("pen-shaped sponge, e = 0.5 (no drag)", *rectangle(0.14, 0.01), 0.5),
        ("pen-shaped sponge, e = 0.3 (no drag)", *rectangle(0.14, 0.01), 0.3),
        ("sponge block 10 x 3 cm, e = 0.5", *rectangle(0.10, 0.03), 0.5),
        ("sponge block 10 x 3 cm, e = 0.3", *rectangle(0.10, 0.03), 0.3),
    ]
    # Time-grid check for the impact finder, same initial conditions.
    a = [polygon_drop(*rectangle(0.14, 0.01), 0.7, 0.5, t0, 0.0, n_grid=1200)[0] for t0 in np.linspace(0.1, 3.0, 60)]
    b = [polygon_drop(*rectangle(0.14, 0.01), 0.7, 0.5, t0, 0.0, n_grid=6000)[0] for t0 in np.linspace(0.1, 3.0, 60)]
    print(f"     impact-finder grid check, 60 drops: median rise {np.median(a):.4f} m vs {np.median(b):.4f} m on a 5x finer grid")
    check("Monte Carlo median not sensitive to the impact-finder grid", abs(np.median(a) - np.median(b)) < 0.01)
    print("     case                                      median COM rise  P(rise > 10% of drop)  P(first hit sends COM up)  median rotation share")
    table = {}
    h0 = 0.5
    for name, corners, Iom, e in cases:
        peaks, ups, rots = [], [], []
        for _ in range(N):
            th0 = rng.uniform(0, 2 * np.pi)
            om0 = rng.normal(0, 3.0)
            p, ratio, rs, _ = polygon_drop(corners, Iom, e, h0, th0, om0)
            peaks.append(p / h0); ups.append(ratio > 0); rots.append(rs)
        peaks = np.array(peaks)
        table[name] = (np.median(peaks), np.mean(peaks > 0.1))
        print(f"     {name:40s}  {100 * np.median(peaks):8.1f} %       {100 * np.mean(peaks > 0.1):7.1f} %            {100 * np.mean(ups):7.1f} %              {100 * np.median(rots):5.1f} %")
    print(f"     ball with e = 0.7 for comparison: COM rise always {100 * 0.7**2:.0f} % of the drop height,"
          f" rotation share 0 % (the normal passes through the centre)")
    check("pen with e = 0.7 rises less than a ball of the same e",
          table["pen, 14 x 1 cm, e_tip = 0.7"][0] < 0.7**2)
    check("the COM-rise chance grows with tip e",
          table["pen, 14 x 1 cm, e_tip = 0.7"][1] > table["pen-shaped sponge, e = 0.3 (no drag)"][1])

    print()
    if FAILS:
        print("FAILED:", FAILS)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
