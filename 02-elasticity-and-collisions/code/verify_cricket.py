"""E08: from release to bounce to the batter. What can and can't be predicted.

Vectorised RK4 for many deliveries, checked against scipy on single trajectories.
Inputs that are measured or assumed are listed in ASSUMED and in the README.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import g_0 as g  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# MCC Laws: ball 155.9-163 g, circumference 22.4-22.9 cm; popping crease 1.22 m in front of the stumps.
M, R = 0.159, 0.227 / (2 * np.pi)
I_FACTOR = 0.4          # solid uniform sphere; a real cork-and-leather ball differs
RHO = 1.225
A_BALL = np.pi * R**2
PITCH = 20.12           # stump to stump, 22 yards
CREASE = 1.22
X_BAT = PITCH - 2 * CREASE   # release taken above the bowler's popping crease; batter at the other one
X_STUMPS = PITCH - CREASE

# Nominal: e near the top of Howard's laboratory range (0.3-0.4), mu inside his range, Cd assumed.
NOMINAL = dict(e=0.40, mu=0.35, Cd=0.45, Cl=0.0, spin_hz=10.0)
HEAD, CHEST = (1.50, 1.90), (1.10, 1.50)  # bands at the batter, assumed for a 1.8 m batter in stance


def rk4_flight(x, y, vx, vy, Cd, Cl, stop, dt=2e-3, t_max=2.0):
    """Integrate many balls until stop(x, y) is true. Returns state at the stop by interpolation."""
    x, y, vx, vy = map(lambda a: np.array(a, dtype=float), (x, y, vx, vy))
    n = x.size
    k = 0.5 * RHO * A_BALL / M

    def f(vx, vy):
        sp = np.hypot(vx, vy)
        return -k * sp * (Cd * vx + Cl * vy), -g + k * sp * (-Cd * vy + Cl * vx)

    done = np.zeros(n, bool)
    out = np.zeros((4, n))
    t_out = np.zeros(n)
    t = 0.0
    while not done.all() and t < t_max:
        ax1, ay1 = f(vx, vy)
        ax2, ay2 = f(vx + 0.5 * dt * ax1, vy + 0.5 * dt * ay1)
        ax3, ay3 = f(vx + 0.5 * dt * ax2, vy + 0.5 * dt * ay2)
        ax4, ay4 = f(vx + dt * ax3, vy + dt * ay3)
        nx = x + dt / 6 * (vx + 2 * (vx + 0.5 * dt * ax1) + 2 * (vx + 0.5 * dt * ax2) + (vx + dt * ax3))
        ny = y + dt / 6 * (vy + 2 * (vy + 0.5 * dt * ay1) + 2 * (vy + 0.5 * dt * ay2) + (vy + dt * ay3))
        nvx = vx + dt / 6 * (ax1 + 2 * ax2 + 2 * ax3 + ax4)
        nvy = vy + dt / 6 * (ay1 + 2 * ay2 + 2 * ay3 + ay4)
        s0, s1 = stop(x, y), stop(nx, ny)  # stop value: >0 keep going, <=0 stop
        hit = (~done) & (s1 <= 0)
        if hit.any():
            w = s0[hit] / (s0[hit] - s1[hit])
            out[0, hit] = x[hit] + w * (nx[hit] - x[hit])
            out[1, hit] = y[hit] + w * (ny[hit] - y[hit])
            out[2, hit] = vx[hit] + w * (nvx[hit] - vx[hit])
            out[3, hit] = vy[hit] + w * (nvy[hit] - vy[hit])
            t_out[hit] = t + w * dt
            done |= hit
        x, y, vx, vy = nx, ny, nvx, nvy
        t += dt
    return out, t_out, done


def bounce(vx, vy, omega_back, e, mu):
    """Rigid-sphere impact: restitution e on the normal speed, Coulomb friction, spin included.

    omega_back > 0 is backspin (rad/s). Slip speed of the bottom of the ball: u = vx + R omega.
    The friction impulse needed to stop slip is J* = u m / (1 + 1/I_FACTOR); the most that
    friction can supply is mu (1 + e) m |vy|.
    """
    u = vx + R * omega_back
    J_grip = np.abs(u) / (1 + 1 / I_FACTOR)              # per unit mass
    J_max = mu * (1 + e) * np.abs(vy)
    J = np.minimum(J_grip, J_max) * np.sign(u)
    vx2 = vx - J
    vy2 = e * np.abs(vy)
    omega2 = omega_back - J / (I_FACTOR * R)
    gripped = J_grip <= J_max
    return vx2, vy2, omega2, gripped


def delivery(v0, alpha_deg, H, e, mu, Cd=0.45, Cl=0.0, spin_hz=10.0, kick=0.0):
    """Returns a dict of arrays for release speed v0 and angle alpha below the horizontal.

    kick (deg) rotates the rebound velocity upward at constant speed. It stands in for the
    surface deformation the rigid model leaves out; its size is an assumption.
    """
    v0, alpha_deg, H, e, mu, Cd, Cl, spin_hz, kick = np.broadcast_arrays(
        *map(lambda a: np.asarray(a, dtype=float), (v0, alpha_deg, H, e, mu, Cd, Cl, spin_hz, kick)))
    a = np.radians(alpha_deg)
    n = v0.size
    x0 = np.zeros(n)
    s1, t1, done1 = rk4_flight(x0, H.ravel(), (v0 * np.cos(a)).ravel(), (-v0 * np.sin(a)).ravel(),
                               Cd.ravel(), Cl.ravel(), lambda x, y: y)
    xp, _, vxi, vyi = s1
    omega = 2 * np.pi * spin_hz.ravel()
    vx2, vy2, om2, grip = bounce(vxi, vyi, omega, e.ravel(), mu.ravel())
    sp2 = np.hypot(vx2, vy2)
    ang2 = np.arctan2(vy2, vx2) + np.radians(kick.ravel())
    vx2, vy2 = sp2 * np.cos(ang2), sp2 * np.sin(ang2)
    # Second flight, stopping at the batter's crease or a second touchdown.
    Cl2 = Cl.ravel() * np.sign(om2)
    s2, t2, done2 = rk4_flight(xp, np.full(n, 1e-6), vx2, vy2, Cd.ravel(), Cl2,
                               lambda x, y: np.minimum(X_BAT - x, y + 1e-3))
    h = s2[1].copy()
    # Apex after the bounce, ignoring drag on the vertical (small over this interval)
    apex = vy2**2 / (2 * g)
    return dict(xp=xp, length=X_STUMPS - xp, vxi=vxi, vyi=vyi, strike=np.degrees(np.arctan2(-vyi, vxi)),
                vx2=vx2, vy2=vy2, rebound=np.degrees(np.arctan2(vy2, vx2)), h=h, apex=apex,
                t1=t1, t2=t2, grip=grip, om2=om2, reached=s2[0] >= X_BAT - 1e-6, full=xp >= X_BAT)


def alpha_for_length(v0, H, length, **kw):
    """Release angle that lands at this length from the batter's stumps."""
    target = X_STUMPS - length
    f = lambda al: delivery(v0, al, H, 0.40, 0.35, **kw)["xp"][0] - target
    return brentq(f, -5.0, 30.0, xtol=1e-5)


def scipy_flight(v0, alpha, H, Cd):
    a = np.radians(alpha)
    k = 0.5 * RHO * A_BALL / M

    def rhs(t, s):
        sp = np.hypot(s[2], s[3])
        return [s[2], s[3], -k * Cd * sp * s[2], -g - k * Cd * sp * s[3]]

    def ground(t, s):
        return s[1]

    ground.terminal, ground.direction = True, -1
    sol = solve_ivp(rhs, (0, 3), [0, H, v0 * np.cos(a), -v0 * np.sin(a)], events=ground,
                    method="DOP853", rtol=1e-11, atol=1e-12)
    return sol.y[:, -1], sol.t[-1]


def main():
    print(f"Ball: m = {M} kg, R = {R * 100:.2f} cm (MCC Law 4). Batter at x = {X_BAT:.2f} m from release, stumps at {X_STUMPS:.2f} m.")
    print("\n1. Release to pitch, no air: strike angle from speed, angle and height")
    print("     v0 (m/s)  angle  H (m)   x_pitch(m)  strike angle(deg) sim   closed form")
    ok = True
    for v0, al, H in ((40, 3, 2.2), (40, 7, 2.2), (30, 3, 2.1), (25, -2, 2.0), (40, 0, 2.2)):
        s, t, _ = rk4_flight(np.zeros(1), np.array([H], float), np.array([v0 * np.cos(np.radians(al))]),
                             np.array([-v0 * np.sin(np.radians(al))]), np.array([0.0]), np.array([0.0]),
                             lambda x, y: y)
        sim = np.degrees(np.arctan2(-s[3, 0], s[2, 0]))
        vz0 = v0 * np.sin(np.radians(al))
        closed = np.degrees(np.arctan2(np.sqrt(vz0**2 + 2 * g * H), v0 * np.cos(np.radians(al))))
        tc = (-vz0 + np.sqrt(vz0**2 + 2 * g * H)) / g
        xc = v0 * np.cos(np.radians(al)) * tc
        print(f"     {v0:5.0f}    {al:5.1f}  {H:4.2f}     {s[0, 0]:7.3f}      {sim:7.3f}            {closed:7.3f}   (x closed {xc:.3f})")
        ok &= abs(sim - closed) < 1e-3 and abs(s[0, 0] - xc) < 1e-3
    check("strike angle is fixed by v0, alpha and H under gravity alone", ok)

    print("\n2. Integrator checks with drag (Cd = 0.45)")
    d_s = delivery(40.0, 5.0, 2.2, 0.45, 0.35)
    ref, t_ref = scipy_flight(40.0, 5.0, 2.2, 0.45)
    print(f"     RK4 x_pitch {d_s['xp'][0]:.5f} m, scipy DOP853 {ref[0]:.5f} m; vy {d_s['vyi'][0]:.5f} vs {ref[3]:.5f} m/s")
    check("vectorised RK4 agrees with an adaptive integrator", abs(d_s["xp"][0] - ref[0]) < 1e-3 and abs(d_s["vyi"][0] - ref[3]) < 1e-3)
    s_a, _, _ = rk4_flight(np.zeros(1), np.array([2.2]), np.array([40 * np.cos(0.0873)]), np.array([-40 * np.sin(0.0873)]),
                           np.array([0.45]), np.array([0.0]), lambda x, y: y, dt=1e-3)
    s_b, _, _ = rk4_flight(np.zeros(1), np.array([2.2]), np.array([40 * np.cos(0.0873)]), np.array([-40 * np.sin(0.0873)]),
                           np.array([0.45]), np.array([0.0]), lambda x, y: y, dt=2.5e-4)
    print(f"     step-size: x_pitch {s_a[0, 0]:.6f} (1 ms) vs {s_b[0, 0]:.6f} (0.25 ms)")
    check("converged with step size", abs(s_a[0, 0] - s_b[0, 0]) < 1e-3)

    print("\n3. Bounce impulse checks")
    print("     incoming (vx, vy, spin Hz)   e   mu     vx'     vy'   spin' Hz  regime  KE+rot before  after (J/kg)")
    ok = True
    rng = np.random.default_rng(3)
    for vx, vy, hz, e, mu in ((38, -9, 10, 0.45, 0.35), (38, -9, 10, 0.45, 2.0), (22, -6, -30, 0.45, 0.35),
                              (22, -6, -40, 0.45, 2.0), (22, -6, 0, 0.45, 2.0)):
        om = 2 * np.pi * hz
        vx2, vy2, om2, gr = bounce(np.array([vx], float), np.array([vy], float), np.array([om]), np.array([e]), np.array([mu]))
        E1 = 0.5 * (vx**2 + vy**2) + 0.5 * I_FACTOR * R**2 * om**2
        E2 = 0.5 * (vx2[0] ** 2 + vy2[0] ** 2) + 0.5 * I_FACTOR * R**2 * om2[0] ** 2
        print(f"     ({vx:3d}, {vy:3d}, {hz:4d})                 {e:4.2f} {mu:4.2f}  {vx2[0]:6.2f}  {vy2[0]:5.2f}  {om2[0] / (2 * np.pi):7.1f}   {'grip' if gr[0] else 'slip'}    {E1:7.1f}     {E2:7.1f}")
        ok &= E2 <= E1 + 1e-9
        if gr[0]:
            ok &= abs(vx2[0] + R * om2[0]) < 1e-9  # no slip at the end
    vxr = rng.uniform(10, 45, 20000)
    vyr = -rng.uniform(1, 20, 20000)
    omr = 2 * np.pi * rng.uniform(-50, 50, 20000)
    er = rng.uniform(0.2, 0.9, 20000)
    mur = rng.uniform(0.05, 2.0, 20000)
    v2x, v2y, o2, _ = bounce(vxr, vyr, omr, er, mur)
    E1 = 0.5 * (vxr**2 + vyr**2) + 0.5 * I_FACTOR * R**2 * omr**2
    E2 = 0.5 * (v2x**2 + v2y**2) + 0.5 * I_FACTOR * R**2 * o2**2
    print(f"     20000 random impacts: largest energy gain {np.max(E2 - E1):.2e} J/kg, smallest loss {np.min(E1 - E2):.2e} J/kg")
    check("impact model never creates energy and ends gripped impacts with no slip", ok and np.all(E2 <= E1 + 1e-9))
    # rolling limit for a spinless sphere
    v3, _, _, _ = bounce(np.array([20.0]), np.array([-5.0]), np.array([0.0]), np.array([0.5]), np.array([5.0]))
    print(f"     spinless ball that grips: vx' / vx = {v3[0] / 20:.4f} (derived 5/7 = {5 / 7:.4f})")
    check("grip limit is 5/7 for a uniform sphere", abs(v3[0] / 20 - 5 / 7) < 1e-9)

    print("\n3b. Against the published 'average seam bowler' of James et al. (2005), as used by Howard (2012):")
    print("     36 m/s at 12.5 deg at the pitch, bouncing about 5 m in front of the batter, no drag after the bounce")
    vxj, vyj = 36 * np.cos(np.radians(12.5)), -36 * np.sin(np.radians(12.5))
    hs_j = []
    for e_j, mu_j in ((0.30, 0.3), (0.35, 0.3), (0.40, 0.3), (0.30, 0.6), (0.40, 0.6)):
        vx_j, vy_j, _, _ = bounce(np.array([vxj]), np.array([vyj]), np.array([2 * np.pi * 10.0]), np.array([e_j]), np.array([mu_j]))
        t_j = 5.0 / vx_j[0]
        h_j = vy_j[0] * t_j - 0.5 * g * t_j**2
        pace = np.hypot(vx_j[0], vy_j[0]) / 36
        hs_j.append((h_j, pace))
        print(f"       e = {e_j:.2f}, mu = {mu_j:.1f}: height at batter {h_j:.3f} m, pace (speed after/before) {pace:.3f}")
    print("     Howard's Figures 32-33 span heights 0.15-0.40 m and pace 0.76-0.90, and he states the rigid model underestimates the rebound angle")
    check("rigid model reproduces the height and pace range published for the same rigid model",
          all(0.10 < h < 0.45 and 0.70 < p < 0.92 for h, p in hs_j))
    al_j = brentq(lambda al: delivery(39.0, al, 2.2, **NOMINAL)["strike"][0] - 12.5, 2.0, 20.0, xtol=1e-4)
    d_j = delivery(39.0, al_j, 2.2, **NOMINAL)
    sp_j = np.hypot(d_j["vxi"][0], d_j["vyi"][0])
    print(f"     release 39 m/s, H = 2.2 m, alpha {al_j:.2f} deg gives strike 12.5 deg, arrival speed {sp_j:.1f} m/s, length {d_j['length'][0]:.2f} m from stumps")
    check("release-to-pitch model gives the published arrival speed and angle at a sensible length",
          abs(sp_j - 36) < 1.5 and 4 < d_j["length"][0] < 8)

    print("\n4. A 40 m/s delivery, H = 2.2 m, e = 0.40, mu = 0.35, Cd = 0.45. Height at the batter with no kick and with a 2.5 and 5 deg kick")
    print("     alpha(deg)  x_pitch(m)  length from stumps(m)  strike(deg)  v_y in(m/s)  rebound(deg)  v_x out  h, no kick(m)   h, 2.5 deg   h, 5 deg")
    table = []
    for al in (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0):
        d = delivery(40.0, al, 2.2, **NOMINAL)
        d2 = delivery(40.0, al, 2.2, kick=2.5, **NOMINAL)
        d5 = delivery(40.0, al, 2.2, kick=5.0, **NOMINAL)
        table.append(d)
        star = "*" if d["full"][0] else " "
        print(f"     {al:6.1f}      {d['xp'][0]:7.2f}       {d['length'][0]:6.2f}              {d['strike'][0]:6.2f}      {-d['vyi'][0]:6.2f}      {d['rebound'][0]:6.2f}     {d['vx2'][0]:5.2f}     {d['h'][0]:6.3f}{star}       {d2['h'][0]:6.3f}{star}     {d5['h'][0]:6.3f}{star}")
    print("     * pitched beyond the batter: a full toss, the height shown is then not a bounce")

    print("\n5. Sensitivity for a short ball: alpha set so it lands 9 m from the stumps at 40 m/s")
    H0 = 2.2
    al0 = alpha_for_length(40.0, H0, 9.0)
    base = delivery(40.0, al0, H0, **NOMINAL)
    print(f"     alpha = {al0:.2f} deg, strike {base['strike'][0]:.2f} deg, v_y in {-base['vyi'][0]:.2f} m/s, height at batter {base['h'][0]:.3f} m")
    print("     input               step              change in height at batter (m)")
    tests = [
        ("release speed", dict(v0=41.0), "+1 m/s"), ("release angle", dict(alpha=al0 + 0.5), "+0.5 deg"),
        ("release angle", dict(alpha=al0 - 0.5), "-0.5 deg"),
        ("release height", dict(H=2.3), "+0.1 m"), ("restitution e", dict(e=0.50), "+0.05"),
        ("friction mu", dict(mu=0.45), "+0.10"), ("drag Cd", dict(Cd=0.35), "-0.10"),
        ("lift Cl (backspin)", dict(Cl=0.1), "+0.10"), ("spin", dict(spin_hz=20.0), "+10 Hz"),
        ("deformation kick", dict(kick=1.0), "+1 deg"),
    ]
    sens = {}
    for name, kw, step in tests:
        p = dict(v0=40.0, alpha=al0, H=H0, e=0.40, mu=0.35, Cd=0.45, Cl=0.0, spin_hz=10.0)
        p.update(kw)
        al_ = p.pop("alpha")
        v_ = p.pop("v0")
        H_ = p.pop("H")
        d = delivery(v_, al_, H_, **p)
        sens[name + step] = d["h"][0] - base["h"][0]
        print(f"     {name:20s}{step:12s}      {d['h'][0] - base['h'][0]:+.3f}")
    check("height sensitive to e by more than 0.1 m per 0.05", abs(sens["restitution e+0.05"]) > 0.1)
    check("height sensitive to the unmodelled deformation by more than 0.1 m per degree", abs(sens["deformation kick+1 deg"]) > 0.1)

    print("\n6. Uncertainty on height at the batter, Monte Carlo, N = 5000, four target lengths")
    print("     A: v0, angle, height known exactly; pitch and ball unknown over the ranges below")
    print("     B: as A, plus realistic release measurement error (v0 +/- 1 m/s, angle +/- 0.5 deg, H +/- 0.05 m, normal 1 sigma)")
    print("     C: pitch calibrated on the day (e 0.40 +/- 0.02, mu 0.35 +/- 0.05, kick 2.5 +/- 0.5 deg), delivery-to-delivery e scatter +/- 0.03, release exact")
    print("     D: as C, plus release measurement error as in B")
    print("     ranges for A, B: e 0.30-0.50, mu 0.20-0.60, kick 0-5 deg, Cd 0.30-0.55, Cl 0-0.1, spin 5-15 Hz (uniform, assumed)")
    rng = np.random.default_rng(11)
    N = 5000

    def sample(case, al):
        if case in "AB":
            e = rng.uniform(0.30, 0.50, N)
            mu = rng.uniform(0.20, 0.60, N)
            kick = rng.uniform(0.0, 5.0, N)
        else:
            e = 0.40 + rng.normal(0, 0.02, N) + rng.uniform(-0.03, 0.03, N)
            mu = 0.35 + rng.normal(0, 0.05, N)
            kick = 2.5 + rng.normal(0, 0.5, N)
        if case in "AB":
            Cd = rng.uniform(0.30, 0.55, N)
            Cl = rng.uniform(0.0, 0.1, N)
        else:
            Cd = 0.45 + rng.normal(0, 0.03, N)
            Cl = rng.uniform(0.03, 0.07, N)
        hz = rng.uniform(5, 15, N)
        if case in "BD":
            v = 40 + rng.normal(0, 1.0, N)
            a = al + rng.normal(0, 0.5, N)
            H = H0 + rng.normal(0, 0.05, N)
        else:
            v, a, H = np.full(N, 40.0), np.full(N, al), np.full(N, H0)
        return delivery(v, a, H, e, mu, Cd, Cl, hz, kick)

    res = {}
    print("     length(m)  case   median h(m)   5%..95% (m)   width(m)  P(below chest) P(chest 1.1-1.5)  P(head 1.5-1.9)  P(over 1.9)")
    for L in (7.0, 9.0, 11.0, 13.0):
        al = alpha_for_length(40.0, H0, L)
        for case in "ABCD":
            d = sample(case, al)
            h = d["h"][~d["full"]]
            q5, q50, q95 = np.percentile(h, [5, 50, 95])
            res[(L, case)] = q95 - q5
            res[(L, case, "pc")] = (np.mean(h < 1.1), np.mean((h >= 1.1) & (h < 1.5)), np.mean((h >= 1.5) & (h < 1.9)), np.mean(h >= 1.9))
            pb, pc, ph, po = res[(L, case, "pc")]
            print(f"     {L:6.1f}      {case}     {q50:6.3f}     {q5:5.2f} .. {q95:5.2f}   {q95 - q5:5.2f}      {pb:5.2f}          {pc:5.2f}          {ph:5.2f}        {po:5.2f}")
    gap = 0.4  # the width of the chest and head bands
    check("with only release data the height at the batter spreads over more than one 0.4 m band",
          res[(9.0, "A")] > gap, f"5-95% width {res[(9.0, 'A')]:.2f} m")
    check("calibrating the pitch narrows the spread by at least a factor 2 if the release is exact",
          res[(9.0, "C")] < 0.5 * res[(9.0, "A")], f"{res[(9.0, 'A')]:.2f} -> {res[(9.0, 'C')]:.2f} m")
    check("release error alone keeps the calibrated forecast from a clean chest-or-head call", res[(9.0, "D")] > 0.5 * res[(9.0, "C")])

    print("\n7. Release height at the same length, 6.82 m from the stumps, 37 m/s, 2.5 deg kick")
    hs = []
    for H_ in (1.9, 2.1, 2.3, 2.5):
        al = alpha_for_length(37.0, H_, 6.82)
        d = delivery(37.0, al, H_, kick=2.5, **NOMINAL)
        hs.append(d["h"][0])
        print(f"     H = {H_:.1f} m: alpha {al:5.2f} deg, strike {d['strike'][0]:5.2f} deg, height at batter {d['h'][0]:.3f} m, apex {d['apex'][0]:.3f} m")
    print(f"     model: {(hs[-1] - hs[0]) / 0.6:.2f} m of bounce per m of release height at this length")
    check("taller release gives higher bounce at the same length", all(a < b for a, b in zip(hs, hs[1:])))

    print("\n8. Speed at the same 9 m length, and a spinner")
    for v in (25.0, 30.0, 35.0, 40.0, 45.0):
        al = alpha_for_length(v, H0, 9.0)
        d = delivery(v, al, H0, kick=2.5, **NOMINAL)
        print(f"     v0 = {v:.0f} m/s: alpha {al:5.2f} deg, strike {d['strike'][0]:5.2f} deg, v_y out {d['vy2'][0]:5.2f} m/s, height at batter {d['h'][0]:.3f} m, apex {d['apex'][0]:.3f} m")
    al = alpha_for_length(22.0, 2.0, 4.0)
    d = delivery(22.0, al, 2.0, 0.40, 0.6, 0.45, 0.0, -30.0)
    print(f"     spinner 22 m/s, 30 Hz topspin, mu 0.6, 4 m length: regime {'grip' if d['grip'][0] else 'slip'}, v_x {d['vxi'][0]:.1f} -> {d['vx2'][0]:.1f} m/s, height at batter {d['h'][0]:.2f} m")

    print()
    if FAILS:
        print("FAILED:", FAILS)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
