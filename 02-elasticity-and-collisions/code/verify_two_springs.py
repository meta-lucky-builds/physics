"""E07: a ball between two opposed springs. Does it go on forever?

Units: m = k = L0 = 1 unless stated.
  1. Massless springs, one-sided contact: touching, gapped and preloaded cases, each against an
     independently derived period, with a symplectic integrator.
  2. Springs with their own mass: discretised chains. Total energy stays, the ball's does not.
  3. Damping: spring dashpots, dry friction and air drag give three different decay laws.
  4. Sideways stability: static buckling threshold and parametric growth while oscillating,
     Floquet analysis against the full two-dimensional simulation.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_spring import helix_fe  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def verlet_1d(force, potential, x0, T, dt=5e-3):
    n = int(T / dt)
    x, v = x0, 0.0
    a = force(x)
    crossings = []
    xmax = abs(x0)
    E0 = 0.5 * v**2 + potential(x)
    Emax = 0.0
    for i in range(n):
        v += 0.5 * dt * a
        x_new = x + dt * v
        a = force(x_new)
        v += 0.5 * dt * a
        if x > 0 >= x_new:  # downward zero crossing, interpolated
            crossings.append((i + x / (x - x_new)) * dt)
        x = x_new
        if i % 50 == 0:
            Emax = max(Emax, abs(0.5 * v**2 + potential(x) - E0))
    periods = np.diff(crossings)
    return periods, Emax / E0, x, v


def main():
    print("1. Massless springs, one-sided contact (no damping)")
    # a. Touching, no gap, no preload: exactly one spring acts on either side of centre.
    A = 0.3
    F = lambda x: -x
    V = lambda x: 0.5 * x**2
    per, drift, x_end, v_end = verlet_1d(F, V, A, 200 * np.pi)
    print(f"   a. touching: mean period {per.mean():.6f} (2 pi = {2 * np.pi:.6f}), energy drift {drift:.1e} over 100 periods")
    check("touching springs are exactly a linear oscillator", abs(per.mean() / (2 * np.pi) - 1) < 1e-4)
    check("no amplitude loss over 100 periods", drift < 1e-5)

    # b. Gap: springs do not reach the ball at centre.
    print("   b. gap of half-width g: T = 2 pi + 4 g / (A - g) (A = turning point), simulated:")
    ok = True
    for g_, A_ in ((0.05, 0.3), (0.1, 0.3), (0.1, 0.15)):
        Fg = lambda x, g_=g_: (-(x - g_) if x > g_ else (-(x + g_) if x < -g_ else 0.0))
        Vg = lambda x, g_=g_: (0.5 * (abs(x) - g_) ** 2 if abs(x) > g_ else 0.0)
        per, drift, *_ = verlet_1d(Fg, Vg, A_, 200)
        pred = 2 * np.pi + 4 * g_ / (A_ - g_)
        ok &= abs(per.mean() / pred - 1) < 2e-3
        print(f"      g = {g_:.2f}, A = {A_:.2f}: period {per.mean():.4f}, predicted {pred:.4f}, energy drift {drift:.1e}")
    check("gap makes the period depend on amplitude, as derived", ok)

    # c. Preload p: both springs act for |x| < p.
    print("   c. preload p: stiffness 2k inside |x| < p, k outside. Period from the energy integral:")
    ok = True
    for p_, A_ in ((0.1, 0.05), (0.1, 0.2), (0.1, 0.4)):
        Fp = lambda x, p_=p_: (-2 * x if abs(x) < p_ else -(x + np.sign(x) * p_))
        Vp = lambda x, p_=p_: (x**2 if abs(x) < p_ else 0.5 * (abs(x) + p_) ** 2 - p_**2)
        E = Vp(A_)
        # substitute x = A sin(u) to remove the end-point singularity
        integrand = lambda u: A_ * np.cos(u) / np.sqrt(2 * max(E - Vp(A_ * np.sin(u)), 1e-300))
        T_pred = 4 * quad(integrand, 0, np.pi / 2, limit=400, points=[np.arcsin(min(p_ / A_, 1.0))] if p_ < A_ else None)[0]
        per, drift, *_ = verlet_1d(Fp, Vp, A_, 200)
        ok &= abs(per.mean() / T_pred - 1) < 2e-3
        print(f"      p = {p_:.2f}, A = {A_:.2f}: period {per.mean():.4f}, energy integral {T_pred:.4f}")
    check("preloaded system is nonlinear: period from simulation equals the energy integral", ok)

    print("\n2. Springs with their own mass (chains of 12 masses; contact is one-sided)")
    print("   spring mass/ball mass  first period  rule 2 pi sqrt(1 + ms/3)  energy drift   ball's energy share: cycle 1  late mean  late min")

    def chains(ms, T=60.0, N=12, kc=60.0, A=0.3, cfl=0.12):
        mn = ms / N
        kl = float(N)
        wmax = max(2 * np.sqrt(kl / mn), np.sqrt((kc + 2 * kl) / mn))
        dt = min(cfl / wmax, 0.01)
        # chain displacements measured from rest; index 0 next to the anchor, N-1 is the free end
        up = np.linspace(1, N, N) / N * A   # upper chain compressed with the ball at x = A
        lo = np.zeros(N)
        vu, vl = np.zeros(N), np.zeros(N)
        x, v = A, 0.0

        def acc(up, lo, x):
            # upper spring: positive displacement = compression toward the anchor
            du = np.diff(np.r_[0.0, up, 0.0])
            au = kl * (np.r_[du[1:-1], 0.0] - du[:-1])
            dl = np.diff(np.r_[0.0, lo, 0.0])
            al = kl * (np.r_[dl[1:-1], 0.0] - dl[:-1])
            Fu = kc * max(0.0, x - up[-1])
            Fl = kc * max(0.0, -x - lo[-1])
            au[-1] += Fu
            al[-1] += Fl
            return au / mn, al / mn, (Fl - Fu)

        def energy():
            du = np.diff(np.r_[0.0, up])
            dl = np.diff(np.r_[0.0, lo])
            Es = (0.5 * mn * (vu @ vu + vl @ vl) + 0.5 * kl * (du @ du + dl @ dl))
            Ec = 0.5 * kc * (max(0.0, x - up[-1]) ** 2 + max(0.0, -x - lo[-1]) ** 2)
            return 0.5 * v**2 + Es + Ec

        # let the initial chain settle against the ball statically: linear compress is its equilibrium
        au, al, a = acc(up, lo, x)
        n = int(T / dt)
        E0 = energy()
        zero_cross = []
        share, drift = [], 0.0
        win_vmax, vmax_list = 0.0, []
        window = int(2 * np.pi / dt)
        for i in range(n):
            vu += 0.5 * dt * au; vl += 0.5 * dt * al; v += 0.5 * dt * a
            xold = x
            up += dt * vu; lo += dt * vl; x += dt * v
            au, al, a = acc(up, lo, x)
            vu += 0.5 * dt * au; vl += 0.5 * dt * al; v += 0.5 * dt * a
            if xold > 0 >= x:
                zero_cross.append((i + xold / (xold - x)) * dt)
            win_vmax = max(win_vmax, abs(v))
            if (i + 1) % window == 0:
                vmax_list.append(win_vmax)
                win_vmax = 0.0
                E = energy()
                drift = max(drift, abs(E - E0) / E0)
                share.append(0.5 * vmax_list[-1] ** 2 / E0)
        share = np.array(share)
        # first period: time of the second downward crossing minus the first
        first = zero_cross[1] - zero_cross[0] if len(zero_cross) > 1 else np.nan
        return first, drift, share

    res = {}
    for ms, Nm, cfl in ((0.01, 6, 0.12), (0.1, 12, 0.12), (1.0, 12, 0.06)):
        first, drift, share = chains(ms, N=Nm, cfl=cfl)
        res[ms] = (first, drift, share)
        late = share[len(share) // 2:]
        print(f"      {ms:5.2f}                 {first:7.3f}       {2 * np.pi * np.sqrt(1 + ms / 3):7.3f}                     {drift:.1e}        {share[0]:5.3f}          {late.mean():5.3f}   {late.min():5.3f}")
    check("total energy conserved in the chain model (drift below 0.5%)", all(r[1] < 5e-3 for r in res.values()))
    check("first period follows the one-third-of-spring-mass rule within 5% for ms <= 0.1",
          all(abs(res[ms][0] / (2 * np.pi * np.sqrt(1 + ms / 3)) - 1) < 0.05 for ms in (0.01, 0.1)))
    print("      mesh check, share after the first cycle and late-time mean:")
    f_n = []
    for ms, Nm, cfl in ((0.1, 24, 0.12), (1.0, 24, 0.06)):
        _, dr_n, s_n = chains(ms, N=Nm, cfl=cfl)
        f_n.append((ms, s_n[0], s_n[len(s_n) // 2:].mean()))
        print(f"        ms = {ms}, N = {Nm}: first cycle {s_n[0]:.3f} (N = 12: {res[ms][2][0]:.3f}), late mean {f_n[-1][2]:.3f}, drift {dr_n:.1e}")
    check("loss in the first cycle is mesh-independent to 0.03",
          all(abs(q[1] - res[q[0]][2][0]) < 0.03 for q in f_n))
    check("light spring keeps nearly all the energy in the ball", res[0.01][2][len(res[0.01][2]) // 2:].mean() > 0.9)
    check("ms = 0.1: the ball slowly leaks energy into the spring (late share 0.6-0.95, not converged in detail)",
          0.6 < res[0.1][2][len(res[0.1][2]) // 2:].mean() < 0.95)
    check("ms = 1: most energy has left the ball within a few cycles, with none lost overall",
          res[1.0][2][len(res[1.0][2]) // 2:].mean() < 0.25 and f_n[1][2] < 0.25)

    print("\n3. Three kinds of damping, three decay laws (A = 0.3, touching springs, one-sided)")
    # a. Dashpot inside each spring. The spring end can lag behind the ball when it relaxes.
    def dashpot(zeta, T=150.0, dt=1e-3, A=0.3):
        """Each spring is k in parallel with a dashpot, massless, and can only push.

        In contact its end follows the ball. It lets go when its force reaches zero, then its
        end relaxes on its own, du' = -du/c, burning all the energy still stored in it. It is
        picked up again when the end and the ball meet.
        """
        c = 2 * zeta
        x, v = A, 0.0
        comp = [A, 0.0]             # compression of upper, lower
        touching = [True, False]
        heat, E0 = 0.0, 0.5 * A**2
        peaks, t_peaks, separated, prev_v = [], [], 0, 0.0
        for i in range(int(T / dt)):
            F = 0.0
            for j, sgn in ((0, 1.0), (1, -1.0)):
                reach = sgn * x      # how far the ball presses toward this spring's anchor
                rate = sgn * v
                if touching[j]:
                    if reach > 0 and reach + c * rate > 0:
                        comp[j] = reach
                        F -= sgn * (reach + c * rate)
                        heat += c * rate**2 * dt
                        continue
                    touching[j] = False
                relaxed = comp[j] * np.exp(-dt / c)
                if reach >= relaxed and reach > 0 and reach + c * rate > 0:
                    heat += 0.5 * (comp[j] ** 2 - reach**2) if comp[j] > reach else 0.0
                    comp[j] = reach
                    touching[j] = True
                else:
                    heat += 0.5 * (comp[j] ** 2 - relaxed**2)
                    comp[j] = relaxed
                    if relaxed > 1e-6:
                        separated += 1
            v += dt * F
            x += dt * v
            if prev_v * v < 0:
                peaks.append(abs(x)); t_peaks.append(i * dt)
            prev_v = v
        E_end = 0.5 * v**2 + 0.5 * comp[0] ** 2 + 0.5 * comp[1] ** 2
        return np.array(t_peaks), np.array(peaks), separated * dt / T, (E_end + heat) / E0, heat / E0

    print("   a. dashpot in the spring: amplitude should fall as exp(-zeta t) if contact never broke")
    print("      zeta   fitted decay rate  rate/zeta  fraction of time a compressed spring is out of contact  energy audit")
    ok = True
    for z in (0.005, 0.02, 0.05):
        t, p, sep, audit, _ = dashpot(z)
        good = p > 1e-5
        rate = -np.polyfit(t[good], np.log(p[good]), 1)[0]
        print(f"      {z:5.3f}      {rate:8.5f}       {rate / z:5.3f}          {sep:6.3f}                                          {audit:.4f}")
        ok &= 0.7 < rate / z < 1.5 and abs(audit - 1) < 0.02
    check("exponential decay close to the linear rate, energy audit closes", ok)

    # b. Dry friction
    def dry(mu, T=300.0, dt=1e-3, A=0.3):
        x, v = A, 0.0
        peaks = [A]
        for i in range(int(T / dt)):
            if v == 0.0 and abs(x) <= mu:
                return np.array(peaks), i * dt   # static friction holds it
            a = -x - mu * (np.sign(v) if v != 0.0 else -np.sign(x))
            v_new = v + dt * a
            if v != 0.0 and v * v_new <= 0:      # it stops inside this step
                v_new = 0.0
                peaks.append(abs(x))
            x += dt * v_new
            v = v_new
        return np.array(peaks), None

    mu = 0.01
    p, t_stop = dry(mu)
    drop = -np.diff(p[:10]).mean()
    n_half = int(np.ceil((0.3 - mu) / (2 * mu)))
    print(f"   b. dry friction F = {mu} k: amplitude drop per half cycle {drop:.4f} (derived 2F/k = {2 * mu:.4f});"
          f" stops at t = {t_stop:.1f}, derived about {n_half} half cycles = {n_half * np.pi:.1f}")
    check("dry friction: linear decay and a finite stopping time",
          abs(drop / (2 * mu) - 1) < 0.05 and t_stop is not None and abs(t_stop / (n_half * np.pi) - 1) < 0.1)

    # c. Air drag
    def air(b, T=600.0, dt=2e-3, A=0.3):
        def rhs(t, y):
            return [y[1], -y[0] - b * y[1] * abs(y[1])]

        def turn(t, y):
            return y[1]

        sol = solve_ivp(rhs, (0, T), [A, 0.0], method="DOP853", rtol=1e-9, atol=1e-11, events=turn)
        return sol.t_events[0], np.abs(sol.y_events[0][:, 0])

    b = 0.5
    t, p = air(b)
    p = p[1:]; t = t[1:]
    slope = np.polyfit(t, 1 / p, 1)[0]
    print(f"   c. air drag F = b v|v|, b = {b}: slope of 1/A vs t = {slope:.5f} (derived 4 b / (3 pi) = {4 * b / (3 * np.pi):.5f})")
    check("air drag: amplitude falls as 1/t, never exactly zero", abs(slope / (4 * b / (3 * np.pi)) - 1) < 0.05)

    print("\n4. Sideways stability")
    # a. Static: preloaded springs hinged at the anchors, lateral stiffness k_s each.
    def sim_2d(k_s, L0, Lanchor, A, y0, T, dt=2e-3, zeta_y=0.0):
        z, y = A, y0
        vz, vy = 0.0, 0.0

        def acc(z, y, vy):
            az, ay = 0.0, -2 * zeta_y * vy
            ru = np.hypot(Lanchor - z, y)
            cu = L0 - ru
            if cu > 0:
                az += -cu * (Lanchor - z) / ru
                ay += cu * y / ru - k_s * y
            rl = np.hypot(Lanchor + z, y)
            cl = L0 - rl
            if cl > 0:
                az += cl * (Lanchor + z) / rl
                ay += cl * y / rl - k_s * y
            return az, ay

        n = int(T / dt)
        rec_every = max(1, int(0.1 / dt))
        ts, ys = [], []
        az, ay = acc(z, y, vy)
        for i in range(n):
            vz += 0.5 * dt * az; vy += 0.5 * dt * ay
            z += dt * vz; y += dt * vy
            az, ay = acc(z, y, vy)
            vz += 0.5 * dt * az; vy += 0.5 * dt * ay
            if i % rec_every == 0:
                ts.append(i * dt); ys.append(abs(y))
            if abs(y) > 0.2:
                break
        return np.array(ts), np.array(ys)

    L0, La = 1.0, 0.9
    P = (L0 - La)
    lam_pred = np.sqrt(2 * P / La)
    t, y = sim_2d(0.0, L0, La, 0.0, 1e-6, 40.0, dt=4e-3)
    sel = (y > 1e-5) & (y < 1e-3)
    lam = np.polyfit(t[sel], np.log(y[sel]), 1)[0]
    print(f"   a. preload strain {P / La:.3f}, no lateral stiffness: growth rate {lam:.4f}, derived sqrt(2 P / (L m)) = {lam_pred:.4f}")
    check("ball exactly on the line is unstable under preload if the springs have no lateral stiffness",
          abs(lam / lam_pred - 1) < 0.02)
    ks_crit = P / La
    t1, y1 = sim_2d(0.9 * ks_crit, L0, La, 0.0, 1e-6, 150.0, dt=4e-3)
    t2, y2 = sim_2d(1.1 * ks_crit, L0, La, 0.0, 1e-6, 150.0, dt=4e-3)
    print(f"      threshold k_s/k = P/(k L) = {ks_crit:.4f}: at 0.9 of it y reaches {y1.max():.1e}, at 1.1 of it y stays at {y2.max():.1e}")
    check("buckling threshold k_s = P/L", y1.max() > 1e-3 and y2.max() < 2e-6)

    print("   b. What is k_s for a real coil? Finite-element helix (music wire, d 2 mm, D 20 mm, pitch 6 mm), ends kept parallel")
    d, D, G, pitch = 2e-3, 20e-3, 72e9, 6e-3
    print("      coils  length/D   k_axial (N/m)  k_lateral (N/m)  lateral/axial")
    ratios = []
    for n in (4, 7, 10, 15, 25):
        k_a, _ = helix_fe(d, D, pitch, n, G, per_turn=24)
        k_l = helix_fe(d, D, pitch, n, G, per_turn=24, mode="lateral")
        ratios.append(k_l / k_a)
        print(f"      {n:4d}   {n * pitch / D:6.2f}    {k_a:10.0f}      {k_l:10.0f}       {k_l / k_a:6.3f}")
    check("slender springs have much less relative lateral stiffness (they buckle earlier)", ratios[-1] < 0.3 * ratios[0])

    # c. Floquet, touching springs with no preload: lateral stiffness oscillates at twice the bounce rate.
    print("   c. Touching, no preload, A/L0 = 0.3: lateral motion y'' + [k_s - |z|/(L0-|z|)] y = 0, z = A cos t")
    A = 0.3

    def floquet(k_s, zeta_y=0.0):
        def rhs(t, s):
            z = abs(A * np.cos(t))
            q = k_s - z / (1 - z)
            return [s[1], -q * s[0] - 2 * zeta_y * s[1], s[3], -q * s[2] - 2 * zeta_y * s[3]]

        sol = solve_ivp(rhs, (0, np.pi), [1, 0, 0, 1], method="DOP853", rtol=1e-9, atol=1e-11)
        M = sol.y[:, -1].reshape(2, 2, order="F")
        rho = np.abs(np.linalg.eigvals(M)).max()
        return np.log(rho) / np.pi

    grid = np.array([0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0, 3.5, 3.9, 4.1, 4.5, 6.0])
    rates = np.array([floquet(k) for k in grid])
    print("      k_s/k:     " + " ".join(f"{k:6.2f}" for k in grid[:12]))
    print("      growth/s:  " + " ".join(f"{r:6.3f}" for r in rates[:12]))
    print("      k_s/k:     " + " ".join(f"{k:6.2f}" for k in grid[12:]))
    print("      growth/s:  " + " ".join(f"{r:6.3f}" for r in rates[12:]))
    fine = np.linspace(0.02, 6.0, 120)
    fr = np.array([floquet(k) for k in fine])
    unstable = fr > 1e-3
    bands, start = [], None
    for k_, u in zip(fine, unstable):
        if u and start is None:
            start = k_
        if not u and start is not None:
            bands.append((start, k_)); start = None
    if start is not None:
        bands.append((start, fine[-1]))
    print("      unstable bands of k_s/k, undamped: " + ", ".join(f"{a:.2f}-{b:.2f}" for a, b in bands))
    fr_d = np.array([floquet(k, 0.02) for k in fine])
    print(f"      with lateral damping ratio 0.02 (relative to bounce frequency): {100 * np.mean(fr_d > 1e-3):.0f} % of the scanned range stays unstable"
          f" (undamped {100 * np.mean(unstable):.0f} %)")
    i_peak = int(np.argmax(fr * (fine > 0.3)))
    k_tongue = fine[i_peak]
    print(f"      strongest band above k_s/k = 0.3: k_s/k = {k_tongue:.2f}, growth {fr[i_peak]:.3f}")
    check("resonant band near k_s = k (lateral frequency = bounce frequency)", 0.5 < k_tongue < 1.3)
    k_stable = next(k for k, r in zip(fine, fr) if k > 1.5 and r < 1e-4)
    t_a, y_a = sim_2d(k_tongue, 1.0, 1.0, A, 1e-7, 100.0, dt=4e-3)
    t_b, y_b = sim_2d(k_stable, 1.0, 1.0, A, 1e-7, 100.0, dt=4e-3)
    sel = (y_a > 3e-7) & (y_a < 1e-3)
    if sel.sum() > 10:
        env_t, env_y = t_a[sel], y_a[sel]
        # fit the upper envelope with a window maximum
        tb = np.arange(env_t[0], env_t[-1], 2 * np.pi)
        m = np.array([env_y[(env_t >= a) & (env_t < a + 2 * np.pi)].max() for a in tb if ((env_t >= a) & (env_t < a + 2 * np.pi)).any()])
        tbm = np.array([a for a in tb if ((env_t >= a) & (env_t < a + 2 * np.pi)).any()])
        rate_sim = np.polyfit(tbm, np.log(m), 1)[0]
    else:
        rate_sim = float("nan")
    print(f"      full 2D simulation, 1e-7 sideways offset: at k_s/k = {k_tongue:.2f} growth {rate_sim:.3f} (Floquet {fr[i_peak]:.3f});"
          f" at k_s/k = {k_stable:.2f} y stays below {y_b.max():.1e}")
    check("full simulation confirms Floquet growth rate within 25%", abs(rate_sim / fr[i_peak] - 1) < 0.25)
    check("away from the bands the on-line motion stays stable", y_b.max() < 3e-6)

    print()
    if FAILS:
        print("FAILED:", FAILS)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
