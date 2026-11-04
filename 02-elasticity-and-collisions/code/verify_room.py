"""E09: bounces in a room. Geometric series, the Zeno point, and what really happens there.

  1. Event-driven ball under gravity against the series for bounce times and heights.
  2. The same ball with a finite-stiffness contact: a finite number of bounces, no Zeno.
  3. Where the ideal model breaks: flight time vs contact time, height vs atomic size.
  4. Gravity-free box, e = 1 and e < 1, with and without air drag.
  5. The 3D room with gravity, by axes.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import g_0 as g, k_B  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def event_driven(h0, e, n_max=400, v_stop=1e-12):
    """Exact flight parabolas, instantaneous bounces. No series formula is used."""
    t, times, heights = 0.0, [], []
    # first fall
    v = np.sqrt(2 * g * h0)
    t = v / g
    times.append(t)
    for _ in range(n_max):
        v *= e
        if v < v_stop:
            break
        heights.append(v**2 / (2 * g))
        t += 2 * v / g  # time up and back
        times.append(t)
    return np.array(times), np.array(heights)


def soft_contact(h0, e_target, m, k, c_zeta, T_cap=30.0):
    """Finite-stiffness contact, spring-damper that cannot pull, integrated phase by phase."""
    w = np.sqrt(k / m)
    c = 2 * c_zeta * m * w

    def flight(t, y):
        return [y[1], -g]

    def touch(t, y):
        return y[0]

    touch.terminal, touch.direction = True, -1

    def contact(t, y):
        z, v = y
        return [v, -g + max(0.0, -k * z - c * v) / m]

    def lift(t, y):
        return -k * y[0] - c * y[1]  # force falling to zero

    lift.terminal, lift.direction = True, -1
    t, z, v = 0.0, h0, 0.0
    bounce_times, v_in, v_out, contact_durations, peak = [], [], [], [], []
    tau = np.pi / w
    while t < T_cap:
        tf = max(2 * abs(v) / g, 2 * np.sqrt(2 * max(z, 0) / g)) + 1e-3
        sol = solve_ivp(flight, (t, t + 4 * tf + 1), [z, v], events=touch, method="DOP853",
                        rtol=1e-12, atol=1e-15)
        t = sol.t[-1]
        v_i = sol.y[1, -1]
        bounce_times.append(t)
        v_in.append(-v_i)
        sol = solve_ivp(contact, (t, t + 200 * tau), [0.0, v_i], events=lift, method="DOP853",
                        rtol=1e-11, atol=1e-14, first_step=tau * 1e-4)
        if sol.t_events[0].size == 0:
            # Force never reached zero in 200 contact times: it is resting on the floor.
            return dict(n=len(bounce_times) - 1, t_rest=t, times=np.array(bounce_times), v_in=np.array(v_in),
                        v_out=np.array(v_out), z_final=sol.y[0, -1], tau=tau, dur=np.array(contact_durations))
        contact_durations.append(sol.t[-1] - t)
        t, z, v = sol.t[-1], sol.y[0, -1], sol.y[1, -1]
        v_out.append(v)
        # after lift-off with z < 0 it is still in the spring's zone, but the force is zero,
        # and the flight phase carries it to z = 0 only if v is big enough
        if z < 0:
            # height it can reach above z
            if z + v**2 / (2 * g) <= 0:
                # re-engages within the zone; continue as contact from here
                sol = solve_ivp(contact, (t, t + 200 * tau), [z, v], method="DOP853",
                                rtol=1e-11, atol=1e-14, first_step=tau * 1e-4)
                return dict(n=len(bounce_times), t_rest=t, times=np.array(bounce_times), v_in=np.array(v_in),
                            v_out=np.array(v_out), z_final=sol.y[0, -1], tau=tau, dur=np.array(contact_durations))
            dt_up = (v - np.sqrt(v**2 + 2 * g * z)) / g
            t += dt_up
            v = v - g * dt_up
            z = 0.0
    return dict(n=len(bounce_times), t_rest=t, times=np.array(bounce_times), v_in=np.array(v_in),
                v_out=np.array(v_out), z_final=z, tau=tau, dur=np.array(contact_durations))


def e_of_zeta(z):
    s = np.sqrt(1 - z**2)
    phi = np.arctan2(2 * z * s, 2 * z**2 - 1)
    t = phi / s
    return -np.exp(-z * t) * (np.cos(phi) - z / s * np.sin(phi))


def zeta_of_e(e):
    from scipy.optimize import brentq
    return brentq(lambda z: e_of_zeta(z) - e, 1e-6, 0.99)


def axis_events(x0, v0, L, e, T, a=0.0, n_max=200000, v_stop=1e-12):
    """One coordinate in [0, L], acceleration a (<= 0), walls with restitution e.
    Returns the list of (time, wall) up to time T. Exact flights."""
    t, x, v = 0.0, x0, v0
    out = []
    while t < T and len(out) < n_max:
        times = []
        if a == 0.0:
            if v > 0:
                times.append(((L - x) / v, "high"))
            elif v < 0:
                times.append((x / -v, "low"))
        else:
            # low wall
            disc = v**2 - 2 * a * x
            root = (-v - np.sqrt(disc)) / a
            if root > 1e-15:
                times.append((root, "low"))
            disc_h = v**2 + 2 * a * (L - x)
            if disc_h >= 0:
                r = (-v + np.sqrt(disc_h)) / a
                if r > 1e-15:
                    times.append((r, "high"))
        if not times:
            break
        dt, wall = min(times, key=lambda q: q[0])
        if t + dt > T:
            break
        x = (L if wall == "high" else 0.0)
        v = -e * (v + a * dt)
        t += dt
        out.append((t, wall))
        if abs(v) < v_stop:
            out.append((t, "stopped"))
            break
    return out


def main():
    h0, e = 1.0, 0.8
    t0 = np.sqrt(2 * h0 / g)
    print(f"1. Point ball under gravity, h0 = {h0} m, e = {e}")
    times, heights = event_driven(h0, e)
    T_series = t0 * (1 + 2 * e / (1 - e))
    n = np.arange(1, 11)
    t_series_n = t0 * (1 + 2 * e * (1 - e ** (n - 1)) / (1 - e))
    print(f"   t0 = {t0:.4f} s; series limit T = t0 (1+e)/(1-e) = {T_series:.6f} s;"
          f" event-driven time at {len(times)} bounces {times[-1]:.6f} s")
    print("   bounce n    time (event-driven)   time (series)    height after bounce (m)   e^(2n) h0")
    for k in (1, 2, 3, 5, 10):
        print(f"     {k:3d}       {times[k - 1]:12.6f}      {t_series_n[k - 1]:12.6f}      {heights[k - 1]:.6e}          {e ** (2 * k) * h0:.6e}")
    check("event-driven bounce times equal the series", np.allclose(times[:10], t_series_n, rtol=1e-12))
    check("heights equal e^(2n) h0", np.allclose(heights[:50], e ** (2 * np.arange(1, 51)) * h0, rtol=1e-10))
    check("simulation accumulates at the series limit (infinitely many bounces, finite time)",
          abs(times[-1] / T_series - 1) < 1e-9 and len(times) > 100, f"{len(times)} bounces before v < 1e-12 m/s")
    d_total = h0 * (1 + e**2) / (1 - e**2)
    print(f"   total path length also finite: h0 (1+e^2)/(1-e^2) = {d_total:.4f} m")
    times_e1, _ = event_driven(h0, 1.0, n_max=1000)
    gaps = np.diff(times_e1)
    print(f"   e = 1: gap between bounces constant at {gaps.mean():.6f} s (2 t0 = {2 * t0:.6f}), {len(times_e1)} bounces take {times_e1[-1]:.2f} s")
    check("e = 1: bounce time grows without limit", np.ptp(gaps) < 1e-9 and times_e1[-1] > 1000 * t0)

    print("\n2. Same ball with a real contact: spring-damper, tension-free, gravity on during contact")
    print("     m (kg)  k (N/m)    contact time (ms)   e_target   bounces   rest at t (s)   Zeno T (s)   sag mg/k (mm), final depth (mm)")
    runs = {}
    for m, k_c, e_t in ((0.05, 2e5, 0.8), (0.05, 2e4, 0.8), (0.05, 2e5, 0.9), (0.05, 2e5, 0.6), (0.005, 2e6, 0.95)):
        z = zeta_of_e(e_t)
        r = soft_contact(h0, e_t, m, k_c, z)
        r["e"] = e_t
        runs[(m, k_c, e_t)] = r
        T_z = t0 * (1 + e_t) / (1 - e_t)
        print(f"     {m:5.3f}  {k_c:9.0f}   {1e3 * r['tau']:8.3f}          {e_t:4.2f}     {r['n']:5d}      {r['t_rest']:8.4f}     {T_z:8.4f}     {1e3 * m * g / k_c:.3f}, {-1e3 * r['z_final']:.3f}")
    r = runs[(0.05, 2e5, 0.8)]
    pred_n = np.log(2 * np.sqrt(2 * g * h0) / (g * r["tau"])) / np.log(1 / 0.8)
    print(f"   crossover when flight time = contact time: n* = ln(2 v0 / (g tau)) / ln(1/e) = {pred_n:.1f}, simulated bounces {r['n']}")
    check("finite number of bounces with a real contact", all(q["n"] < 500 for q in runs.values()))
    check("bounce count is close to the flight-time = contact-time estimate", abs(r["n"] - pred_n) < 0.3 * pred_n + 3)
    check("it comes to rest near the Zeno time", all(abs(q["t_rest"] / (t0 * (1 + q["e"]) / (1 - q["e"])) - 1) < 0.05 for q in runs.values()))
    effective = r["v_out"] / r["v_in"][: len(r["v_out"])]
    print(f"   per-bounce restitution: first {effective[0]:.4f}, middle {effective[len(effective) // 2]:.4f}, last three {np.round(effective[-3:], 3)}")
    check("effective e constant for most bounces, falls in the last ones (gravity during contact)",
          abs(effective[len(effective) // 2] - 0.8) < 0.01 and effective[-1] < 0.75)

    print("\n3. Other limits on the series (e = 0.8, h0 = 1 m, m = 50 g)")
    n_atom = np.log(h0 / 1e-10) / (2 * np.log(1 / e))
    m = 0.05
    n_thermal = np.log(m * g * h0 / (k_B * 300)) / (2 * np.log(1 / e))
    print(f"   height below one atomic size (1e-10 m) after n = {n_atom:.0f} bounces")
    print(f"   bounce energy below k_B T at 300 K after n = {n_thermal:.0f} bounces; that hop height would be {k_B * 300 / (m * g):.1e} m, far below an atom, so the point is never reached")
    D = k_B * 300 / (6 * np.pi * 1.8e-5 * 0.025)
    print(f"   50 g, 2.5 cm ball in still air, no gravity: Stokes diffusion D = {D:.1e} m^2/s, time to wander 1 m = {1 / (2 * D) / 3.156e7:.1e} years")

    print("\n4. No gravity: a ball between two walls L = 4 m apart, v0 = 5 m/s")
    L, v0 = 4.0, 5.0
    ev1 = axis_events(0.0, v0, L, 1.0, 1000.0)
    print(f"   e = 1: {len(ev1)} collisions in 1000 s; rate {len(ev1) / 1000:.3f}/s, v0/L = {v0 / L:.3f}/s")
    check("e = 1: collision count grows linearly forever", abs(len(ev1) / 1000 / (v0 / L) - 1) < 0.01)
    e2 = 0.8
    ev2 = axis_events(0.0, v0, L, e2, 1e15, v_stop=1e-300, n_max=120)
    t2 = np.array([q[0] for q in ev2])
    n2 = np.arange(1, len(t2) + 1)
    pred = (L / v0) * (e2 ** (-n2) - 1) / (1 / e2 - 1) * 1.0
    print(f"   e = {e2}: n-th collision time (L/v0) (e^-n - 1)/(1/e - 1). Simulated vs series at n = 10, 50, 100:")
    for q in (10, 50, 100):
        print(f"      n = {q:3d}: {t2[q - 1]:.4e} s, {pred[q - 1]:.4e} s")
    check("e < 1 without gravity: gaps grow geometrically, collisions go on with no finite-time limit",
          np.allclose(t2[:100], pred[:100], rtol=1e-9) and t2[100] / t2[99] > 1 / e2 * 0.99)
    v_at = v0 * e2 ** n2
    print(f"   speed against time: v t (1-e)/(e L) at n = 100 is {v_at[99] * t2[99] * (1 - e2) / (e2 * L):.4f} (derived value 1 as n grows)")
    check("speed falls as 1/t", abs(v_at[99] * t2[99] * (1 - e2) / (e2 * L) - 1) < 0.01)

    # Air drag: separate results for linear and quadratic
    print("   adding air: linear drag b, and the distance it allows")
    b = 1e-3  # kg/s, mass 0.05 kg gives m/b = 50 s
    m_ = 0.05
    stop_dist = m_ * v0 / b

    def box_linear(v0, b_over_m, L, e, T_end=5000.0):
        # exact: between walls x(t) = (v/(b/m)) (1 - exp(-(b/m) t))
        x, v, t, count = 0.0, v0, 0.0, 0
        while t < T_end:
            room = L  # one crossing each time after a hit
            reach = v / b_over_m
            if reach <= room:
                return count, t, False
            dt = -np.log(1 - room * b_over_m / v) / b_over_m
            v = (v - room * b_over_m) * e
            t += dt
            count += 1
        return count, t, True

    count_l, t_l, still = box_linear(v0, b / m_, L, 1.0, 1e5)
    print(f"      linear drag, m/b = {m_ / b:.0f} s, e = 1: total reach v0 m/b = {stop_dist:.0f} m, so {count_l} collisions, derived floor({stop_dist:.0f}/{L:.0f}) = {int(stop_dist // L)}")
    check("linear drag: finite collisions even with e = 1", count_l == int(stop_dist // L) and not still)

    print("\n5. Room with gravity, 4 x 3 x 2.5 m, start (1, 1, 1) m with v = (3, 2, 6) m/s, frictionless, axes uncoupled")
    Lx, Ly, Lz = 4.0, 3.0, 2.5
    x0, y0, z0 = 1.0, 1.0, 1.0
    vx, vy, vz = 3.0, 2.0, 6.0
    for e_r, T in ((1.0, 100.0), (0.8, 100.0)):
        ex = axis_events(x0, vx, Lx, e_r, T)
        ey = axis_events(y0, vy, Ly, e_r, T)
        ez = axis_events(z0, vz, Lz, e_r, T, a=-g)
        floor = [q for q in ez if q[1] == "low"]
        ceil = [q for q in ez if q[1] == "high"]
        print(f"   e = {e_r}: by t = {T:.0f} s, side walls {len(ex) + len(ey)}, ceiling {len(ceil)}, floor {len(floor)}")
        if e_r == 1.0:
            ez2 = axis_events(z0, vz, Lz, e_r, 2 * T, a=-g)
            ratio = len(ez2) / max(1, len(ez))
            check("e = 1: collisions keep accumulating at a steady rate", abs(ratio - 2) < 0.1, f"count ratio over double time {ratio:.2f}")
            E0 = 0.5 * vz**2 + g * z0
            print(f"      vertical energy per kg conserved by construction {E0:.2f} J/kg; ceiling reached since {z0 + vz**2 / (2 * g):.2f} m > {Lz} m")
        else:
            T_collapse = floor[-1][0]
            print(f"      vertical motion stops at T = {T_collapse:.4f} s after {len(floor)} floor hits (speed cut-off 1e-12 m/s)")
            wall_before = sum(1 for q in ex + ey if q[0] <= T_collapse)
            print(f"      side-wall hits before that time {wall_before}; floor hits are infinitely many in the ideal limit, side-wall hits are not")
            late_x = axis_events(x0, vx, Lx, e_r, 1e8, v_stop=1e-300, n_max=60)
            print(f"      with a frictionless floor the ball then slides: side-wall hit no. 10, 30, 50 at t = {late_x[9][0]:.3g}, {late_x[29][0]:.3g}, {late_x[49][0]:.3g} s")
            check("vertical collapse happens in finite time", T_collapse < 30 and len(floor) > 50)
            check("side-wall hits stay finite before the vertical collapse", wall_before < 40)

    print()
    if FAILS:
        print("FAILED:", FAILS)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
