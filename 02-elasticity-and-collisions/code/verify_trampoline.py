"""E05, E06: person and ball on a trampoline.

Part A (E05): person + bed as coupled masses with a one-sided foot contact and a leg whose
length the person controls. Energy audit: mechanical energy change = leg work - heat.
Part B (E06): ball dropped on an axisymmetric membrane solved as a wave equation, so the
"effective bed mass" comes out of the simulation rather than being chosen.

All bed, ball and body numbers are labelled assumptions except where the README cites them.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import g_0 as g  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# ================================================================ Part A
M_P, M_B, K_B = 60.0, 5.0, 5000.0   # kg, kg, N/m (K_B as in Srinivasan et al. 2013)
L_MAX, L_MIN = 1.0, 0.7             # leg length range, m (COM above foot)
K_C = 2e4                           # compliance of a braced leg and foot, N/m
C_C = 2 * np.sqrt(K_C * M_B)        # critically damped against the bed mass, so no chatter


def simulate_person(T, h0, zeta_bed=0.0, mode="rigid", delay=1.0, v_leg=1.5, F_max=np.inf,
                    dt=5e-4, c_contact=C_C, M_b=M_B, stroke=None):
    """Velocity-Verlet with semi-implicit damping.

    mode "rigid": constant leg length.
    mode "pump":  leg retracts in flight (no force, so no work) and extends in contact, starting
                  delay * (quarter bed period) after touchdown. It cannot extend while the
                  force exceeds F_max.
    mode "kill":  leg starts long at touchdown and is shortened while the bed pushes up.
    """
    l_min = L_MIN if stroke is None else L_MAX - stroke
    w_q = 0.5 * np.pi / np.sqrt(K_B / M_P)  # quarter period
    z_b0 = -(M_b + 0) * g / K_B
    c_b = 2 * zeta_bed * np.sqrt(K_B * M_P)
    l = L_MAX if mode != "pump" else l_min
    z_p, v_p = h0 + L_MAX, 0.0
    z_b, v_b = z_b0, 0.0
    t, W, D = 0.0, 0.0, 0.0
    apexes, contact_times, bottoms = [], [], []
    in_contact, t_land = False, None
    prev_v, prev_vb = v_p, v_b
    F = 0.0
    ldot = 0.0
    a_p, a_b = -g, -g - K_B * z_b / M_b

    def energy():
        gap = (z_p - l) - z_b
        return (0.5 * M_P * v_p**2 + 0.5 * M_b * v_b**2 + M_P * g * z_p + M_b * g * z_b
                + 0.5 * K_B * z_b**2 + (0.5 * K_C * gap**2 if gap < 0 else 0.0))

    E0 = energy()
    last_touch = -1.0
    F_prev = 0.0
    for _ in range(int(T / dt)):
        v_p_half = v_p + 0.5 * dt * a_p
        v_b_half = v_b + 0.5 * dt * a_b
        z_p += dt * v_p_half
        z_b += dt * v_b_half
        ldot = 0.0
        if mode == "pump":
            if not in_contact:
                if l > l_min:
                    ldot = -3.0
            elif t - t_land >= delay * w_q and l < L_MAX and F_prev < F_max:
                ldot = v_leg
        elif mode == "kill":
            if not in_contact:
                if l < L_MAX:
                    ldot = 3.0
            elif t - t_land >= w_q and l > L_MIN:
                ldot = -v_leg
        l = min(L_MAX, max(l_min, l + ldot * dt))
        gap = (z_p - l) - z_b
        rate = v_p_half - ldot - v_b_half
        F = max(0.0, -K_C * gap - c_contact * rate) if gap < 0 else 0.0
        a_p = -g + F / M_P
        v_p = v_p_half + 0.5 * dt * a_p
        a_b = -g - F / M_b - K_B * z_b / M_b
        v_b = (v_b_half + 0.5 * dt * a_b) / (1 + 0.5 * dt * c_b / M_b)
        a_b -= c_b * v_b / M_b
        W += F * ldot * dt
        D += c_b * v_b**2 * dt
        if F > 0:
            D += F * 0.0
            D += c_contact * rate**2 * dt if -K_C * gap - c_contact * rate > 0 else 0.0
            last_touch = t
            if not in_contact:
                in_contact, t_land = True, t
        elif in_contact and t - last_touch > 0.03:  # 30 ms gap: really airborne
            in_contact = False
            contact_times.append(last_touch - t_land)
        if in_contact and prev_v < 0 <= v_p:
            bottoms.append(t - t_land)
        if prev_v > 0 >= v_p and not in_contact:
            apexes.append(z_p - L_MAX)
        prev_v = v_p
        F_prev = F
        t += dt
    return dict(apexes=np.array(apexes), contacts=np.array(contact_times), W=W, D=D, E0=E0, E1=energy(),
                bottoms=np.array(bottoms), w_q=w_q)


# ================================================================ Part B
def build_membrane(a, dr, T_m, sigma, r_patch):
    N = int(round(a / dr))
    r = (np.arange(N) + 0.5) * dr
    area = 2 * np.pi * r * dr
    rf = np.arange(1, N + 1) * dr  # outer face of each annulus; the last one is the rim
    g_face = T_m * 2 * np.pi * rf / dr
    in_patch = r <= r_patch
    weights = np.where(in_patch, area, 0.0)
    weights /= weights.sum()
    return dict(N=N, r=r, dr=dr, area=area, g_face=g_face, weights=weights, sigma=sigma,
                in_patch=in_patch, c=np.sqrt(T_m / sigma))


def static_stiffness(mem):
    gf = mem["g_face"]
    diag = gf.copy()
    diag[1:] += gf[:-1]
    A = sp.diags([diag, -gf[:-1], -gf[:-1]], [0, 1, -1], format="csc")
    u = spla.spsolve(A, mem["weights"])  # unit force spread uniformly over the patch
    return 1.0 / (mem["weights"] @ u)


def fundamental_freq(mem):
    gf, area, s = mem["g_face"], mem["area"], mem["sigma"]
    diag = gf.copy()
    diag[1:] += gf[:-1]
    d = diag / (s * area)
    off = -gf[:-1] / (s * np.sqrt(area[:-1] * area[1:]))
    A = np.diag(d) + np.diag(off, 1) + np.diag(off, -1)
    return np.sqrt(np.linalg.eigvalsh(A)[0])


def e_to_zeta(e0):
    def f(z):
        s = np.sqrt(1 - z**2)
        phi = np.arctan2(2 * z * s, 2 * z**2 - 1)
        t = phi / s
        return -np.exp(-z * t) * (np.cos(phi) - z / s * np.sin(phi)) - e0
    return brentq(f, 1e-6, 0.99)


def drop_ball(mem, m, k, zeta, h0, rigid=False, gamma=2.0, T_max=2.0, gravity=True):
    """Ball (mass m, own contact stiffness k, own damping) dropped on the membrane or a floor.

    Leapfrog in time. The ball is free to leave and re-strike the patch.
    e is from the peak height above the undisturbed bed; e_v is lift-off speed / impact speed.
    """
    N, area, gf, wts, s = mem["N"], mem["area"], mem["g_face"], mem["weights"], mem["sigma"]
    c_ball = 2 * zeta * np.sqrt(k * m)
    mi = s * area
    patch_mass = s * area[mem["in_patch"]].sum()
    dt = min(0.4 * mem["dr"] / mem["c"], 0.25 / np.sqrt(k / min(m, patch_mass)))
    g_eff = g if gravity else 0.0
    u, v = np.zeros(N), np.zeros(N)
    v_in = -np.sqrt(2 * g * h0)
    z, vz = 0.0, v_in      # ball's lowest point, starting at the bed surface
    t, heat_ball = 0.0, 0.0
    peak_z, v_lift, contact_total, bed_share = -np.inf, 0.0, 0.0, 0.0
    in_contact, rose = False, False
    E0 = 0.5 * m * v_in**2
    flux = np.zeros(N)
    a_b = np.zeros(N)
    F, az = 0.0, -g_eff
    for _ in range(int(T_max / dt)):
        if not rigid:
            v += 0.5 * dt * a_b
            u += dt * v
        vz += 0.5 * dt * az
        z += dt * vz
        if rigid:
            up, vp = 0.0, 0.0
        else:
            up, vp = wts @ u, wts @ v
        pen = up - z
        rate = vp - vz
        F = max(0.0, k * pen + c_ball * rate) if pen > 0 else 0.0
        if not rigid:
            flux[:-1] = gf[:-1] * (u[1:] - u[:-1])
            flux[-1] = -gf[-1] * u[-1]
            a_b[0] = flux[0]
            a_b[1:] = flux[1:] - flux[:-1]
            a_b /= mi
            a_b -= (F / mi) * wts
            a_b -= gamma * v
            v += 0.5 * dt * a_b
        az = -g_eff + F / m
        vz += 0.5 * dt * az
        t += dt
        if F > 0:
            contact_total += dt
            heat_ball += c_ball * rate**2 * dt
            in_contact = True
        else:
            if in_contact and vz > 0:
                v_lift = vz
                rose = True
                if not rigid:
                    bed_share = (0.5 * np.sum(mi * v**2) + 0.5 * np.sum(gf[:-1] * (u[1:] - u[:-1]) ** 2)
                                 + 0.5 * gf[-1] * u[-1] ** 2) / E0
            in_contact = False
        if rose and F == 0:
            if not gravity and z > 0.05:
                break
            if gravity and vz > 0 and z - up > 0.05 and z > 0:  # clear of the bed and rising: ballistic from here
                peak_z = z + vz**2 / (2 * g_eff)
                break
            peak_z = max(peak_z, z)
            if vz < 0 and z > 0:
                break
    E_end = 0.5 * m * vz**2 + 0.5 * k * max(pen, 0.0) ** 2
    if not rigid:
        E_end += (0.5 * np.sum(mi * v**2) + 0.5 * np.sum(gf[:-1] * (u[1:] - u[:-1]) ** 2)
                  + 0.5 * gf[-1] * u[-1] ** 2)
    e_eff = np.sqrt(max(peak_z, 0.0) / h0) if rose else 0.0
    return dict(e=e_eff, e_v=v_lift / -v_in, contact=contact_total, bed_share=bed_share,
                heat_share=heat_ball / E0, E_ratio=E_end / E0, rose=rose)


def main():
    print(f"Part A: person on a bed (M_p = {M_P:.0f} kg, M_b = {M_B:.0f} kg, K = {K_B:.0f} N/m)")
    print("A1. Passive person, stiff legs, no bed damping: energy and contact time")
    r = simulate_person(20, 3.0, mode="rigid")
    v_land = np.sqrt(2 * g * 3.0)

    def t_contact(K):  # massless bed (derived in README): t_c = (pi + 2 atan(g/(w v))) / w
        w = np.sqrt(K / M_P)
        return (np.pi + 2 * np.arctan(g / (w * v_land))) / w

    K_series = K_B * K_C / (K_B + K_C)
    tc_bed, tc_series = t_contact(K_B), t_contact(K_series)
    print(f"     apexes (m): {np.round(r['apexes'][:6], 3)}")
    print(f"     first contact time {r['contacts'][0]:.3f} s; theory {tc_bed:.3f} s for the bed spring alone,"
          f" {tc_series:.3f} s for bed and leg springs in series")
    check("contact time matches bed and leg springs in series, to 5%",
          abs(r["contacts"][0] / tc_series - 1) < 0.05, f"ratio {r['contacts'][0] / tc_series:.3f}")
    audit0 = (r["E1"] - r["E0"] + r["D"]) / r["D"]
    print(f"     audit: (Delta E + heat)/heat = {audit0:+.1e}")
    check("energy audit closes for the passive run", abs(audit0) < 0.02)
    cons = simulate_person(20, 3.0, mode="rigid", c_contact=0.0, dt=1e-4)
    drift = abs(cons["E1"] - cons["E0"]) / (M_P * g * 3.0)
    print(f"     no damping of any kind, dt = 0.1 ms: total energy drift {drift:.1e} of mgh")
    check("total energy conserved by the integrator with no damping", drift < 2e-3)
    r1 = simulate_person(20, 3.0, mode="rigid", dt=2.5e-4)
    print(f"     step-size check: first apex {r1['apexes'][0]:.4f} m (dt/2) vs {r['apexes'][0]:.4f} m")
    check("first apex converged with step size", abs(r1["apexes"][0] - r["apexes"][0]) < 0.01)

    print("A2. What a passive person loses per bounce")
    print("     M_b (kg)  bed zeta   apexes (m) from 3 m")
    for Mb, z in ((0.5, 0.0), (5.0, 0.0), (5.0, 0.01), (5.0, 0.03)):
        rr = simulate_person(15, 3.0, zeta_bed=z, mode="rigid", M_b=Mb, c_contact=2 * np.sqrt(K_C * Mb))
        print(f"     {Mb:6.1f}    {z:5.2f}     {np.round(rr['apexes'][:6], 2)}")
    r_d = simulate_person(25, 3.0, zeta_bed=0.03, mode="rigid")
    check("passive bounces die away with damping", r_d["apexes"][-1] < 0.5 * r_d["apexes"][0],
          f"{r_d['apexes'][0]:.2f} -> {r_d['apexes'][-1]:.2f} m")

    print("A3. Pumping with the legs, bed zeta = 0.03, 0.3 m stroke, no force limit")
    r_p = simulate_person(40, 0.5, zeta_bed=0.03, mode="pump", delay=0.5)
    ap = r_p["apexes"]
    print(f"     apexes from a 0.5 m start: {np.round(ap[:10], 2)} ... {np.round(ap[-3:], 2)}")
    print(f"     leg work {r_p['W'] / 1e3:.1f} kJ in {len(ap)} bounces, mean {r_p['W'] / len(ap):.0f} J per bounce")
    check("pumping gains height", ap[5] > 1.5 * ap[0], f"{ap[0]:.2f} -> {ap[5]:.2f} m")
    r_c = simulate_person(20, 1.0, zeta_bed=0.03, mode="pump", delay=0.5)
    audit = (r_c["E1"] - r_c["E0"] - r_c["W"] + r_c["D"]) / r_c["W"]
    print(f"     audit: (Delta E - leg work + heat)/leg work = {audit:+.1e}")
    check("energy gained equals leg work minus heat", abs(audit) < 0.02, f"{audit:+.1e}")
    print("     Height after 60 s, where leg work per bounce has come near the loss per bounce:")
    print("     bed zeta  stroke (m)  force limit (body weights)  apex (m)")
    steady = {}
    for z, st, cap in ((0.03, 0.3, np.inf), (0.03, 0.3, 4.5), (0.01, 0.3, np.inf), (0.01, 0.3, 4.5),
                       (0.003, 0.3, 4.5), (0.01, 0.15, np.inf)):
        rr = simulate_person(60, 0.5, zeta_bed=z, mode="pump", delay=0.5, F_max=cap * M_P * g, stroke=st)
        steady[(z, st, cap)] = rr["apexes"][-1]
        print(f"     {z:6.3f}     {st:4.2f}        {cap:6.1f}                  {rr['apexes'][-1]:7.2f}")
    F_peak = lambda h: M_P * g + np.sqrt((M_P * g) ** 2 + 2 * K_B * M_P * g * h)
    print(f"     peak bed force on a person, massless linear bed: {F_peak(2.0) / (M_P * g):.1f} body weights from 2 m,"
          f" {F_peak(8.0) / (M_P * g):.1f} from 8 m")
    check("height falls as bed loss grows, all else equal",
          steady[(0.003, 0.3, 4.5)] > steady[(0.01, 0.3, 4.5)] > steady[(0.03, 0.3, 4.5)])
    check("force limit lowers the steady height", steady[(0.01, 0.3, 4.5)] < steady[(0.01, 0.3, np.inf)])

    print("A4. Timing: apex after 5 bounces from 1 m, by when the leg extension starts (bed zeta = 0.03)")
    print("     start (quarter periods after touchdown)  stroke midpoint  apex after 5 bounces (m)  leg work per bounce (J)")
    best = []
    stroke_time = 0.3 / 1.5
    rq = simulate_person(5, 1.0, zeta_bed=0.03, mode="rigid")
    t_bottom = rq["bottoms"][0]
    w_q = rq["w_q"]
    print(f"     bed bottom with stiff legs is {t_bottom / w_q:.2f} quarter periods after touchdown; the stroke lasts {stroke_time / w_q:.2f}")
    for delay in (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0):
        rr = simulate_person(11, 1.0, zeta_bed=0.03, mode="pump", delay=delay)
        n = min(5, len(rr["apexes"]) - 1)
        mid = delay + 0.5 * stroke_time / w_q
        best.append((delay, mid, rr["apexes"][n]))
        print(f"       {delay:4.2f}                                  {mid:5.2f}           {rr['apexes'][n]:6.2f}                  {rr['W'] / max(1, len(rr['apexes'])):7.0f}")
    bd = max(best, key=lambda x: x[2])
    worst = min(best, key=lambda x: x[2])
    check("best stroke straddles the bottom of the bed (largest force)",
          abs(bd[1] - t_bottom / w_q) < 0.7, f"best midpoint {bd[1]:.2f}, bottom {t_bottom / w_q:.2f}")
    check("timing changes the outcome by a large factor", bd[2] > 2 * max(worst[2], 0.05),
          f"{worst[2]:.2f} to {bd[2]:.2f} m")
    r_k = simulate_person(12, 3.0, zeta_bed=0.03, mode="kill")
    r_n = simulate_person(12, 3.0, zeta_bed=0.03, mode="rigid")
    print(f"     'killing' the bounce (shorten the leg as the bed pushes up), from 3 m:"
          f" apexes {np.round(r_k['apexes'][:4], 2)}, leg work {r_k['W'] / 1e3:.2f} kJ")
    print(f"     for comparison, stiff legs: apexes {np.round(r_n['apexes'][:4], 2)}")
    check("shortening the leg under load removes energy",
          r_k["W"] < 0 and r_k["apexes"][0] < r_n["apexes"][0] - 0.3)

    # ------------------------------------------------------------------ Part B
    print("\nPart B: ball on a circular membrane, solved as a wave equation")
    a, dr, sigma = 1.5, 0.01, 1.0
    r_foot = 0.15
    T_m = 5000.0 * (np.log(a / r_foot) + 0.25) / (2 * np.pi)  # sets K = 5000 N/m under a foot patch
    mem_foot = build_membrane(a, dr, T_m, sigma, r_foot)
    K_num = static_stiffness(mem_foot)
    print(f"     tension {T_m:.0f} N/m, areal density {sigma} kg/m^2, radius {a} m (assumed),"
          f" wave speed {mem_foot['c']:.1f} m/s")
    print(f"     static stiffness under a 15 cm foot patch {K_num:.0f} N/m (formula 2 pi T / (ln(a/r) + 1/4) = 5000)")
    check("numerical static stiffness matches the formula", abs(K_num / 5000 - 1) < 0.01)
    w1 = fundamental_freq(mem_foot)
    w1_exact = 2.40483 * mem_foot["c"] / a
    print(f"     lowest mode {w1:.2f} rad/s, Bessel-zero value {w1_exact:.2f} rad/s")
    check("membrane's lowest mode matches j01 c/a", abs(w1 / w1_exact - 1) < 0.005)
    mem = build_membrane(a, dr, T_m, sigma, 0.03)
    t_e = drop_ball(mem, 0.15, 2e4, 0.0, 1.0, gamma=0.0, gravity=False, T_max=0.08)
    print(f"     lossless, gravity-free run: final total energy / initial = {t_e['E_ratio']:.4f}")
    check("wave solver conserves energy with no damping", abs(t_e["E_ratio"] - 1) < 0.01)

    print("     Named balls from 1 m. 'Floor' uses the same code with a rigid floor.")
    print("     ball                  m(kg)  k_c(kN/m)  e_floor  e_bed  contact: floor(ms)  bed(ms)  bed mass reached/m  KE left in bed at lift-off")
    balls = [  # name, m, k, e0, patch radius
        ("superball", 0.05, 490e3, 0.85, 0.02),
        ("cricket ball", 0.156, 680e3, 0.55, 0.02),
        ("tennis ball", 0.058, 23e3, 0.75, 0.03),
        ("basketball", 0.60, 18e3, 0.80, 0.08),
        ("dead foam ball", 0.02, 1e3, 0.30, 0.03),
        ("medicine ball", 5.0, 50e3, 0.30, 0.10),
        ("person, lossy legs", 60.0, 2e4, 0.30, 0.15),
    ]
    named, floor_ok, budget = {}, True, []
    h0 = 1.0
    for name, m, k, e0, rp in balls:
        mem = build_membrane(a, dr, T_m, sigma, rp)
        z = e_to_zeta(e0)
        fl0 = drop_ball(mem, m, k, z, h0, rigid=True, gravity=False, T_max=0.5)
        floor_ok &= abs(fl0["e_v"] - e0) < 0.01
        fl = drop_ball(mem, m, k, z, h0, rigid=True)
        bd = drop_ball(mem, m, k, z, h0)
        reach = sigma * np.pi * (mem["c"] * fl["contact"]) ** 2
        named[name] = (fl["e"], bd["e"], bd["bed_share"], m)
        budget.append(bd["e"] ** 2 + bd["bed_share"] + bd["heat_share"])
        print(f"     {name:20s} {m:6.3f}   {k / 1e3:6.0f}    {fl['e']:6.3f}  {bd['e']:6.3f}   {1e3 * fl['contact']:8.2f}       {1e3 * bd['contact']:7.1f}   {reach / m:8.2f}          {100 * bd['bed_share']:5.1f} %")
    check("gravity-free floor runs recover each ball's intended e to 0.01", floor_ok)
    print(f"     energy budget, rebound e^2 + bed share + heat in ball, each ball: {np.round(budget, 3)}")
    check("rebound, bed and ball shares add to about 1 (bed damping and gravity in contact aside)",
          all(0.8 < b < 1.05 for b in budget))
    check("balls do bounce on the bed", all(named[n][1] > 0.1 for n in ("tennis ball", "basketball", "cricket ball")))
    check("a person with lossy legs bounces better on the bed than on a floor",
          named["person, lossy legs"][1] > named["person, lossy legs"][0] + 0.2)
    check("a tennis ball and a superball bounce less on the bed than on a floor",
          named["tennis ball"][1] < named["tennis ball"][0] and named["superball"][1] < named["superball"][0])

    print("     Mass scan: k_c = 20 kN/m, e_floor = 0.75, 5 cm patch")
    print("     m (kg)    e_floor   e_bed    contact (ms)   KE left in bed at lift-off")
    mem = build_membrane(a, dr, T_m, sigma, 0.05)
    z = e_to_zeta(0.75)
    scan = []
    for m in (0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 60.0):
        bd = drop_ball(mem, m, 2e4, z, h0)
        fl = drop_ball(mem, m, 2e4, z, h0, rigid=True)
        scan.append((m, bd["e"]))
        print(f"     {m:6.2f}    {fl['e']:6.3f}   {bd['e']:6.3f}   {1e3 * bd['contact']:8.1f}       {100 * bd['bed_share']:5.1f} %")
    es = np.array([s[1] for s in scan])
    i_min = int(np.argmin(es))
    print(f"     lowest bounce at m = {scan[i_min][0]} kg (e_bed = {es[i_min]:.3f}); highest at m = {scan[int(np.argmax(es))][0]} kg (e_bed = {es.max():.3f})")
    check("in the heavy range, more mass does not lose the bed's rebound advantage",
          es[-1] > 0.75 or es[-2] > 0.75)
    check("rebound rises from light masses to about 10 kg", es[6] > es[0] + 0.3)

    print("     Is the jump set by the echo from the rim? Mass where e_bed crosses 0.75, three bed radii.")
    print("     Prediction: contact time pi sqrt(m/K_s) equals the wave return time 2a/c")
    print("     radius (m)  return time (ms)  predicted m* (kg)  simulated m* (kg)")
    ok_echo = True
    z = e_to_zeta(0.75)
    for a_i in (1.0, 1.5, 2.25):
        mem = build_membrane(a_i, dr, T_m, sigma, 0.05)
        Kp = static_stiffness(mem)
        Ks = Kp * 2e4 / (Kp + 2e4)
        t_ret = 2 * a_i / mem["c"]
        m_pred = Ks * (t_ret / np.pi) ** 2
        lo, hi = 0.05 * m_pred, 20 * m_pred
        for _ in range(10):
            mid = np.sqrt(lo * hi)
            if drop_ball(mem, mid, 2e4, z, h0)["e"] > 0.75:
                hi = mid
            else:
                lo = mid
        m_sim = np.sqrt(lo * hi)
        ok_echo &= 0.3 < m_sim / m_pred < 3
        print(f"     {a_i:6.2f}       {1e3 * t_ret:7.1f}          {m_pred:8.2f}          {m_sim:8.2f}")
        if a_i == 1.0:
            m_first = m_sim
        if a_i == 2.25:
            m_last = m_sim
    check("threshold mass within a factor 3 of the echo prediction at all radii", ok_echo)
    check("threshold mass grows with bed radius, as an echo requires", m_last > 2 * m_first,
          f"{m_first:.2f} -> {m_last:.2f} kg")

    print("     Sensitivity of the tennis ball to the assumed bed density:")
    for sig_i in (0.5, 1.0, 2.0):
        e_i = drop_ball(build_membrane(a, dr, T_m, sig_i, 0.03), 0.058, 23e3, e_to_zeta(0.75), h0)["e"]
        print(f"       areal density {sig_i} kg/m^2: e_bed = {e_i:.3f}")

    print("     Grid check on the tennis ball, dr = 1 cm vs 0.5 cm:")
    m, k, e0, rp = 0.058, 23e3, 0.75, 0.03
    z = e_to_zeta(e0)
    e_coarse = drop_ball(build_membrane(a, 0.01, T_m, sigma, rp), m, k, z, h0)["e"]
    e_fine = drop_ball(build_membrane(a, 0.005, T_m, sigma, rp), m, k, z, h0)["e"]
    print(f"     e_bed {e_coarse:.3f} vs {e_fine:.3f}")
    check("tennis-ball result stable under grid refinement", abs(e_coarse - e_fine) < 0.03)

    print()
    if FAILS:
        print("FAILED:", FAILS)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
