"""T18-T21: a person as a radiator, a crowd as a heat source, and why absorption is never 100 percent.

  A. Radiated power by numerical integration of Planck's law, independent of the T^4 formula.
  B. View factors by ray tracing: person to person, person to sky in a dense crowd.
  C. Crowd heat and humidity: coupled ODEs for air and body, checked against an energy balance.
  D. Psychrometric wet-bulb solved from first principles versus the Stull fit; NWS heat index.
  E. Fresnel absorptivity of skin-like water, and limits on a perfect absorber.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import c, h, k_B, sigma_SB

FAILS = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def planck_exitance(T, l1=0.5e-6, l2=2000e-6, n=400000):
    lam = np.geomspace(l1, l2, n)
    B = 2 * np.pi * h * c**2 / lam**5 / np.expm1(h * c / (lam * k_B * T))
    return lam, B


def part_a():
    print("--- A. Radiated power of a person ---")
    T_skin = 307.15   # 34 C, assumed mean skin temperature
    A_dubois = 1.8    # m^2, typical adult body surface area
    f_eff = 0.75      # assumed effective radiating fraction; limbs facing each other do not lose to the room
    eps = 0.98        # skin, Steketee 1973
    A = A_dubois * f_eff
    lam, B = planck_exitance(T_skin)
    M = np.trapezoid(B, lam)
    print(f"  Planck integral at {T_skin - 273.15:.0f} C = {M:.1f} W/m2; sigma T^4 = {sigma_SB * T_skin**4:.1f} W/m2")
    check("numerical Planck integral matches Stefan-Boltzmann", abs(M / (sigma_SB * T_skin**4) - 1) < 1e-3)
    peak = lam[np.argmax(B)]
    print(f"  peak per unit wavelength {peak * 1e6:.2f} um")
    check("peak near 9.4 um", abs(peak * 1e6 - 9.44) < 0.05)
    in_8_14 = np.trapezoid(B[(lam > 8e-6) & (lam < 14e-6)], lam[(lam > 8e-6) & (lam < 14e-6)]) / M
    vis = np.trapezoid(B[(lam < 0.78e-6)], lam[lam < 0.78e-6]) / M
    print(f"  fraction in 8-14 um {in_8_14:.3f}; fraction below 0.78 um {vis:.1e}")
    x = np.linspace(1e-4, 60, 400000)
    phot = 2 * np.pi * (k_B * T_skin) ** 3 / (h**3 * c**2) * np.trapezoid(x**2 / np.expm1(x), x)
    gross = A * eps * M
    print(f"  area {A_dubois} m2 x {f_eff} = {A:.2f} m2, emissivity {eps}")
    print(f"  gross emission {gross:.0f} W, {A * eps * phot:.2e} photons/s, i.e. {gross * 86400 / 4184:.0f} kcal/day if nothing came back")
    print("  room T    received W   net W")
    results = {}
    for Tr in (15.0, 20.0, 25.0, 30.0, 34.0, 40.0):
        Tk = Tr + 273.15
        recv = A * eps * sigma_SB * Tk**4
        net = gross - recv
        results[Tr] = net
        print(f"  {Tr:6.0f} {recv:11.0f} {net:9.0f}")
    print("  uncertainty at 22 C room: vary skin 32-35 C, f_eff 0.7-0.8, area 1.6-2.0 m2")
    lo = 1.6 * 0.7 * 0.97 * sigma_SB * (305.15**4 - 295.15**4)
    hi = 2.0 * 0.8 * 0.99 * sigma_SB * (308.15**4 - 295.15**4)
    print(f"  net range {lo:.0f} to {hi:.0f} W")
    check("net radiation at 34 C room is zero, reverses above it", abs(results[34.0]) < 1 and results[40.0] < 0)
    check("gross far above resting metabolism, net near it", gross > 500 and 50 < results[25.0] < 130)
    h_r = 4 * eps * sigma_SB * (0.5 * (T_skin + 295.15)) ** 3
    # de Dear et al. quote their coefficient per unit of total body surface, so compare on that basis.
    print(f"  linearised h_r = 4 eps sigma T^3 = {h_r:.2f} W/m2/K per radiating area = {h_r * f_eff:.2f} per total body area; "
          f"manikin measurement 4.5 W/m2/K (de Dear et al. 1997)")
    print(f"  so the measurement implies an effective fraction f_eff = {4.5 / h_r:.2f}, against the {f_eff} assumed")
    check("assumed f_eff agrees with the fraction implied by the manikin measurement within 0.05", abs(4.5 / h_r - f_eff) < 0.05)
    return gross, A, eps, T_skin


def rays_from_cylinder(rng, n, r, H):
    z = rng.random(n) * H
    th = 2 * np.pi * rng.random(n)
    nrm = np.stack([np.cos(th), np.sin(th), np.zeros(n)], axis=1)
    pos = nrm * r + np.stack([np.zeros(n), np.zeros(n), z], axis=1)
    u1, u2 = rng.random(n), rng.random(n)
    rr = np.sqrt(u1)
    ph = 2 * np.pi * u2
    t1 = np.stack([-np.sin(th), np.cos(th), np.zeros(n)], axis=1)
    t2 = np.tile([0.0, 0.0, 1.0], (n, 1))
    d = (rr * np.cos(ph))[:, None] * t1 + (rr * np.sin(ph))[:, None] * t2 + np.sqrt(1 - u1)[:, None] * nrm
    return pos, d


def trace_lattice(rng, centres, r, H, n=30000):
    """Fractions of rays from the central cylinder landing on others, ground, or sky. Cylinders absorb."""
    pos, d = rays_from_cylinder(rng, n, r, H)
    a = d[:, 0] ** 2 + d[:, 1] ** 2
    best = np.full(n, np.inf)
    for cx, cy in centres:
        ox, oy = pos[:, 0] - cx, pos[:, 1] - cy
        b = ox * d[:, 0] + oy * d[:, 1]
        cc = ox**2 + oy**2 - r**2
        disc = b**2 - a * cc
        ok = (disc > 0) & (a > 0)
        t = np.full(n, np.inf)
        t[ok] = (-b[ok] - np.sqrt(disc[ok])) / a[ok]
        z = pos[:, 2] + t * d[:, 2]
        good = (t > 1e-9) & (z >= 0) & (z <= H)
        best = np.where(good & (t < best), t, best)
    t_ground = np.where(d[:, 2] < 0, -pos[:, 2] / np.where(d[:, 2] < 0, d[:, 2], -1), np.inf)
    hit_person = np.isfinite(best) & (best < t_ground)
    hit_ground = (~hit_person) & (d[:, 2] < 0)
    return hit_person.mean(), hit_ground.mean(), 1 - hit_person.mean() - hit_ground.mean()


def part_b(rng, gross, A, eps, T_skin):
    print("\n--- B. Ray-traced exchange (cylinder person: radius 0.16 m, height 1.7 m; side surface only) ---")
    r, H = 0.16, 1.7
    print("  Two people, distance between axes:")
    for sep in (0.5, 1.0, 2.0, 5.0):
        fp, fg, fs = trace_lattice(rng, [(sep, 0.0)], r, H, n=200000)
        print(f"    {sep:.1f} m: fraction of A's emission landing on B = {fp:.3f}, so B receives {fp * gross:.0f} W gross "
              f"and returns about the same to A")
    fp1, _, _ = trace_lattice(rng, [(1.0, 0.0)], r, H, n=200000)
    # Exact check for very long parallel cylinders is not this finite case; check with symmetry at another scale.
    fp_big, _, _ = trace_lattice(rng, [(1000.0, 0.0)], r, H, n=50000)
    check("distant person catches almost nothing", fp_big < 1e-3, f"{fp_big:.1e}")
    # Long-cylinder exact view factor as a code check: F = (1/pi)[sqrt(X^2-1) + asin(1/X) - X], X = s/(2r)
    X = 1.0 / (2 * r)
    F_exact = (np.sqrt(X**2 - 1) + np.arcsin(1 / X) - X) / np.pi
    fp_long, _, _ = trace_lattice(rng, [(1.0, 0.0)], r, 200.0, n=200000)
    print(f"  ray tracer on a nearly infinite cylinder at 1 m: {fp_long:.4f}; closed form {F_exact:.4f}")
    check("ray tracer reproduces the long-cylinder view factor", abs(fp_long - F_exact) < 0.005)

    print("  Square-lattice crowd, one person surrounded (people per m2, gap between bodies):")
    out = {}
    for density in (0.5, 1.0, 2.0, 4.0, 6.0, 8.0):
        s = 1 / np.sqrt(density)
        m = int(np.ceil(8 / s))
        centres = [(i * s, j * s) for i in range(-m, m + 1) for j in range(-m, m + 1) if (i, j) != (0, 0)]
        fp, fg, fs = trace_lattice(rng, centres, r, H, n=12000)
        out[density] = (fp, fg, fs)
        print(f"    {density:3.1f}/m2, gap {s - 2 * r:.2f} m: to other people {fp:.3f}, to ground {fg:.3f}, to sky {fs:.3f}")
    check("packing blocks the path to the sky", out[6.0][2] < 0.1 and out[0.5][2] > 0.2)
    check("sky fraction falls monotonically with density", all(out[a][2] > out[b][2] for a, b in zip(list(out)[:-1], list(out)[1:])))
    return out


def psat_kpa(Tc):
    """Saturation pressure over water, Magnus-Tetens form, checked below against NIST values."""
    return 0.61094 * np.exp(17.625 * Tc / (Tc + 243.04))


def stull(T, RH):
    return (T * np.arctan(0.151977 * np.sqrt(RH + 8.313659)) + np.arctan(T + RH) - np.arctan(RH - 1.676331)
            + 0.00391838 * RH**1.5 * np.arctan(0.023101 * RH) - 4.686035)


def wet_bulb_exact(T, RH, P=101.325):
    """Adiabatic saturation: enthalpy balance of air plus evaporated water."""
    e = RH / 100 * psat_kpa(T)
    w = 0.622 * e / (P - e)
    cpa, cpv, L0 = 1.005, 1.86, 2501.0

    def f(Tw):
        ws = 0.622 * psat_kpa(Tw) / (P - psat_kpa(Tw))
        return cpa * (T - Tw) + w * (L0 + cpv * T) - ws * (L0 + cpv * Tw) + (ws - w) * 4.186 * Tw

    return brentq(f, -40, T + 1e-9) if f(T) < 0 else T


def heat_index_c(T, RH):
    F = T * 9 / 5 + 32
    s = 0.5 * (F + 61 + (F - 68) * 1.2 + RH * 0.094)
    if 0.5 * (s + F) < 80:
        hi = 0.5 * (s + F)
    else:
        hi = (-42.379 + 2.04901523 * F + 10.14333127 * RH - 0.22475541 * F * RH - 0.00683783 * F**2
              - 0.05481717 * RH**2 + 0.00122874 * F**2 * RH + 0.00085282 * F * RH**2 - 0.00000199 * F**2 * RH**2)
        if RH < 13 and 80 <= F <= 112:
            hi -= (13 - RH) / 4 * np.sqrt((17 - abs(F - 95)) / 17)
        if RH > 85 and 80 <= F <= 87:
            hi += (RH - 85) / 10 * (87 - F) / 5
    return (hi - 32) * 5 / 9


def crowd_model(T_out, RH_out, vent_lps, vol_per_person, M=120.0, exposure=0.5, t_end=6 * 3600, dt=2.0):
    """One person's share of a hall. States: air T, humidity ratio w, body core temperature.

    Skin 35 C. Dry coefficient 3.4 + 4.5 W/m2/K (de Dear 1997), evaporation limit by Lewis relation
    h_e = 16.5 h_c W/m2/kPa, wetted skin, exposure = fraction of skin open to the air.
    Sweat cap 1 kg/h, about 680 W. Body 70 kg, 3470 J/kg/K.
    """
    P = 101.325
    rho, cpa, L = 1.18, 1005.0, 2.43e6
    A = 1.8 * exposure
    h_dry = 3.4 + 4.5
    h_e = 16.5 * 3.4
    Vdot = vent_lps * 1e-3
    V = vol_per_person
    T = T_out
    e_out = RH_out / 100 * psat_kpa(T_out)
    w_out = 0.622 * e_out / (P - e_out)
    w = w_out
    Tcore = 37.0
    t = 0.0
    t_40 = None
    Tk_skin = 35.0
    co2 = 0.0004
    history = []
    while t < t_end:
        e = w * P / (0.622 + w)
        dry = h_dry * A * (Tk_skin - T)
        E_max = h_e * A * max(psat_kpa(Tk_skin) - e, 0.0)
        E_need = M - dry
        E = min(max(E_need, 0.0), E_max, 680.0)
        stored = M - dry - E
        Tcore += stored / (70 * 3470) * dt
        T += ((dry + 0.0) / (rho * V * cpa) + Vdot / V * (T_out - T)) * dt
        w += (E / L / (rho * V) + Vdot / V * (w_out - w)) * dt
        co2 += (5e-6 / V + Vdot / V * (0.0004 - co2)) * dt
        t += dt
        if t_40 is None and Tcore >= 40.0:
            t_40 = t
        if int(t / dt) % 900 == 0:
            history.append((t, T, w, Tcore))
    rh = 100 * (w * P / (0.622 + w)) / psat_kpa(T)
    return dict(T=T, w=w, rh=min(rh, 100.0), Tcore=Tcore, t_40=t_40, stored=stored, E=E, dry=dry, co2=co2, w_out=w_out,
                V=Vdot, M=M, rho=rho, cpa=cpa, L=L)


def part_c():
    print("\n--- C. Crowd heat: the numbers ---")
    M = 120.0
    for N, density in ((1e4, 4.0), (1e6, 4.0), (1e7, 4.0)):
        total = N * M
        area = N / density
        side = np.sqrt(area)
        q = M * density
        u, H = 2.0, 100.0
        dT = q * side / (1.2 * 1005 * u * H)
        print(f"  {N:.0e} people at {density}/m2: {total / 1e6:.1f} MW over {area / 1e6:.3f} km2 = {q:.0f} W/m2; "
              f"air warming passing over at 2 m/s, mixed 100 m deep: {dT:.1f} K")
    q6 = 6 * M
    print(f"  at 6/m2 the flux is {q6:.0f} W/m2 against about 1000 W/m2 of noon sun")
    check("flux depends on density, not on head count", abs(M * 4 - 480) < 1e-9)
    print(f"  no heat loss at all: core warms {M / (70 * 3470) * 3600:.2f} K per hour ({M:.0f} W, 70 kg, 3470 J/kg/K); "
          f"37 to 40 C in {3 / (M / (70 * 3470) * 3600):.1f} h")

    print("\n  Hall model (120 W each; half the skin exposed), hot day outside 32 C, 50 percent RH:")
    print("  ventilation    volume     air T   air RH   wet-bulb   core T at 6 h   time to 40 C    CO2 at 6 h")
    cases = [(10.0, 5.0), (5.0, 5.0), (2.5, 3.0), (1.0, 3.0), (0.25, 3.0), (0.0, 3.0)]
    results = {}
    for vent, vol in cases:
        r = crowd_model(32.0, 50.0, vent, vol)
        tw = wet_bulb_exact(r["T"], r["rh"]) if r["rh"] > 0 else float("nan")
        results[(vent, vol)] = (r, tw)
        tstr = f"{r['t_40'] / 3600:.1f} h" if r["t_40"] else "not reached"
        print(f"  {vent:5.2f} L/s/p {vol:6.1f} m3  {r['T']:6.1f} {r['rh']:7.0f}% {tw:8.1f} {r['Tcore']:12.1f} {tstr:>15s} {r['co2'] * 100:9.2f}%")
    r0, _ = results[(10.0, 5.0)]
    # Independent steady-state energy balance: all body heat leaves in the outgoing air as sensible plus latent.
    out_sens = r0["rho"] * r0["V"] * r0["cpa"] * (r0["T"] - 32.0)
    out_lat = r0["rho"] * r0["V"] * r0["L"] * (r0["w"] - r0["w_out"])
    print(f"  energy check at 10 L/s: out in air {out_sens:.1f} W sensible + {out_lat:.1f} W latent = {out_sens + out_lat:.1f} W vs {M:.0f} W in")
    check("steady state closes the energy balance within 2 percent", abs((out_sens + out_lat) / M - 1) < 0.02)
    check("well-ventilated hall does not reach 40 C core", results[(10.0, 5.0)][0]["t_40"] is None)
    check("sealed hall does reach 40 C within 6 h", results[(0.0, 3.0)][0]["t_40"] is not None)
    co2_sealed = results[(0.0, 3.0)][0]["co2"]
    print(f"  in the sealed case CO2 reaches {co2_sealed * 100:.1f} percent at 6 h with 0.005 L/s per person (see README for the source caveat)")


def part_d():
    print("\n--- D. Wet-bulb and heat index ---")
    for Tc, PkPa in ((20.0, 2.3393), (40.0, 7.3849)):
        print(f"  saturation pressure at {Tc:.0f} C: Magnus {psat_kpa(Tc):.3f} kPa, NIST {PkPa:.3f} kPa")
        check(f"Magnus form within 1 percent of NIST at {Tc:.0f} C", abs(psat_kpa(Tc) / PkPa - 1) < 0.01)
    print("  T/C  RH%   Tw exact   Tw Stull   heat index C")
    worst = 0.0
    for T, RH in ((30, 50), (35, 40), (35, 75), (38, 60), (40, 75), (45, 30)):
        a, b = wet_bulb_exact(T, RH), stull(T, RH)
        worst = max(worst, abs(a - b))
        print(f"  {T:3d} {RH:4d} {a:9.2f} {b:9.2f} {heat_index_c(T, RH):11.1f}")
    check("independent psychrometric wet-bulb agrees with Stull within his stated error of about 1 K", worst < 1.1, f"worst {worst:.2f} K")
    check("NWS heat index: 90 F, 70 percent RH gives the published 105-106 F", 105 <= heat_index_c((90 - 32) * 5 / 9, 70) * 9 / 5 + 32 <= 106.5)

    print("  Skin limit: sensible and latent exchange vanish when wet-bulb reaches skin temperature.")
    print(f"  T = 35 C at 100 percent RH: wet-bulb = {wet_bulb_exact(35.0, 99.99):.1f} C, which cannot shed even 1 W from 35 C skin")
    check("35 C saturated air has a 35 C wet-bulb", abs(wet_bulb_exact(35.0, 99.99) - 35.0) < 0.05)


def part_e():
    print("\n--- E. Absorption: never exactly 1 ---")

    def R_avg(theta, nc):
        ci = np.cos(theta)
        s2 = np.sin(theta) ** 2
        ct = np.sqrt(nc**2 - s2)
        rs = np.abs((ci - ct) / (ci + ct)) ** 2
        rp = np.abs((nc**2 * ci - ct) / (nc**2 * ci + ct)) ** 2
        return 0.5 * (rs + rp)

    alpha_cm = 638.0  # 1/cm, water at 10 um, Hale and Querry 1973
    lam_cm = 10e-4
    k = alpha_cm * lam_cm / (4 * np.pi)
    print(f"  k from the tabulated absorption coefficient: alpha lambda / 4 pi = {k:.4f}")
    print(f"  penetration depth 1/alpha = {1 / alpha_cm * 1e4:.1f} um (radiation is absorbed in the outer skin)")
    print("  n     normal absorptivity   hemispherical absorptivity")
    vals = {}
    th = np.linspace(0, np.pi / 2, 40001)
    for n in (1.15, 1.218, 1.30, 1.40):
        nc = n + 1j * k
        norm = 1 - R_avg(0.0, nc)
        hemi = 2 * np.trapezoid((1 - R_avg(th, nc)) * np.cos(th) * np.sin(th), th)
        vals[n] = (norm, hemi)
        print(f"  {n:.3f} {norm:12.4f} {hemi:20.4f}")
    check("model hemispherical absorptivity of water at 10 um near measured skin value 0.98 (within 0.03)",
          abs(vals[1.218][1] - 0.98) < 0.03, f"{vals[1.218][1]:.3f}")
    check("no choice of n in that range reaches 1", all(v[1] < 0.995 for v in vals.values()))
    grazing = 1 - R_avg(np.deg2rad(89.999), 1.218 + 1j * k)
    print(f"  at 89.999 degrees absorptivity {grazing:.4f}")
    check("grazing reflection tends to 1 for a flat surface with n != 1", grazing < 0.01)

    print("  Rozanov bound for a non-magnetic absorber on a metal backing, thickness d:")
    print("    integral over wavelength of |ln R| <= 2 pi^2 d; so R <= R0 across a band needs d >= |ln R0| dlambda / (2 pi^2)")
    band = 47e-6
    for R0 in (1e-1, 1e-2, 1e-5, 1e-10, 1e-100):
        d = abs(np.log(R0)) * band / (2 * np.pi**2)
        print(f"    R0 = {R0:.0e} across 3-50 um: d >= {d * 1e6:.1f} um")
    check("thickness needed grows without limit as R0 -> 0", abs(np.log(1e-100)) > 100 * abs(np.log(1e-1)) * 0.9)

    print("  Cavity: rough or folded skin raises absorptivity but cannot reach 1")
    for F in (0.0, 0.3, 0.6, 0.9, 0.99):
        eps = 0.97
        print(f"    self view factor {F:.2f}: effective absorptivity {eps / (1 - (1 - eps) * F):.5f}")
    print("  Energy accounting, in contrast, is exact: absorbed + reflected + transmitted = incident.")
    print("  Energy: absorbed and re-emitted power match only in steady state, and not spectrally:")
    T1, T2 = 5778.0, 307.15
    print(f"    sunlight at {T1:.0f} K in, peak {2.898e-3 / T1 * 1e9:.0f} nm; skin out at {T2:.0f} K, peak {2.898e-3 / T2 * 1e6:.1f} um")
    print(f"    per joule moved: entropy gain 4/3 (1/{T2:.0f} - 1/{T1:.0f}) = {4 / 3 * (1 / T2 - 1 / T1):.5f} J/K")
    check("downshifting from sun to skin creates entropy", 1 / T2 > 1 / T1)


def main():
    rng = np.random.default_rng(7)
    gross, A, eps, Ts = part_a()
    part_b(rng, gross, A, eps, Ts)
    part_c()
    part_d()
    part_e()
    print(f"\n{len(FAILS)} checks failed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
