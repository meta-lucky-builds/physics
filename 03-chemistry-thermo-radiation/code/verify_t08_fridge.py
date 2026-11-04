"""T08: vapour-compression cycle from tabulated refrigerant properties (NIST Chemistry WebBook fluid data).

All properties below are copied from NIST tables; reference state IIR, H = 200 kJ/kg and S = 1 kJ/kg/K
for saturated liquid at 0 C. Only differences are used.
R-134a equation of state: Tillner-Roth and Baehr (1994). Isobutane data from the same NIST service.
"""
import sys

import numpy as np

FAILS = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


# (P/MPa): (T_sat/C, h_l, s_l, h_v, s_v)
R134A_SAT = {
    0.10: (-26.361, 165.44, 0.86756, 382.60, 1.7475),
    0.20: (-10.076, 186.60, 0.95027, 392.62, 1.7334),
    0.30: (0.6721, 200.90, 1.0033, 398.99, 1.7267),
    0.40: (8.9306, 212.11, 1.0433, 403.72, 1.7226),
    0.50: (15.735, 221.50, 1.0759, 407.47, 1.7197),
    0.60: (21.572, 229.68, 1.1037, 410.57, 1.7175),
    0.70: (26.713, 236.99, 1.1280, 413.20, 1.7156),
    0.80: (31.328, 243.65, 1.1497, 415.46, 1.7140),
    0.90: (35.526, 249.78, 1.1695, 417.43, 1.7126),
    1.00: (39.388, 255.50, 1.1876, 419.16, 1.7113),
    1.10: (42.969, 260.87, 1.2044, 420.69, 1.7100),
    1.20: (46.315, 265.95, 1.2201, 422.04, 1.7087),
}
# superheated vapour at 1.0 MPa: (T/C, h, s)
R134A_1MPA = [(39.388, 419.16, 1.7113), (40.0, 419.86, 1.7135), (45.0, 425.44, 1.7312), (50.0, 430.88, 1.7482),
              (55.0, 436.24, 1.7646), (60.0, 441.53, 1.7806)]

# Isobutane (R-600a)
ISOB_LOW = (0.10, -12.085, 172.75, 0.89833, 538.15, 2.2980)
ISOB_HIGH_SAT = (0.55, 41.329, 299.66, 1.3369, 609.55, 2.3223)
ISOB_055 = [(41.329, 609.55, 2.3223), (45.0, 616.65, 2.3448), (50.0, 626.36, 2.3751), (55.0, 636.11, 2.4050), (60.0, 645.92, 2.4347)]


def interp_h(table, s):
    ss = [row[2] for row in table]
    hh = [row[1] for row in table]
    tt = [row[0] for row in table]
    if s < ss[0] or s > ss[-1]:
        raise ValueError("s outside table")
    return float(np.interp(s, ss, hh)), float(np.interp(s, ss, tt))


def cycle(name, T_e, h1, s1, s_l_e, h_l_e, T_c, h3, s3, h2s, T2s, eta):
    q_l = h1 - h3
    w_s = h2s - h1
    w_real = w_s / eta
    h2 = h1 + w_real
    q_h = h2 - h3
    cop_ideal = q_l / w_s
    cop_real = q_l / w_real
    TL, TH = T_e + 273.15, T_c + 273.15
    cop_carnot = TL / (TH - TL)
    x4 = (h3 - h_l_e) / (h1 - h_l_e)
    s4 = s_l_e + x4 * (s1 - s_l_e)
    print(f"  {name}: evaporate {T_e:.1f} C, condense {T_c:.1f} C")
    print(f"    q_L = {q_l:.1f} kJ/kg, w(isentropic) = {w_s:.1f} kJ/kg, T after compression {T2s:.1f} C")
    print(f"    flash fraction after expansion valve x4 = {x4:.3f}, entropy gain in valve {s4 - s3:.4f} kJ/kg/K")
    print(f"    COP Carnot {cop_carnot:.2f}, COP ideal cycle {cop_ideal:.2f}, COP with eta = {eta} compressor {cop_real:.2f}")
    return dict(q_l=q_l, w_s=w_s, w=w_real, q_h=q_h, cop_c=cop_carnot, cop_i=cop_ideal, cop_r=cop_real, ds_valve=s4 - s3,
                TL=TL, TH=TH)


def main():
    print("--- Boiling point claim ---")
    print(f"  R-134a boils at {R134A_SAT[0.10][0]:.1f} C at 0.10 MPa (about atmospheric). Water boils at 99.6 C at that pressure.")
    check("R-134a boils far below cabinet temperature at atmospheric pressure", R134A_SAT[0.10][0] < -20)

    print("\n--- Consistency of the copied table ---")
    print("  At a phase change Gibbs energy is equal in both phases, so h_v - h_l = T (s_v - s_l) exactly.")
    Ps = sorted(R134A_SAT)
    worst = 0.0
    for P in Ps:
        T, hl, sl, hv, sv = R134A_SAT[P]
        lhs, rhs = hv - hl, (T + 273.15) * (sv - sl)
        worst = max(worst, abs(lhs / rhs - 1))
    print(f"  worst mismatch over {len(Ps)} R-134a rows: {100 * worst:.3f} percent")
    check("latent heat equals T ds in every R-134a row within 0.1 percent", worst < 1e-3)
    for label, row in (("isobutane 0.10 MPa", ISOB_LOW), ("isobutane 0.55 MPa", ISOB_HIGH_SAT)):
        _, T, hl, sl, hv, sv = row
        m = abs((hv - hl) / ((T + 273.15) * (sv - sl)) - 1)
        check(f"{label}: latent heat equals T ds within 0.1 percent", m < 1e-3, f"{100 * m:.3f} percent")
    print("  Ideal-vapour Clausius-Clapeyron, an approximation, for comparison (M = 102.03 g/mol):")
    errs = []
    for P1, P2 in ((0.1, 0.2), (0.5, 0.6), (1.0, 1.1)):
        T1, T2 = R134A_SAT[P1][0] + 273.15, R134A_SAT[P2][0] + 273.15
        L_cc = np.log(P2 / P1) * 8.314 / (1 / T1 - 1 / T2) / 102.03
        L_table = 0.5 * (R134A_SAT[P1][3] - R134A_SAT[P1][1] + R134A_SAT[P2][3] - R134A_SAT[P2][1])
        errs.append(L_cc / L_table - 1)
        print(f"    {P1:.1f}-{P2:.1f} MPa: L = {L_cc:.0f} kJ/kg, table {L_table:.0f} kJ/kg, overestimate {100 * errs[-1]:.0f} percent")
    check("ideal-vapour estimate worsens with pressure, as a real vapour nearing its critical point should", errs[0] < errs[1] < errs[2])

    print("\n--- Cycle calculations ---")
    res = {}
    eta = 0.70  # assumed isentropic efficiency, typical small hermetic compressor order of magnitude
    for label, P_e in (("R-134a fridge", 0.20), ("R-134a freezer", 0.10)):
        Te, hle, sle, hve, sve = R134A_SAT[P_e]
        Tc, h3, s3, _, _ = R134A_SAT[1.00]
        h2s, T2s = interp_h(R134A_1MPA, sve)
        res[label] = cycle(label, Te, hve, sve, sle, hle, Tc, h3, s3, h2s, T2s, eta)

    # Isobutane: isentropic compression from saturated vapour lands in the wet region at this pressure.
    _, Te, hle, sle, h1, s1 = ISOB_LOW
    _, Tc, h3, s3, hv3, sv3 = ISOB_HIGH_SAT
    if s1 < sv3:
        x2 = (s1 - s3) / (sv3 - s3)
        h2s = h3 + x2 * (hv3 - h3)
        T2s = Tc
        print(f"  note: isobutane from saturated vapour ends at quality {x2:.3f}; real compressors run on superheated suction gas")
    else:
        h2s, T2s = interp_h(ISOB_055, s1)
    res["R-600a"] = cycle("R-600a fridge", Te, h1, s1, sle, hle, Tc, h3, s3, h2s, T2s, eta)

    print("\n--- Checks on each cycle ---")
    for label, r in res.items():
        check(f"{label}: first law q_H = q_L + w", abs(r["q_h"] - r["q_l"] - r["w"]) < 1e-9)
        check(f"{label}: COP real < COP ideal cycle < COP Carnot", r["cop_r"] < r["cop_i"] < r["cop_c"],
              f"{r['cop_r']:.2f} < {r['cop_i']:.2f} < {r['cop_c']:.2f}")
        check(f"{label}: expansion valve generates entropy", r["ds_valve"] > 0)
        # Entropy of reservoirs per kg refrigerant at cabinet and room temperature, cabinet 4 C and room 25 C
        Tcab, Troom = 277.15, 298.15
        ds_universe = r["q_h"] / Troom - r["q_l"] / Tcab
        check(f"{label}: entropy of the surroundings rises over a cycle", ds_universe > 0, f"{ds_universe:.3f} kJ/kg/K per kg circulated")

    print("\n--- What if there were no compressor? ---")
    r = res["R-134a fridge"]
    ds_no_work = r["q_l"] / 298.15 - r["q_l"] / 277.15
    print(f"  Moving q_L = {r['q_l']:.0f} kJ/kg from a 4 C cabinet to a 25 C room with no work: dS = {ds_no_work:+.4f} kJ/kg/K")
    check("zero work requires negative entropy change, forbidden", ds_no_work < 0)
    w_min = r["q_l"] * (298.15 / 277.15 - 1)
    print(f"  smallest work that makes dS = 0: {w_min:.1f} kJ/kg (reversible machine between 4 C and 25 C), COP {r['q_l'] / w_min:.1f}")
    T_sat_low = R134A_SAT[0.20][0]
    print(f"  Vapour leaving the evaporator is at {T_sat_low:.1f} C and 0.20 MPa, colder than the room, so it could never give heat to the room.")
    Ps_arr = np.array(Ps)
    Ts_arr = np.array([R134A_SAT[p][0] for p in Ps])
    P_min = float(np.interp(25.0, Ts_arr, Ps_arr))
    print(f"  It only condenses against 25 C air if its saturation temperature exceeds 25 C: pressure above about {P_min:.2f} MPa")
    check("compressor must raise pressure by more than 3 times", P_min / 0.20 > 3, f"factor {P_min / 0.20:.1f}")

    print("\n--- The apple: 0.20 kg cooled from 30 C to 4 C ---")
    m_apple, c_apple = 0.20, 3550.0  # 85 percent water times 4180 J/kg/K
    Q = m_apple * c_apple * 26
    mr = Q / 1000 / r["q_l"]
    print(f"  heat to remove {Q / 1000:.1f} kJ; refrigerant circulated {mr * 1000:.0f} g; compressor work {mr * r['w']:.1f} kJ; "
          f"heat dumped to the room {mr * r['q_h']:.1f} kJ")
    check("room receives more heat than the apple lost", mr * r["q_h"] > Q / 1000)

    print(f"\n{len(FAILS)} checks failed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
