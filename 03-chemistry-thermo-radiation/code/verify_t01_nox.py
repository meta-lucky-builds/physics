"""T01: enthalpy and equilibrium of N2 + 2 O2 -> 2 NO2, and the real route through NO.

Inputs: NIST-JANAF Shomate parameters as listed in the NIST Chemistry WebBook
(Chase 1998), and textbook mean bond energies. Everything else is computed.

The equilibrium is solved twice by unrelated routes:
  (a) from equilibrium constants K(T) = exp(-dG/RT) with a root finder;
  (b) by directly minimising the total Gibbs energy of the mixture over the
      composition, with no K anywhere.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve, minimize

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import R

FAILS = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


# (Tmax, A, B, C, D, E, F, G) with t = T/1000. H(T)-H(298) + dfH(298) = ... + F
SHOMATE = {
    "N2": [(500, 28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 226.4168),
           (2000, 19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 212.3900),
           (6000, 35.51872, 1.128728, -0.196103, 0.014662, -4.553760, -18.97091, 224.9810)],
    "O2": [(700, 31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945),
           (2000, 30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663),
           (6000, 20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 237.6185)],
    "NO": [(1200, 23.83491, 12.58878, -1.139011, -1.497459, 0.214194, 83.35783, 237.1219),
           (6000, 35.99169, 0.957170, -0.148032, 0.009974, -3.004088, 73.10787, 246.1619)],
    "NO2": [(1200, 16.10857, 75.89525, -54.38740, 14.30777, 0.239423, 26.17464, 240.5386),
            (6000, 56.82541, 0.738053, -0.144721, 0.009777, -5.459911, 2.846456, 290.5056)],
}


def thermo(species, T):
    """Return H (kJ/mol, formation basis at 298 K) and S (J/mol/K) at T."""
    for row in SHOMATE[species]:
        if T <= row[0]:
            _, A, B, C, D, E, F, G = row
            break
    t = T / 1000.0
    H = A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F
    S = A * np.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G
    return H, S


def G0(species, T):
    H, S = thermo(species, T)
    return H * 1000.0 - T * S  # J/mol


def dG_NO(T):   # N2 + O2 -> 2 NO
    return 2 * G0("NO", T) - G0("N2", T) - G0("O2", T)


def dG_NO2(T):  # 2 NO + O2 -> 2 NO2
    return 2 * G0("NO2", T) - 2 * G0("NO", T) - G0("O2", T)


def dH_net(T):
    return 2 * thermo("NO2", T)[0] - thermo("N2", T)[0] - 2 * thermo("O2", T)[0]


def dH_NO(T):
    return 2 * thermo("NO", T)[0] - thermo("N2", T)[0] - thermo("O2", T)[0]


def equilibrium_from_K(T, n0=(0.79, 0.21)):
    """Mole fractions from K1, K2 at P = 1 bar, solved for ln(n_NO), ln(n_NO2)."""
    lnK1 = -dG_NO(T) / (R * T)
    lnK2 = -dG_NO2(T) / (R * T)
    N2_0, O2_0 = n0

    def amounts(p):
        u, v = np.exp(p)  # moles NO, NO2
        return np.array([N2_0 - (u + v) / 2, O2_0 - (u + v) / 2 - v / 2, u, v])

    def eqs(p):
        n = amounts(p)
        if np.any(n <= 0):
            return [1e6, 1e6]
        tot = n.sum()
        return [2 * np.log(n[2]) - np.log(n[0]) - np.log(n[1]) - lnK1,
                2 * np.log(n[3]) - 2 * np.log(n[2]) - np.log(n[1]) + np.log(tot) - lnK2]

    lnu0 = 0.5 * (lnK1 + np.log(N2_0 * O2_0))
    lnv0 = lnu0 + 0.5 * (lnK2 + np.log(O2_0))
    p, _, ier, msg = fsolve(eqs, [lnu0, lnv0], full_output=True, xtol=1e-12)
    if ier != 1:
        raise RuntimeError(f"no convergence at T={T}: {msg}")
    n = amounts(p)
    return n / n.sum()


def equilibrium_by_gibbs(T, n0=(0.79, 0.21)):
    """Minimise G_total(n) directly over two free composition variables."""
    g = np.array([G0(s, T) for s in ("N2", "O2", "NO", "NO2")]) / (R * T)
    N2_0, O2_0 = n0

    def G_total(p):
        u, v = np.exp(p[0]), np.exp(p[1])
        n = np.array([N2_0 - (u + v) / 2, O2_0 - (u + v) / 2 - v / 2, u, v])
        if np.any(n <= 0):
            return 1e30
        return float(np.sum(n * (g + np.log(n / n.sum()))))

    # Deliberately poor start: a factor e away from the K solution in each amount.
    x0 = equilibrium_from_K(T, n0)
    res = minimize(G_total, [np.log(x0[2]) + 1.0, np.log(x0[3]) - 1.0], method="Nelder-Mead",
                   options={"xatol": 1e-7, "fatol": 1e-15, "maxiter": 8000})
    u, v = np.exp(res.x)
    n = np.array([N2_0 - (u + v) / 2, O2_0 - (u + v) / 2 - v / 2, u, v])
    return n / n.sum()


def main():
    print("--- Standard enthalpies at 298.15 K ---")
    dfH = {"NO": 90.29, "NO2": 33.10, "N": 472.68, "O": 249.18}
    dH_formation = 2 * dfH["NO2"]
    dH_step1 = 2 * dfH["NO"]
    dH_step2 = 2 * dfH["NO2"] - 2 * dfH["NO"]
    print(f"N2 + 2 O2 -> 2 NO2 : dH = 2 x 33.10 = {dH_formation:+.1f} kJ per mole N2")
    print(f"step N2 + O2 -> 2 NO   : {dH_step1:+.1f} kJ")
    print(f"step 2 NO + O2 -> 2 NO2: {dH_step2:+.1f} kJ")
    check("steps sum to the net enthalpy (Hess)", abs(dH_step1 + dH_step2 - dH_formation) < 1e-9)
    check("Shomate fit reproduces dfH at 298.15 K", abs(dH_net(298.15) - dH_formation) < 0.1,
          f"{dH_net(298.15):+.2f} kJ")

    print("\n--- Bond energies from atomisation enthalpies (molecule-specific) ---")
    D_N2 = 2 * dfH["N"]
    D_O2 = 2 * dfH["O"]
    D_NO = dfH["N"] + dfH["O"] - dfH["NO"]
    atom_NO2 = dfH["N"] + 2 * dfH["O"] - dfH["NO2"]
    D_ONO = atom_NO2 - D_NO
    print(f"D(N2) = {D_N2:.1f}, D(O2) = {D_O2:.1f}, D(N-O in NO) = {D_NO:.1f} kJ/mol")
    print(f"NO2 atomisation = {atom_NO2:.1f} kJ/mol: first N-O bond (ONO->NO+O) {D_ONO:.1f}, "
          f"second {D_NO:.1f}, mean {atom_NO2 / 2:.1f}")
    broken = D_N2 + 2 * D_O2
    formed = 2 * atom_NO2
    print(f"bonds broken {broken:.1f} kJ, bonds formed {formed:.1f} kJ, net {broken - formed:+.1f} kJ")
    check("specific bond energies return the formation route", abs(broken - formed - dH_formation) < 0.01,
          "this is Hess's law, an identity, not an independent measurement")

    print("\n--- Same reaction from mean bond tables (N#N 945, O=O 498, N=O 607, N-O 201) ---")
    table_double = 945 + 2 * 498 - 4 * 607
    table_lewis = 945 + 2 * 498 - 2 * (607 + 201)
    print(f"treating NO2 as O=N=O  : {table_double:+d} kJ (wrong sign)")
    print(f"treating NO2 as O=N-O  : {table_lewis:+d} kJ (wrong size)")
    print(f"thermochemical value   : {dH_formation:+.1f} kJ")
    check("both Lewis-structure estimates from mean bond energies miss by more than 200 kJ",
          abs(table_double - dH_formation) > 200 and abs(table_lewis - dH_formation) > 200)

    print("\n--- Temperature dependence, per mole N2 ---")
    print("   T/K   dH(net)  dG(net)   dG(N2+O2->2NO)  dG(2NO+O2->2NO2)")
    temps = [298.15, 500, 1000, 1500, 2000, 2500, 3000]
    for T in temps:
        print(f"{T:7.1f} {dH_net(T):8.1f} {(dG_NO(T) + dG_NO2(T)) / 1000:8.1f} "
              f"{dG_NO(T) / 1000:12.1f} {dG_NO2(T) / 1000:16.1f}")
    T_scan = np.linspace(298.15, 5000, 800)
    net = np.array([dG_NO(T) + dG_NO2(T) for T in T_scan])
    check("net dG for N2 + 2 O2 -> 2 NO2 is positive at every T from 298 to 5000 K", np.all(net > 0),
          f"minimum {net.min() / 1000:.1f} kJ")
    no = np.array([dG_NO(T) for T in T_scan])
    check("dG(N2+O2->2NO) falls monotonically with T (entropy gain, weak)", np.all(np.diff(no) < 0))
    no2 = np.array([dG_NO2(T) for T in T_scan])
    T_cross = T_scan[np.argmax(no2 > 0)]
    print(f"2 NO + O2 -> 2 NO2 changes sign from favourable to unfavourable near {T_cross:.0f} K")
    check("NO -> NO2 favoured only when cool", no2[0] < 0 and no2[-1] > 0)

    # van't Hoff, numerically, as a check that G(T) and H(T) are mutually consistent.
    for T in (1000.0, 2000.0):
        eps = 1.0
        slope = (dG_NO(T + eps) / (T + eps) - dG_NO(T - eps) / (T - eps)) / (1 / (T + eps) - 1 / (T - eps))
        check(f"Gibbs-Helmholtz at {T:.0f} K", abs(slope / 1000 - dH_NO(T)) < 0.3,
              f"d(dG/T)/d(1/T) = {slope / 1000:.1f} vs dH = {dH_NO(T):.1f} kJ")

    print("\n--- Equilibrium composition of air (79 N2 : 21 O2, 1 bar) ---")
    print("   T/K    x(NO)       x(NO2)      x(NO2)/x(NO)")
    for T in [298.15, 500, 800, 1000, 1500, 1800, 2000, 2500, 3000]:
        x = equilibrium_from_K(T)
        print(f"{T:7.1f} {x[2]:10.3e} {x[3]:10.3e} {x[3] / x[2]:12.3e}")
    for T in (1000.0, 2000.0, 3000.0):
        a = equilibrium_from_K(T)
        b = equilibrium_by_gibbs(T)
        check(f"K route and direct Gibbs minimisation agree at {T:.0f} K",
              abs(a[2] - b[2]) / a[2] < 1e-3 and abs(a[3] - b[3]) / a[3] < 5e-3,
              f"x_NO {a[2]:.4e} / {b[2]:.4e}, x_NO2 {a[3]:.4e} / {b[3]:.4e}")
    x_series = [equilibrium_from_K(T)[2] for T in (500, 1000, 1500, 2000, 2500)]
    check("equilibrium NO rises strongly with T", np.all(np.diff(x_series) > 0))

    print("\n--- Why it is not burning: energy the reaction takes from the gas ---")
    for T in (2000.0, 3000.0):
        x = equilibrium_from_K(T)
        dH_per_mol_gas = x[2] * dH_NO(T) / 2 * 1000  # J per mol of mixture
        cp = 33.0  # J/mol/K, hot air, from the Shomate fits to within 10 percent
        print(f"T = {T:.0f} K: forming x_NO = {x[2]:.4f} absorbs {dH_per_mol_gas:.0f} J/mol of air, "
              f"cooling it by about {dH_per_mol_gas / cp:.0f} K")
    check("forming NO cools the gas (endothermic at flame temperature)", dH_NO(2000.0) > 0,
          f"dH(N2+O2->2NO, 2000 K) = {dH_NO(2000.0):+.1f} kJ")

    print("\n--- Kinetics: Zeldovich rate-limiting step O + N2 -> NO + N ---")
    dH_zeld = dfH["NO"] + dfH["N"] - dfH["O"]
    Ta_min = dH_zeld * 1000 / R
    Ta_used = 38370.0
    print(f"enthalpy of O + N2 -> NO + N = {dH_zeld:+.1f} kJ/mol, so activation temperature >= {Ta_min:.0f} K")
    check("published activation temperature (38370 K) respects the endothermic floor", Ta_used >= Ta_min,
          f"38370 K vs floor {Ta_min:.0f} K")
    print("  T/K   x_NO(eq)   time to reach x_NO(eq) at initial rate")
    times = {}
    for T in (1200.0, 1500.0, 1800.0, 2000.0, 2200.0, 2500.0):
        k = 1.8e8 * np.exp(-Ta_used / T)  # m^3 mol^-1 s^-1
        c = 1e5 / (R * T)  # mol/m^3
        x = equilibrium_from_K(T)
        # partial equilibrium 1/2 O2 <-> O using 298 K values dfH = 249.18, dS = 161.06 - 205.15/2
        dG_O = 249180.0 - T * (161.06 - 205.15 / 2)
        x_O = np.exp(-dG_O / (R * T)) * np.sqrt(x[1])
        rate = 2 * k * x_O * x[0] * c  # d x_NO / dt, s^-1
        tau = x[2] / rate
        times[T] = tau
        print(f"{T:6.0f} {x[2]:10.3e} {tau:12.3e} s")
    check("NO formation is frozen below flame temperatures: time scale falls by >1e5 from 1500 to 2500 K",
          times[1500.0] / times[2500.0] > 1e5, f"{times[1500.0]:.2e} s -> {times[2500.0]:.2e} s")

    print(f"\n{len(FAILS)} checks failed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
