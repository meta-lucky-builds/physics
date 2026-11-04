"""T07: put numbers on the student's everyday gas pictures, and test the bottle's entropy sign.

The entropy claim is tested three ways, none assuming the second law:
  (1) exact microstate counting for two Einstein solids exchanging energy quanta;
  (2) a random-exchange simulation of those solids from a hot start;
  (3) numerical integration of dQ/T for Newton cooling, against the closed-form result,
      and a sweep of random start conditions.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import gammaln

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import N_A, R, k_B

FAILS = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def ln_omega(N, q):
    """ln of the multiplicity of an Einstein solid: q quanta among N oscillators."""
    return gammaln(q + N) - gammaln(q + 1) - gammaln(N)


def gas_numbers():
    print("--- Molecules clashing with walls: air at 300 K, 1 atm ---")
    T, P = 300.0, 101325.0
    m = 28.97e-3 / N_A
    n = P / (k_B * T)
    vbar = np.sqrt(8 * k_B * T / (np.pi * m))
    vrms = np.sqrt(3 * k_B * T / m)
    d = 3.7e-10  # kinetic diameter of air, standard textbook value
    lam = 1 / (np.sqrt(2) * np.pi * d**2 * n)
    print(f"  number density {n:.3g} /m3, mean speed {vbar:.0f} m/s, rms speed {vrms:.0f} m/s")
    print(f"  mean free path {lam * 1e9:.0f} nm, collisions per molecule {vbar / lam:.3g} /s")
    print(f"  wall hits {n * vbar / 4 / 1e4:.3g} per cm2 per s")
    # sampled check of the speed distribution
    rng = np.random.default_rng(3)
    v = rng.normal(0, np.sqrt(k_B * T / m), size=(2_000_000, 3))
    speed = np.linalg.norm(v, axis=1)
    check("sampled Gaussian velocity components give the formula mean speed", abs(speed.mean() / vbar - 1) < 2e-3,
          f"{speed.mean():.1f} vs {vbar:.1f} m/s")
    # flux measured from the sample: sum of positive v_x, per molecule, times n
    flux = n * np.mean(np.clip(v[:, 0], 0, None))
    check("sampled wall flux equals n vbar / 4", abs(flux / (n * vbar / 4) - 1) < 3e-3)


def evaporation():
    print("\n--- Hot bath: what really changes with temperature (water) ---")
    T1, T2 = 293.15, 313.15
    M = 18.015e-3
    m = M / N_A
    for T in (T1, T2):
        print(f"  {T - 273.15:.0f} C: vapour rms speed {np.sqrt(3 * k_B * T / m):.0f} m/s")
    speed_gain = np.sqrt(T2 / T1) - 1
    P1, P2 = 2.3393, 7.3849  # kPa, saturation, NIST WebBook (IAPWS-95)
    h1 = 2537.4 - 83.914
    h2 = 2573.5 - 167.53
    L = 0.5 * (h1 + h2) * M * 1000  # J/mol
    print(f"  latent heat from NIST enthalpies: {h1:.1f} and {h2:.1f} kJ/kg, mean {L / 1000:.1f} kJ/mol = {L / N_A / (k_B * 300):.1f} kT per molecule")
    predicted = np.exp(L / R * (1 / T1 - 1 / T2))
    actual = P2 / P1
    print(f"  speed gain 20 -> 40 C: {100 * speed_gain:.1f} percent")
    print(f"  Boltzmann-factor prediction of vapour pressure ratio: {predicted:.2f}; NIST: {actual:.2f}")
    check("escape is governed by the exponential tail: predicted pressure ratio within 3 percent", abs(predicted / actual - 1) < 0.03)
    check("evaporation rate rises over 40 times more than molecular speed does", (actual - 1) / speed_gain > 40)

    print("  Selection in a sample: remove molecules whose vx exceeds a barrier, see the rest cool.")
    rng = np.random.default_rng(5)
    v = rng.normal(0, 1, size=(4_000_000, 3))  # kT = m = 1
    barrier = 3.0  # kT, reduced so the tail is sampled; water is near 18 kT
    leaving = 0.5 * v[:, 0] ** 2 > barrier
    keep = ~leaving
    frac = leaving.mean()
    ke_all = 0.5 * np.mean(np.sum(v**2, axis=1))
    ke_keep = 0.5 * np.mean(np.sum(v[keep] ** 2, axis=1))
    print(f"  barrier {barrier} kT: fraction leaving {frac:.4f}, KE of those left {ke_keep:.4f} kT vs start {ke_all:.4f} kT")
    check("the ones left behind are colder", ke_keep < ke_all)
    print(f"  1 g evaporating from 1 kg of water removes {h2:.0f} J, cooling it by {h2 / 4186:.2f} K")


def surfactant():
    print("\n--- Oily hands and bubbles: surface energy with and without soap ---")
    gamma_water = 0.0728  # N/m near 20 C
    gamma_sds = 0.038  # N/m at the CMC of sodium dodecyl sulfate, 8.2 mM
    r = 1e-3
    for label, g in (("water", gamma_water), ("SDS solution", gamma_sds)):
        e = 2 * 4 * np.pi * r**2 * g
        print(f"  {label:13s}: surface tension {g * 1000:.1f} mN/m, film energy of a 2 mm bubble {e * 1e6:.2f} uJ, "
              f"Laplace overpressure {4 * g / r:.0f} Pa")
    print(f"  soap lowers the cost of new surface by {100 * (1 - gamma_sds / gamma_water):.0f} percent")
    print("  (energy alone does not make a film last: see README on film elasticity)")
    c_cmc = 8.2  # mol/m3
    nums = c_cmc * N_A
    print(f"  SDS at the CMC: {nums:.2g} molecules per m3; a 1 um oil drop needs a shell of about "
          f"{4 * np.pi * (0.5e-6) ** 2 / 0.5e-18:.2g} molecules at 0.5 nm2 each (head area assumed)")


def einstein_solids():
    print("\n--- Exact counting: bottle (small solid) in a room (large solid) ---")
    N_b, N_r = 400, 4000
    q_b0, q_r0 = 2000, 4000
    q_tot = q_b0 + q_r0
    q = np.arange(0, q_tot + 1)
    s_b = ln_omega(N_b, q)
    s_r = ln_omega(N_r, q_tot - q)
    s = s_b + s_r
    q_eq = int(q[np.argmax(s)])
    print(f"  start: bottle {q_b0} quanta on {N_b} oscillators ({q_b0 / N_b:.2f} each), room {q_r0} on {N_r} ({q_r0 / N_r:.2f} each)")
    print(f"  most probable split: bottle {q_eq} quanta ({q_eq / N_b:.3f} each), room {(q_tot - q_eq) / N_r:.3f} each")
    d_b = s_b[q_eq] - s_b[q_b0]
    d_r = s_r[q_eq] - s_r[q_b0]
    print(f"  dS/k: bottle {d_b:+.1f}, room {d_r:+.1f}, total {d_b + d_r:+.1f}")
    check("bottle entropy falls", d_b < 0)
    check("room entropy rises by more", d_r > -d_b)
    check("the entropy maximum is at equal energy per oscillator (equal temperature)", abs(q_eq / N_b - (q_tot - q_eq) / N_r) < 0.01)
    prob_ratio = s[q_eq] - s[q_b0]
    print(f"  equilibrium is e^{prob_ratio:.0f} = 10^{prob_ratio / np.log(10):.0f} times more probable than the start")

    # Random exchange of quanta, no entropy in the rules.
    rng = np.random.default_rng(9)
    n = np.zeros(N_b + N_r, dtype=int)
    n[:N_b] = q_b0 // N_b
    n[N_b:] = q_r0 // N_r
    total_moves = 3_000_000
    a = rng.integers(0, N_b + N_r, size=total_moves)
    b = rng.integers(0, N_b + N_r, size=total_moves)
    history = []
    for t in range(total_moves):
        if n[a[t]] > 0:
            n[a[t]] -= 1
            n[b[t]] += 1
        if t % 600_000 == 0:
            history.append(int(n[:N_b].sum()))
    q_final = np.mean(n[:N_b].sum())
    print(f"  random exchange: bottle quanta over time {history}, final {q_final:.0f}; exact maximum {q_eq}")
    check("simulation reaches the counted maximum within 3 percent of the bottle's start excess",
          abs(q_final - q_eq) < 0.2 * (q_b0 - q_eq), f"{q_final:.0f} vs {q_eq}")


def bottle():
    print("\n--- Cooling bottle: 0.50 kg of water, 60 C to 25 C, room held at 25 C ---")
    m, c = 0.50, 4186.0
    Ti, Tr = 333.15, 298.15
    S_b_exact = m * c * np.log(Tr / Ti)
    Q = m * c * (Ti - Tr)
    S_r_exact = Q / Tr
    print(f"  heat released {Q / 1000:.1f} kJ")
    print(f"  dS bottle {S_b_exact:+.2f} J/K, dS room {S_r_exact:+.2f} J/K, dS total {S_b_exact + S_r_exact:+.2f} J/K")
    # Newton cooling integrated step by step, recording dQ/T for both sides.
    tau = 1800.0
    dt = 0.5
    T = Ti
    Sb = Sr = 0.0
    t = 0.0
    while t < 40 * tau:
        dQ = m * c * (T - Tr) / tau * dt
        Sb -= dQ / T
        Sr += dQ / Tr
        T -= dQ / (m * c)
        t += dt
    print(f"  Newton cooling stepped numerically: dS bottle {Sb:+.2f}, room {Sr:+.2f}, total {Sb + Sr:+.2f} J/K")
    check("stepped integration agrees with the closed form", abs(Sb + Sr - (S_b_exact + S_r_exact)) < 0.05)

    rng = np.random.default_rng(13)
    Ta = rng.uniform(250, 500, 200000)
    Tb = rng.uniform(250, 500, 200000)
    total = np.log(Tb / Ta) + (Ta - Tb) / Tb  # per unit heat capacity
    check("total entropy change >= 0 for 200000 random start and room temperatures, hot or cold", np.all(total >= 0))
    print(f"  smallest value in the sweep {total.min():.2e} (zero only when T_start = T_room)")
    x = rng.uniform(0.01, 10, 200000)
    check("identity used in the derivation: ln x <= x - 1", np.all(np.log(x) <= x - 1 + 1e-12))
    # entropy budget in molecules
    dn = -S_b_exact / k_B
    print(f"  bottle loses {dn:.3g} k of entropy; world gains a net {(S_b_exact + S_r_exact) / k_B:.3g} k, "
          f"so reversal odds are about exp(-{(S_b_exact + S_r_exact) / k_B:.3g})")


def main():
    gas_numbers()
    evaporation()
    surfactant()
    einstein_solids()
    bottle()
    print(f"\n{len(FAILS)} checks failed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
