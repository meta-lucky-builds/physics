"""A03, A07, A12-A18: the full cycle, event by event, with a strict energy ledger.

Part 1 runs a stochastic cycle. A photon pair is emitted, bounces with reflectivity R,
does work on the piston, is blueshifted back to the regeneration energy, and either
converts or is lost. Energy is audited per event and the net work is measured.
Part 2 tests the cascade idea with conserved charges tracked as integers.
Part 3 checks the thermodynamic statement with reaction equilibrium.
"""
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from common.constants import c, e, h, hbar, k_B, m_e, m_p  # noqa: E402
from scipy.constants import physical_constants  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


r_e = physical_constants["classical electron radius"][0]
sigma_T = 8 * math.pi * r_e**2 / 3
MeV = 1e6 * e
M = m_e * c**2 / e / 1e3     # keV
rng = np.random.default_rng(1234)


def sigma_bw(beta):
    return (3 / 16) * sigma_T * (1 - beta**2) * (
        (3 - beta**4) * math.log((1 + beta) / (1 - beta)) - 2 * beta * (2 - beta**2))


def cycle(R, beta_pair, recover, n_cycles, E_f=255.5, recoil_eV=0.0, K_loss_extra=0.0):
    """One-pair-at-a-time cycle, keV. Returns (mean net work per pair consumed, stats).

    Each attempt starts from one stored pair (2M). Expansion extracts work down to E_f per
    photon. Compression returns the photons to E_c per photon, where E_c = gamma M sets the
    pair speed beta_pair in the CM. Photons then cross once per round trip with chance P,
    and each round trip can lose them with probability 1 - R^2. A made pair carries
    kinetic energy, of which the fraction `recover` returns as work. A lost photon pair
    is replaced by external energy 2M, the only way to keep the stock of fuel.
    """
    gamma = 1 / math.sqrt(1 - beta_pair**2)
    E_c = gamma * M
    lam = h * c / (E_c * 1e3 * e)
    P = sigma_bw(beta_pair) / (lam**2 / math.pi)
    ledger_error = 0.0
    W_total, external, made, crossings, lost = 0.0, 0.0, 0, 0, 0
    for _ in range(n_cycles):
        W = 2 * (M - E_f)                       # expansion work
        heat = 0.0
        W -= 2 * (E_c - E_f)                    # compression work paid back
        energy_in_photons = 2 * E_c
        # Bounce-by-bounce conversion is a geometric race between P and loss.
        survive = R * R
        p_convert = P / (P + (1 - P) * (1 - survive)) if P > 0 else 0.0
        k_cross = rng.geometric(P + (1 - P) * (1 - survive)) if (P + (1 - P) * (1 - survive)) > 0 else 10**12
        crossings += k_cross
        if rng.uniform() < p_convert:
            kinetic = energy_in_photons - 2 * M
            kinetic_after_recoil = kinetic - recoil_eV / 1e3 - K_loss_extra
            W += recover * kinetic_after_recoil
            heat += (1 - recover) * kinetic_after_recoil + recoil_eV / 1e3 + K_loss_extra
            stored = 2 * M
            made += 1
        else:
            heat += energy_in_photons           # absorbed in the mirror
            stored = 0.0
            external += 2 * M                   # fresh pair bought from outside
            lost += 1
        # First law for this event. In: stored pair (2M) and any replacement energy.
        energy_in = 2 * M
        energy_out = W + heat + stored
        ledger_error = max(ledger_error, abs(energy_in - energy_out))
        W_total += W
    net_per_cycle = (W_total - external) / n_cycles
    return net_per_cycle, dict(P=P, made=made, lost=lost, mean_crossings=crossings / n_cycles,
                               ledger=ledger_error, E_c=E_c)


print("1. Stochastic cycle, per stored pair. Net = work out - work in - replacement energy (keV)")
print("   R          beta_pair  recover   P per cross  made/total   net [keV]    ideal bound 0?")
rows = [
    (0.29, 0.70, 0.5),
    (0.99, 0.70, 0.5),
    (0.99, 0.70, 1.0),
    (0.999999, 0.70, 1.0),
    (1 - 1e-9, 0.70, 1.0),
    (1 - 1e-9, 0.10, 1.0),
    (1 - 1e-11, 0.01, 1.0),
    (1.0, 0.01, 0.99),
    (1.0, 0.01, 1.0),
]
for R, beta, rec in rows:
    n = 200000 if R < 1 else 20000
    net, st = cycle(R, beta, rec, n)
    print(f"   {R:<10.12g} {beta:8.3f} {rec:7.2f} {st['P']:13.2e} {st['made'] / n:10.2e} {net:12.4f}"
          f"     crossings/cycle {st['mean_crossings']:.1e}")
    check(f"first law holds per event at R={R:.12g}, beta={beta}", st["ledger"] < 1e-9)
    check(f"net work <= 0 at R={R:.12g}, beta={beta}, recover={rec}", net <= 1e-9)
    if R < 1 or rec < 1:
        check(f"net work strictly negative at R={R:.12g}, beta={beta}", net < 0, f"{net:.3f} keV")
    # Compare with the analytic expectation of the same race.
    KE_c = 2 * (st["E_c"] - M)
    p_loss = (1 - st["P"]) * (1 - R * R)
    p_conv = st["P"] / (st["P"] + p_loss)
    expected = -(1 - p_conv) * 2 * st["E_c"] - p_conv * (1 - rec) * KE_c
    spread = 2 * st["E_c"] * math.sqrt(p_conv * (1 - p_conv) / n) + 1e-9
    check(f"Monte Carlo net matches analytic expectation at R={R:.12g}, beta={beta}",
          abs(net - expected) < 5 * spread, f"expected {expected:.4f} keV")

print("\n   Sweep of expansion depth E_f with R=0.99, beta=0.7, recover=1: net must not depend on E_f")
values = []
for E_f in (500.0, 255.5, 100.0, 10.0):
    rng = np.random.default_rng(99)
    net, _ = cycle(0.99, 0.7, 1.0, 100000, E_f=E_f)
    values.append(net)
    print(f"   E_f = {E_f:6.1f} keV: net = {net:.4f} keV")
check("work taken by the piston is exactly repaid: net is independent of expansion depth",
      max(values) - min(values) < 1e-6)

print("\n2. Rate price of the zero-loss limit: conversion P shrinks with the pair speed")
print("   Assumes cavity length 1 m. Time for one conversion even with R = 1.")
print("   beta_pair   KE per pair [keV]   P per crossing   crossings needed   mean cycle time [s]")
L = 1.0
for beta in (0.7, 0.1, 0.01, 1e-3, 1e-4, 1e-6):
    gamma = 1 / math.sqrt(1 - beta**2)
    lam = h * c / (gamma * M * 1e3 * e)
    P = sigma_bw(beta) / (lam**2 / math.pi)
    KE = 2 * (gamma - 1) * M
    print(f"   {beta:9.1e} {KE:14.4e} {P:18.2e} {1 / P:15.2e} {1 / P * 2 * L / c:18.2e}")
check("P tends to 0 as the pair kinetic energy tends to 0",
      sigma_bw(1e-6) / (h * c / (M * 1e3 * e)) ** 2 < sigma_bw(0.1) / (h * c / (M * 1e3 * e)) ** 2)

print("\n3. Nuclear recoil in A16: gamma + nucleus. Budget per pair at threshold (eV)")
for name, A in [("proton", 1.00728), ("Au", 196.967)]:
    mass_keV = A * 931494.0
    recoil = 2 * M * M / mass_keV * 1e3
    print(f"   {name:7s}: recoil excess = {recoil:.3f} eV of {2 * M * 1e3:.0f} eV")
    # Worst case: the nucleus keeps this energy and it ends as heat in its lattice.
    net, st = cycle(1.0, 0.1, 1.0, 20000, recoil_eV=recoil)
    check(f"recoil on {name} makes the cycle strictly negative when not recovered",
          abs(net + recoil / 1e3) < 1e-9, f"net {net * 1e3:.3f} eV")

print("\n4. A07: cascade with charges tracked as integers")
masses_MeV = {"e": 0.51100, "mu": 105.658, "pi": 139.570, "p": 938.272}
print("   X    2 m c^2 [MeV]   least cost on a proton target [MeV]     gain     net [MeV]")
target = 938.272
for X, m in masses_MeV.items():
    cost = 6 * target if X == "p" else 2 * m * (1 + m / target)
    gain = 2 * m
    print(f"   {X:3s} {2 * m:12.3f} {cost:22.3f} {gain:20.3f} {gain - cost:12.3f}")
    check(f"cascade step net is not positive for {X}", gain - cost <= 0)


def run_cascade(E0, m, generations, eps, n_debris=1000):
    """Debris X meets X-bar, and X-bar can only be made in an X X-bar pair. eps = production efficiency."""
    energy = E0
    n_x, n_xbar = n_debris, 0
    log = []
    for g in range(generations):
        pairs = min(int(energy * eps / (2 * m)), n_x)
        energy -= pairs * 2 * m / eps       # production bill
        n_x += pairs                        # pair creation adds one of each
        n_xbar += pairs
        burned = min(n_x, n_xbar)           # every X-bar annihilates with some X
        n_x -= burned
        n_xbar -= burned
        energy += burned * 2 * m            # annihilation returns 2 m c^2 each
        log.append((g, pairs, energy, n_x, n_xbar))
    return log


print("\n   Recycling all output into new antiparticles, electrons as debris, 1e9 MeV to start")
for eps in (1.0, 0.5, 1e-8):
    log = run_cascade(1e9, 0.511, 8, eps)
    print(f"   eps = {eps:<7g} pairs per generation {[row[1] for row in log][:4]}...,"
          f" final energy {log[-1][2]:.6e} MeV, debris X = {log[-1][3]}, leftover X-bar = {log[-1][4]}")
    check(f"energy never grows at eps={eps}", all(row[2] <= 1e9 * (1 + 1e-12) for row in log))
    check(f"net lepton number X - Xbar is unchanged at eps={eps}", log[-1][3] - log[-1][4] == 1000)
    check(f"debris count is still 1000 after 8 generations at eps={eps}", log[-1][3] == 1000)

print("\n   Lifetimes compared with a typical AD cycle of ~100 s")
for name, tau in [("pi0", 8.4e-17), ("pi+-", 2.6033e-8), ("mu", 2.197e-6)]:
    print(f"   {name:5s} tau = {tau:.2e} s; a 100 s production cycle is {100 / tau:.1e} lifetimes")
check("a charged pion is gone 1e9 lifetimes before a new antipion could be made",
      100 / 2.6033e-8 > 1e9)
G_F = 1.1663787e-5              # GeV^-2
s_GeV2 = (1.022e-3) ** 2
sigma_nu = G_F**2 * s_GeV2 / math.pi * 0.3894e-27   # cm^2, order of magnitude
n_e_pb = 11.34 * 6.022e23 * 82 / 207.2
mfp_cm = 1 / (n_e_pb * sigma_nu)
print(f"   neutrino at 1 MeV CM energy: sigma ~ G_F^2 s/pi = {sigma_nu:.1e} cm^2;"
      f" mean free path in lead ~ {mfp_cm / 9.46e17:.0f} light years")
check("a ~MeV neutrino passes through light years of lead", mfp_cm / 9.46e17 > 1)

print("\n5. A03: reaction equilibrium, Kelvin and Carnot")
print("   Reversible means zero affinity: n+ n- = n_eq^2 with n_eq = 2 (2 pi m k T/h^2)^(3/2) exp(-mc^2/kT)")
print("    T [K]        kT/mc^2     log10(n_eq [1/m^3])")
for T in (38e-12, 2.725, 300.0, 1500.0, 1e8, 1e9, 6e9, 1e10):
    kT_mc2 = k_B * T / (m_e * c**2)
    log10n = (math.log10(2) + 1.5 * math.log10(2 * math.pi * m_e * k_B * T / h**2)
              - (m_e * c**2 / (k_B * T)) / math.log(10))
    print(f"   {T:9.2e} {kT_mc2:12.3e} {log10n:18.4e}")
log_1500 = (1.5 * math.log10(2 * math.pi * m_e * k_B * 1500 / h**2)
            - (m_e * c**2 / (k_B * 1500)) / math.log(10))
check("at 1500 K pair equilibrium density is below 10^-1e6 per m^3", log_1500 < -1e6, f"10^{log_1500:.3e}")
n_ph_6e9 = 0.2436 * (k_B * 6e9 / (hbar * c)) ** 3
n_eq_6e9 = 2 * (2 * math.pi * m_e * k_B * 6e9 / h**2) ** 1.5 * math.exp(-m_e * c**2 / (k_B * 6e9))
print(f"   at 6e9 K (non-relativistic form, rough): n_eq/n_gamma ~ {n_eq_6e9 / n_ph_6e9:.2f}: pairs appear freely")
check("thermal pair creation only matters near 6e9 K", n_eq_6e9 / n_ph_6e9 > 0.05)

print("\n   Entropy of a photon in a cavity, compared with its energy (A21 D)")
V = 1.0
dE_E = 2e-3
lam_511 = h * c / (M * 1e3 * e)
modes = 8 * math.pi * (V / lam_511**3) * dE_E
ln_modes = math.log(modes)
T_amb = 300.0
cost_eV = ln_modes * k_B * T_amb / e
print(f"   modes in 1 m^3 within dE/E=0.2%: {modes:.2e}; S = k ln(modes) = {ln_modes:.0f} k = {ln_modes / math.log(2):.0f} bits")
print(f"   T S at 300 K = {cost_eV:.2f} eV per photon = {cost_eV / (M * 1e3):.1e} of its energy")
landauer_eV = k_B * T_amb * math.log(2) / e
print(f"   Landauer cost at 300 K = {landauer_eV:.4f} eV per bit")
check("the entropy price of knowing a 511 keV photon's mode is below 1e-5 of its energy",
      cost_eV / (M * 1e3) < 1e-5)

print("\n   Mixing entropy of one e+ and one e-: k ln 2 each at most")
print(f"   2 k T ln 2 at 300 K = {2 * landauer_eV:.3f} eV of 1.022 MeV = {2 * landauer_eV / (2 * M * 1e3):.1e}")

print("\n   Kelvin test with a real heat bath: W <= Q_h (1 - T_c/T_h), with no pair term at all")
for T_h, T_c in [(1500.0, 300.0), (1500.0, 2.725), (1500.0, 38e-12)]:
    print(f"   T_h={T_h:.0f} K, T_c={T_c:.3g} K: eta_C = 1 - {T_c / T_h:.3e}")
print("   Over a closed cycle the pair stock returns to start, so the annihilation term is zero.")
print("   Any positive W must be ordinary heat-engine output paid for by Q_h from the bath.")

print("\n6. Three-photon decay of ortho-positronium cannot be undone by two-photon collision")
n = 300000
# Random massless 3-body decays with total invariant mass 2m: use energy simplex.
x = rng.dirichlet([1, 1, 1], n) * 2       # x1+x2+x3 = 2 in units of m, requires each x <= 1
mask = (x <= 1).all(axis=1)
x = x[mask]
s12 = 4 * (1 - x[:, 2])                    # s12 = M^2 - 2 M E3 (units m^2) with M=2m
print(f"   {len(x)} sampled decays: largest s of any two photons / (2m)^2 = {s12.max() / 4:.4f}")
check("for three photons, any two are below the pair threshold", s12.max() < 4)
print("   Ratio of lifetimes: 142 ns / 125 ps = "
      f"{142e-9 / 125e-12:.0f}, so with 3:1 spin weight a spin-averaged free 2-gamma:3-gamma ratio"
      f" near {142e-9 / 125e-12 / 3:.0f}:1 (standard value about 372:1).")

print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILED: ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
