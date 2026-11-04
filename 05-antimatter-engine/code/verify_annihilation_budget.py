"""A19, A20: where the energy of proton-antiproton annihilation at rest actually goes.

A Monte Carlo generates pions by n-body phase space with the Klempt et al. mean
multiplicities, then decays every pion and muon with explicit four-vectors. It is a model:
real annihilation goes partly through resonances and 5% of events contain kaons. Running
two different multiplicity assumptions shows how far the answer depends on the model.
"""
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from common.constants import c, e, m_e, m_p  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


MeV = 1e6 * e
M_P = m_p * c**2 / MeV
M_E = m_e * c**2 / MeV
M_PI = 139.570        # MeV, PDG
M_PI0 = 134.977       # MeV, PDG
M_MU = 105.658        # MeV, PDG
W = 2 * M_P
rng = np.random.default_rng(19)


def two_body_p(M, m1, m2):
    return math.sqrt(max(0.0, (M * M - (m1 + m2) ** 2) * (M * M - (m1 - m2) ** 2))) / (2 * M)


def isotropic():
    cos = rng.uniform(-1, 1)
    phi = rng.uniform(0, 2 * math.pi)
    sin = math.sqrt(1 - cos * cos)
    return np.array([sin * math.cos(phi), sin * math.sin(phi), cos])


def boost(p4, beta):
    """Boost four-vector p4=(E,px,py,pz) by velocity beta (3-vector)."""
    b2 = beta @ beta
    if b2 == 0:
        return p4.copy()
    g = 1 / math.sqrt(1 - b2)
    bp = beta @ p4[1:]
    E = g * (p4[0] + bp)
    p = p4[1:] + ((g - 1) * bp / b2 + g * p4[0]) * beta
    return np.concatenate(([E], p))


def genbod(masses):
    """Raubold-Lynch n-body phase space. Returns (four-vectors, weight)."""
    n = len(masses)
    if n == 2:
        p = two_body_p(W, masses[0], masses[1])
        u = isotropic()
        return [np.concatenate(([math.hypot(p, masses[0])], p * u)),
                np.concatenate(([math.hypot(p, masses[1])], -p * u))], p
    T = W - sum(masses)
    r = np.sort(rng.uniform(size=n - 2))
    M_inv = [masses[0]]
    for i in range(1, n - 1):
        M_inv.append(sum(masses[: i + 1]) + r[i - 1] * T)
    M_inv.append(W)
    weight = 1.0
    momenta = [np.concatenate(([masses[0]], np.zeros(3)))]
    for i in range(1, n):
        p = two_body_p(M_inv[i], M_inv[i - 1], masses[i])
        weight *= p
        u = isotropic()
        new = np.concatenate(([math.hypot(p, masses[i])], p * u))
        beta = -p * u / math.hypot(p, M_inv[i - 1])
        momenta = [boost(q, beta) for q in momenta]
        momenta.append(new)
    return momenta, weight


def pick_final_state(mode):
    """Number of (charged, neutral) pions under two assumptions."""
    if mode == "A":
        # Charged prong weights chosen to give <n+->=3.05; neutral: Poisson mean tuned below.
        n_ch = rng.choice([0, 2, 4, 6], p=[0.04, 0.44, 0.46, 0.06])
        n_0 = rng.poisson(1.93)
        while n_ch + n_0 < 2:
            n_0 = rng.poisson(1.93)
    else:
        # Fixed mixture with the same mean charged and neutral counts
        n_ch, n_0 = (2, 3) if rng.uniform() < 0.5 else (4, 1)
    return n_ch, n_0


def decay_chain(p4_pi, charged):
    """Returns energies by final class: gamma, neutrino, electron/positron."""
    gamma_E = nu_E = e_E = 0.0
    beta_pi = p4_pi[1:] / p4_pi[0]
    if not charged:
        u = isotropic()
        for sign in (1, -1):
            p = boost(np.concatenate(([M_PI0 / 2], sign * M_PI0 / 2 * u)), beta_pi)
            gamma_E += p[0]
        return gamma_E, nu_E, e_E
    p_mu = two_body_p(M_PI, M_MU, 0.0)
    u = isotropic()
    nu = boost(np.concatenate(([p_mu], -p_mu * u)), beta_pi)
    mu = boost(np.concatenate(([math.hypot(p_mu, M_MU)], p_mu * u)), beta_pi)
    nu_E += nu[0]
    # Muon decay: Michel spectrum x^2(3-2x) for the electron energy, isotropic, polarisation ignored.
    while True:
        x = rng.uniform()
        if rng.uniform() * 1.0 <= x * x * (3 - 2 * x):
            break
    E_e = max(x * M_MU / 2, M_E)
    p_e = math.sqrt(E_e * E_e - M_E * M_E)
    u_e = isotropic()
    ele = boost(np.concatenate(([E_e], p_e * u_e)), mu[1:] / mu[0])
    e_E += ele[0]
    nu_E += mu[0] - ele[0]        # the two muon-decay neutrinos carry the rest
    return gamma_E, nu_E, e_E


def run(mode, n_events=8000):
    wsum = 0.0
    tot = dict(gamma=0.0, nu=0.0, ele=0.0, E_ch=0.0, E_0=0.0, n_ch=0.0, n_0=0.0, ph=0.0)
    worst = 0.0
    for _ in range(n_events):
        n_ch, n_0 = pick_final_state(mode)
        masses = [M_PI] * n_ch + [M_PI0] * n_0
        if sum(masses) >= W:
            continue
        momenta, w = genbod(masses)
        total = sum(momenta)
        worst = max(worst, abs(total[0] - W), np.abs(total[1:]).max())
        # The phase-space weight applies inside a multiplicity class. Dividing by the class's
        # mean weight keeps each class at the probability assumed in pick_final_state.
        w_class = w / class_scale(n_ch, n_0, masses)
        g_e, n_e, e_e = 0.0, 0.0, 0.0
        for i, p in enumerate(momenta):
            g_, nu_, el_ = decay_chain(p, i < n_ch)
            g_e += g_
            n_e += nu_
            e_e += el_
        E_ch = sum(p[0] for p in momenta[:n_ch])
        E_0 = sum(p[0] for p in momenta[n_ch:])
        wsum += w_class
        tot["gamma"] += w_class * g_e
        tot["nu"] += w_class * n_e
        tot["ele"] += w_class * e_e
        tot["E_ch"] += w_class * E_ch
        tot["E_0"] += w_class * E_0
        tot["n_ch"] += w_class * n_ch
        tot["n_0"] += w_class * n_0
        tot["ph"] += w_class * 2 * n_0
    return {k: v / wsum for k, v in tot.items()}, worst


_scale_cache = {}


def class_scale(n_ch, n_0, masses):
    """Mean phase-space weight of a class, so classes keep their chosen probabilities."""
    key = (n_ch, n_0)
    if key not in _scale_cache:
        _scale_cache[key] = max(np.mean([genbod(masses)[1] for _ in range(100)]), 1e-300)
    return _scale_cache[key]


print("1. Monte Carlo of p pbar at rest, 1876.5 MeV")
print("   mode A: charged 0/2/4/6 prongs (4/44/46/6 %), Poisson neutrals; mode B: 50% (2,3), 50% (4,1).")
print("   mode  <n+->  <n0>   charged pions  neutral pions   gamma    nu     e+-    sum    photons/ann  <E_gamma>")
results = {}
for mode in "AB":
    res, worst = run(mode)
    results[mode] = res
    fracs = {k: res[k] / W for k in ("gamma", "nu", "ele", "E_ch", "E_0")}
    sum_f = fracs["gamma"] + fracs["nu"] + fracs["ele"]
    mean_gamma = res["gamma"] / res["ph"]
    print(f"   {mode:4s} {res['n_ch']:6.2f} {res['n_0']:5.2f} {fracs['E_ch']:13.3f} {fracs['E_0']:14.3f}"
          f" {fracs['gamma']:8.3f} {fracs['nu']:6.3f} {fracs['ele']:6.3f} {sum_f:7.4f}"
          f" {res['ph']:9.2f} {mean_gamma:10.1f} MeV")
    check(f"mode {mode}: four-momentum conserved per event", worst < 1e-6, f"worst {worst:.1e} MeV")
    check(f"mode {mode}: final classes sum to the initial energy (per event)", abs(sum_f - 1) < 1e-3)
    check(f"mode {mode}: <n+-> within 10% of 3.05", abs(res["n_ch"] - 3.05) / 3.05 < 0.10, f"{res['n_ch']:.2f}")
    check(f"mode {mode}: <n0> within 20% of 1.93", abs(res["n_0"] - 1.93) / 1.93 < 0.20, f"{res['n_0']:.2f}")

gamma_f = np.mean([results[m]["gamma"] / W for m in "AB"])
nu_f = np.mean([results[m]["nu"] / W for m in "AB"])
ele_f = np.mean([results[m]["ele"] / W for m in "AB"])
print(f"\n   Rounded result: photons {gamma_f:.0%}, neutrinos {nu_f:.0%}, electrons and positrons {ele_f:.0%}.")
check("neutrinos take roughly 40-50% of the energy", 0.38 < nu_f < 0.52)
check("photons take well under the 1-2 GeV stated in A19", gamma_f * W < 1000, f"{gamma_f * W:.0f} MeV")

print("\n   Single decay check: one charged pion at rest (exact two-body plus mean Michel spectrum)")
p_mu = two_body_p(M_PI, M_MU, 0)
E_nu1 = p_mu
E_mu = math.hypot(p_mu, M_MU)
mean_x = 0.7
E_e_mean = mean_x * M_MU / 2 * E_mu / M_MU
nu_share = (M_PI - E_e_mean) / M_PI
print(f"   E(nu_mu) = {E_nu1:.2f} MeV, <E(e)> = {E_e_mean:.1f} MeV: neutrino share of a charged pion = {nu_share:.1%}")
check("neutrino share of a charged pion is about 72%", abs(nu_share - 0.72) < 0.01)

print("\n2. A19 ledger per annihilation, using the Monte Carlo means (MeV)")
E_gamma, E_nu, E_ele = gamma_f * W, nu_f * W, ele_f * W
n_pairs = np.mean([results[m]["n_ch"] for m in "AB"]) / 2
regen = 2 * M_E
W_max = E_gamma + E_ele     # generous: every photon and e+- joule extracted as work
print(f"   fuel rest energy consumed                 {W:9.1f}")
print(f"   lost to neutrinos                         {E_nu:9.1f}")
print(f"   in photons                                {E_gamma:9.1f}")
print(f"   in e+- (rest + kinetic)                   {E_ele:9.1f}")
print(f"   upper bound on work                       {W_max:9.1f}  ({W_max / W:.1%} of fuel)")
print(f"   e+e- pairs the pion decays leave behind   {n_pairs:9.2f}  per annihilation")
print(f"   cost to remake one e+e- pair              {regen:9.3f}  ({regen / W:.3%} of fuel)")
net_closed = W_max - W
print(f"\n   Close the p pbar loop too, at the bare threshold: net = {W_max:.1f} - {W:.1f} = {net_closed:.1f} MeV")
check("closing the proton loop is strictly negative by at least the neutrino energy",
      net_closed < 0 and abs(net_closed + E_nu) < 1)
print("   Do not close it: work comes from burning one antiproton per cycle. That is fuel.")

print("\n3. A20: the claimed formula against a real energy accounting")
print("   E_gamma [MeV]   claimed 1 - m_e c^2/E    meaning")
for E_val, label in [(0.511, "e+e- at rest"), (1.0, "1 MeV photon"), (E_gamma / max(results["A"]["ph"], 1), "mean pi0 photon")]:
    print(f"   {E_val:10.3f} {1 - M_E / E_val:18.5f}    {label}")
check("the claimed formula gives 0 for e+e- at rest, where 100% is not excluded by conservation",
      abs(1 - M_E / 0.511) < 2e-3)
print("   Derived open-loop accounting: eta = f_usable * (1 - E_final/E_initial) * capture")
print("     e+e-    : f_usable = 1, so eta -> 1 as E_final -> 0 (no threshold matters if fuel is not remade)")
print(f"     p pbar  : f_usable <= {W_max / W:.3f}, so eta <= {W_max / W:.3f} for any E_gamma")
check("open-loop p pbar ceiling is below the 1 the formula implies", W_max / W < 0.7)

print("\n4. Real production cost (sourced inputs, order of magnitude)")
p_per_pulse = 1.5e13            # Torregrosa Martin et al. 2021
p_momentum = 26e3               # MeV/c
KE_p = math.hypot(p_momentum, M_P) - M_P
pbar_per_cycle = (3e7, 4e7)     # Bartmann/ELENA paper figures, per AD cycle
wall_plug = 0.05                # LaPointe & Reddy 2001 assumption
print(f"   beam energy per pulse: {p_per_pulse:.1e} x {KE_p / 1000:.2f} GeV = {p_per_pulse * KE_p * MeV:.1e} J")
for N_pbar in pbar_per_cycle:
    beam_eff = N_pbar * 2 * M_P / (p_per_pulse * KE_p)
    print(f"   {N_pbar:.0e} pbar/cycle: annihilation energy / beam energy = {beam_eff:.1e};"
          f" with 5% wall plug = {beam_eff * wall_plug:.1e}")
eff_nominal = 3.5e7 * 2 * M_P / (p_per_pulse * KE_p) * wall_plug
eff_capture = eff_nominal * 0.01
print(f"   LaPointe & Reddy figure: 1 pbar per 2.5e6 protons, beam efficiency"
      f" {2 * M_P / (2.5e6 * 26000):.1e} (their 3e-8)")
print(f"   wall-plug to delivered beam ~ {eff_nominal:.0e}; foil capture below 1% makes it below {eff_capture:.0e}")
check("delivered-beam wall-plug efficiency is in the 1e-9 to 1e-7 range", 1e-9 < eff_nominal < 1e-7,
      f"{eff_nominal:.1e}")

print("\n5. Round-trip energy storage: efficiency of put-in energy back as work")
for eps in (1e-9, 1e-8):
    print(f"   production {eps:.0e} x usable {W_max / W:.2f} (ideal conversion) = {eps * W_max / W:.1e}")
print("   Fuel energy density: 2 x 9e16 J per kg of antimatter = 1.8e17 J/kg")
E_gram = 1.8e14
print(f"   making 1 g of antiprotons: rest energy of the pair mass {E_gram:.1e} J; wall-plug at 1e-8:"
      f" {E_gram / 1e-8:.1e} J")
rate = 3.5e7 / 100          # antiprotons/s through the AD
N_gram = 1e-3 / m_p
print(f"   AD rate {rate:.1e} pbar/s means 1 g takes {N_gram / rate / (3.156e7):.1e} years")
check("one gram at the AD rate takes longer than the age of the universe",
      N_gram / rate / 3.156e7 > 1.38e10)

print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILED: ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
