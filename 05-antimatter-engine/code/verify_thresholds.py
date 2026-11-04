"""A01, A02, A03, A16: production thresholds, tested by searching final states directly.

The invariant-mass derivation in the README predicts each threshold. This script does
not use that formula for the main tests. It treats "what is the least beam energy that
conserves energy and momentum?" as a constrained minimisation over final-state momenta.
"""
import pathlib
import sys

import numpy as np
from scipy.optimize import brentq, minimize

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from common.constants import c, e, m_e, m_p  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


MeV = 1e6 * e
mec2 = m_e * c**2 / MeV
mpc2 = m_p * c**2 / MeV
rng = np.random.default_rng(5)

print("1. Fixed-target threshold from direct minimisation (units of the particle mass)")


def min_beam_energy(n_final, seed):
    """Least E_beam with beam(m=1) + target(m=1, at rest) -> n_final particles of m=1."""
    local = np.random.default_rng(seed)

    def beam_energy(p):
        total = p.reshape(n_final, 3).sum(axis=0)
        return np.sqrt(total @ total + 1.0)

    def energy_balance(p):
        energies = np.sqrt((p.reshape(n_final, 3) ** 2).sum(axis=1) + 1.0)
        return beam_energy(p) + 1.0 - energies.sum()

    start = local.normal(0, 3, 3 * n_final)
    res = minimize(beam_energy, start, constraints={"type": "eq", "fun": energy_balance},
                   method="SLSQP", options={"maxiter": 2000, "ftol": 1e-13})
    return res.fun, res.x.reshape(n_final, 3), res.success


best, momenta = np.inf, None
for s in range(5):
    value, p, ok = min_beam_energy(4, s)
    if ok and value < best:
        best, momenta = value, p
kinetic = best - 1.0
energies = np.sqrt((momenta**2).sum(axis=1) + 1.0)
speeds = np.linalg.norm(momenta, axis=1) / energies
print(f"   beam energy {best:.6f} m c^2, beam kinetic {kinetic:.6f} m c^2,"
      f" total lab energy {best + 1:.6f} m c^2")
print(f"   final-state speeds at the minimum: {np.round(speeds, 5)}")
check("p p -> p p p pbar threshold kinetic energy is 6 m_p c^2", abs(kinetic - 6) < 1e-4,
      f"{kinetic * mpc2 / 1000:.3f} GeV")
check("total lab energy is 8 m c^2 with beam 7 m c^2 + target 1 m c^2", abs(best - 7) < 1e-4)
check("all four final particles share beta = sqrt(3)/2 at threshold",
      np.all(np.abs(speeds - np.sqrt(3) / 2) < 1e-3), f"beta={speeds.mean():.5f}")

print("\n2. Photon on a nucleus: solve the conservation laws by root finding (no s formula)")
for name, M_MeV in [("proton", mpc2), ("Au-197 nucleus", 196.967 * 931.494 - 79 * mec2)]:
    # At threshold everything moves together with speed b. Require consistency between
    # k = momentum carried off and k + M = energy carried off.
    total_mass = M_MeV + 2 * mec2

    def mismatch(b):
        g = 1 / np.sqrt(1 - b * b)
        return (g * total_mass - M_MeV) - g * b * total_mass

    b = brentq(mismatch, 1e-12, 0.9, xtol=1e-16)
    k = 1 / np.sqrt(1 - b * b) * b * total_mass
    excess_eV = (k - 2 * mec2) * 1e6
    predicted_eV = 2 * mec2**2 / M_MeV * 1e6
    print(f"   {name:15s}: k = {k:.9f} MeV, excess above 2 m_e c^2 = {excess_eV:.3f} eV"
          f" (2m^2/M = {predicted_eV:.3f} eV)")
    check(f"recoil excess on {name} matches 2 m^2 c^2 / M",
          abs(excess_eV - predicted_eV) / predicted_eV < 2e-3)

print("\n3. A single free photon cannot become a free pair (A03)")
n = 1_000_000
p1 = rng.normal(0, 1, (n, 3)) * 10 ** rng.uniform(-3, 4, (n, 1))
p2 = rng.normal(0, 1, (n, 3)) * 10 ** rng.uniform(-3, 4, (n, 1))
E_sum = np.sqrt((p1**2).sum(1) + 1) + np.sqrt((p2**2).sum(1) + 1)
P_sum = np.linalg.norm(p1 + p2, axis=1)
gap = (E_sum - P_sum) * (E_sum + P_sum)  # E^2 - P^2, numerically stable
print(f"   {n} random pairs: min (E^2-P^2)/(mc^2)^2 = {gap.min():.6f} (a photon needs 0)")
check("every pair has invariant mass >= 2m, a photon has 0", gap.min() >= 4 - 1e-9)
print("   min (E - P c)/(m c^2) over samples:"
      f" {(E_sum - P_sum).min():.3e}  (always > 0, so E = Pc never holds)")
check("E - Pc stays strictly positive", (E_sum - P_sum).min() > 0)

print("\n4. Two photons: threshold on s means total energy >= 2 m c^2 in every lab frame (A16, A18)")
E1 = 10 ** rng.uniform(-3, 3, 2_000_000)
E2 = 10 ** rng.uniform(-3, 3, 2_000_000)
cos = rng.uniform(-1, 1, 2_000_000)
s_over_4m2 = E1 * E2 * (1 - cos) / 2
over = s_over_4m2 >= 1
print(f"   {over.sum()} of {len(E1)} random photon pairs are above threshold;"
      f" least E1+E2 among them = {(E1 + E2)[over].min():.5f} m c^2")
check("no above-threshold photon pair has E1 + E2 < 2 m c^2", (E1 + E2)[over].min() >= 2)
check("sub-threshold redshifted head-on pair fails",
      (0.9 * 0.9 * (1 - (-1)) / 2) < 1, "E1=E2=0.9 m c^2")
check("asymmetric head-on pair passes when E1 E2 >= m^2 c^4",
      (0.5 * 2.0 * 2 / 2) >= 1, "E1=0.5, E2=2.0 m c^2, total 2.5 m c^2")

print("\n5. Kinetic energy adds to the bill (A02): pair made moving with Lorentz factor g")
print("     beta      gamma    energy/2mc^2   forward photon/mc^2   backward photon/mc^2")
for beta in (0, 0.5, 0.866, 0.99, 0.999):
    g = 1 / np.sqrt(1 - beta**2)
    forward, backward = g * (1 + beta), g * (1 - beta)
    print(f"   {beta:7.3f} {g:9.3f} {g:12.3f} {forward:18.3f} {backward:21.4f}")
    check(f"in-flight annihilation returns exactly 2 g m c^2 at beta={beta}",
          abs(forward + backward - 2 * g) < 1e-12)

print("\n6. How hard a fast antiproton is to hold (A02)")
print("   p [MeV/c]    KE [MeV]     beta    B rho [T m]   KE / 10 kV well")
for p_MeV in (3500.0, 100.0, 13.7):
    E_total = np.hypot(p_MeV, mpc2)
    KE = E_total - mpc2
    beta = p_MeV / E_total
    Brho = p_MeV * 1e6 / c
    print(f"   {p_MeV:9.1f} {KE:11.3f} {beta:9.4f} {Brho:12.3f} {KE * 1e6 / 1e4:14.3e}")
KE_100 = np.hypot(100.0, mpc2) - mpc2
check("100 MeV/c antiproton has the AD's quoted 5.3 MeV", abs(KE_100 - 5.3) < 0.05,
      f"{KE_100:.3f} MeV")

print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILED: ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
