"""A10, A14, A15, A18: a photon doing work on moving mirrors, simulated bounce by bounce.

Each reflection is done by Lorentz-transforming into the mirror frame, reversing the
photon there and transforming back. The Doppler formula in the README is not used in
the simulation; it is compared with what comes out. Energy and momentum are then
audited in the lab frame, which the transform does not assume.
"""
import math
import pathlib
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from common.constants import c, e, h, m_e  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


def reflect_clamped(E, direction, beta):
    """Photon (energy E, moving +1 or -1) reflects off a mirror held at velocity beta."""
    g = 1 / math.sqrt(1 - beta * beta)
    E_mirror = g * E * (1 - direction * beta)      # into the mirror frame
    return g * E_mirror * (1 - direction * beta), -direction   # reversed, back to the lab


def reflect_free(E, direction, beta, M):
    """Photon reflects off a free mirror of rest energy M. Returns E', dir', beta'."""
    g = 1 / math.sqrt(1 - beta * beta)
    E_m = g * E * (1 - direction * beta)           # photon energy in mirror rest frame
    E_m_after = E_m / (1 + 2 * E_m / M)            # two-body recoil in that frame
    p_recoil = direction * (E_m + E_m_after)       # mirror momentum in that frame (c=1)
    E_recoil = math.sqrt(M * M + p_recoil * p_recoil)
    beta_rel = p_recoil / E_recoil
    beta_new = (beta_rel + beta) / (1 + beta_rel * beta)
    E_after = g * E_m_after * (1 - direction * beta)
    return E_after, -direction, beta_new


print("1. Single reflection off a mirror receding at fixed speed (the load holds its speed)")
print("    beta     E'/E (sim)   (1-b)/(1+b)    work/E")
for b in (1e-6, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999999):
    E_after, _ = reflect_clamped(1.0, +1, b)
    formula = (1 - b) / (1 + b)
    print(f"   {b:8.6f} {E_after:12.8f} {formula:12.8f} {1 - E_after:10.8f}")
    check(f"Lorentz reflection gives (1-b)/(1+b) at beta={b}", abs(E_after - formula) < 1e-12)
E_up, _ = reflect_clamped(1.0, +1, -0.1)
check("an approaching mirror blueshifts, so work is paid in", E_up > 1, f"E'/E={E_up:.4f}")
check("no finite beta < 1 removes all the energy", reflect_clamped(1.0, 1, 1 - 1e-12)[0] > 0)

print("\n2. Audit of one free-mirror reflection against lab-frame conservation laws")
E0, M = 2.0, 3.0
E1, d1, b1 = reflect_free(E0, +1, 0.2, M)
g0, g1 = 1 / math.sqrt(1 - 0.2**2), 1 / math.sqrt(1 - b1**2)
print(f"   energy before {E0 + g0 * M:.12f}, after {E1 + g1 * M:.12f}")
print(f"   momentum before {E0 + g0 * M * 0.2:.12f}, after {d1 * E1 + g1 * M * b1:.12f}")
check("energy conserved", abs(E0 + g0 * M - E1 - g1 * M) < 1e-12)
check("momentum conserved", abs(E0 + g0 * M * 0.2 - d1 * E1 - g1 * M * b1) < 1e-12)


def lab_solve(E, beta, M):
    """Solve the same collision directly from conservation in the lab, independently."""
    g = 1 / math.sqrt(1 - beta * beta)

    def momentum_gap(b2):
        g2 = 1 / math.sqrt(1 - b2 * b2)
        E_photon = E + g * M - g2 * M           # energy conservation
        return E + g * M * beta + E_photon - g2 * M * b2   # momentum conservation

    b2 = brentq(momentum_gap, beta, 1 - 1e-12, xtol=1e-15)
    return E + g * M - M / math.sqrt(1 - b2 * b2), b2


Ep_direct, b_direct = lab_solve(E0, 0.2, M)
check("direct lab-frame solve agrees with the transform", abs(Ep_direct - E1) < 1e-9,
      f"E'={Ep_direct:.9f}")

print("\n3. One bounce off a free mirror at rest: E' = E/(1 + 2E/Mc^2), and Compton as a check")
E_gamma = m_e * c**2 / e / 1e3     # keV
mass_e = E_gamma
E_back, _, b_e = reflect_free(E_gamma, +1, 0.0, mass_e)
lam = h * c / (E_gamma * 1e3 * e)
lam_c = h / (m_e * c)
E_compton = h * c / (lam + 2 * lam_c) / (1e3 * e)
print(f"   mirror = one free electron: E' = {E_back:.3f} keV; Compton 180 deg formula {E_compton:.3f} keV")
print(f"   the electron leaves with beta = {b_e:.4f}, carrying {E_gamma - E_back:.1f} keV")
check("a free electron used as a mirror reproduces Compton backscatter", abs(E_back - E_compton) < 1e-6)
check("backscatter leaves one third of 511 keV", abs(E_back / E_gamma - 1 / 3) < 1e-9)
for eps in (1e-6, 1e-2, 1, 1e2, 1e6):
    Ep, _, bp = reflect_free(1.0, +1, 0.0, 1 / eps)
    print(f"   E/Mc^2 = {eps:8.0e}: fraction extracted in one bounce = {1 - Ep:.6f}, mirror beta = {bp:.8f}")
    check(f"mirror stays below c at E/Mc^2={eps:.0e}", bp < 1)

print("\n4. Optical piston with a clamped mirror receding at fixed beta: event-driven run")
print("   (left mirror fixed at x=0, right mirror starts at L0=1, c=1)")
print("    beta   bounces  E_final/E0   W_total/E0   sum-1      E L/(first E L)  time")
for beta in (0.001, 0.01, 0.1, 0.5):
    E, x, t, L, direction = 1.0, 0.0, 0.0, 1.0, +1
    W, bounces, product, first = 0.0, 0, None, None
    while E > 0.01 and bounces < 2_000_000:
        if direction > 0:
            dt = (L - x) / (1 - beta)
            t += dt
            L += beta * dt
            x = L
            E_new, direction = reflect_clamped(E, +1, beta)
            W += E - E_new
            E = E_new
            bounces += 1
            product = E * L
            first = first or product
        else:
            t += x
            L += beta * x
            x = 0.0
            direction = +1
    total = W + E
    print(f"   {beta:5.3f} {bounces:8d} {E:11.5f} {W:12.6f} {total - 1:9.1e} {product / first:12.8f} {t:10.2f}")
    check(f"W + E_photon = E0 at beta={beta}", abs(total - 1) < 1e-9)
    check(f"E L is exactly invariant at bounce instants at beta={beta}", abs(product / first - 1) < 1e-9)

print("\n5. Two free mirrors with no load: photon energy flows into mirror kinetic energy")
print("   E0/Mc^2    bounces(to L=100 L0)  E_f/E0     sum check   beta_mirror   sqrt(E0/M)   N*sqrt(eps)")
results = []
for eps in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
    M = 1.0 / eps
    E, xA, xB, bA, bB, x, t, direction = 1.0, 0.0, 1.0, 0.0, 0.0, 0.5, 0.0, +1
    bounces = 0
    while (xB - xA) < 100 and bounces < 5_000_000:
        if direction > 0:
            dt = (xB - x) / (1 - bB)
            t += dt
            xA += bA * dt
            xB += bB * dt
            x = xB
            E, direction, bB = reflect_free(E, +1, bB, M)
        else:
            dt = (x - xA) / (1 + bA)
            t += dt
            xA += bA * dt
            xB += bB * dt
            x = xA
            E, direction, bA = reflect_free(E, -1, bA, M)
        bounces += 1
    # Kinetic energies computed from the small-beta series to avoid cancellation error.
    ke = 0.0
    for b in (bA, bB):
        g = 1 / math.sqrt(1 - b * b)
        ke += M * (g * g * b * b) / (g + 1)
    total = E + ke
    v_expected = math.sqrt(eps)
    product = E * (xB - xA)
    print(f"   {eps:8.0e} {bounces:18d} {E:12.6f} {total - 1:10.1e} {abs(bB):13.4e} {v_expected:11.4e}"
          f" {bounces * math.sqrt(eps):9.3f}")
    results.append((eps, bounces, E, total, abs(bB), v_expected, product))
for eps, bounces, E, total, b, v, product in results:
    check(f"energy ledger closes at eps={eps:.0e}", abs(total - 1) < 1e-8, f"err {abs(total - 1):.1e}")
    check(f"adiabatic invariant E L holds at eps={eps:.0e}", abs(product - 1) < 0.03, f"E L={product:.4f}")
    check(f"mirror speed tends to sqrt(E0/M) at eps={eps:.0e}", abs(b - v) / v < 0.02)
slope = np.polyfit(np.log([r[0] for r in results]), np.log([r[1] for r in results]), 1)[0]
print(f"   fitted N ~ eps^{slope:.3f}, expected eps^-0.5")
check("bounce count scales as sqrt(M/E0)", abs(slope + 0.5) < 0.02)

print("\n6. Extrapolating to a real pair of mirrors, one 511 keV photon (c restored)")
E_ph = m_e * c**2
N_coeff = np.mean([r[1] * math.sqrt(r[0]) for r in results])
print(f"   measured coefficient: N = {N_coeff:.2f} * sqrt(Mc^2/E0) to reach L = 100 L0")
print("   mirror mass    bounces N      needed 1-R (for survival 1/e)   final mirror speed")
for mass in (1e-3, 1e-9, 1e-15):
    eps = E_ph / (mass * c**2)
    N_real = N_coeff / math.sqrt(eps)
    print(f"   {mass:9.0e} kg {N_real:11.2e} {1 / N_real:21.2e} {math.sqrt(eps) * c:21.2e} m/s")
N_1g = N_coeff / math.sqrt(E_ph / (1e-3 * c**2))
check("a 1 g mirror pair needs over 1e13 lossless reflections per photon", N_1g > 1e13, f"N={N_1g:.1e}")

print("\n7. Closed cycle in the piston alone: expand, then recompress to threshold (A18)")
print("   E_final/mc^2 per photon: W_out and the W_in needed to climb back (per photon, units mc^2)")
for E_f_target in (0.99, 0.5, 0.1):
    beta_ex, beta_in = 0.01, -0.01
    E, W_out = 1.0, 0.0
    while E > E_f_target:
        E_new, _ = reflect_clamped(E, +1, beta_ex)
        W_out += E - E_new
        E = E_new
    W_in = 0.0
    while E < 1.0:
        E_new, _ = reflect_clamped(E, +1, beta_in)
        W_in += E_new - E
        E = E_new
    net = W_out - W_in
    print(f"   target {E_f_target:4.2f}: W_out = {W_out:.6f}, W_in = {W_in:.6f}, net = {net:+.6f}")
    check(f"returning to threshold costs at least the work taken (target {E_f_target})", net <= 1e-12)

print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILED: ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
