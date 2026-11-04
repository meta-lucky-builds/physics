"""A03, A16, A18, A21(C): the gamma-gamma cross section and how often a pair is remade.

The closed-form Breit-Wheeler and Dirac cross sections are compared with a numerical
angular integral of the tree-level QED squared amplitude (Peskin & Schroeder eq. 5.105,
crossed), so the three expressions test one another. Detailed balance is then checked.
"""
import math
import pathlib
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from common.constants import c, e, h, hbar, m_e  # noqa: E402
from scipy.constants import alpha, physical_constants  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


r_e = physical_constants["classical electron radius"][0]
sigma_T = 8 * math.pi * r_e**2 / 3
lambda_bar = hbar / (m_e * c)
barn = 1e-28


def amp2(beta, cos_t):
    """(1/4) sum |M|^2 / e^4 for gamma gamma <-> e+ e- in the CM frame, in units m=1."""
    g2 = 1 / (1 - beta**2)
    E2 = g2
    pk1 = E2 * (1 - beta * cos_t)
    pk2 = E2 * (1 + beta * cos_t)
    inv = 1 / pk1 + 1 / pk2
    return 2 * (pk2 / pk1 + pk1 / pk2 + 2 * inv - inv**2)


def sigma_from_amplitude(beta, channel):
    """Numerical total cross section in m^2. channel: 'pair' or 'annihilation'."""
    g2 = 1 / (1 - beta**2)
    s = 4 * g2                 # m=1
    e4 = (4 * math.pi * alpha) ** 2
    integral = quad(lambda x: amp2(beta, x) * 2 * math.pi, -1, 1, limit=200,
                    points=[0.0] if beta < 0.999 else None)[0]
    if channel == "pair":
        prefactor = e4 * beta / (64 * math.pi**2 * s)          # final speed beta over initial 1
    else:
        prefactor = 0.5 * e4 / beta / (64 * math.pi**2 * s)    # identical final photons: 1/2
    return prefactor * integral * lambda_bar**2


def sigma_bw(beta):
    """Breit-Wheeler total cross section, closed form (Gould & Schreder form)."""
    return (3 / 16) * sigma_T * (1 - beta**2) * (
        (3 - beta**4) * math.log((1 + beta) / (1 - beta)) - 2 * beta * (2 - beta**2))


def sigma_dirac(beta):
    """Dirac e+ e- -> 2 gamma, written with gamma of the positron in the electron rest frame."""
    g_cm2 = 1 / (1 - beta**2)
    g = 2 * g_cm2 - 1
    root = math.sqrt(g * g - 1)
    return math.pi * r_e**2 / (g + 1) * (
        (g * g + 4 * g + 1) / (g * g - 1) * math.log(g + root) - (g + 3) / root)


print("1. Three independent routes to the cross sections")
print("    beta     sigma_BW [b]   amplitude [b]   sigma_ann [b]   amplitude [b]   ratio/2beta^2")
for b in (0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99):
    bw, num = sigma_bw(b), sigma_from_amplitude(b, "pair")
    dirac, num_a = sigma_dirac(b), sigma_from_amplitude(b, "annihilation")
    ratio = bw / dirac / (2 * b * b)
    print(f"   {b:5.2f} {bw / barn:13.5f} {num / barn:14.5f} {dirac / barn:14.4f} {num_a / barn:14.4f} {ratio:12.8f}")
    check(f"closed BW formula equals amplitude integral at beta={b}", abs(bw - num) / num < 1e-6)
    check(f"Dirac formula equals amplitude integral at beta={b}", abs(dirac - num_a) / num_a < 1e-6)
    check(f"detailed balance sigma_pair/sigma_ann = 2 beta^2 at beta={b}", abs(ratio - 1) < 1e-6)

print("\n2. Threshold behaviour")
for b in (1e-1, 1e-2, 1e-3, 1e-4):
    print(f"   beta={b:.0e}: sigma_BW/(pi r_e^2 beta) = {sigma_bw(b) / (math.pi * r_e**2 * b):.6f}")
check("sigma_BW -> pi r_e^2 beta, which is zero at threshold",
      abs(sigma_bw(1e-3) / (math.pi * r_e**2 * 1e-3) - 1) < 1e-4)
check("sigma_BW vanishes at beta=0", sigma_bw(1e-12) < 1e-40)
v_sigma = 2 * 1e-3 * c * sigma_dirac(1e-3)
print(f"   slow annihilation rate coefficient sigma v = {v_sigma:.3e} m^3/s,"
      f" pi r_e^2 c = {math.pi * r_e**2 * c:.3e} m^3/s (no Coulomb focusing)")
check("sigma v -> pi r_e^2 c at low speed", abs(v_sigma / (math.pi * r_e**2 * c) - 1) < 1e-2)

peak = minimize_scalar(lambda b: -sigma_bw(b), bounds=(0.1, 0.99), method="bounded", options={"xatol": 1e-6})
b_peak, s_peak = peak.x, sigma_bw(peak.x)
g_peak = 1 / math.sqrt(1 - b_peak**2)
print(f"   peak: beta = {b_peak:.4f}, sigma = {s_peak:.3e} m^2 = {s_peak / sigma_T:.4f} sigma_T,"
      f" head-on photons of {g_peak * 511:.0f} keV each")
check("peak is about 0.26 sigma_T near beta 0.7", abs(s_peak / sigma_T - 0.26) < 0.01 and abs(b_peak - 0.7) < 0.02)
print("   E_gamma head-on [keV]    beta      sigma [b]")
for E_keV in (511.0, 512, 520, 600, 722, 1000, 3000, 10000):
    gamma = E_keV / (m_e * c**2 / e / 1e3)
    if gamma <= 1:
        print(f"   {E_keV:12.0f} {0.0:13.4f} {0.0:12.4f}")
        continue
    b = math.sqrt(1 - 1 / gamma**2)
    print(f"   {E_keV:12.0f} {b:13.4f} {sigma_bw(b) / barn:12.4f}")

print("\n3. Alignment: must the photons be perfectly head-on? (A21 C)")
E_each = 1.5    # units of m c^2
print(f"   each photon {E_each} m c^2 = {E_each * 511:.0f} keV")
print("   crossing angle [deg]    beta_cm    sigma/sigma(head-on)")
b_head = math.sqrt(1 - 1 / (E_each**2))
for angle in (180, 150, 120, 90, 80):
    s_ratio = E_each**2 * (1 - math.cos(math.radians(angle))) / 2
    ratio = sigma_bw(math.sqrt(1 - 1 / s_ratio)) / sigma_bw(b_head) if s_ratio > 1 else 0.0
    print(f"   {angle:12d} {math.sqrt(max(0, 1 - 1 / s_ratio)):16.4f} {ratio:16.3f}")
s_90 = E_each**2 * (1 - math.cos(math.radians(90))) / 2
check("90 degree crossing still works for 1.5 m c^2 photons", s_90 > 1)

print("\n4. Chance that a single photon pair converts per crossing, best possible focus")
print("   (assumes Gaussian waist w0 = lambda/pi, area A = lambda^2/pi; beam geometry is an estimate)")
print("    E_gamma[keV]   beta      lambda[pm]   P = sigma/A")
P_best = 0
for E_keV in (512, 520, 600, 722, 1000, 3000):
    gamma = E_keV / (m_e * c**2 / e / 1e3)
    b = math.sqrt(1 - 1 / gamma**2)
    lam = h * c / (E_keV * 1e3 * e)
    P = sigma_bw(b) / (lam**2 / math.pi)
    P_best = max(P, P_best)
    print(f"   {E_keV:10.0f} {b:8.4f} {lam * 1e12:12.3f} {P:12.2e}")
print(f"   largest P in this list {P_best:.2e}; alpha^2 = {alpha**2:.2e}")
check("single-pair conversion is below 1e-3 per crossing even at best focus", P_best < 1e-3)

print("\n5. Conversion before the photons are lost: eps = P/(P + 1 - R^2) for a cavity")
print("   (P at 722 keV, the cross-section peak; R is reflectivity per bounce, assumed)")
P_722 = sigma_bw(math.sqrt(1 - (511 / 722) ** 2)) / ((h * c / (722e3 * e)) ** 2 / math.pi)
print("         R         eps       photons lost per pair made")
for R in (0.29, 0.99, 0.999999, 1 - 1e-9):
    eps = P_722 / (P_722 + 1 - R**2)
    print(f"   {R:12.9f} {eps:10.2e} {(1 - eps) / eps:15.3e}")
check("even R=0.999999 converts only a part of the photons",
      P_722 / (P_722 + 1 - 0.999999**2) < 1)

print("\n6. Bulk photon beams: optical depth tau = N sigma / A")
A_spot = math.pi * (1e-6) ** 2
N_tau1 = A_spot / s_peak
E_tau1 = N_tau1 * 722e3 * e
print(f"   1 um radius spot, sigma at peak: tau = 1 needs N = {N_tau1:.2e} photons per pulse"
      f" = {E_tau1:.2e} J per pulse of 722 keV gammas")
check("optical depth 1 needs >1 kJ of MeV gammas in a micron spot", E_tau1 > 1e3)

print("\n7. Static fields supply momentum, not energy; a pair from the field alone needs Schwinger's scale")
E_S = m_e**2 * c**3 / (e * hbar)
print(f"   E_S = m^2 c^3/(e hbar) = {E_S:.2e} V/m; work over one reduced Compton wavelength = "
      f"{e * E_S * lambda_bar / e / 1e3:.0f} keV")
check("work over a reduced Compton wavelength at E_S equals m c^2",
      abs(e * E_S * lambda_bar - m_e * c**2) / (m_e * c**2) < 1e-12)

print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILED: ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
