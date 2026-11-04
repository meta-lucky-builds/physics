"""A14, A15, A21(A): how well can matter reflect a 511 keV photon?

Inputs from outside: densities and masses of the elements, and NIST mass attenuation
coefficients for Au and Pb. Everything else is computed. Fresnel reflectivity is computed
from the full complex index, and its small-angle loss is compared with an analytic limit.
"""
import cmath
import math
import pathlib
import sys

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from common.constants import N_A, c, e, h, m_e  # noqa: E402
from scipy.constants import physical_constants  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


r_e = physical_constants["classical electron radius"][0]
E_keV = m_e * c**2 / e / 1e3
lam = h * c / (E_keV * 1e3 * e)
print(f"511 keV photon: lambda = {lam * 1e12:.4f} pm")


def log_interp(E, E1, v1, E2, v2):
    return math.exp(math.log(v1) + (math.log(E) - math.log(E1)) / (math.log(E2) - math.log(E1))
                    * (math.log(v2) - math.log(v1)))


mu_rho_au = log_interp(0.511, 0.5, 0.1530, 0.6, 0.1194)  # cm^2/g, NIST XCOM table
mu_rho_pb = log_interp(0.511, 0.5, 0.1614, 0.6, 0.1248)
print(f"NIST-interpolated mu/rho at 511 keV: Au {mu_rho_au:.4f}, Pb {mu_rho_pb:.4f} cm^2/g")

print("\n1. Refractive decrement and critical angle (free-electron scattering)")
print("   element   n_e [1/m^3]    delta       theta_c [mrad]   theta_c [arcmin]")
materials = {"C": (6, 12.011, 3510.0), "Si": (14, 28.086, 2330.0), "W": (74, 183.84, 19300.0),
             "Au": (79, 196.967, 19320.0), "Pb": (82, 207.2, 11340.0)}
theta_c = {}
delta = {}
for name, (Z, A, rho) in materials.items():
    n_e = rho * 1e3 * N_A * Z / A
    d = r_e * lam**2 * n_e / (2 * math.pi)
    delta[name] = d
    theta_c[name] = math.sqrt(2 * d)
    print(f"   {name:6s} {n_e:12.3e} {d:12.3e} {theta_c[name] * 1e3:12.4f} {math.degrees(theta_c[name]) * 60:15.3f}")
check("gold critical angle at 511 keV is below 0.2 mrad", theta_c["Au"] < 2e-4,
      f"{theta_c['Au'] * 1e3:.3f} mrad")

print("\n2. Fresnel reflectivity for gold, exact complex index")
mu_au = mu_rho_au * 19.32 * 100         # 1/m
beta_au = mu_au * lam / (4 * math.pi)
n_au = complex(1 - delta["Au"], beta_au)
print(f"   delta = {delta['Au']:.3e}, beta = {beta_au:.3e}, beta/delta = {beta_au / delta['Au']:.3e}")


def fresnel(theta, n=n_au):
    """Reflectivity at grazing angle theta (s-polarisation form, valid at small angles too)."""
    kz1 = math.sin(theta)
    kz2 = cmath.sqrt(n * n - math.cos(theta) ** 2)
    return abs((kz1 - kz2) / (kz1 + kz2)) ** 2


def roughness(theta, sigma):
    kz1 = 2 * math.pi / lam * math.sin(theta)
    kz2 = 2 * math.pi / lam * cmath.sqrt(n_au * n_au - math.cos(theta) ** 2)
    return math.exp(-4 * (kz1 * kz2).real * sigma**2)


print("   theta [mrad]   R (smooth)   R (sigma = 0.3 nm, Nevot-Croce)")
thc = theta_c["Au"]
for frac in (0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 3.0, 10.0):
    th = frac * thc
    R = fresnel(th)
    print(f"   {th * 1e3:10.4f} {R:12.6f} {R * roughness(th, 0.3e-9):12.6f}")
R_normal = fresnel(math.pi / 2)
print(f"   normal incidence: R = {R_normal:.3e}, (delta/2)^2 = {(delta['Au'] / 2) ** 2:.3e}")
check("normal-incidence reflectivity follows (delta/2)^2", abs(R_normal / (delta["Au"] / 2) ** 2 - 1) < 0.01)
check("no flat metal mirror reflects 511 keV at normal incidence above 1e-15", R_normal < 1e-15)
check("above 10 theta_c gold reflects less than 1e-4", fresnel(10 * thc) < 1e-4)
print("   Note: at this energy the real and imaginary parts come from incoherent processes"
      " too, so this is the coherent specular part only.")

print("\n3. A closed sphere cannot use grazing incidence (A14)")
for name in ("C", "Au"):
    frac = (1 - math.cos(theta_c[name])) / 2
    print(f"   {name}: photons in a cone of half-angle theta_c around one axis: {frac:.2e} of 4 pi")
check("grazing-capable solid angle is below 1e-7 of the sphere for Au",
      (1 - math.cos(theta_c["Au"])) / 2 < 1e-7)
print("   Mirror placed around the source: radial photons meet a concentric sphere at normal incidence.")

print("\n4. Photon storage ring from many grazing bounces (A15)")
print("   Bending by 2 pi needs N = pi/theta bounces at grazing angle theta.")
print("    theta/theta_c      N       -ln(survival per lap)   predicted 4 pi beta/theta_c^3")
pred = 4 * math.pi * beta_au / thc**3
results = []
for frac in (0.5, 0.1, 0.01, 0.001, 1e-4):
    th = frac * thc
    N = math.pi / th
    R = fresnel(th)
    loss = -N * math.log(R)
    results.append(loss)
    print(f"   {frac:12.4f} {N:10.2e} {loss:16.2f} {pred:20.2f}")
print(f"   Survival per lap in the smooth limit: exp(-{results[-1]:.1f}) = 10^{-results[-1] / math.log(10):.1f}")
check("per-lap loss tends to 4 pi beta/theta_c^3, independent of ring size",
      abs(results[-1] - pred) / pred < 0.02)
check("even a perfectly smooth gold gallery loses the photon in one lap", results[-1] > 50)

print("\n5. Crystal diffraction")
d_si111 = 5.431e-10 / math.sqrt(3)
theta_b = lam / (2 * d_si111)
F_si = 4 * math.sqrt(2) * 10.5        # |F_111| with f = 10.5 near sin(theta)/lambda = 0.16 per angstrom (assumed)
V_si = (5.431e-10) ** 3
width = 2 * r_e * lam**2 * F_si / (math.pi * V_si * math.sin(2 * theta_b))
print(f"   Si(111), d = {d_si111 * 1e10:.3f} A: Bragg angle {theta_b * 1e3:.3f} mrad = {math.degrees(theta_b):.4f} deg")
print(f"   Darwin width estimate {width * 1e6:.2f} urad = {math.degrees(width) * 3600:.3f} arcsec, deflection only 2 theta_B")
check("Si(111) Bragg angle at 511 keV is under 0.25 degree", math.degrees(theta_b) < 0.25)
check("Si(111) accepts under 1 arcsec", math.degrees(width) * 3600 < 1)

d_needed = lam / 2
q = 4 * math.pi / lam
u_rms = 7.5e-12               # assumed Si room-temperature thermal displacement, 0.075 A
a_1s = 5.29177e-11 / 14       # hydrogen-like 1s radius in Si, assumed form of the core electron
form = 1 / (1 + (q * a_1s / 2) ** 2) ** 2
print(f"   true backscatter needs d = lambda/2 = {d_needed * 1e12:.2f} pm, far below any lattice spacing")
print(f"   backscatter momentum q = {q:.3e} 1/m; 1s form factor per core electron {form:.2e}")
print(f"   thermal amplitude factor exp(-q^2 u^2 / 2) = 10^{-(q * u_rms) ** 2 / 2 / math.log(10):.0f}")
check("backscatter Debye-Waller suppression is below 1e-100", (q * u_rms) ** 2 / 2 > 230)

theta_au = lam / (2 * 2.355e-10)
N_laue = 2 * math.pi / (2 * theta_au)
print(f"   Measured Laue peak reflectivity: Au(111) 29% at 588 keV (Barriere et al. 2009)")
print(f"   Au(111) at 511 keV: theta_B = {theta_au * 1e3:.3f} mrad, so {N_laue:.0f} diffractions per 2 pi turn;"
      f" 0.29^{N_laue:.0f} = 10^{N_laue * math.log10(0.29):.0f}")
check("a Laue-crystal ring would lose the photon in a few steps", 0.29 ** N_laue < 1e-100)

print("\n6. Compton scattering is what the photon actually does in matter")


def kn_ds(cos_t, x):
    ratio = 1 / (1 + x * (1 - cos_t))
    return 2 * math.pi * 0.5 * r_e**2 * ratio**2 * (ratio + 1 / ratio - (1 - cos_t**2))


x = E_keV * 1e3 * e / (m_e * c**2)
sigma_kn = quad(lambda ct: kn_ds(ct, x), -1, 1)[0]
mean_loss = quad(lambda ct: kn_ds(ct, x) * (1 - 1 / (1 + x * (1 - ct))), -1, 1)[0] / sigma_kn
analytic = 2 * math.pi * r_e**2 * (
    (1 + x) / x**2 * (2 * (1 + x) / (1 + 2 * x) - math.log(1 + 2 * x) / x)
    + math.log(1 + 2 * x) / (2 * x) - (1 + 3 * x) / (1 + 2 * x) ** 2)
print(f"   Klein-Nishina at 511 keV: sigma = {sigma_kn / 1e-28:.4f} b per electron"
      f" (closed form {analytic / 1e-28:.4f} b), mean energy lost per scatter {mean_loss:.1%}")
check("numerical Klein-Nishina integral matches the closed form", abs(sigma_kn - analytic) / analytic < 1e-6)
mass_compton = sigma_kn * 1e4 * N_A * 79 / 196.967
print(f"   Compton share for Au: {mass_compton:.4f} cm^2/g of NIST total {mu_rho_au:.4f} cm^2/g"
      f" ({mass_compton / mu_rho_au:.0%}); the rest is mainly photoabsorption")

print("\n7. Absorbing is easy: thickness of lead needed to stop 511 keV")
mu_pb = mu_rho_pb * 11.34  # 1/cm
for frac in (0.9, 0.99, 0.999):
    print(f"   stop {frac:.1%}: {-math.log(1 - frac) / mu_pb:.2f} cm of Pb (attenuation, first interaction)")
check("3 cm of Pb attenuates > 99% of 511 keV", math.exp(-mu_pb * 3) < 0.01)

print("\n8. What electron density would make a true metal mirror at 511 keV?")
n_crit = math.pi / (r_e * lam**2)
rho_crit = n_crit * 2 * 1.6726e-27
print(f"   plasma frequency = gamma frequency at n_e = {n_crit:.2e} 1/m^3 = {n_crit / 4.66e30 * 1e-0:.1e} times solid gold,"
      f" mass density ~ {rho_crit:.1e} kg/m^3 for Z/A = 1/2")
check("critical density is at white-dwarf-core scale, above 1e11 kg/m^3", rho_crit > 1e11)
print("   (At that density electrons are relativistic and degenerate, so this is only a scale.)")

print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILED: ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
