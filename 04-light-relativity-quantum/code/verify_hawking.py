"""L10: Hawking temperature by dimensional analysis and Euclidean regularity.

- Solve the unit equations for the only combination of hbar, c, G, M that has
  the units of energy.
- Numerically find the period of Euclidean time that removes the conical singularity
  at the Schwarzschild horizon. This sets the numerical factor 1/(8 pi).
- Numbers: temperatures, wavelength against the hole's size, and a light-based probe
  of quantum gravity (GRB 090510).
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad

from checklib import K, check, finish, section

section("Dimensional analysis")
a, b, d, e = sp.symbols("a b d e")
# exponents of [kg, m, s] in hbar^a c^b G^d M^e, set equal to those of energy (1, 2, -2)
eqs = [
    a - d + e - 1,            # kg: hbar=kg m^2/s, G=m^3/(kg s^2), M=kg
    2 * a + b + 3 * d - 2,    # m
    -a - b - 2 * d + 2,       # s
]
sol = sp.solve(eqs, [a, b, d], dict=True)[0]
print(f"  hbar^a c^b G^d M^e with units of energy: a = {sol[a]}, b = {sol[b]}, d = {sol[d]}, e free")
# Require T to depend on quantum mechanics linearly: a = 1
e_val = sp.solve(sol[a] - 1, e)[0]
result = {s: sol[s].subs(e, e_val) for s in sol}
result[e] = e_val
print(f"  With energy linear in hbar: {result}")
check("kT ~ hbar c^3 / (G M)", result[a] == 1 and result[b] == 3 and result[d] == -1 and result[e] == -1)
print("  Note: units alone leave one free exponent, because (hbar c^5/G)^(1/2) is itself an energy, the")
print("  Planck energy. The step 'linear in hbar' is an extra assumption, and units cannot give the number in front.")

section("Euclidean regularity at the horizon (units c = 1, r_s = 1)")
rs = 1.0


def proper_radial(eps):
    return quad(lambda r: 1 / np.sqrt(1 - rs / r), rs, rs + eps, limit=500)[0]


for eps in [1e-2, 1e-4, 1e-6]:
    rho = proper_radial(eps)
    circ_per_beta = np.sqrt(1 - rs / (rs + eps))
    beta_needed = 2 * np.pi * rho / circ_per_beta
    print(f"  r = r_s + {eps:g}: Euclidean period required for a smooth origin = {beta_needed:.6f} r_s/c")
check("the period is 4 pi r_s / c", abs(beta_needed / (4 * np.pi) - 1) < 1e-4, f"found {beta_needed / np.pi:.5f} pi")
print("  T = hbar/(k beta) with beta = 4 pi r_s/c and r_s = 2GM/c^2 gives T = hbar c^3/(8 pi G M k)")


def T_H(M):
    return K.hbar * K.c**3 / (8 * np.pi * K.G * M * K.k_B)


section("Numbers")
T_sun = T_H(K.M_sun)
print(f"  Solar mass: T = {T_sun:.3e} K")
check("solar-mass hole is near 6e-8 K", 5.9e-8 < T_sun < 6.3e-8)
M_cmb = K.hbar * K.c**3 / (8 * np.pi * K.G * K.k_B * 2.725)
print(f"  Hole as warm as the CMB (2.725 K): M = {M_cmb:.2e} kg, about {M_cmb / 7.35e22:.2f} Moon masses")
print(f"  So every stellar-mass hole now absorbs more CMB than it radiates: T_sun/T_CMB = {T_sun / 2.725:.1e}")
for name, M in [("10 M_sun", 10 * K.M_sun), ("Sgr A* (4.3e6 M_sun)", 4.3e6 * K.M_sun), ("1e12 kg", 1e12)]:
    print(f"  {name}: T = {T_H(M):.2e} K")
lam_over_rs = 2 * np.pi * K.hbar * K.c / (K.k_B * T_sun) / (2 * K.G * K.M_sun / K.c**2)
check("thermal wavelength hc/(kT) = 8 pi^2 r_s, far larger than the hole", abs(lam_over_rs / (8 * np.pi**2) - 1) < 1e-10, f"{lam_over_rs:.1f} r_s")
L_pl = np.sqrt(K.hbar * K.G / K.c**3)
E_pl = np.sqrt(K.hbar * K.c**5 / K.G)
print(f"  Planck length {L_pl:.3e} m, Planck energy {E_pl / K.eV / 1e9:.3e} GeV")

section("Light as a probe: a Planck-suppressed delay from GRB 090510")
H0 = 70e3 / 3.0857e22  # 1/s, round value
Om, OL = 0.3, 0.7
z = 0.903  # redshift of GRB 090510 (general knowledge)
Kz = quad(lambda zz: (1 + zz) / np.sqrt(Om * (1 + zz) ** 3 + OL), 0, z)[0]
dE = 31.0  # GeV, highest-energy photon in the first second
E_QG = E_pl / K.eV / 1e9
dt = dE / E_QG * Kz / H0
print(f"  Linear Lorentz violation at the Planck energy would delay a {dE} GeV photon by {dt:.2f} s")
print(f"  That delay is what makes the burst's sub-second structure a Planck-scale test")
check("delay at the Planck scale is order 1 s, so the test really reaches that scale", 0.3 < dt < 3)

finish()
