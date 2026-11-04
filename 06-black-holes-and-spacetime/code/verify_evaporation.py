"""Hawking temperature, evaporation ODE, lifetimes, CMB balance, Page curve toy, nested holes."""
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import digamma

from bh_tools import C, M_M87, M_SGRA, Report, T_hawking, l_P, m_P, r_g, t_evap, t_g

R = Report("verify_evaporation: Hawking temperature, lifetimes, final stages, Page curve")

# ---------------------------------------------------------------- temperature
R.section("1. Temperature from Euclidean regularity at the horizon (M=1 numerics)")
# Near r=2M with tE = i t: ds^2 = rho^2 (dtE/4M)^2 + drho^2 + ..., rho = proper distance.
# No conical singularity requires period 8 pi M. Check rho(r) numerically: f ~ (rho/4M)^2.
for x in (1e-3, 1e-5):
    rr = 2 + x
    rho = quad(lambda q: 1 / np.sqrt(1 - 2 / q), 2, rr)[0]
    print(f"  r-2M = {x:.0e}: sqrt(f)/(rho/4M) = {np.sqrt(1 - 2 / rr) / (rho / 4):.6f}")
    last = np.sqrt(1 - 2 / rr) / (rho / 4)
R.check("Near-horizon metric is flat polar form with angular variable tE/4M", abs(last - 1) < 1e-4)
print("  so tE has period 8 pi GM/c^3, and T = hbar/(k * period) = hbar c^3/(8 pi G M k)")
print(f"  T_H(1 M_sun) = {T_hawking(C.M_sun):.3e} K")
d = 1.0
tol = C.hbar * C.c / (2 * np.pi * C.k_B * d)
print(f"  Local check: hover at proper distance d=1 m: Tolman T_H/sqrt(f) = Unruh hbar a/(2 pi c k), a = c^2/d = {tol:.3e} K")
Mt = 1e6 * C.M_sun
rs = 2 * r_g(Mt)
# exact radius at proper distance 1 m: d ~ 2 sqrt(rs x) for x<<rs
x = d**2 / (4 * rs)
sf = np.sqrt(x / rs)
R.check("Tolman-redshifted T_H equals Unruh temperature for the hovering observer", abs(T_hawking(Mt) / sf / tol - 1) < 1e-6)

# ---------------------------------------------------------------- ODE
R.section("2. Mass-loss ODE integrated numerically, with power from a numerically integrated Planck law")
planck = quad(lambda x: x**3 / np.expm1(x), 0, 600, limit=200)[0]   # pi^4/15
sb = C.k_B**4 / (4 * np.pi**2 * C.c**2 * C.hbar**3) * planck
R.check("Planck integral gives Stefan-Boltzmann constant", abs(sb / C.sigma_SB - 1) < 1e-6)


def power(M):
    A = 4 * np.pi * (2 * r_g(M)) ** 2
    return A * sb * T_hawking(M) ** 4


M0 = 1e11  # kg
tc = t_evap(M0)
s = solve_ivp(lambda t, y: [-power(y[0]) / C.c**2], [0, tc * (1 - 1e-6)], [M0], rtol=1e-11, atol=1e3, method="DOP853",
              dense_output=True)
ok = True
for fr in (0.5, 0.9, 0.99, 0.999):
    mn = s.sol(fr * tc)[0]
    mc = M0 * (1 - fr) ** (1 / 3)
    print(f"  t/t_evap = {fr:5.3f}: M/M0 numeric {mn / M0:.6f}, closed form (1-t/t_evap)^(1/3) {mc / M0:.6f}")
    ok &= abs(mn / mc - 1) < 1e-5
R.check("ODE agrees with t = 5120 pi G^2 M^3/(hbar c^4)", ok)
R.check("Power formula P = hbar c^6/(15360 pi G^2 M^2)", abs(power(M0) / (C.hbar * C.c**6 / (15360 * np.pi * C.G**2 * M0**2)) - 1) < 1e-6)

alpha_page = 2e-4
print(f"  Page 1976 total power 2e-4 hbar c^6/(G^2 M^2), vs blackbody photon-only {1 / (15360 * np.pi):.3e}: ratio {alpha_page * 15360 * np.pi:.1f}")
print(f"  (massless nu, photon, graviton, with greybody factors).  Geometric-optics area 27/16 alone would give 1.69.")

print("\n  Object     M/M_sun     T_H(K)       t_evap blackbody (yr)    Page-corrected (yr)   T_H/T_CMB")
Tcmb = 2.7255
rows = (("1 M_sun", C.M_sun), ("10 M_sun", 10 * C.M_sun), ("Sgr A*", M_SGRA), ("M87*", M_M87), ("2e12 M_sun", 2.04e12 * C.M_sun))
for nm, M in rows:
    tb = t_evap(M) / C.year
    tp = C.G**2 * M**3 / (3 * alpha_page * C.hbar * C.c**4) / C.year
    print(f"  {nm:10s} {M / C.M_sun:9.3g}  {T_hawking(M):10.3e}  {tb:18.3e}  {tp:18.3e}  {T_hawking(M) / Tcmb:12.2e}")
tH = 13.8e9 * C.year
M_today_bb = (tH * C.hbar * C.c**4 / (5120 * np.pi * C.G**2)) ** (1 / 3)
M_today_pg = (tH * 3 * alpha_page * C.hbar * C.c**4 / C.G**2) ** (1 / 3)
print(f"\n  Mass that evaporates in 13.8 Gyr: blackbody {M_today_bb * 1e3:.2e} g; Page's power {M_today_pg * 1e3:.2e} g; Page's published (5+-1)e14 g")
R.check("Age-of-universe evaporation mass is ~1e14-1e15 g", 1e14 < M_today_pg * 1e3 < 1e15)
R.check("Known holes: lifetimes exceed the age of the universe by > 1e50", t_evap(3 * C.M_sun) / tH > 1e50)

# ---------------------------------------------------------------- CMB
R.section("3. Net energy balance today: absorbers, not evaporators")
M_eq = C.hbar * C.c**3 / (8 * np.pi * C.G * C.k_B * Tcmb)
print(f"  T_H = T_CMB at M = {M_eq:.2e} kg (Moon = 7.3e22 kg)")
R.check("Every stellar hole is colder than the CMB", T_hawking(3 * C.M_sun) < Tcmb)
M = 10 * C.M_sun
A = 27 * np.pi * r_g(M) ** 2
Pabs = A * C.sigma_SB * Tcmb**4
Pem = power(M)
print(f"  10 M_sun: absorbs ~{Pabs:.2e} W from CMB (capture area 27 pi M^2), emits ~{Pem:.2e} W; ratio {Pabs / Pem:.1e}")
h, Om = 0.6736, 0.3153
H0 = h * 100e3 / 3.0857e22
OL = 1 - Om


def t_of_a(a):
    return quad(lambda q: 1 / (q * H0 * np.sqrt(Om / q**3 + OL)), 1, a, limit=200)[0] / C.year / 1e9


print("  Time from now until the CMB cools to T_H (Planck 2018 LCDM, no future surprises):")
for nm, M in (("10 M_sun", 10 * C.M_sun), ("M87*", M_M87)):
    a_ = Tcmb / T_hawking(M)
    print(f"   {nm}: scale factor {a_:.2e}, t = {t_of_a(a_):.0f} Gyr")
HL = H0 * np.sqrt(OL)
T_dS = C.hbar * HL / (2 * np.pi * C.k_B)
M_dS = C.hbar * C.c**3 / (8 * np.pi * C.G * C.k_B * T_dS) / C.M_sun
print(f"  de Sitter horizon temperature (Gibbons-Hawking) hbar H/(2 pi k) = {T_dS:.2e} K; T_H falls below this only above {M_dS:.1e} M_sun")

# ---------------------------------------------------------------- final stage
R.section("4. Last moments (blackbody formula, photon-only; species raise the power)")
print("  time left    M left (kg)    E left (J)     T (K)         kT (GeV)")
coef = 5120 * np.pi * C.G**2 / (C.hbar * C.c**4)
for tl in (3.156e7, 1.0, 1e-3, 1e-6):
    Ml = (tl / coef) ** (1 / 3)
    Tl = T_hawking(Ml)
    print(f"  {tl:9.3e} s  {Ml:11.3e}  {Ml * C.c**2:11.3e}  {Tl:11.3e}  {C.k_B * Tl / (C.eV * 1e9):11.3e}")
M1s = (1.0 / coef) ** (1 / 3)
print(f"  Energy released in the last second: {M1s * C.c**2:.2e} J = {M1s * C.c**2 / 4.184e15:.1e} Mt TNT")
R.check("Last-second energy is about 2e22 J in this approximation", 1.5e22 < M1s * C.c**2 < 2.5e22)
print(f"  Semiclassical control: mean quantum ~2.7 kT over Mc^2 = 2.7/(8 pi) (m_P/M)^2; reaches 1 at M = {np.sqrt(2.7 / (8 * np.pi)) * m_P:.2e} kg = {np.sqrt(2.7 / (8 * np.pi)):.2f} m_P")
t_at_mp = coef * m_P**3
print(f"  Blackbody time left at M = m_P: {t_at_mp:.1e} s; T there = {T_hawking(m_P):.1e} K")

# ---------------------------------------------------------------- nested
R.section("5. A small hole inside a big one")
print("  Outer hole    pi*M (s)    inner mass that can evaporate inside (kg)    its r_s (m)")
for nm, M in (("Sgr A*", M_SGRA), ("M87*", M_M87)):
    tm = np.pi * t_g(M)
    mi = (tm / coef) ** (1 / 3)
    print(f"  {nm:10s}  {tm:.3e}   {mi:.3e}                              {2 * r_g(mi):.3e}")
    R.check(f"{nm}: such an inner hole is far smaller than the outer curvature radius", 2 * r_g(mi) < 1e-6 * r_g(M))

# ---------------------------------------------------------------- Page curve
R.section("6. Page curve: toy model with a Haar-random pure state of 14 qubits")
rng = np.random.default_rng(14)
n = 14
samples = 3
print("  k radiation qubits   <S>/ln2 numeric   Page formula   'Hawking' thermal count k")
vals = []
for k in range(0, n + 1):
    Sk, f = 0.0, 0.0
    dA, dB = 2**k, 2 ** (n - k)
    m_, n_ = min(dA, dB), max(dA, dB)
    page = digamma(m_ * n_ + 1) - digamma(n_ + 1) - (m_ - 1) / (2 * n_)
    for _ in range(samples):
        psi = rng.normal(size=2**n) + 1j * rng.normal(size=2**n)
        psi /= np.linalg.norm(psi)
        sv = np.linalg.svd(psi.reshape(dA, dB), compute_uv=False)
        p = sv**2
        p = p[p > 1e-16]
        Sk += -(p * np.log(p)).sum() / samples
    vals.append((k, Sk))
    print(f"   {k:2d}                  {Sk / np.log(2):8.3f}        {page / np.log(2):8.3f}         {k:2d}")
    if 0 < k < n:
        R.check(f"Page average matches at k={k}", abs(Sk - page) < 0.02)
R.check("Entropy returns to zero when all qubits are in the radiation", abs(vals[-1][1]) < 1e-9)
R.check("Maximum at the halfway point", max(vals, key=lambda q: q[1])[0] == n // 2)
print("  Page time for a hole, with radiation entropy = beta x hole entropy lost:")
for beta in (1.0, 4 / 3):
    mp_ = np.sqrt(beta / (1 + beta))
    print(f"   beta = {beta:.2f}: S_rad = S_BH at M/M0 = {mp_:.3f}, t/t_evap = {1 - mp_**3:.3f}")
print("  Scrambling time (r_s/c) ln(S/k), Hayden-Preskill scale:")
for nm, M in (("10 M_sun", 10 * C.M_sun), ("Sgr A*", M_SGRA)):
    Sb = 4 * np.pi * (2 * r_g(M)) ** 2 / (4 * l_P**2)
    print(f"   {nm}: {2 * t_g(M) * np.log(Sb):.2e} s")

# ---------------------------------------------------------------- information has no separate stress energy
R.section("7. Energy per bit swallowed")
for nm, M in (("10 M_sun", 10 * C.M_sun), ("M87*", M_M87)):
    dE = C.k_B * T_hawking(M) * np.log(2)
    print(f"  {nm}: adding one bit of horizon entropy needs dM c^2 = k T_H ln 2 = {dE:.2e} J, a quantum of wavelength ~ {C.h * C.c / dE:.2e} m ({C.h * C.c / dE / (2 * r_g(M)):.1f} r_s)")

R.close()
