"""Entropy: definitions, not-fluctuations, black hole entropy from T, scaling, potentials, Sakharov toy."""
import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln, zeta

from bh_tools import C, M_M87, M_SGRA, Report, T_hawking, l_P, r_g

R = Report("verify_entropy: definitions, scaling, potentials")


def S_vN(rho):
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-14]
    return float(-(w * np.log(w)).sum())


# ---------------------------------------------------------------- entropy is not fluctuation
R.section("1. Entropy and fluctuations are different quantities (k=1, hbar=m=omega=1)")
N = 60
a = np.diag(np.sqrt(np.arange(1, N)), 1)
xop = (a + a.T) / np.sqrt(2)
H = np.diag(np.arange(N) + 0.5)
ground = np.zeros((N, N))
ground[0, 0] = 1.0
varx = np.trace(ground @ xop @ xop) - np.trace(ground @ xop) ** 2
print(f"  Oscillator ground state: S = {S_vN(ground):.6f}, Var(x) = {varx:.6f}, Var(p) = 0.5 too. Pure state, real fluctuations, zero entropy.")
R.check("Vacuum-like state: S = 0 and Var(x) = 1/2", abs(S_vN(ground)) < 1e-12 and abs(varx - 0.5) < 1e-12)
g = 16
micro = np.eye(g) / g
Hm = np.eye(g) * 3.7
varE = np.trace(micro @ Hm @ Hm) - np.trace(micro @ Hm) ** 2
print(f"  Equal mixture of {g} degenerate states: S = {S_vN(micro):.4f} = ln {g}, Var(E) = {varE:.2e}. Entropy with no energy fluctuation.")
R.check("Degenerate mixture: S = ln 16 and Var(E) = 0", abs(S_vN(micro) - np.log(16)) < 1e-12 and abs(varE) < 1e-12)
Ts = np.array([0.2, 0.42, 0.8, 2, 10, 100])
print("  Two-level system, gap 1:  T      S(T)     Var(E)")
Ss, Vs = [], []
for T_ in Ts:
    p = 1 / (1 + np.exp(1 / T_))
    Ss.append(-(p * np.log(p) + (1 - p) * np.log(1 - p)))
    Vs.append(p * (1 - p))
    print(f"                       {T_:6.2f}   {Ss[-1]:.4f}   {Vs[-1]:.4f}")
R.check("S rises monotonically with T", all(np.diff(Ss) > 0))
Tf = np.linspace(0.05, 200, 200000)
pf = 1 / (1 + np.exp(1 / Tf))
Vf = pf * (1 - pf)
Tcv = Tf[np.argmax(Vf / Tf**2)]
print(f"  Heat capacity C = Var(E)/T^2 peaks at T = {Tcv:.3f} (Schottky peak) and then falls, while S keeps rising: no function S(fluctuation).")
R.check("C(T) has an interior maximum, S has none", 0.3 < Tcv < 0.5)

# ---------------------------------------------------------------- BH entropy from first law
R.section("2. Black hole entropy obtained by integrating dS = c^2 dM / T_H")
print("  M/M_sun     S numeric/k       k A/(4 l_P^2)/k    S/(k ln 2) bits")
ok = True
for name, M in (("1", C.M_sun), ("10", 10 * C.M_sun), ("Sgr A*", M_SGRA), ("M87*", M_M87)):
    Sn = quad(lambda m: C.c**2 / T_hawking(m * M) * M, 0, 1)[0] / C.k_B
    A = 4 * np.pi * (2 * r_g(M)) ** 2
    Sa = A / (4 * l_P**2)
    ok &= abs(Sn / Sa - 1) < 1e-9
    print(f"  {name:8s} {Sn:14.4e} {Sa:16.4e} {Sa / np.log(2):16.4e}")
R.check("First law with Hawking T reproduces S = k A/(4 l_P^2)", ok)
S_BH_sun = 4 * np.pi * (2 * r_g(C.M_sun)) ** 2 / (4 * l_P**2)
rb = 2 * r_g(C.M_sun)
S_bek = 2 * np.pi * rb * C.M_sun * C.c**2 / (C.hbar * C.c)
R.check("Schwarzschild saturates Bekenstein's bound S <= 2 pi k R E/(hbar c)", abs(S_bek / S_BH_sun - 1) < 1e-12)
Ai = 36**2 + 29**2
print(f"  GW150914 (central values, no initial spin): S_final/S_initial = {0.5 * 62**2 * (1 + np.sqrt(1 - 0.67**2)) / Ai:.2f}")


# Sackur-Tetrode for a fully ionised solar-mass gas at the Sun's mean density
def sackur(N_, m, T_, V):
    nQ = (m * C.k_B * T_ / (2 * np.pi * C.hbar**2)) ** 1.5
    return N_ * (np.log(V * nQ / N_) + 2.5)


R_sun = 6.96e8
V_sun = 4 / 3 * np.pi * R_sun**3
Np = C.M_sun / C.m_p
Tmean = 4.5e6  # K, rough mass-weighted figure: an estimate
S_sun = sackur(Np, C.m_p, Tmean, V_sun) + sackur(Np, C.m_e, Tmean, V_sun)
print(f"  Sun as ideal ionised H at T ~ 4.5e6 K (estimate): S/k = {S_sun:.2e}, i.e. {S_sun / Np:.0f} per baryon")
print(f"  Same mass as a black hole: S/k = {S_BH_sun:.2e}; ratio {S_BH_sun / S_sun:.1e}")
R.check("Collapsing a star to a hole raises entropy by ~ 1e19", 1e18 < S_BH_sun / S_sun < 1e20)

# ---------------------------------------------------------------- scaling
R.section("3. S against mass: extensive matter vs black hole")
Ms = np.geomspace(1, 1e6, 7) * C.M_sun
S_gas = np.array([sackur(m / C.m_p, C.m_p, Tmean, V_sun * m / C.M_sun) for m in Ms])
S_bh = np.array([4 * np.pi * (2 * r_g(m)) ** 2 / (4 * l_P**2) for m in Ms])
sl_g = np.polyfit(np.log(Ms), np.log(S_gas), 1)[0]
sl_b = np.polyfit(np.log(Ms), np.log(S_bh), 1)[0]
R.check("Gas at fixed T and density: S ~ M^1", abs(sl_g - 1) < 1e-6, f"slope {sl_g:.4f}")
R.check("Black hole: S ~ M^2", abs(sl_b - 2) < 1e-9, f"slope {sl_b:.4f}")
S1 = 4 * np.pi * (2 * r_g(C.M_sun)) ** 2 / (4 * l_P**2)
S2 = 4 * np.pi * (2 * r_g(2 * C.M_sun)) ** 2 / (4 * l_P**2)
R.check("Black hole entropy is not extensive: S(2M) = 2 x [S(M)+S(M)]", abs(S2 / (2 * S1) - 2) < 1e-12)

# ---------------------------------------------------------------- local decrease
R.section("4. Local entropy decrease with export: freezing water at 263 K")
dH, Tm, Ts_ = 6010.0, 273.15, 263.15
dSw, dSs = -dH / Tm, dH / Ts_
print(f"  Water (latent heat 6.01 kJ/mol, heat-capacity terms neglected): dS_water = {dSw:.2f}, dS_surroundings = {dSs:.2f}, total = {dSw + dSs:+.2f} J/(mol K)")
R.check("Crystal entropy falls, total rises", dSw < 0 and dSw + dSs > 0)
print(f"  Landauer: erasing 1 bit at 300 K costs at least kT ln 2 = {C.k_B * 300 * np.log(2):.2e} J")

# ---------------------------------------------------------------- radiation entropy
R.section("5. Entropy carried by blackbody photons (Zurek's 4/3)")
u_int = quad(lambda x: x**3 / np.expm1(x), 0, 600, limit=200)[0]
s_int = quad(lambda x: x**2 * (x / np.expm1(x) - np.log(-np.expm1(-x))), 0, 600, limit=200)[0]
print(f"  From the Planck distribution: S T / U = {s_int / u_int:.6f}")
R.check("Blackbody radiation carries S = (4/3) U/T, more than the hole's dS = dU/T", abs(s_int / u_int - 4 / 3) < 1e-6)

# ---------------------------------------------------------------- entropic force, Massieu
R.section("6. 'Entropy potential' uses that do exist")
Nchain, bl = 10000, 1.0
xs = np.arange(0, 400, 2)
Sx = gammaln(Nchain + 1) - gammaln((Nchain + xs) / 2 + 1) - gammaln((Nchain - xs) / 2 + 1)
F_num = np.gradient(Sx, xs)          # f/(kT) = dS/dx
F_ana = -xs / (Nchain * bl**2)
err = np.abs(F_num[5:-5] - F_ana[5:-5]).max()
R.check("Polymer entropic force from counting chains: f = -kT x/(N b^2), a real force with no energy term", err < 2e-4, f"max error {err:.1e}")
T_ = 1.3
Z = lambda b: 1 + np.exp(-b)
beta = 1 / T_
p1 = np.exp(-beta) / Z(beta)
E_ = p1
S_ = -(p1 * np.log(p1) + (1 - p1) * np.log(1 - p1))
Mass = S_ - E_ / T_
R.check("Massieu function S - E/T equals ln Z", abs(Mass - np.log(Z(beta))) < 1e-12)
h = 1e-5
dM = (np.log(Z(beta + h)) - np.log(Z(beta - h))) / (2 * h)
R.check("d(Massieu)/d(1/T) = -E: it really is an entropy-based potential", abs(dM + E_) < 1e-8)

# ---------------------------------------------------------------- Sakharov
R.section("7. Baryon asymmetry numbers and the equilibrium argument")
Ob_h2 = 0.02237
H100 = 100e3 / 3.0857e22
rho_c1 = 3 * H100**2 / (8 * np.pi * C.G)
nb = Ob_h2 * rho_c1 / C.m_p
Tcmb = 2.7255
ng = (2 * zeta(3) / np.pi**2) * (C.k_B * Tcmb / (C.hbar * C.c)) ** 3
print(f"  n_b = {nb:.3f} /m^3, n_gamma = {ng:.3e} /m^3, eta = n_b/n_gamma = {nb / ng:.2e}")
R.check("eta about 6e-10", 5.5e-10 < nb / ng < 6.6e-10)
E_b = np.array([1.0, 1.0, 2.5, 2.5])
Bq = np.array([1, -1, 1, -1])
for d, lab in ((0.0, "CPT-symmetric spectrum"), (0.05, "particle/antiparticle masses split by hand")):
    En = E_b + d * Bq
    w = np.exp(-En / 0.7)
    print(f"  {lab}: thermal <B> = {(w * Bq).sum() / w.sum():+.4f}")
w = np.exp(-E_b / 0.7)
R.check("In equilibrium with CPT, <B> = 0 exactly", abs((w * Bq).sum()) < 1e-15)

# ---------------------------------------------------------------- budget
R.section("8. Comparison with the universe's entropy budget")
print(f"  10 M_sun hole: {S_BH_sun * 100:.1e} k.  One Sgr A*: {4 * np.pi * (2 * r_g(M_SGRA)) ** 2 / (4 * l_P**2):.1e} k.  M87*: {4 * np.pi * (2 * r_g(M_M87)) ** 2 / (4 * l_P**2):.1e} k.")

R.close()
