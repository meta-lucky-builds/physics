"""C01-C03, C27: Friedmann dynamics, age, horizons, reachability in flat LCDM.

Parameters: Planck 2018 VI (arXiv:1807.06209) H0 = 67.4 km/s/Mpc, Omega_m = 0.315.
Radiation: T_CMB = 2.72548 K (Fixsen 2009) with N_eff = 3.046.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from _harness import Checks
from common.constants import G, c, year, ly, k_B, hbar

C = Checks("Cosmology: Friedmann, age, horizons, reachability")
Mpc = 3.0856775814913673e22
Gyr = 1e9 * year
Gly = 1e9 * ly


class LCDM:
    def __init__(self, H0_kms=67.4, Om=0.315, T=2.72548, Neff=3.046):
        self.H0 = H0_kms * 1e3 / Mpc
        rho_crit = 3 * self.H0**2 / (8 * np.pi * G)
        rho_gamma = (np.pi**2 / 15) * (k_B * T)**4 / (hbar**3 * c**5)
        self.Or = rho_gamma / rho_crit * (1 + 0.2271 * Neff)
        self.Om = Om
        self.OL = 1 - Om - self.Or

    def E(self, a):
        return np.sqrt(self.Or / a**4 + self.Om / a**3 + self.OL)

    def age(self, a=1.0):
        # Integrate in ln(a): the integrand is smooth over many decades.
        return quad(lambda u: 1 / (self.H0 * self.E(np.exp(u))), -60, np.log(a),
                    limit=800, epsabs=0, epsrel=1e-10)[0]

    def chi(self, a1, a2):
        """Comoving distance light covers between two scale factors, in metres."""
        return quad(lambda u: c * np.exp(-u) / (self.H0 * self.E(np.exp(u))),
                    np.log(a1), np.log(a2), limit=800, epsabs=0, epsrel=1e-10)[0]

    def chi_event(self, a=1.0):
        return self.chi(a, a * 1e12)

    def chi_particle(self, a=1.0):
        return self.chi(1e-20, a)


m = LCDM()

C.section("Newtonian shell versus Friedmann equation (matter + Lambda)")
# Integrate the shell's acceleration (Newton's inverse square inside a uniform ball, plus a
# repulsive Lambda term) without using the first Friedmann equation, then test whether that
# equation emerges as a constant of the motion.
Om, OL = 0.315, 0.685
H0 = m.H0


def rhs(t, y):
    a, v = y
    return [v, -0.5 * Om * H0**2 / a**2 + OL * H0**2 * a]


sol = solve_ivp(rhs, [0, 40 * Gyr], [1.0, H0], rtol=1e-11, atol=1e-13, method="DOP853",
                dense_output=True)
ts = np.linspace(0, 40 * Gyr, 200)
a_t, v_t = sol.sol(ts)
K = (v_t**2 - H0**2 * (Om / a_t + OL * a_t**2)) / H0**2
print(f"  a after 40 Gyr = {a_t[-1]:.3f}; max |k c^2 / H0^2| over the run = {np.abs(K).max():.2e}")
C.check("a'^2 - (8 pi G rho/3) a^2 - Lambda a^2/3 stays constant (0 for flat)",
        np.abs(K).max() < 1e-8)

C.section("Age and horizons (Planck 2018 parameters)")
t0 = m.age() / Gyr
print(f"  Omega_r = {m.Or:.3e}, Omega_Lambda = {m.OL:.4f}, 1/H0 = {1/m.H0/Gyr:.2f} Gyr")
print(f"  age t0 = {t0:.3f} Gyr (Planck base-LCDM table: 13.797 +/- 0.023)")
C.check("age agrees with Planck to 0.3%", abs(t0 - 13.797) / 13.797 < 0.003)
# Independent analytic check: without radiation the flat-LCDM age has a closed form.
m0 = LCDM(T=1e-6)
ana = 2 / (3 * m0.H0 * np.sqrt(m0.OL)) * np.arcsinh(np.sqrt(m0.OL / m0.Om)) / Gyr
C.check(f"numeric age without radiation matches the analytic asinh formula ({ana:.4f} Gyr)",
        abs(m0.age() / Gyr - ana) < 1e-4)


def march_photon(model, n=400000):
    """Fixed-grid trapezoid in ln a, a different integration route from adaptive quadrature."""
    u = np.linspace(np.log(1e-14), 0, n)
    f = np.exp(-u) / (model.H0 * model.E(np.exp(u)))
    return c * np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(u))


chi_p = m.chi_particle()
chi_m = march_photon(m)
print(f"  particle horizon = {chi_p/Gly:.2f} Gly; diameter = {2*chi_p/Gly:.1f} Gly; c*t0 = {t0:.2f} Gly")
C.check("particle horizon: adaptive quadrature and fixed-grid march agree to 0.01%",
        abs(chi_p - chi_m) / chi_p < 1e-4)
C.check("diameter lies in 90-95 Gly and exceeds 2 c t0",
        90 < 2 * chi_p / Gly < 95 and chi_p / Gly > t0)

zrec = 1090
chi_rec = m.chi(1 / (1 + zrec), 1)
print(f"  CMB light: source now at {chi_rec/Gly:.1f} Gly; proper distance at emission "
      f"{chi_rec/(1+zrec)/(1e6*ly):.1f} Mly")
vrec = m.H0 * m.E(1 / (1 + zrec)) * chi_rec / (1 + zrec) / c
v0 = m.H0 * chi_rec / c
print(f"  its recession rate: {vrec:.1f} c at emission, {v0:.2f} c now (coordinate rates, not local speeds)")
C.check("that source recedes faster than c, yet its light reached us", v0 > 1 and vrec > 1)

# Local speed limit: for a photon moving toward us, d(a chi)/dt = H (a chi) - c.
# Evaluate the left side by finite difference along an actual inbound ray.
a_probe = 0.5
d_before = a_probe * m.chi(a_probe, 1)
da = a_probe * 1e-6
dt = da / (a_probe * m.H0 * m.E(a_probe))
d_after = (a_probe + da) * m.chi(a_probe + da, 1)
rate = (d_after - d_before) / dt
expect = m.H0 * m.E(a_probe) * d_before - c
print(f"  at a = 0.5: d(proper distance)/dt = {rate/c:+.4f} c; H d - c = {expect/c:+.4f} c")
C.check("proper-distance rate = H d - c: the photon moves at c through local space",
        abs(rate - expect) / c < 1e-4)

C.section("Hubble tension sensitivity")
m_sh = LCDM(H0_kms=73.04)
print(f"  c/H0 = {c/m.H0/Gly:.2f} Gly (Planck), {c/m_sh.H0/Gly:.2f} Gly (SH0ES)")
print(f"  age with H0 = 73.04, Om = 0.315 held fixed = {m_sh.age()/Gyr:.2f} Gyr (illustrative, not a fit)")
sig = (73.04 - 67.4) / np.hypot(1.04, 0.5)
print(f"  tension = {73.04-67.4:.2f} km/s/Mpc = {sig:.1f} sigma from the quoted errors")
C.check("the two H0 values differ by close to the 5 sigma Riess et al. quote", 4.5 < sig < 5.5)

C.section("Lookback, Hubble radius, event horizon, the 94% figure")
for z in [1, 3, 10, 14.44]:
    print(f"  z = {z:5}: lookback {(m.age() - m.age(1/(1+z)))/Gyr:6.2f} Gyr, "
          f"comoving distance {m.chi(1/(1+z),1)/Gly:5.1f} Gly")
chi_e = m.chi_event()
RH = c / m.H0
print(f"  Hubble radius c/H0 = {RH/Gly:.2f} Gly; comoving event horizon = {chi_e/Gly:.2f} Gly")
f_all = (chi_e / chi_p)**3
print(f"  unreachable fraction of the observable VOLUME = 1-(chi_e/chi_p)^3 = {100*(1-f_all):.1f}%")
chi_g = m.chi(1 / (1 + 14.44), 1)
f_gal = (chi_e / chi_g)**3
print(f"  counting only out to the most distant confirmed galaxy, z = 14.44: {100*(1-f_gal):.1f}%")
C.check("student's 94% lies between the two counting assumptions",
        (1 - f_gal) * 100 < 94 < (1 - f_all) * 100)
# Direct test of the event horizon: launch a photon now and march it with a fixed grid.
u = np.linspace(0, np.log(1e8), 2000001)
f = np.exp(-u) / (m.H0 * m.E(np.exp(u)))
reach = c * np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(u))
print(f"  a photon launched now and followed to a = 1e8 covers {reach/Gly:.2f} Gly of comoving distance")
C.check("marched photon saturates at the computed event horizon (agree to 0.1%)",
        reach < chi_e and abs(reach - chi_e) / chi_e < 1e-3)

C.section("Future reachability of today's observable volume")
print(f"  {'cosmic age (Gyr)':>17} {'a':>6} {'event horizon (Gly)':>20} {'reachable % of volume':>22}")
for a in [1, 2, 4, 8, 16, 64]:
    chi_ea = m.chi_event(a)
    print(f"  {m.age(a)/Gyr:17.1f} {a:6d} {chi_ea/Gly:20.2f} {100*(chi_ea/chi_p)**3:22.4f}")


def time_when(frac):
    aa = brentq(lambda a: (m.chi_event(a) / chi_p)**3 - frac, 1, 1e6)
    return m.age(aa) / Gyr


t99, t9999 = time_when(0.01), time_when(1e-4)
print(f"  99% unreachable at cosmic age {t99:.0f} Gyr; 99.99% unreachable at age {t9999:.0f} Gyr")
print(f"  late-time e-folding of the scale factor = {1/(m.H0*np.sqrt(m.OL))/Gyr:.1f} Gyr")
C.check("99.99% unreachable arrives in tens of Gyr, not trillions of years", 20 < t9999 < 200)

C.section("Andromeda: light delay, expansion and binding")
D = 2.5e6 * ly
print(f"  light delay = {D/c/year:.3g} yr; Hubble flow at 2.5 Mly would be {m.H0*D/1e3:.0f} km/s recession")
print("  measured radial velocity = -109.3 km/s (van der Marel et al. 2012), the opposite sign")
Msun = 1.98841e30
Mneed = (109.3e3)**2 * D / (2 * G) / Msun
print(f"  mass that makes -109 km/s bound at this separation = {Mneed:.2e} Msun; "
      "measured M_LG = 3.17e12 Msun")
C.check("Local Group mass exceeds the binding mass, so the Hubble flow does not apply here",
        3.17e12 > Mneed)
print(f"  fractional stretch of M31's light if expansion did apply, H0 D/c = {m.H0*D/c:.1e}")
print(f"  fractional change of the CMB temperature over 2000 yr = {m.H0*2000*year:.1e}")

C.section("Backward extrapolation: scale factor reaches 0 at finite time")
print(f"  age at a = 1e-12 is {m.age(1e-12)/year:.2e} yr")
C.check("the integral of da/(aH) converges as a -> 0: finite age within this model",
        m.age(1e-12) < 1e-6 * m.age())
C.finish()
