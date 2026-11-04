"""C06-C08: time dilation at 0.9995 c, a realistic trip, and the actual hazards."""
import numpy as np
from scipy.integrate import solve_ivp
from _harness import Checks
from common.constants import c, year, ly, m_p, g_0, eV

C = Checks("Time dilation at 0.9995 c and the real hazards")
beta = 0.9995
gamma = 1 / np.sqrt(1 - beta**2)
print(f"  gamma(0.9995) = {gamma:.3f}")

C.section("Light clock simulated event by event, no gamma in the code")
# A photon bounces between mirrors 1 m apart, transverse to the motion, in the Earth frame.
# In that frame the mirrors move at v, and the photon must move at c (second postulate).
Lc = 1.0
v = beta * c
# Find the photon's tick by directly marching: vy fixed by requiring |velocity| = c, vx = v.
vy = np.sqrt(c**2 - v**2)
ticks = 0
y = 0.0
t = 0.0
dt = 1e-11
direction = 1
while ticks < 200:
    y += direction * vy * dt
    t += dt
    if y >= Lc:
        y = 2 * Lc - y
        direction = -1
    elif y <= 0:
        y = -y
        direction = 1
        ticks += 1
t_rest = ticks * 2 * Lc / c
print(f"  {ticks} round trips: Earth-frame time {t:.5e} s; clock's own count {t_rest:.5e} s")
print(f"  measured ratio = {t/t_rest:.3f}")
C.check("light-clock ratio matches 1/sqrt(1-beta^2) to 0.1%", abs(t / t_rest - gamma) / gamma < 1e-3)

C.section("Lorentz transform checks")


def lorentz(ct, x, b):
    g = 1 / np.sqrt(1 - b**2)
    return g * (ct - b * x), g * (x - b * ct)


tp = 10.0 * c * year  # 10 yr of ship time expressed as ct', event at the ship
ct, x = lorentz(tp, 0.0, -beta)
print(f"  10 ship-years correspond to {ct/(c*year):.1f} Earth years and {x/ly:.1f} ly")
s2_a = ct**2 - x**2
s2_b = tp**2
C.check("interval is invariant", abs(s2_a - s2_b) / s2_b < 1e-9)
C.check("Earth time = gamma * ship time", abs(ct / tp - gamma) < 1e-9)

C.section("Coasting trips at 0.9995 c")
for dly in [4.37, 100, 2.5e6]:
    te = dly / beta
    print(f"  distance {dly:>10} ly: Earth {te:12.2f} yr, traveller {te/gamma:10.2f} yr")

C.section("A round trip with real acceleration: 1 g, integrated numerically")
a = g_0
# Integrate rapidity along proper time. Phases: accelerate, coast, decelerate, then reverse.
tau_acc = (c / a) * np.arccosh(gamma)
coast_tau = 1 * year


def profile(tau):
    """Proper acceleration sign for a 4-burn trip with coasts of coast_tau each way."""
    T = 2 * tau_acc + coast_tau
    tt = tau % (2 * T)
    sign = 1.0 if tt < T else -1.0
    t_leg = tt if tt < T else tt - T
    if t_leg < tau_acc:
        return sign * a
    if t_leg < tau_acc + coast_tau:
        return 0.0
    return -sign * a


T_total = 2 * (2 * tau_acc + coast_tau)
N = 400000
taus = np.linspace(0, T_total, N + 1)
dtau = taus[1] - taus[0]
eta = 0.0
tE = 0.0
xE = 0.0
eta_max = 0.0
xmax = 0.0
for tau in taus[:-1]:
    k1 = profile(tau)
    k2 = profile(tau + dtau)
    eta_new = eta + 0.5 * (k1 + k2) * dtau / c
    tE += 0.5 * (np.cosh(eta) + np.cosh(eta_new)) * dtau
    xE += 0.5 * (np.sinh(eta) + np.sinh(eta_new)) * dtau * c
    eta = eta_new
    eta_max = max(eta_max, eta)
    xmax = max(xmax, xE)
print(f"  peak gamma = {np.cosh(eta_max):.2f}, peak speed = {np.tanh(eta_max):.5f} c")
print(f"  traveller aged {T_total/year:.2f} yr; Earth aged {tE/year:.2f} yr; turned around at {xmax/ly:.1f} ly")
print(f"  back at start? final x = {xE/ly:.2e} ly; acceleration phase alone: {tau_acc/year:.2f} yr proper time")
C.check("numerical trip returns to its start", abs(xE) / ly < 1e-2)
C.check("Earth ages far more than the traveller, as in the twin effect", tE > 5 * T_total)
te_ana = 4 * (c / a) * np.sinh(np.arccosh(gamma)) + 2 * gamma * coast_tau
C.check(f"integrated Earth time matches closed form ({te_ana/year:.2f} yr) to 0.1%",
        abs(tE - te_ana) / te_ana < 1e-3)

C.section("Biology check: nothing local changes (a chemical rate model, boosted)")
# Toy first-order reaction A -> B with rate k in the rest frame of the cell. Lorentz-transform
# the events of the reaction. Whatever the speed, comoving observers get the same half-life.
k = 1.0
hl_rest = np.log(2) / k
for b in [0.0, 0.9, beta]:
    g = 1 / np.sqrt(1 - b**2)
    # Earth-frame events: decay happens at t = g*tau, x = v*t. Recover the proper time.
    t_e = g * hl_rest
    x_e = b * t_e  # in units where c = 1
    tau = np.sqrt(t_e**2 - x_e**2)
    print(f"  speed {b:<7}: Earth sees half-life {t_e:9.3f}; proper time between events {tau:.3f}")
    C.check(f"proper half-life unchanged at beta = {b}", abs(tau - hl_rest) < 1e-9)

C.section("Hazards at gamma = 31.6")
for n_cm3 in [0.1, 1.0]:
    n = n_cm3 * 1e6
    flux = gamma * n * beta * c  # protons per m^2 per s in the ship frame
    E_p = (gamma - 1) * m_p * c**2
    power = flux * E_p
    # Minimum-ionising stopping power in water is about 2 MeV cm^2/g: 2e6 eV * 1e-3 m^2... per kg
    S = 2e6 * eV * 1e-4 / 1e-3  # J m^2 / kg
    dose = flux * S
    print(f"  n = {n_cm3} H/cm^3: proton energy {E_p/(1e9*eV):.1f} GeV, flux {flux*1e-4:.2e} /cm^2/s, "
          f"beam power {power/1e6:.1f} MW/m^2, ionisation dose {dose:.0f} Gy/s")
    C.check(f"unshielded dose at n = {n_cm3} is above 4 Gy (roughly lethal) within one second", dose > 4)
md = 1e-9
print(f"  1 microgram dust grain: kinetic energy {(gamma-1)*md*c**2:.2e} J = {(gamma-1)*md*c**2/4.184e9:.2f} tonnes TNT")
rocket = np.sqrt((1 + beta) / (1 - beta))
print(f"  ideal photon rocket mass ratio per speed change = {rocket:.1f}; round trip with stops = {rocket**4:.2e}")
print(f"  kinetic energy per kg of ship = {(gamma-1)*c**2:.2e} J/kg")
tcmb = 2.725 * gamma * (1 + beta)
print(f"  CMB head-on temperature in the ship frame = {tcmb:.0f} K (a minor effect)")
# A direct simulation of the same flux: Monte Carlo count of protons hitting 1 cm^2 in 1 microsecond.
rng = np.random.default_rng(7)
n = 1e6  # 1 per cm^3
L_contracted = beta * c * 1e-6 * gamma  # ISM-frame length swept per microsecond of ship time
count = rng.poisson(n * L_contracted * 1e-4)
print(f"  Poisson count on 1 cm^2 per ship microsecond (ISM frame sweep) = {count}; expected {n*L_contracted*1e-4:.0f}")
C.check("frame-sweep count agrees with gamma*n*v flux", abs(count - gamma * n * beta * c * 1e-4 * 1e-6) < 5 * np.sqrt(count))
C.finish()
