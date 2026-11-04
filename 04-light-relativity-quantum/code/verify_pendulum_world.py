"""L14, L15: a world keeping time only with cesium pendulums, and a pendulum in an ultracold bath.

- Integrate the pendulum equation to measure rate against temperature and amplitude.
- A toy world: only ratios of clocks are observable; a common drift is invisible;
  independent processes expose it.
- Debye model and Grueneisen relation, tested against copper's published contraction.
- Energy budget for cooling a metal bob with an ultracold gas.
"""
import numpy as np
from scipy.integrate import quad, solve_ivp

from checklib import K, check, finish, section

rng = np.random.default_rng(23)
g0 = K.g_0
day = 86400.0


def period_sim(L, theta0, g=g0, n=40):
    """Period from direct integration (mean over n cycles, located by zero crossings)."""
    w = np.sqrt(g / L)

    def rhs(t, y):
        return [y[1], -w**2 * np.sin(y[0])]

    def ev(t, y):
        return y[0]

    ev.direction = 1
    T0 = 2 * np.pi / w
    sol = solve_ivp(rhs, [0, T0 * (n + 2)], [theta0, 0.0], events=ev, method="DOP853",
                    rtol=1e-12, atol=1e-13)
    te = sol.t_events[0]
    return (te[-1] - te[0]) / (len(te) - 1)


section("Pendulum rate")
L0 = g0 * (2 / (2 * np.pi)) ** 2  # seconds pendulum, period 2 s
print(f"  Length for a 2 s period: L = {L0:.4f} m")
T_small = period_sim(L0, 1e-3)
check("small-swing period = 2 pi sqrt(L/g)", abs(T_small / 2.0 - 1) < 1e-6, f"{T_small:.8f} s")
alpha_cs = 97e-6
print("  Metal rod: temperature coefficient, seconds/day per K")
for name, alpha in [("cesium", 97e-6), ("brass", 19e-6), ("steel", 12e-6), ("invar", 1.2e-6)]:
    print(f"    {name:7s} alpha = {alpha * 1e6:5.1f} ppm/K -> {alpha / 2 * day:.3f} s/day/K")
dTs = np.array([-20.0, -10.0, 0.0, 10.0, 20.0])
Ps = np.array([period_sim(L0 * (1 + alpha_cs * dT), 1e-3) for dT in dTs])
slope = np.polyfit(dTs, Ps / 2.0, 1)[0]
check("simulated dT/T per kelvin = alpha/2", abs(slope / (alpha_cs / 2) - 1) < 2e-3, f"{slope * 1e6:.2f} ppm/K against {alpha_cs / 2 * 1e6:.2f} ppm/K")
print(f"  So a cesium pendulum changes rate by {slope * day:.2f} s/day per kelvin; a 10 K swing is {slope * day * 10:.0f} s/day")
print(f"  A compound (uniform rod) pendulum has T = 2 pi sqrt(2L/(3g)): same scale law, same alpha/2")
T5 = period_sim(L0, np.radians(5))
amp = T5 / T_small - 1
check("amplitude matters: T grows by theta0^2/16", abs(amp / (np.radians(5) ** 2 / 16) - 1) < 2e-3, f"{amp * day:.1f} s/day at 5 degrees")
g_eq, g_pole = 9.7803, 9.8322  # m/s^2 (general knowledge)
gdiff = 0.5 * (g_pole - g_eq) / g0
print(f"  Moving between equator and pole: {gdiff * day:.0f} s/day, equal to {gdiff / (alpha_cs / 2):.0f} K of cesium")

section("Toy world: only ratios can be measured")
n_days = 2000
t_days = np.arange(n_days)
T_common = 15 + 10 * np.sin(2 * np.pi * t_days / 365.25) + rng.normal(0, 3, n_days)
rate = lambda dT: 1 - alpha_cs / 2 * dT  # pendulum ticks per unit true time
clocks = np.array([np.cumsum(rate(T_common - 15)) for _ in range(5)])
ratio_spread = np.ptp(clocks, axis=0).max()
check("identical clocks at the same temperature never disagree, despite drift", ratio_spread == 0,
      f"yet each drifts {np.ptp(clocks[0] - t_days - 1) * day:.0f} s against model time")
# Two rooms with different temperature: pendulums disagree
room = np.cumsum(rate(T_common - 15 + 8))
diff = (room - clocks[0]) * day
check("a room 8 K warmer shows up immediately", abs(diff[-1] / n_days + alpha_cs / 2 * 8 * day) < 1e-6, f"{diff[-1] / n_days:.2f} s/day")

# Independent processes that keep their own time in the model: Earth spin, Moon, Mercury.
periods = {"Earth spin": 1.0, "Moon orbit": 27.32, "Mercury orbit": 87.97}
slopes = {}
print("  Apparent rate in pendulum time against daily temperature, 3 ppm measurement noise on each process:")
for name, Pp in periods.items():
    # apparent rate of process, per pendulum day
    app = (1 / Pp) * 1 / rate(T_common - 15)
    meas = app * (1 + rng.normal(0, 3e-6, n_days))
    s = np.polyfit(T_common - 15, meas * Pp, 1)[0]
    slopes[name] = s
    print(f"    {name:14s} fractional slope {s * 1e6:7.3f} ppm/K")
vals = np.array(list(slopes.values()))
check("all processes appear to speed up on hot days by the same fraction, alpha/2",
      np.all(np.abs(vals / (alpha_cs / 2) - 1) < 0.02), f"alpha/2 = {alpha_cs / 2 * 1e6:.2f} ppm/K")
print("  Same fractional slope in unrelated processes is the signature: blame the clock, not the Sun, Moon and Mercury")

section("Debye model and Grueneisen relation")
R = K.R


def debye_C(T, th):
    """Lattice heat capacity per mole."""
    x = th / T
    if x > 200:
        return 12 * np.pi**4 / 5 * R * (T / th) ** 3
    I = quad(lambda s: s**4 * np.exp(s) / np.expm1(s) ** 2, 0, x, limit=200)[0]
    return 9 * R * (T / th) ** 3 * I


def debye_U(T, th):
    return quad(lambda t: debye_C(t, th), 1e-3, T, limit=200, points=[min(th / 10, T)])[0]


th_cu = 343.0  # K, Debye temperature of copper (Kittel, not re-verified online)
gam_cu = 1.96  # Grueneisen parameter, copper
B_cu = 140e9  # Pa (general knowledge)
M_cu = 63.546e-3
rho_cu = 8960.0
Vm_cu = M_cu / rho_cu
ge_cu = 0.695e-3  # J/mol/K^2 electronic coefficient (general knowledge)
check("Debye C at 1000 K approaches Dulong-Petit 3R", abs(debye_C(1000, th_cu) / (3 * R) - 1) < 0.01)
check("Debye C at 5 K follows the T^3 law", abs(debye_C(5, th_cu) / (12 * np.pi**4 / 5 * R * (5 / th_cu) ** 3) - 1) < 1e-3)


def alpha_cu(T):
    return gam_cu * (debye_C(T, th_cu) + ge_cu * T) / (3 * B_cu * Vm_cu)


contract, _ = quad(alpha_cu, 4, 293, limit=300, points=[30, 100])
pub = 0.324e-2
check("Grueneisen integral reproduces copper's 293 K to 4 K contraction within 10%", abs(contract / pub - 1) < 0.10,
      f"model {contract * 100:.3f}%, published {pub * 100:.3f}%")
print(f"  alpha(293 K) model = {alpha_cu(293) * 1e6:.1f} ppm/K (measured about 16.5 ppm/K, general knowledge)")
print("  Why cesium expands so much: alpha = gamma C/(3 B V). Bulk modulus is the main factor")
Cs_B, Cs_V = 1.6e9, 132.905e-3 / 1930.0
implied = 97e-6 * 3 * Cs_B * Cs_V / (3 * R)
print(f"    Cu: 3BV = {3 * B_cu * Vm_cu:.2e} J/mol; Cs: 3BV = {3 * Cs_B * Cs_V:.2e} J/mol, {B_cu * Vm_cu / (Cs_B * Cs_V):.1f}x softer")
print(f"    Cs with alpha = 97 ppm/K and C = 3R implies a Grueneisen parameter of {implied:.2f}, an ordinary value")
check("cesium's implied Grueneisen parameter is ordinary (0.5 to 3)", 0.5 < implied < 3)

section("L15: cold pendulum")
print("  T (K)    alpha_L (1/K)   rate change per K (s/day/K)   contraction left between T and 0")
Ts = [293, 77, 4, 1, 1e-3, 1e-6, 1e-9]
alpha_vals = []
for T in Ts:
    a = alpha_cu(T)
    alpha_vals.append(a)
    left = quad(alpha_cu, 1e-12, T, limit=200, points=[min(30, T)])[0] if T > 1 else gam_cu * (ge_cu * T**2 / 2) / (3 * B_cu * Vm_cu)
    print(f"  {T:<7g}  {a:.3e}       {a / 2 * day:.3e}                 {left:.3e}")
check("alpha falls toward zero (third law) from 293 K to 1 nK", all(np.diff(alpha_vals) < 0) and alpha_vals[-1] < 1e-18, f"alpha(1 nK) = {alpha_vals[-1]:.1e} /K")
L_rate = (contract + quad(alpha_cu, 1e-6, 4, limit=200)[0]) / 2 * day
print(f"  Total rate gain for a copper pendulum 293 K -> 0 K: {L_rate:.0f} s/day; almost all of it occurs above 30 K")
frac_above = quad(alpha_cu, 30, 293, limit=300)[0] / contract
check("almost all contraction occurs above 30 K", frac_above > 0.99, f"{frac_above * 100:.2f}%")
# "Classical only" counterfactual
alpha_classical = gam_cu * 3 * R / (3 * B_cu * Vm_cu)
print(f"  Classical (Dulong-Petit down to 0) would keep alpha = {alpha_classical:.2e} /K at all T,")
print(f"  which would give a contraction of {alpha_classical * 293 * 100:.2f}% and no freeze-out: the freeze-out IS the quantum effect")

section("L15: energy budget for a 100 g copper bob")
m_bob = 0.100
n_bob = m_bob / M_cu
r_bob = (3 * m_bob / (4 * np.pi * rho_cu)) ** (1 / 3)
A_bob = 4 * np.pi * r_bob**2
U_bob = n_bob * (debye_U(293, th_cu) + ge_cu * 293**2 / 2)
C_bob = n_bob * (debye_C(293, th_cu) + ge_cu * 293)
print(f"  Bob radius {r_bob * 100:.2f} cm, area {A_bob * 1e4:.0f} cm^2, heat capacity {C_bob:.1f} J/K, heat to remove {U_bob:.2e} J")
N_at = 1e6
m_Rb = 87 * 1.66054e-27
T_gas = 100e-9
n_gas = 1e20  # atoms/m^3, 1e14 cm^-3 peak density typical of a BEC
print(f"  Ultracold cloud: {N_at:.0e} Rb atoms, total mass {N_at * m_Rb:.1e} kg, density {n_gas * m_Rb:.1e} kg/m^3 (air 1.2 kg/m^3)")
U_cloud = 1.5 * N_at * K.k_B * T_gas
E_photon = K.k_B * 293
print(f"  Thermal energy of the whole cloud at {T_gas * 1e9:.0f} nK: {U_cloud:.1e} J; one room-temperature thermal photon: {E_photon:.1e} J")
check("one 293 K thermal photon carries more energy than the entire cold cloud", E_photon > U_cloud)
E_sink = 1.5 * N_at * K.k_B * 293
dT_bob = E_sink / C_bob
print(f"  If every atom warms to 293 K it absorbs {E_sink:.1e} J, cooling the bob by {dT_bob:.1e} K")
check("cloud's maximum heat sink is 1e17+ times too small", U_bob / E_sink > 1e17, f"ratio {U_bob / E_sink:.1e}")
E_ads = 0.5 * K.eV * N_at
print(f"  Atoms that stick release roughly 0.5 eV each (order of magnitude): {E_ads:.1e} J, more than they could remove")
check("adsorption heat exceeds the most the atoms could carry off", E_ads > E_sink)
C_cloud = 1.5 * N_at * K.k_B  # classical gas value; a condensate has less
for T in [1e-3, 1e-6, 1e-9]:
    C_b = n_bob * (debye_C(T, th_cu) + ge_cu * T)
    print(f"  Bob heat capacity at {T:g} K: {C_b:.2e} J/K; cloud {C_cloud:.1e} J/K at most")
check("even at 1 nK the bob's heat capacity exceeds the cloud's", n_bob * ge_cu * 1e-9 > C_cloud)
v_th = np.sqrt(8 * K.k_B * T_gas / (np.pi * m_Rb))
flux = n_gas * v_th / 4 * A_bob
print(f"  Collision rate on the bob at that density: {flux:.1e} atoms/s, cloud lost in {N_at / flux * 1e9:.1f} ns if it enveloped the bob")
rate_load = 1e10  # atoms/s fed by a strong laser-cooling source (order of magnitude)
P_cool = rate_load * 1.5 * K.k_B * 293
eps = 0.02  # polished copper emissivity (general knowledge)
P_rad = eps * K.sigma_SB * (293**4) * A_bob
k_wire, A_wire, L_wire = 15.0, np.pi * (0.1e-3) ** 2, 1.0
P_wire = k_wire * A_wire * 293 / L_wire
print(f"  Sustained cooling with fresh atoms at {rate_load:.0e}/s: at most {P_cool:.1e} W")
print(f"  Heat coming in: radiation from 293 K walls {P_rad:.2e} W, 0.2 mm steel wire to room temperature {P_wire:.2e} W")
check("heat load beats cold-atom cooling by more than 1e8", (P_rad + P_wire) / P_cool > 1e8, f"ratio {(P_rad + P_wire) / P_cool:.1e}")

section("L15: when does the swing itself stop being classical?")
omega = np.sqrt(g0 / 1.0)
for T in [293, 1e-3, 1e-9, 38e-12]:
    nbar = 1 / np.expm1(K.hbar * omega / (K.k_B * T))
    xth = np.sqrt(K.hbar / (2 * m_bob * omega) * (2 * nbar + 1))
    print(f"  T = {T:<8g} K: mean quanta in the swing mode {nbar:.3g}, rms position {xth:.3g} m")
T_q = K.hbar * omega / K.k_B
print(f"  hbar omega / k_B = {T_q:.2e} K. Below this, the swing mode is in its quantum ground state")
print(f"  Zero-point rms = {np.sqrt(K.hbar / (2 * m_bob * omega)):.2e} m")
check("the swing mode nears its ground state only near 20 pK", 1e-11 < T_q < 5e-11)

finish()
