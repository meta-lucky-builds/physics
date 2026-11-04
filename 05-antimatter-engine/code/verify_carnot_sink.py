"""A04, A05, A06, A11: what a 38 pK cold reservoir really does for a 1500 K engine.

The stepwise engine integrates small reversible cycles against a sink whose temperature
rises. The closed forms in the README are compared with that integration. The kinetic
heat-transfer rate is checked with a Monte Carlo gas of atoms hitting a wall.
"""
import math
import pathlib
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from common.constants import c, e, k_B, m_e, sigma_SB  # noqa: E402

failures = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        failures.append(name)


T_h = 1500.0
T_0 = 38e-12          # Deppner et al. 2021: 3/2 k_B * 38 pK of internal kinetic energy
T_a = 300.0           # assumed ambient
T_cmb = 2.725
N = 1e5               # Deppner et al. 2021: about 1e5 Rb-87 atoms
m_Rb = 86.909 * 1.66054e-27
E_pair = 2 * m_e * c**2

print("1. The Carnot number (A06, A11)")
ratio = T_0 / T_h
print(f"   T_c/T_h = {ratio:.4e}, so eta_C = 1 - {ratio:.3e}")
check("1 - eta_C is 2.5e-14, about 0.025 pico rather than 1 pico", abs(ratio - 2.533e-14) < 1e-17)
for name, T in [("38 pK", T_0), ("CMB 2.725 K", T_cmb), ("300 K ambient", T_a)]:
    print(f"   eta_C with T_c = {name:13s}: 1 - {T / T_h:.3e} = {1 - T / T_h:.6f}")

print("\n2. Size of the reservoir (ideal-gas value, generous for a condensate)")
C = 1.5 * N * k_B
E_th = C * T_0
print(f"   C = 3/2 N k_B = {C:.3e} J/K, thermal content C T0 = {E_th:.3e} J")
print(f"   pair energy 2 m c^2 = {E_pair:.3e} J = {E_pair / E_th:.2e} times the thermal content")
Q_c_ideal = E_pair * ratio
print(f"   Carnot waste per pair at fixed T_c: {Q_c_ideal:.3e} J"
      f" = {Q_c_ideal / E_th:.1f} times the thermal content")
check("the waste heat of one pair exceeds the sink's whole thermal content", Q_c_ideal > E_th)
v_1d = math.sqrt(k_B * T_0 / m_Rb)
print(f"   one-axis rms speed at 38 pK: {v_1d * 1e6:.0f} um/s (the APS Viewpoint quotes about 60 um/s)")
check("38 pK corresponds to the reported 60 um/s expansion", abs(v_1d - 60e-6) < 5e-6)
column = N / (math.pi * (150e-6) ** 2)
sigma_atom = 37 * 0.2865e-28
print(f"   Compton optical depth of the cloud for 511 keV (37 electrons per atom, radius 150 um): "
      f"{column * sigma_atom:.1e}")
check("a 511 keV photon goes straight through the cloud", column * sigma_atom < 1e-10)


def run_engine(Q_total, C_sink, T_start=T_0, T_hot=T_h, steps=100_000):
    """March reversible cycles. Stops when the sink reaches T_hot. Midpoint-rule steps."""
    dQ = Q_total / steps
    T, Q, W, Q_c = T_start, 0.0, 0.0, 0.0
    for _ in range(steps):
        if T >= T_hot:
            break
        T_mid = T * (1 + 0.5 * dQ / (C_sink * T_hot))
        rejected = dQ * T_mid / T_hot
        T += rejected / C_sink
        Q += dQ
        W += dQ - rejected
        Q_c += rejected
    return Q, W, Q_c, T


print("\n3. Engine stepped against the warming sink, N = 1e5 atoms, offered Q_h = 2 m c^2")
Q_used, W, Q_c, T_end = run_engine(E_pair, C)
eta_bar = W / Q_used
print(f"   heat accepted {Q_used:.3e} J ({Q_used / E_pair:.1%} of the pair energy) before T_sink = T_h")
print(f"   work {W:.3e} J, rejected {Q_c:.3e} J, final sink temperature {T_end:.1f} K")
print(f"   mean efficiency over the stroke {eta_bar:.4f}")
Q_limit = C * T_h * math.log(T_h / T_0)
eta_limit = 1 - (1 - T_0 / T_h) / math.log(T_h / T_0)
print(f"   closed form: Q_max = C T_h ln(T_h/T_0) = {Q_limit:.3e} J, mean eta = {eta_limit:.4f}")
check("stepped heat limit matches C T_h ln(T_h/T_0)", abs(Q_used - Q_limit) / Q_limit < 2e-3)
check("stepped mean efficiency matches 1 - 1/ln(T_h/T_0)", abs(eta_bar - eta_limit) < 1e-3)
check("first law holds in the stepped run", abs(W + Q_c - Q_used) / Q_used < 1e-9)
check("one pair cannot be fully absorbed by this sink", Q_limit < E_pair,
      f"{Q_limit / E_pair:.2f} of the pair")
print("   instantaneous eta falls from 1 - 2.5e-14 to 0 during the stroke.")

print("\n   Efficiency against heat handled, same sink")
print("     Q_h [J]      T_final [K]     mean eta      1 - mean eta")
for Q in (1e-22, 1e-18, 1e-16, 1e-15, 1e-14, 5e-14):
    T_f = T_0 * math.exp(Q / (C * T_h))
    loss = C * (T_f - T_0) / Q
    print(f"   {Q:9.1e} {T_f:14.3e} {1 - loss:12.6f} {loss:14.3e}")
Q_small = 1e-18
_, W_small, Q_c_small, _ = run_engine(Q_small, C, steps=20000)
check("small heat loads recover the textbook Carnot loss",
      abs(Q_c_small / Q_small - ratio) / ratio < 0.01,
      f"loss {Q_c_small / Q_small:.3e} at Q_h=1e-18 J")

atoms_double = Q_c_ideal / (1.5 * k_B * T_0)
kg_per_kg = c**2 * ratio / (1.5 * k_B * T_0) * m_Rb
print(f"\n   atoms needed so one pair only doubles T_c: {atoms_double:.2e}")
print(f"   Rb needed at 38 pK to sink the waste of 1 kg of annihilated mass the same way: {kg_per_kg:.2e} kg")

print("\n4. Paying for the cold sink: T_c cancels exactly (A06)")
Qh, Th, Tc, Ta = sp.symbols("Q_h T_h T_c T_a", positive=True)
Qc = Qh * Tc / Th
W_engine = Qh - Qc
W_fridge = Qc * (Ta - Tc) / Tc
net = sp.simplify((W_engine - W_fridge) / Qh)
print(f"   sympy: (W_engine - W_fridge)/Q_h = {net}")
check("net efficiency simplifies to 1 - T_a/T_h with T_c absent",
      sp.simplify(net - (1 - Ta / Th)) == 0 and Tc not in net.free_symbols)
f = sp.symbols("f", positive=True)
f_breakeven = sp.solve(sp.Eq(W_engine, W_fridge / f), f)[0]
print(f"   fridge must exceed this fraction of its Carnot COP to net anything: {sp.simplify(f_breakeven)}")
f_num = float(f_breakeven.subs({Ta: T_a, Th: T_h, Tc: T_0, Qh: 1}))
print(f"   numerically {f_num:.4f} for T_a = 300 K, T_h = 1500 K")
check("break-even fridge quality is 20% of Carnot", abs(f_num - 0.2) < 1e-6)
cop_carnot = T_0 / (T_a - T_0)
print(f"   Carnot COP from 38 pK to 300 K = {cop_carnot:.3e}")

print("\n5. Whole-life ledger with a finite sink: prepare, run, return to ambient")


def exergy(C_sink, T):
    """Work needed to bring the sink from T_a to T, equal to the work recoverable going back."""
    return C_sink * (T_a * math.log(T_a / T) - (T_a - T))


W_prep = exergy(C, T_0)
print(f"   ideal least work to cool the sink from 300 K to 38 pK: {W_prep:.3e} J"
      f" = {W_prep / e / 1e3:.0f} keV")
print("     Q_h [J]    W_engine [J]   W_prep [J]   W_recovered [J]   net/Q_h   1-T_a/T_h")
for Q in (1e-16, 1e-15, 1e-14, 5e-14, 9e-14):
    _, W_e, _, T_f = run_engine(Q, C)
    W_back = exergy(C, T_f)
    net_eff = (W_e - W_prep + W_back) / Q
    print(f"   {Q:9.1e} {W_e:12.4e} {W_prep:12.4e} {W_back:15.4e} {net_eff:11.5f} {1 - T_a / T_h:9.3f}")
    check(f"lifetime net efficiency is 1 - T_a/T_h at Q_h={Q:.0e}", abs(net_eff - (1 - T_a / T_h)) < 1e-3)
real_input = 1e3 * 10    # assumed 1 kW for 10 s per shot; order of magnitude only
print(f"   assumed real apparatus input 1 kW x 10 s = {real_input:.0e} J per shot,"
      f" which is {real_input / W_prep:.0e} times the ideal {W_prep:.1e} J")

print("\n6. Heat-rejection rate into a dilute gas: Monte Carlo against kinetic theory")
r_cloud = 150e-6      # assumed cloud radius, m
T_wall = 2 * T_0      # wall warmer than the gas by one T_0
v_bar = math.sqrt(8 * k_B * T_0 / (math.pi * m_Rb))
sigma_g = math.sqrt(k_B * T_0 / m_Rb)
sigma_w = math.sqrt(k_B * T_wall / m_Rb)
P_kinetic = N * (3 * v_bar / (4 * r_cloud)) * 2 * k_B * (T_wall - T_0)
print(f"   mean speed {v_bar * 1e6:.0f} um/s at 38 pK; cloud radius assumed {r_cloud * 1e6:.0f} um")
print(f"   kinetic theory: P_c = N (3 v/4r) 2 k (T_wall - T_gas) = {P_kinetic:.2e} W")

rng = np.random.default_rng(11)
n_mc = 200_000
pos = rng.normal(size=(n_mc, 3))
pos *= (r_cloud * rng.uniform(size=(n_mc, 1)) ** (1 / 3)) / np.linalg.norm(pos, axis=1, keepdims=True)
vel = rng.normal(0, sigma_g, (n_mc, 3))
ke0 = 0.5 * m_Rb * (vel**2).sum(axis=1).mean()
dt = 0.002 * r_cloud / v_bar
gain = 0.0
hits_total = 0
steps = 40
for _ in range(steps):
    pos += vel * dt
    out = (pos**2).sum(axis=1) > r_cloud**2
    hits = int(out.sum())
    if hits:
        normal = pos[out] / np.linalg.norm(pos[out], axis=1, keepdims=True)
        before = 0.5 * m_Rb * (vel[out] ** 2).sum(axis=1)
        # Diffuse re-emission from a wall at T_wall: Rayleigh normal speed, Gaussian tangents.
        v_n = sigma_w * np.sqrt(-2 * np.log(rng.uniform(size=hits)))
        v_new = rng.normal(0, sigma_w, (hits, 3))
        v_new -= (v_new * normal).sum(axis=1, keepdims=True) * normal
        v_new -= v_n[:, None] * normal
        vel[out] = v_new
        pos[out] = normal * r_cloud * (1 - 1e-9)
        gain += (0.5 * m_Rb * (v_new**2).sum(axis=1) - before).sum()
        hits_total += hits
P_mc = gain / (steps * dt) * (N / n_mc)
print(f"   Monte Carlo ({n_mc} atoms, {hits_total} wall hits): P_c = {P_mc:.2e} W")
check("Monte Carlo heat rate agrees with kinetic theory within 10%",
      abs(P_mc - P_kinetic) / P_kinetic < 0.10, f"ratio {P_mc / P_kinetic:.3f}")
P_w = P_kinetic * (T_h / T_0 - 1)
print(f"   work this allows with an ideal engine: P_c (T_h/T_c - 1) = {P_w:.2e} W")
t_pair = Q_c_ideal / P_kinetic
print(f"   time to dump one pair's Carnot waste: {t_pair:.0f} s; drop-tower free fall is 4.74 s")
check("the waste of one pair takes longer than the drop", t_pair > 4.74)
print(f"   gas total thermal energy / rate = {E_th / P_kinetic:.2f} s, the sink's own turnover time")

print("\n7. Radiative cooling of the 38 pK side")
flux = sigma_SB * T_0**4
wien = 2.897771955e-3 / T_0
area = 4 * math.pi * r_cloud**2
print(f"   blackbody flux at 38 pK: {flux:.2e} W/m^2; over the cloud {flux * area:.2e} W")
print(f"   Wien peak wavelength: {wien:.2e} m, compared with Earth's radius 6.4e6 m")
check("blackbody wavelength at 38 pK is larger than the Earth", wien > 6.4e6)

print(f"\n{'ALL CHECKS PASSED' if not failures else 'FAILED: ' + ', '.join(failures)}")
sys.exit(1 if failures else 0)
