"""C40: neural energy and spike budget, plus a demonstration that a signal is a process.

None of this touches the hard problem of consciousness; it says what thoughts cost and
what a moving pattern is, not why any process feels like something.
"""
import numpy as np
from _harness import Checks
from common.constants import k_B, N_A

C = Checks("Brain numbers and signal-as-process")

C.section("Energy budget")
kcal_day = 516.0  # Herculano-Houzel 2012, for 86 billion neurons
P = kcal_day * 4184 / 86400
N = 86e9
print(f"  brain power = {kcal_day} kcal/day = {P:.1f} W; per neuron = {P/N*1e9:.2f} nW")
C.check("brain power is in the familiar 20-25 W range", 20 <= P <= 26)
G_ATP = 50e3 / N_A  # J per ATP, general-knowledge round value (cellular range about 50-60 kJ/mol)
rate_atp = P / G_ATP
print(f"  ATP turnover at 50 kJ/mol = {rate_atp:.1e} ATP/s; mass of ATP cycled per day = {rate_atp*86400/N_A*0.507:.0f} kg")
spike_atp = 2.4e9  # Lennie 2003
spike_rate = 0.16  # spikes/s/neuron, Lennie 2003 estimate for cortex
E_spike = spike_atp * G_ATP
print(f"  one spike at {spike_atp:.1e} ATP = {E_spike:.1e} J")
print(f"  cortex neurons (19% of 86 billion) = {0.19*N:.1e}")
spikes_cx = 0.19 * N * spike_rate
print(f"  cortical spikes/s at {spike_rate}/s: {spikes_cx:.1e}; power for those spikes = {spikes_cx*E_spike:.2f} W")
C.check("spikes alone use a modest part of the total, consistent with budgets where other costs matter",
        spikes_cx * E_spike < P)
for r in [0.16, 1, 5]:
    print(f"    if all 86e9 neurons fired {r}/s at that cost: {r*N*E_spike:.1f} W")
T = 310.0
landauer = k_B * T * np.log(2)
print(f"  Landauer minimum at 310 K = {landauer:.2e} J/bit; a spike costs {E_spike/landauer:.1e} times that")
C.check("a real spike is far above the thermodynamic floor", E_spike / landauer > 1e6)
atoms = 1.4 / 0.018 * N_A * 3
print(f"  atoms in a 1.4 kg water-like brain = {atoms:.1e}; a thought adds none of its own")

C.section("A pulse travels; the atoms do not: mass-spring chain, Newton's laws only")
# 2000 unit masses joined by unit springs, nearest neighbours. Kick the first few, track the
# position of peak strain and the largest displacement of any mass.
M = 2000
x = np.zeros(M)
v = np.zeros(M)
sigma = 8.0
centre = 60.0
idx = np.arange(M)
# A right-moving pulse for the continuum wave equation: x = f(i - t), v = -f'(i - t).
x = 0.01 * np.exp(-(idx - centre)**2 / (2 * sigma**2))
v = 0.01 * (idx - centre) / sigma**2 * np.exp(-(idx - centre)**2 / (2 * sigma**2))
dt = 0.2


def force(x_):
    f = np.zeros_like(x_)
    f[1:-1] = x_[2:] - 2 * x_[1:-1] + x_[:-2]
    f[0] = x_[1] - x_[0]
    f[-1] = x_[-2] - x_[-1]
    return f


peak0 = np.argmax(np.abs(x))
a = force(x)
n_steps = 5000
disp_max = np.abs(x).max()
start = x.copy()
for k in range(n_steps):
    v_half = v + 0.5 * dt * a
    x = x + dt * v_half
    a = force(x)
    v = v_half + 0.5 * dt * a
    disp_max = max(disp_max, np.abs(x).max())
peak1 = np.argmax(np.abs(x))
speed = (peak1 - peak0) / (n_steps * dt)
print(f"  pulse moved {peak1 - peak0} lattice sites in time {n_steps*dt:.0f}; speed = {speed:.3f} (continuum prediction 1)")
print(f"  largest displacement of any mass in the whole run = {disp_max:.4f} lattice spacings")
print(f"  ratio of distance signal went to distance any mass went = {(peak1-peak0)/disp_max:.1e}")
C.check("pulse speed matches the wave speed to 1%", abs(speed - 1) < 0.01)
C.check("signal crosses about 1000 sites while no mass moves one site", peak1 - peak0 > 900 and disp_max < 0.1)
C.check("number of masses is unchanged: the signal is a pattern, not a new object", x.size == M)
C.finish()
