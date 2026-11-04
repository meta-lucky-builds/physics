"""C31-C34: fluctuation estimates and a direct test of P ~ exp(-Delta S) in a toy system.

The numbers below are order-of-magnitude. The toy model tests the statistical logic, not
whether de Sitter space actually produces such fluctuations (that part is disputed).
"""
import numpy as np
from math import lgamma, log
from _harness import Checks
from common.constants import G, c, hbar, k_B, year, N_A

C = Checks("Boltzmann brain estimates")
Mpc = 3.0856775814913673e22
H0 = 67.4e3 / Mpc
OL = 0.685
H_inf = H0 * np.sqrt(OL)

C.section("Late de Sitter phase")
L = c / H_inf
lp2 = G * hbar / c**3
S_dS = np.pi * L**2 / lp2  # entropy in units of k_B = horizon area / 4 l_p^2
T_dS = hbar * H_inf / (2 * np.pi * k_B)
print(f"  H_inf = {H_inf:.3e} /s; horizon radius = {L:.3e} m; T_dS = {T_dS:.2e} K")
print(f"  de Sitter entropy S_dS/k = pi L^2/l_P^2 = {S_dS:.2e}")
C.check("S_dS is about 3e122 for the measured Lambda", 2e122 < S_dS < 4e122)

C.section("Brain estimates (inputs stated, order of magnitude only)")
M = 1.4  # kg, adult brain, general knowledge
Mw = 0.018
n_mol = M / Mw
S_molar_water = 70.0  # J/mol/K, standard molar entropy of liquid water (general knowledge)
S_brain = n_mol * S_molar_water / k_B
print(f"  thermal entropy of a 1.4 kg brain as water: {S_brain:.2e} k")
E = M * c**2
dS_horizon = E / (k_B * T_dS)
print(f"  energy {E:.2e} J at T_dS removes {dS_horizon:.2e} k from the horizon")
print(f"  P(one brain in the bath per trial) ~ exp(-{dS_horizon:.1e}): log10 P = {-dS_horizon/np.log(10):.2e}")
C.check("horizon cost dwarfs the brain's own entropy", dS_horizon > 1e30 * S_brain)
Natoms = n_mol * N_A * 3
print(f"  atoms in brain ~ {Natoms:.1e}")
# Alternative lower estimate: only the information in the brain's arrangement.
S_low = 1e16 * np.log(2)  # one bit per ~1e15-1e16 synapse-scale state (rough)
print(f"  low-end 'information only' cost, 1e16 bits: exp(-{S_low:.1e}) (rough, shows the range)")

C.section("Comparison with recurrence time")
log10_t_rec = S_dS / np.log(10)
log10_t_bb = dS_horizon / np.log(10)
print(f"  recurrence time ~ exp(S_dS): log10 t_rec ~ {log10_t_rec:.2e}")
print(f"  brain waiting time per horizon volume: log10 t_BB ~ {log10_t_bb:.2e}")
print(f"  recurrence/brain exponent ratio = {S_dS/dS_horizon:.1e}")
print(f"  brains per recurrence: exp({S_dS:.2e} - {dS_horizon:.2e}) = exp({S_dS:.2e}), which is huge")
C.check("if the phase lasts a recurrence time, freak brains vastly outnumber ordinary observers",
        S_dS > 1e50 * dS_horizon)
# Vacuum decay: if the vacuum decays faster than the brain rate, the problem disappears.
print(f"  brains avoided only if the vacuum decays in time < exp({dS_horizon:.1e}) (any unit), far longer than 13.8 Gyr")

C.section("Direct test of the statistics: Ehrenfest urn")
# N two-state units, one flipped at random each step. Measure how often the system is found in
# a low-entropy macrostate, and compare with exp(S(n) - S_max). Also compare a small fluctuation
# with a larger one: smaller ones are more common, which is the heart of the argument.
N = 40
steps = 3_000_000
rng = np.random.default_rng(42)
flips = rng.integers(0, N, size=steps)
state = rng.integers(0, 2, size=N)
n = int(state.sum())
hist = np.zeros(N + 1, dtype=np.int64)
s = state.tolist()
for f in flips:
    if s[f]:
        s[f] = 0
        n -= 1
    else:
        s[f] = 1
        n += 1
    hist[n] += 1
S = np.array([lgamma(N + 1) - lgamma(k + 1) - lgamma(N - k + 1) for k in range(N + 1)])
pred = np.exp(S - N * log(2))
print(f"  N = {N}, steps = {steps}")
print(f"  {'n up':>5} {'S/k':>7} {'measured P':>12} {'exp(S-N ln2)':>14}")
good = True
for k in [20, 15, 12, 10, 8]:
    frac = (hist[k] + hist[N - k]) / (2 * steps) if k != N - k else hist[k] / steps
    print(f"  {k:5d} {S[k]:7.2f} {frac:12.3e} {pred[k]:14.3e}")
    good &= abs(frac - pred[k]) / pred[k] < 0.25
C.check("measured macrostate frequency matches exp(Delta S) over 3+ decades (within 25%)", bool(good))
fr_small = (hist[12] + hist[28]) / 2
fr_big = (hist[8] + hist[32]) / 2
print(f"  smaller dip (n=12) is {fr_small/fr_big:.0f} times as common as larger dip (n=8); predicted {pred[12]/pred[8]:.0f}")
C.check("smaller fluctuations are far more common than larger ones", fr_small > 50 * fr_big)

C.section("The undefeatable part, as a Bayes check")
# If hypothesis 'BB' assigns our exact present experience the same likelihood as 'ordinary',
# no experience can move the odds.
prior = np.array([1e-3, 0.5, 0.999])
post = prior * 1.0 / (prior * 1.0 + (1 - prior) * 1.0)
print(f"  priors {prior.tolist()} -> posteriors {post.tolist()} when likelihoods are equal")
C.check("equal likelihood means no experiment can shift the odds", np.allclose(prior, post))
C.finish()
