"""L02: rework the historical measurements of c and test the one-way convention.

- Synthetic Roemer eclipse timings, fitted for c.
- Synthetic Fizeau wheel: find the first dark speed by simulating chopped light.
- Number checks for Bradley, Michelson, Evenson.
- Reichenbach convention: rerun a round trip and a Roemer-style observation in
  coordinates with anisotropic one-way speeds, using proper times only, and confirm
  that nothing observable changes.
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from checklib import K, check, finish, section

c = K.c
AU = 1.495978707e11  # m, IAU 2012 exact

section("Roemer: Io eclipse timing")
P_io = 42.459 * 3600  # s, Io period (general knowledge)
print(f"  Roemer's 22 min across the orbit diameter: c = {2 * AU / (22 * 60) / 1e3:,.0f} km/s "
      "using the modern AU; light really takes "
      f"{2 * AU / c / 60:.1f} min")
v_earth = 2 * np.pi * AU / K.year
print(f"  One Io period shifts by at most v/c * P = {v_earth / c * P_io:.1f} s; "
      "the effect accumulates over months, which is why it could be seen")


def roemer_fit(c_true, seed=2, noise=30.0):
    rng = np.random.default_rng(seed)
    n = np.arange(0, 200)
    t = n * P_io  # emission times at Jupiter
    # Earth on a circular orbit, Jupiter at 5.2 AU on a slower orbit; solve light delay
    def dist(tt):
        e = AU * np.array([np.cos(2 * np.pi * tt / K.year), np.sin(2 * np.pi * tt / K.year)])
        wj = 2 * np.pi / (11.86 * K.year)
        j = 5.2 * AU * np.array([np.cos(wj * tt), np.sin(wj * tt)])
        return np.linalg.norm(e - j, axis=0)
    arrive = t.copy()
    for _ in range(6):
        arrive = t + dist(arrive) / c_true
    arrive = arrive + rng.normal(0, noise, n.size)
    # fit arrival = t0 + P n + D/c: unknowns t0, P, 1/c
    D = dist(arrive)
    A = np.column_stack([np.ones(n.size), n, D])
    sol, *_ = np.linalg.lstsq(A, arrive, rcond=None)
    return 1 / sol[2], sol[1]


for c_true in [c, 2.2e8]:
    c_fit, P_fit = roemer_fit(c_true)
    check(f"recovers an injected c = {c_true:.3e} m/s from eclipse timings with 30 s noise",
          abs(c_fit / c_true - 1) < 0.01, f"fit {c_fit:.4e} m/s, period {P_fit / 3600:.3f} h")

section("Bradley: aberration")
alpha = np.radians(20.5 / 3600)
c_b = 2 * np.pi * AU / K.year / np.tan(alpha)
print(f"  20.5 arcsec with v = {v_earth / 1e3:.2f} km/s gives c = {c_b / 1e3:,.0f} km/s "
      f"({(c_b / c - 1) * 100:+.2f}% from modern)")
check("Bradley is within 0.1% given the modern Earth speed", abs(c_b / c - 1) < 1e-3)

section("Fizeau: toothed wheel")
d, N_t = 8633.0, 720


def wheel_open(angle_phase):
    """True if a light ray at this wheel phase (in tooth-periods) goes through a gap."""
    return (angle_phase % 1.0) < 0.5


def throughput(f, c_true, n=20000, seed=3):
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 1.0, n)
    phase_out = N_t * f * t
    phase_back = N_t * f * (t + 2 * d / c_true)
    return np.mean(wheel_open(phase_out) & wheel_open(phase_back)) / 0.5


fs = np.linspace(0.5, 30, 1181)
for c_true in [c, 313.3e6]:
    T = np.array([throughput(f, c_true) for f in fs])
    f_dark = fs[T.argmin()]
    c_rec = 4 * d * N_t * f_dark
    check(f"simulated wheel goes dark where c = 4 d N f gives back c_true = {c_true:.4e}",
          abs(c_rec / c_true - 1) < 2e-3, f"first dark at {f_dark:.2f} rev/s, min throughput {T.min():.3f}")
c_f = 4 * d * N_t * 12.6
print(f"  Fizeau's own 12.6 rev/s: c = {c_f / 1e3:,.0f} km/s ({(c_f / c - 1) * 100:+.1f}% from modern)")
print(f"  Sensitivity: true c would be dark at {c / (4 * d * N_t):.2f} rev/s, so a 0.5 rev/s reading "
      "error matches the 4.5% offset")
check("Fizeau's number reproduces about 313,000 km/s", abs(c_f / 1e3 - 313274) < 1)

section("Foucault and Michelson: rotating mirrors")
D_fo, f_fo = 20.0, 400.0  # m, rev/s: illustrative round numbers
th = 2 * (2 * np.pi * f_fo) * (2 * D_fo / c)
print(f"  Foucault-style setup, {D_fo:.0f} m path, {f_fo:.0f} rev/s: beam swing {th * 1e3:.2f} mrad "
      "(mirror turns while light is away; a reflected ray turns by twice the mirror angle)")
D_m = 35.4e3
f_m = c / (8 * 2 * D_m)
print(f"  Michelson-style octagon, {D_m / 1e3:.1f} km baseline: returns on the next face at "
      f"{f_m:.0f} rev/s, so a 1 km/s error in c is a {1e3 / c:.1e} fractional rate error")
check("octagon null rate is reachable mechanically (hundreds of rev/s)", 400 < f_m < 700, f"{f_m:.0f} rev/s")

section("Evenson 1972: laser frequency times wavelength")
f_l, lam = 88.376181627e12, 3.392231376e-6
c_e = f_l * lam
check("frequency x wavelength reproduces the reported value to within its 1.1 m/s error",
      abs(c_e - 299792456.2) < 1.1, f"{f_l:.9e} Hz x {lam:.9e} m = {c_e:,.1f} m/s")
print(f"  Difference from today's exact value: {c - c_e:.1f} m/s")

section("Reichenbach convention: nothing observable changes")
# Coordinates t' = t + kappa x / c. Light moving +x has speed c/(1+kappa),
# moving -x has speed c/(1-kappa). ds^2 = -c^2 (dt' - kappa dx/c)^2 + dx^2.
L = 1.0e3
for kappa in [0.0, 0.5, -0.7, 0.9]:
    cp, cm = c / (1 + kappa), c / (1 - kappa)
    rt = L / cp + L / cm
    check(f"kappa = {kappa:+.1f}: one-way speeds {cp / c:.3f}c and {cm / c:.3f}c, round trip unchanged",
          abs(rt / (2 * L / c) - 1) < 1e-14)

A = 1.5e11  # m
Omega = 2 * np.pi / K.year
X = 8.0e11  # source position, m
P = 5 * 86400.0


def x_obs(tau):
    return A * np.sin(Omega * tau)


def xdot(tau):
    return A * Omega * np.cos(Omega * tau)


def arrival_taus(kappa, ticks):
    """Observer proper times at which source ticks arrive, in kappa coordinates.

    The observer's motion is specified against their own clock. Their coordinate time
    comes from integrating the metric, so kappa enters through dt' = kappa dx/c +
    sqrt(1 + xdot^2/c^2) dtau. The source sits still, so dt' = dtau there, with its
    t' offset set by the same coordinate definition at x = X."""
    cm = c / (1 - kappa)

    def t_obs(tau):
        I, _ = quad(lambda s: np.sqrt(1 + (xdot(s) / c) ** 2), 0, tau, limit=400, epsabs=1e-5)
        return kappa * x_obs(tau) / c + I

    out = []
    for n in ticks:
        t_emit = n * P + kappa * X / c
        f = lambda tau: t_obs(tau) - t_emit - (X - x_obs(tau)) / cm
        out.append(brentq(f, n * P - 1e5, n * P + 1e5 * 50, xtol=1e-5))
    return np.array(out)


ticks = np.arange(0, 74, 2)
base = arrival_taus(0.0, ticks)
Dd = X - x_obs(base)
# Fit in well-scaled columns; a raw fit with D in metres is badly conditioned.
Dbar, Dscale = Dd.mean(), Dd.std()
Am = np.column_stack([np.ones(ticks.size), ticks - ticks.mean(), (Dd - Dbar) / Dscale])
sol, *_ = np.linalg.lstsq(Am, base, rcond=None)
c_fit = Dscale / sol[2]
check("Roemer-style fit recovers c in standard coordinates", abs(c_fit / c - 1) < 1e-3,
      f"fit {c_fit:.5e} m/s")
for kappa in [0.5, -0.7, 0.9]:
    other = arrival_taus(kappa, ticks)
    worst = np.abs(other - base).max()
    check(f"Roemer arrival times at kappa = {kappa:+.1f} match kappa = 0", worst < 1e-3, f"worst difference {worst:.1e} s")
print("  So a Roemer-style observation also yields the two-way-equivalent c. The moving observer's")
print("  clock term cancels the anisotropy, the same way slow clock transport does.")

finish()
