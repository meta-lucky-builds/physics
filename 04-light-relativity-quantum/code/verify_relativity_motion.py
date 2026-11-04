"""L11: looking in a mirror while moving close to c.

- Simulate the light and mirror in the outside frame, then Lorentz transform the
  events to the traveller's frame.
- Symbolic velocity addition, exact-rational stacking of boosts.
- Constant force integrated with p = gamma m v.
"""
from fractions import Fraction

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from checklib import K, check, finish, section

c = K.c

section("Mirror carried with you: events in the outside frame, transformed to yours")


def lorentz(t, x, beta):
    g = 1 / np.sqrt((1 - beta) * (1 + beta))
    return g * (t - beta * x / c), g * (x - beta * c * t)


d = 1.0  # rest distance from face to mirror, m
for beta in [0.0, 0.5, 0.9, 0.99, 0.9999, 0.999999]:
    g = 1 / np.sqrt((1 - beta) * (1 + beta))
    v = beta * c
    mirror = lambda t: d / g + v * t  # length contracted separation, moving mirror
    face = lambda t: v * t
    # Numerically locate where the photon meets each moving object
    t_max = 4 * g * d / c + 1e-6
    t1 = brentq(lambda t: c * t - mirror(t), 0, t_max, xtol=1e-22, rtol=1e-14)
    x1 = c * t1
    t2 = brentq(lambda t: x1 - c * (t - t1) - face(t), t1, t_max, xtol=1e-22, rtol=1e-14)
    x2 = x1 - c * (t2 - t1)
    tp1, xp1 = lorentz(t1, x1, beta)
    tp2, xp2 = lorentz(t2, x2, beta)
    ok = abs(xp1 - d) < 1e-6 and abs(tp2 - 2 * d / c) < 1e-6 * (2 * d / c) and abs(tp1 - d / c) < 1e-6 * d / c
    check(f"beta = {beta}: you see your image after 2d/c, mirror at d", ok,
          f"your frame {tp2 * 1e9:.4f} ns; outside frame: out {t1 * 1e9:.3g} ns, back {t2 * 1e9 - t1 * 1e9:.3g} ns, total {t2 * 1e9:.3g} ns")


section("Velocity addition from the Lorentz transformation")
u, v_, cc = sp.symbols("u v c", positive=True)
g = 1 / sp.sqrt(1 - v_**2 / cc**2)
s = sp.symbols("s", positive=True)
x_lab = u * s
t_lab = s
xp = g * (x_lab - v_ * t_lab)
tp = g * (t_lab - v_ * x_lab / cc**2)
vel = sp.simplify(xp / tp)
check("u' = (u - v)/(1 - uv/c^2)", sp.simplify(vel - (u - v_) / (1 - u * v_ / cc**2)) == 0)
check("u = c gives u' = c for any v", sp.simplify(vel.subs(u, cc) - cc) == 0)
w = Fraction(0)
for n in range(200):
    w = (w + Fraction(1, 2)) / (1 + w / 2)
check("200 successive boosts of c/2 stay below c (exact arithmetic)", w < 1, f"1 - w = {float(1 - w):.2e}")

section("Force cannot bring a mass to c")
m = 70.0
F = 70.0 * 9.80665  # 1 g
year = K.year


def rhs(t, y):
    p = y[0]
    return [F]


sol = solve_ivp(rhs, [0, 20 * year], [0.0], t_eval=np.array([0.5, 1, 2, 5, 10, 20]) * year, rtol=1e-10)
print("  t (yr)   Newton v/c   true v/c")
for t, p in zip(sol.t, sol.y[0]):
    b = p / np.sqrt((m * c) ** 2 + p**2)
    print(f"  {t / year:5.1f}   {F * t / m / c:9.3f}   {b:.10f}")
check("v stays below c for all force and time tried", np.all(sol.y[0] / np.sqrt((m * c) ** 2 + sol.y[0] ** 2) < 1))
world = 6e20  # J, roughly world annual primary energy (general knowledge, order only)
print("  Energy for a 70 kg person, (gamma-1) m c^2:")
for dev in [1e-1, 1e-2, 1e-6, 1e-12, 1e-24]:
    gam = 1 / np.sqrt(dev * (2 - dev))
    KE = (gam - 1) * m * c**2
    print(f"    1 - v/c = {dev:g}: gamma = {gam:.3e}, energy = {KE:.3e} J = {KE / world:.2e} years of world energy use")

print("  Bertozzi-type check, electrons in MeV kinetic energy:")
mec2 = K.m_e * c**2 / K.eV / 1e6
for T in [0.5, 1.5, 4.5, 15]:
    gam = 1 + T / mec2
    print(f"    {T:5.1f} MeV: v/c = {np.sqrt(1 - 1 / gam**2):.4f}; Newton would say {np.sqrt(2 * T / mec2):.2f}")

section("What an outside observer sees (photon four-momentum)")


def boost(E, p, beta):
    gam = 1 / np.sqrt((1 - beta) * (1 + beta))
    return gam * (E + beta * p), gam * (p + beta * E)  # from traveller frame to outside (c = 1)


for beta in [0.5, 0.99]:
    E_out, _ = boost(1.0, 1.0, beta)  # forward photon
    E_back, _ = boost(1.0, -1.0, beta)  # reflected photon returning
    check(f"beta = {beta}: forward photon blue, returning photon red by the Doppler factor",
          abs(E_out * E_back - 1) < 1e-12 and abs(E_out - np.sqrt((1 + beta) / (1 - beta))) < 1e-10,
          f"outside energies {E_out:.3f} and {E_back:.4f} of your 1.000")

section("Mirror fixed to the ground, you running at it")
for beta in [0.5, 0.99, 0.9999]:
    D = 1000.0
    v = beta * c
    t_m = D / c  # the mirror is reached by light before you get there: at t = D/c you are at beta D
    # Your face at t = D/c is at x = beta D < D, so the photon arrives first.
    t_see = 2 * D / (c * (1 + beta))  # returning photon meets you: D - c (t - D/c) = v t
    t_numeric = brentq(lambda t: D - c * (t - D / c) - v * t, D / c * 0.5, 2 * D / c)
    gam = 1 / np.sqrt((1 - beta) * (1 + beta))
    shift = (1 + beta) / (1 - beta)
    check(f"beta = {beta}: photon still leads you and comes back", abs(t_numeric / t_see - 1) < 1e-10,
          f"you see the image at {t_numeric * 1e6:.3f} us (lab) when you are {D - v * t_numeric:.2f} m from the mirror; light blue-shifted x{shift:.1f}")
print(f"  550 nm reflected at beta = 0.99 arrives at {550 / ((1 + 0.99) / (1 - 0.99)):.2f} nm, in the soft X-ray band")

section("No rest frame at v = c")
print("  gamma = 1/sqrt(1 - v^2/c^2) grows as: " + ", ".join(
    f"{1 / np.sqrt(dev * (2 - dev)):.1e}" for dev in [1e-2, 1e-6, 1e-10]))
print("  At v = c the transformation is undefined. Worse, postulate 2 says any inertial frame measures")
print("  that light at c, so a frame where it stands still contradicts the theory that defines frames.")
finish()
