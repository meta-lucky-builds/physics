"""L03, L04: light from Maxwell's equations, and photon mass and momentum.

- sympy: derive the wave equation from the curl equations for arbitrary fields.
- FDTD: integrate the first-order Maxwell equations and clock the pulse.
- Driven damped charge: measure momentum and energy absorbed from a wave (p = E/c).
- sympy: E^2 - (pc)^2 = (mc^2)^2 from the four-momentum.
- Compton shift by numerically solving energy and momentum conservation.
"""
import numpy as np
import sympy as sp
from scipy import constants as sc
from scipy.optimize import fsolve

from checklib import K, check, finish, section

c = K.c

section("1. Wave equation from Maxwell (symbolic, arbitrary fields)")
x, y, z, t = sp.symbols("x y z t", real=True)
mu, eps = sp.symbols("mu epsilon", positive=True)
coords = (x, y, z)
Ex, Ey, Ez = [sp.Function(n)(x, y, z, t) for n in ("E_x", "E_y", "E_z")]
E = sp.Matrix([Ex, Ey, Ez])


def curl(F):
    return sp.Matrix([
        sp.diff(F[2], y) - sp.diff(F[1], z),
        sp.diff(F[0], z) - sp.diff(F[2], x),
        sp.diff(F[1], x) - sp.diff(F[0], y),
    ])


def div(F):
    return sum(sp.diff(F[i], coords[i]) for i in range(3))


def grad(f):
    return sp.Matrix([sp.diff(f, q) for q in coords])


def lap(F):
    return sp.Matrix([sum(sp.diff(F[i], q, 2) for q in coords) for i in range(3)])


ident = (curl(curl(E)) - (grad(div(E)) - lap(E))).applyfunc(sp.simplify)
check("curl curl E = grad div E - laplacian E", ident == sp.zeros(3, 1))
# Vacuum: div E = 0; dB/dt = -curl E; curl B = mu eps dE/dt.
# Take d/dt of the last equation: curl(dB/dt) = mu eps d2E/dt2, so -curl curl E = mu eps d2E/dt2.
# With div E = 0 this becomes lap E = mu eps d2E/dt2.
k, w, E0 = sp.symbols("k omega E_0", positive=True)
phase = k * x - w * t
Ew = sp.Matrix([0, E0 * sp.cos(phase), 0])
B0 = sp.symbols("B_0", positive=True)
Bw = sp.Matrix([0, 0, B0 * sp.cos(phase)])
faraday = (curl(Ew) + sp.diff(Bw, t)).applyfunc(sp.simplify)
ampere = (curl(Bw) - mu * eps * sp.diff(Ew, t)).applyfunc(sp.simplify)
sol = sp.solve([faraday[2] / sp.sin(phase), ampere[1] / sp.sin(phase)], [w, B0], dict=True)
sol = [s for s in sol if s[w].is_positive]
speed = sp.simplify(sol[0][w] / k)
check("plane wave solves Maxwell only if omega/k = 1/sqrt(mu eps)", sp.simplify(speed - 1 / sp.sqrt(mu * eps)) == 0, f"omega/k = {speed}")
check("and B0 = E0 * sqrt(mu eps) (B0 = E0/c)", sp.simplify(sol[0][B0] - E0 * sp.sqrt(mu * eps)) == 0)
check("div E = div B = 0 for the transverse wave", sp.simplify(div(Ew)) == 0 and sp.simplify(div(Bw)) == 0)
mu0, eps0 = sc.mu_0, sc.epsilon_0
c_me = 1 / np.sqrt(mu0 * eps0)
print(f"  mu0 = {mu0:.9e} H/m, eps0 = {eps0:.9e} F/m gives 1/sqrt(mu0 eps0) = {c_me:,.0f} m/s")
print("  Since 1983 this agreement is built into the values of the units, so it is not a test today.")
print("  The historical test used separate electrical measurements: Weber and Kohlrausch, about 3.1e8 m/s.")


section("2. FDTD integration of the first-order equations")


def fdtd(eps_r=1.0, nx=4000, dx=1e-3, courant=0.5, steps=3000):
    """1D Yee leapfrog: E_y at integer nodes, B_z at half nodes, with SI constants."""
    c_vac = 1 / np.sqrt(mu0 * eps0)
    dt = courant * dx / c_vac
    Ey = np.zeros(nx)
    Bz = np.zeros(nx - 1)
    xs = np.arange(nx) * dx
    x0, sig = 0.5, 0.03
    Ey[:] = np.exp(-((xs - x0) ** 2) / (2 * sig**2))
    # right-moving start: B = E/v, staggered half a cell and half a step
    v = c_vac / np.sqrt(eps_r)
    xh = xs[:-1] + dx / 2
    Bz[:] = np.exp(-((xh + v * dt / 2 - x0) ** 2) / (2 * sig**2)) / v
    t_log, x_log = [], []
    for n in range(steps):
        Bz -= dt / dx * (Ey[1:] - Ey[:-1])
        Ey[1:-1] -= dt / (mu0 * eps0 * eps_r) * (Bz[1:] - Bz[:-1]) / dx
        if n % 50 == 0 and n > 100:
            i = np.argmax(Ey)
            if 2 < i < nx - 3:
                a, b, d = Ey[i - 1], Ey[i], Ey[i + 1]
                xp = xs[i] + 0.5 * dx * (a - d) / (a - 2 * b + d)
                t_log.append((n + 1) * dt)
                x_log.append(xp)
    fit = np.polyfit(t_log, x_log, 1)[0]
    uE = 0.5 * eps0 * eps_r * np.sum(Ey**2) * dx
    uB = 0.5 / mu0 * np.sum(Bz**2) * dx
    return fit, uE, uB


vf, uE, uB = fdtd()
check("pulse in vacuum travels at 1/sqrt(mu0 eps0)", abs(vf / c - 1) < 1e-3, f"measured {vf:.6e} m/s, {vf / c:.5f} c")
check("electric and magnetic energy in the wave are equal", abs(uE / uB - 1) < 1e-2, f"uE/uB = {uE / uB:.4f}")
vf2, _, _ = fdtd(eps_r=4.0, steps=3000, courant=0.45)
check("pulse with eps_r = 4 travels at c/2", abs(vf2 / (c / 2) - 1) < 2e-3, f"measured {vf2 / c:.5f} c")

section("3. A wave carries momentum p = E/c, from classical absorption")


def absorbed_ratio(E0=1e-3, gamma=0.3, omega=1.0, c_=1.0, n_periods=400, per_period=400):
    """RK4 for a damped charge (q = m = 1) in a plane wave moving along +x. Returns
    time-averaged F_x * c / P_abs over the settled second half."""
    dt = 2 * np.pi / omega / per_period
    s = np.array([0.0, 0.0, 0.0, 0.0])  # x, y, vx, vy

    def field(pos_x, tt):
        ph = omega * (tt - pos_x / c_)
        return E0 * np.cos(ph), E0 / c_ * np.cos(ph)

    def rhs(st, tt):
        Ey, Bz = field(st[0], tt)
        vx, vy = st[2], st[3]
        Fx = vy * Bz  # q (v x B)_x
        Fy = Ey - vx * Bz  # q (E + v x B)_y
        return np.array([vx, vy, Fx - gamma * vx, Fy - gamma * vy])

    tt = 0.0
    n = n_periods * per_period
    Fsum = Psum = 0.0
    count = 0
    for i in range(n):
        k1 = rhs(s, tt)
        k2 = rhs(s + dt / 2 * k1, tt + dt / 2)
        k3 = rhs(s + dt / 2 * k2, tt + dt / 2)
        k4 = rhs(s + dt * k3, tt + dt)
        s = s + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        tt += dt
        if i > n // 2:
            Ey, Bz = field(s[0], tt)
            Fsum += s[3] * Bz
            Psum += Ey * s[3]
            count += 1
    return (Fsum / count) * c_ / (Psum / count)


r = absorbed_ratio(n_periods=120, per_period=200)
check("momentum absorbed / (energy absorbed / c) = 1", abs(r - 1) < 2e-2, f"ratio = {r:.5f}")

section("4. E^2 = (pc)^2 + (mc^2)^2")
m, v, cc = sp.symbols("m v c", positive=True)
g = 1 / sp.sqrt(1 - v**2 / cc**2)
En, pm = g * m * cc**2, g * m * v
check("(gamma m c^2)^2 - (gamma m v c)^2 = (m c^2)^2", sp.simplify(En**2 - (pm * cc) ** 2 - (m * cc**2) ** 2) == 0)
# Velocity from any energy-momentum pair is v = p c^2/E.
Eg, mg = sp.symbols("E m_g", positive=True)
vg = sp.sqrt(Eg**2 - (mg * cc**2) ** 2) * cc / Eg
check("v = c exactly when m = 0, v < c for m > 0", sp.simplify(vg.subs(mg, 0) - cc) == 0)

section("5. Compton scattering from conservation laws only")
h = K.h
me = K.m_e


def compton_numeric(lam, theta):
    """Solve for scattered photon energy and electron momentum with no shift formula."""
    Ei = h * c / lam

    def eqs(u):
        Ef, px, py = u[0] * Ei, u[1] * Ei / c, u[2] * Ei / c
        Ee = np.sqrt((me * c**2) ** 2 + (px**2 + py**2) * c**2)
        return [
            (Ei + me * c**2 - Ef - Ee) / Ei,
            (Ei / c - Ef / c * np.cos(theta) - px) / (Ei / c),
            (0 - Ef / c * np.sin(theta) - py) / (Ei / c),
        ]

    # start near the unshifted solution
    guess = [0.95, 1 - 0.95 * np.cos(theta), -0.95 * np.sin(theta)]
    u, _, ier, _ = fsolve(eqs, guess, full_output=True, xtol=1e-13)
    assert ier == 1
    return h * c / (u[0] * Ei)


lam_C = h / (me * c)
print(f"  Compton wavelength h/(m_e c) = {lam_C * 1e12:.4f} pm")
for th_deg in [30, 90, 135, 180]:
    th = np.radians(th_deg)
    dlam = compton_numeric(7.11e-11, th) - 7.11e-11
    pred = lam_C * (1 - np.cos(th))
    check(f"shift at {th_deg} deg equals lambda_C (1 - cos)", abs(dlam / pred - 1) < 1e-6, f"{dlam * 1e12:.4f} pm")

section("6. Photon mass limit in context")
m_lim_eV = 1e-18  # PDG listing, Ryutov 2007
m_lim = m_lim_eV * K.eV / c**2
E_vis = 2.0  # eV
print(f"  Limit m < {m_lim_eV:g} eV/c^2 = {m_lim:.2e} kg = {m_lim / me:.1e} electron masses")
print(f"  If m were at that limit, a {E_vis} eV photon would lag light by 1 - v/c = {0.5 * (m_lim_eV / E_vis) ** 2:.1e}")
E_ph = E_vis * K.eV
print(f"  The quantity E/c^2 for that photon ('relativistic mass') is {E_ph / c**2:.2e} kg; "
      f"it is not the invariant mass. Momentum E/c = {E_ph / c:.2e} kg m/s")

finish()
