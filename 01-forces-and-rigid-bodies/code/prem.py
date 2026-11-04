"""Radial Earth model used by the F09 and F10 scripts (helper, not a check itself).

Density is the isotropic density polynomial of PREM, Dziewonski & Anderson (1981),
Phys. Earth Planet. Inter. 25, 297, with x = r / 6371 km and rho in g/cm^3. The
layer table gives (top radius in km, polynomial coefficients in increasing powers of x).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import G  # noqa: E402

R_MEAN = 6371.0e3  # m, NASA Earth Fact Sheet volumetric mean radius

LAYERS = [
    (1221.5, (13.0885, 0.0, -8.8381, 0.0)),            # inner core
    (3480.0, (12.5815, -1.2638, -3.6426, -5.5281)),    # outer core
    (5701.0, (7.9565, -6.4761, 5.5283, -3.0807)),      # lower mantle
    (5771.0, (5.3197, -1.4836, 0.0, 0.0)),             # transition zone
    (5971.0, (11.2494, -8.0298, 0.0, 0.0)),
    (6151.0, (7.1089, -3.8045, 0.0, 0.0)),
    (6346.6, (2.6910, 0.6924, 0.0, 0.0)),              # LVZ and LID
    (6356.0, (2.9, 0.0, 0.0, 0.0)),                    # lower crust
    (6368.0, (2.6, 0.0, 0.0, 0.0)),                    # upper crust
    (6371.0, (1.02, 0.0, 0.0, 0.0)),                   # ocean
]


def density(r):
    """PREM density in kg/m^3 at radius r in m (array or scalar)."""
    r = np.asarray(r, dtype=float)
    x = r / R_MEAN
    rho = np.zeros_like(x)
    lower = 0.0
    for top_km, c in LAYERS:
        top = top_km * 1e3
        mask = (r >= lower) & (r <= top)
        rho[mask] = 1e3 * (c[0] + c[1] * x[mask] + c[2] * x[mask] ** 2 + c[3] * x[mask] ** 3)
        lower = top
    return rho


def grid(n_per_layer=2000):
    """Radius grid that places nodes on both sides of each density jump."""
    edges = [0.0] + [t * 1e3 for t, _ in LAYERS]
    return [np.linspace(a, b, n_per_layer) for a, b in zip(edges[:-1], edges[1:])]


def profile():
    """Return r, rho, enclosed mass m(r), gravity g(r) and moment of inertia integrand total.

    Integrates layer by layer with Simpson-free cumulative trapezoids on a fine grid.
    Returns a dict of arrays plus totals M and I about a diameter.
    """
    rs, rhos, ms, Is = [], [], [], []
    m_acc = 0.0
    i_acc = 0.0
    edges = [0.0] + [t * 1e3 for t, _ in LAYERS]
    for (a, b), (_, c) in zip(zip(edges[:-1], edges[1:]), LAYERS):
        r = np.linspace(a, b, 2000)
        x = r / R_MEAN
        rho = 1e3 * (c[0] + c[1] * x + c[2] * x ** 2 + c[3] * x ** 3)
        dm = 4 * np.pi * r ** 2 * rho
        di = (8 * np.pi / 3) * r ** 4 * rho
        mm = m_acc + np.concatenate(([0.0], np.cumsum(0.5 * (dm[1:] + dm[:-1]) * np.diff(r))))
        ii = i_acc + np.concatenate(([0.0], np.cumsum(0.5 * (di[1:] + di[:-1]) * np.diff(r))))
        m_acc, i_acc = mm[-1], ii[-1]
        rs.append(r)
        rhos.append(rho)
        ms.append(mm)
        Is.append(ii)
    r = np.concatenate(rs)
    rho = np.concatenate(rhos)
    m = np.concatenate(ms)
    inertia = np.concatenate(Is)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(r > 0, G * m / r ** 2, 0.0)
    return {"r": r, "rho": rho, "m": m, "g": g, "I_enc": inertia, "M": m[-1], "I": inertia[-1]}


def gravity_function():
    """Return a callable g(r) in m/s^2, g(0) = 0, interpolated from the PREM profile."""
    p = profile()
    r, g = p["r"], p["g"]
    # r has repeated nodes at layer boundaries; g is continuous there, so drop repeats.
    keep = np.concatenate(([True], np.diff(r) > 0))
    r, g = r[keep], g[keep]
    return lambda s: np.interp(s, r, g)
