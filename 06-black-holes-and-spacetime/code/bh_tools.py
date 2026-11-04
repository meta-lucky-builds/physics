"""Helpers shared by the verify_*.py scripts of this section (not a check itself)."""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from common import constants as C  # noqa: E402

l_P = (C.hbar * C.G / C.c**3) ** 0.5      # Planck length, m
t_P = l_P / C.c                           # Planck time, s
m_P = (C.hbar * C.c / C.G) ** 0.5         # Planck mass, kg

# Masses used throughout, in kg.  Sources are given in the README.
M_SGRA = 4.15e6 * C.M_sun                 # GRAVITY 2019
M_M87 = 6.5e9 * C.M_sun                   # EHT 2019


def t_g(M):
    """GM/c^3 in seconds: one geometric unit of time for mass M (kg)."""
    return C.G * M / C.c**3


def r_g(M):
    """GM/c^2 in metres: one geometric unit of length for mass M (kg)."""
    return C.G * M / C.c**2


def T_hawking(M):
    return C.hbar * C.c**3 / (8 * math.pi * C.G * M * C.k_B)


def t_evap(M):
    """Photon-only blackbody evaporation time, s (derived in D7)."""
    return 5120 * math.pi * C.G**2 * M**3 / (C.hbar * C.c**4)


class Report:
    def __init__(self, title):
        self.fails = 0
        self.count = 0
        print(title)
        print("=" * len(title))

    def section(self, title):
        print("\n--- " + title)

    def note(self, text):
        print(text)

    def check(self, name, ok, detail=""):
        self.count += 1
        ok = bool(ok)
        if not ok:
            self.fails += 1
        print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""))

    def close(self):
        print(f"\n{self.count - self.fails}/{self.count} checks passed")
        sys.exit(1 if self.fails else 0)


def curvature(g, x):
    """Riemann (all indices down), Ricci and Kretschmann of a metric, symbolically."""
    import sympy as sp

    n = len(x)
    idx = range(n)
    ginv = g.inv().applyfunc(sp.simplify)
    Gam = [[[sp.simplify(sum(ginv[a, d] * (sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b])
                                           - sp.diff(g[b, c], x[d])) for d in idx) / 2)
             for c in idx] for b in idx] for a in idx]

    def r_up(a, b, c, d):
        e = sp.diff(Gam[a][d][b], x[c]) - sp.diff(Gam[a][c][b], x[d])
        e += sum(Gam[a][c][k] * Gam[k][d][b] - Gam[a][d][k] * Gam[k][c][b] for k in idx)
        return e

    Rup = [[[[sp.simplify(r_up(a, b, c, d)) for d in idx] for c in idx] for b in idx] for a in idx]
    Ric = sp.zeros(n, n)
    for b in idx:
        for d in idx:
            Ric[b, d] = sp.simplify(sum(Rup[a][b][a][d] for a in idx))
    Rdn = [[[[sp.simplify(sum(g[a, k] * Rup[k][b][c][d] for k in idx))
              for d in idx] for c in idx] for b in idx] for a in idx]
    K = 0
    nz = [(p, q, s, t) for p in idx for q in idx for s in idx for t in idx if Rdn[p][q][s][t] != 0]
    for a, b, c, d in nz:
        v = sum(ginv[a, p] * ginv[b, q] * ginv[c, s] * ginv[d, t] * Rdn[p][q][s][t] for p, q, s, t in nz)
        K += Rdn[a][b][c][d] * v
    return Rdn, Ric, sp.simplify(K), ginv, Rup
