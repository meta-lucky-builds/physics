"""E01, E04: from the interatomic potential to Hooke's law to a coil spring.

Checks, each done without re-evaluating the formula under test:
  1. Taylor expansion of Lennard-Jones and Morse potentials (symbolic), then a relaxed
     atomic chain whose force-extension curve is measured numerically.
  2. The chain's breaking force, found by trying to solve for equilibrium.
  3. Argon's bulk modulus from an fcc lattice sum of the pair potential, against measurement.
  4. k = G d^4 / (8 D^3 n) against a 3D finite-element beam model of the actual helix.
  5. Rubber: stiffness of a random chain measured from its sampled end-to-end statistics.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.optimize as opt
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import sympy as sy

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import k_B  # noqa: E402

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- 3D beam model of a helix
def _beam_stiffness(L, E, G, A, I, J):
    """12x12 local Euler-Bernoulli beam stiffness (axis = local x)."""
    k = np.zeros((12, 12))
    a, t = E * A / L, G * J / L
    b, c, d, e = 12 * E * I / L**3, 6 * E * I / L**2, 4 * E * I / L, 2 * E * I / L
    k[np.ix_([0, 6], [0, 6])] = [[a, -a], [-a, a]]
    k[np.ix_([3, 9], [3, 9])] = [[t, -t], [-t, t]]
    # bending in local x-y: dof 1, 5, 7, 11
    k[np.ix_([1, 5, 7, 11], [1, 5, 7, 11])] = [
        [b, c, -b, c], [c, d, -c, e], [-b, -c, b, -c], [c, e, -c, d]]
    # bending in local x-z: dof 2, 4, 8, 10
    k[np.ix_([2, 4, 8, 10], [2, 4, 8, 10])] = [
        [b, -c, -b, -c], [-c, d, c, e], [-b, c, b, c], [-c, e, c, d]]
    return k


def helix_nodes(D, pitch, n_turns, per_turn):
    theta = np.linspace(0, 2 * np.pi * n_turns, int(round(n_turns * per_turn)) + 1)
    return np.c_[D / 2 * np.cos(theta), D / 2 * np.sin(theta), pitch * theta / (2 * np.pi)]


def helix_fe(d, D, pitch, n_turns, G, nu=0.29, per_turn=48, mode="axial"):
    """Solve a clamped helix.

    mode="axial": force F along the coil axis, applied through a rigid arm at the free end.
                  Returns k_axial and the shares of strain energy in torsion/bending/stretch.
    mode="lateral": free end kept parallel (rotations fixed), sideways force.
                  Returns k_lateral.
    """
    E = 2 * G * (1 + nu)
    A, I, J = np.pi * d**2 / 4, np.pi * d**4 / 64, np.pi * d**4 / 32
    X = helix_nodes(D, pitch, n_turns, per_turn)
    N = len(X)
    rows, cols, vals, elems = [], [], [], []
    for i in range(N - 1):
        dx = X[i + 1] - X[i]
        L = np.linalg.norm(dx)
        ex = dx / L
        ez = np.cross(ex, [0.0, 0.0, 1.0])
        ez /= np.linalg.norm(ez)
        ey = np.cross(ez, ex)
        R3 = np.vstack([ex, ey, ez])
        T = np.zeros((12, 12))
        for j in range(4):
            T[3 * j:3 * j + 3, 3 * j:3 * j + 3] = R3
        kl = _beam_stiffness(L, E, G, A, I, J)
        ke = T.T @ kl @ T
        dofs = np.r_[6 * i:6 * i + 12]
        rr, cc = np.meshgrid(dofs, dofs, indexing="ij")
        rows.append(rr.ravel()); cols.append(cc.ravel()); vals.append(ke.ravel())
        elems.append((dofs, T, L))
    K = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(6 * N, 6 * N))
    fixed = list(range(6))  # clamp node 0
    f = np.zeros(6 * N)
    tip = 6 * (N - 1)
    F = 1.0

    if mode == "axial":
        Fvec = np.array([0, 0, F])
        arm = np.array([0, 0, X[-1, 2]]) - X[-1]  # from the wire end to the coil axis
        f[tip:tip + 3] = Fvec
        f[tip + 3:tip + 6] = np.cross(arm, Fvec)
    else:
        fixed += list(range(tip + 3, tip + 6))
        results = []
        for direction in (0, 1):
            f[:] = 0
            f[tip + direction] = F
            free = np.setdiff1d(np.arange(6 * N), fixed)
            u = np.zeros(6 * N)
            u[free] = spla.spsolve(K[free][:, free], f[free])
            results.append(F / u[tip + direction])
        return float(np.mean(results))

    free = np.setdiff1d(np.arange(6 * N), fixed)
    u = np.zeros(6 * N)
    u[free] = spla.spsolve(K[free][:, free], f[free])
    deflection = f @ u / F  # the work-conjugate displacement of the load on the axis
    U = np.zeros(3)  # torsion, bending, stretch
    for dofs, T, L in elems:
        ul = T @ u[dofs]
        U[0] += 0.5 * G * J / L * (ul[9] - ul[3]) ** 2
        U[2] += 0.5 * E * A / L * (ul[6] - ul[0]) ** 2
        for idx, sgn in (([1, 5, 7, 11], 1.0), ([2, 4, 8, 10], -1.0)):
            q = ul[idx] * np.array([1, sgn, 1, sgn])
            ke = np.array([[12, 6 * L, -12, 6 * L], [6 * L, 4 * L**2, -6 * L, 2 * L**2],
                           [-12, -6 * L, 12, -6 * L], [6 * L, 2 * L**2, -6 * L, 4 * L**2]])
            U[1] += 0.5 * E * I / L**3 * q @ ke @ q
    return F / deflection, U / U.sum()


def main():
    # -------------------------------------------------------------- 1. Taylor expansion
    print("1. Interatomic potential -> Hooke's law")
    r, eps, sig, D, a, r0 = sy.symbols("r epsilon sigma D a r_0", positive=True)
    lj = 4 * eps * ((sig / r) ** 12 - (sig / r) ** 6)
    rm = 2 ** sy.Rational(1, 6) * sig
    check("LJ force vanishes at r0 = 2^(1/6) sigma",
          sy.simplify(sy.diff(lj, r).subs(r, rm)) == 0)
    k_lj = sy.simplify(sy.diff(lj, r, 2).subs(r, rm))
    c3_lj = sy.simplify(sy.diff(lj, r, 3).subs(r, rm))
    print(f"     LJ: U''(r0) = {sy.nsimplify(k_lj)} = {float(k_lj / eps * sig**2):.3f} eps/sig^2;"
          f" U'''(r0) = {float(c3_lj / eps * sig**3):.2f} eps/sig^3")
    morse = D * (1 - sy.exp(-a * (r - r0))) ** 2
    k_m = sy.simplify(sy.diff(morse, r, 2).subs(r, r0))
    c3_m = sy.simplify(sy.diff(morse, r, 3).subs(r, r0))
    print(f"     Morse: U''(r0) = {k_m}, U'''(r0) = {c3_m}")
    check("Morse k = 2 D a^2", sy.simplify(k_m - 2 * D * a**2) == 0)

    # numerical LJ pair law, reduced units eps = sig = 1
    r_eq = 2 ** (1 / 6)
    F_pair = lambda x: 24 * (2 / x**13 - 1 / x**7)  # repulsion positive
    k_num = 57.146

    # relaxed 1D chain: solve for interior atoms at imposed end displacement
    Nb = 20

    def chain_force(strain):
        L = Nb * r_eq * (1 + strain)
        x0 = np.linspace(0, L, Nb + 1)[1:-1]

        def resid(xi):
            x = np.r_[0.0, xi, L]
            f = F_pair(np.diff(x))
            return f[1:] - f[:-1]

        sol = opt.root(resid, x0, tol=1e-12)
        x = np.r_[0.0, sol.x, L]
        return -F_pair(np.diff(x))[0]  # tension positive

    strains = np.array([-0.02, -0.01, -0.005, -0.002, -0.001, 0.001, 0.002, 0.005, 0.01, 0.02])
    tens = np.array([chain_force(s) for s in strains])
    secant = tens / (strains * r_eq)
    print("     chain secant stiffness per bond vs strain (eps/sig^2):")
    for s, ks in zip(strains, secant):
        print(f"       strain {s:+.3f}: {ks:7.2f}   departure from k: {100 * (ks / k_num - 1):+6.2f} %")
    # Secant = k + (U'''/2) x + O(x^2). The even part isolates k and the odd part isolates U'''.
    s_plus, s_minus = secant[strains == 0.001][0], secant[strains == -0.001][0]
    k_chain = 0.5 * (s_plus + s_minus)
    c3_chain = (s_plus - s_minus) / (0.001 * r_eq)  # = U''' from the two-sided difference
    print(f"     chain, even part of the secant at 0.1% strain: {k_chain:.3f} (U'' = {k_num:.3f});"
          f" odd part gives U''' = {c3_chain:.0f} (exact {float(c3_lj / eps * sig**3):.0f})")
    check("relaxed chain reproduces U''(r0) to 0.1%", abs(k_chain / k_num - 1) < 1e-3)
    check("relaxed chain reproduces U'''(r0) to 2%",
          abs(c3_chain / float(c3_lj / eps * sig**3) - 1) < 0.02)
    s_big = chain_force(0.05) / (0.05 * r_eq)
    c_big = chain_force(-0.05) / (-0.05 * r_eq)
    print(f"     at 5% strain: stretch {s_big:.1f}, compression {c_big:.1f} (eps/sig^2)")
    check("anharmonic: compression is stiffer than stretch at 5% strain", c_big > 1.3 * s_big)

    # -------------------------------------------------------------- 2. breaking force
    print("2. Elastic limit of a perfect chain, by searching for equilibrium under a dead load")
    r_infl = (26 / 7) ** (1 / 6)
    F_max_pred = -F_pair(r_infl)

    def has_equilibrium(load):
        grid = np.linspace(r_eq, 3.0, 200001)
        return np.any(-F_pair(grid) >= load)

    lo, hi = 0.1, 10.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if has_equilibrium(mid) else (lo, mid)
    print(f"     largest load with an equilibrium: {lo:.4f} eps/sig, at inflection predicted"
          f" {F_max_pred:.4f} eps/sig, strain {100 * (r_infl / r_eq - 1):.1f} %")
    check("breaking load equals the force at the potential's inflection point",
          abs(lo / F_max_pred - 1) < 1e-3)

    # -------------------------------------------------------------- 3. argon lattice
    print("3. Argon bulk modulus from a Lennard-Jones fcc lattice sum (static, T = 0)")
    n_shell = 8
    rng = range(-n_shell, n_shell + 1)
    pts = np.array([(i, j, k) for i in rng for j in rng for k in rng if (i + j + k) % 2 == 0
                    and (i, j, k) != (0, 0, 0)], dtype=float)
    # fcc sites: integer triples with even sum, in units of nearest-neighbour distance/sqrt(2)
    dist = np.linalg.norm(pts, axis=1) / np.sqrt(2)  # unit nearest-neighbour distance
    # keep a complete sphere of neighbours only
    dist = dist[dist <= n_shell / np.sqrt(2)]

    def energy_per_atom(rnn):
        rr = rnn * dist
        return 0.5 * np.sum(4 * (rr ** -12 - rr ** -6))

    res = opt.minimize_scalar(energy_per_atom, bounds=(0.9, 1.3), method="bounded",
                              options={"xatol": 1e-9})
    rnn0 = res.x

    def e_of_v(v):  # v = volume per atom in sig^3; fcc: v = rnn^3 / sqrt(2)
        return energy_per_atom((v * np.sqrt(2)) ** (1 / 3))

    v0 = rnn0**3 / np.sqrt(2)
    dv = 0.002 * v0
    B_reduced = v0 * (e_of_v(v0 + dv) - 2 * e_of_v(v0) + e_of_v(v0 - dv)) / dv**2
    print(f"     nearest neighbour {rnn0:.4f} sig, B = {B_reduced:.2f} eps/sig^3"
          f" (independent numerical second derivative)")
    B_measured = 2.88  # GPa, Wilkins 1973 (T=0, P=0), 28.8 +/- 0.4 kbar
    sets = [("eps/k=119.8 K, sig=3.405 A (widely used)", 119.8, 3.405),
            ("eps/k=124 K, sig=3.42 A", 124.0, 3.42),
            ("eps/k=125.7 K, sig=3.345 A (White 1999)", 125.7, 3.345)]
    ok = True
    for name, ek, s in sets:
        B = B_reduced * ek * k_B / (s * 1e-10) ** 3 / 1e9
        print(f"     {name}: B = {B:.2f} GPa, measured 2.88 GPa, ratio {B / B_measured:.2f}")
        ok = ok and 0.8 < B / B_measured < 1.5
    check("gas/solid fitted pair potentials give argon's B within a factor 1.5", ok)

    # -------------------------------------------------------------- 4. the coil
    print("4. Coil spring: finite-element helix vs k = G d^4 / (8 D^3 n)")
    G = 72e9  # music wire, MakeItFrom listing
    nu = 0.29
    d, Dm, n_turns, pitch = 2e-3, 20e-3, 10, 6e-3
    formula = lambda d, Dm, n: G * d**4 / (8 * Dm**3 * n)
    k_prev = None
    print(f"     wire d = {d * 1e3:.1f} mm, coil D = {Dm * 1e3:.0f} mm, n = {n_turns}, pitch = {pitch * 1e3:.0f} mm")
    for per_turn in (12, 24, 48, 96):
        k_fe, shares = helix_fe(d, Dm, pitch, n_turns, G, nu, per_turn)
        print(f"       {per_turn:3d} elements/turn: k_FE = {k_fe:8.2f} N/m,"
              f" k_FE/k_formula = {k_fe / formula(d, Dm, n_turns):.4f}")
        k_prev = k_fe
    k_fe, shares = helix_fe(d, Dm, pitch, n_turns, G, nu, 96)
    print(f"     k_formula = {formula(d, Dm, n_turns):.2f} N/m")
    print(f"     strain-energy shares: torsion {100 * shares[0]:.2f} %,"
          f" bending {100 * shares[1]:.2f} %, stretch {100 * shares[2]:.3f} %")
    check("wire torsion carries over 98% of the energy at small pitch", shares[0] > 0.98)
    check("FE stiffness within 1% of the formula at small pitch",
          abs(k_fe / formula(d, Dm, n_turns) - 1) < 0.01)

    def scale_ratio(**kw):
        p = dict(d=d, D=Dm, n=n_turns)
        k0, _ = helix_fe(p["d"], p["D"], pitch, p["n"], G, nu, 48)
        p2 = {**p, **kw}
        k1, _ = helix_fe(p2["d"], p2["D"], pitch, p2["n"], G, nu, 48)
        return k1 / k0

    r_d, r_D, r_n = scale_ratio(d=2 * d), scale_ratio(D=2 * Dm), scale_ratio(n=2 * n_turns)
    print(f"     FE scaling: 2d -> x{r_d:.2f} (formula 16), 2D -> x{r_D:.4f} (formula 0.125),"
          f" 2n -> x{r_n:.3f} (formula 0.5)")
    check("FE scaling exponents d^4, D^-3, n^-1",
          abs(r_d / 16 - 1) < 0.03 and abs(r_D / 0.125 - 1) < 0.03 and abs(r_n / 0.5 - 1) < 0.03)

    print("     steep coils: exact helix energy result vs FE (formula only holds at small pitch)")
    E = 2 * G * (1 + nu)
    A, I, J = np.pi * d**2 / 4, np.pi * d**4 / 64, np.pi * d**4 / 32
    ok = True
    for alpha_deg in (5, 15, 30, 45):
        al = np.radians(alpha_deg)
        p = np.pi * Dm * np.tan(al)
        R = Dm / 2
        L_wire = n_turns * np.pi * Dm / np.cos(al)
        flex = L_wire * (R**2 * np.cos(al) ** 2 / (G * J) + R**2 * np.sin(al) ** 2 / (E * I)
                         + np.sin(al) ** 2 / (E * A))
        k_exact = 1 / flex
        k_f, sh = helix_fe(d, Dm, p, n_turns, G, nu, 96)
        print(f"       pitch angle {alpha_deg:2d} deg: k_FE {k_f:8.1f}, helix energy {k_exact:8.1f},"
              f" simple formula {formula(d, Dm, n_turns):8.1f} N/m; torsion share {100 * sh[0]:.1f} %")
        ok = ok and abs(k_f / k_exact - 1) < 0.01
    check("FE agrees with the exact small-strain helix energy at all pitch angles", ok)

    # -------------------------------------------------------------- 5. rubber
    print("5. Entropic stiffness of a random chain (rubber), from sampled statistics")
    rng2 = np.random.default_rng(1)
    Nlink, samples = 100, 400000
    # freely jointed chain: x-projection of each random unit link is uniform on [-1, 1]
    x = rng2.uniform(-1, 1, size=(samples, Nlink)).sum(axis=1)
    hist, edges = np.histogram(x, bins=np.linspace(-20, 20, 41))
    centres = 0.5 * (edges[1:] + edges[:-1])
    good = hist > 200
    slope = np.polyfit(centres[good] ** 2, np.log(hist[good]), 1)[0]
    k_ent = -2 * slope  # free energy / kT = -ln P = (k_ent/2) x^2
    print(f"     fitted k = {k_ent:.5f} kT/b^2 for N = {Nlink} links, Gaussian theory 3/N = {3 / Nlink:.5f}")
    check("entropic spring constant 3 kT / (N b^2) within 3%", abs(k_ent * Nlink / 3 - 1) < 0.03)
    T = 300.0
    rho, Mc = 900.0, 5.0  # kg/m^3, kg/mol between crosslinks (assumed, typical soft vulcanised rubber)
    n_c = rho / Mc * 6.022e23
    print(f"     affine network G = n k_B T with n = {n_c:.2e} /m^3 at {T:.0f} K:"
          f" G = {n_c * k_B * T / 1e6:.2f} MPa, E = 3G = {3 * n_c * k_B * T / 1e6:.2f} MPa")

    print()
    if FAILS:
        print("FAILED:", FAILS)
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
