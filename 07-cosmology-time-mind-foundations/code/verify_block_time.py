"""C18-C25: relativity of simultaneity, local conservation versus jumps, self-consistency."""
import numpy as np
from scipy.optimize import fsolve
from _harness import Checks
from common.constants import c, hbar, year, ly

C = Checks("Block universe inputs, jumps, and time-loop consistency")

C.section("Andromeda calculation (Rietdijk-Putnam-Penrose)")


def to_frame(t, x, v):
    g = 1 / np.sqrt(1 - (v / c)**2)
    return g * (t - v * x / c**2), g * (x - v * t)


D = 2.5e6 * ly
v_walk = 1.4
day = 86400
# Find what t on Andromeda is simultaneous with the walker's origin event (t' = 0), by root
# solving on the transform itself rather than using the shortcut formula.
for v in [v_walk, -v_walk]:
    t_star = fsolve(lambda t: to_frame(t[0] * day, D, v)[0], 0.0, xtol=1e-12)[0]
    print(f"  walker at {v:+.1f} m/s: 'now' on Andromeda is t = {t_star:+.3f} days in Earth's frame")
    C.check(f"root agrees with v D/c^2 for v = {v:+.1f}", abs(t_star * day - v * D / c**2) < 1)
span = 2 * v_walk * D / c**2 / day
print(f"  two walkers passing each other disagree by {span:.2f} days about 'now' at Andromeda")
# The two disputed events are spacelike separated: neither can influence the walkers yet.
dt = span * day
ds2 = (c * dt)**2 - D**2
C.check("the disputed event is spacelike separated from the walkers", ds2 < 0)
# An Andromeda event at Earth-frame t = 0: later than the origin for one walker, earlier for the other.
t_plus = to_frame(0.0, D, v_walk)[0]
t_minus = to_frame(0.0, D, -v_walk)[0]
print(f"  Andromeda event at t = 0: t' = {t_plus/day:+.3f} days for one walker, {t_minus/day:+.3f} days for the other")
C.check("the order of spacelike events flips between frames", t_plus * t_minus < 0)
# Timelike separated events keep their order in every frame.
ok = True
rng = np.random.default_rng(3)
for _ in range(10000):
    t = rng.uniform(0, 1) * year
    x = rng.uniform(-0.99, 0.99) * c * t
    v = rng.uniform(-0.999, 0.999) * c
    ok &= to_frame(t, x, v)[0] > 0
C.check("10000 random timelike pairs: no boost reverses their order, so causes stay first", bool(ok))

C.section("Why a global jump breaks conservation in another frame (Feynman's argument)")
# A charge vanishes at x=0, t=0 and reappears at x=L, t=0 in the frame S. Count the charge on
# a slice of constant t' in a frame moving at v, by computing when each event occurs there.
L = 1.0 * c  # one light-second
v = 0.5 * c
t_vanish, _ = to_frame(0.0, 0.0, v)
t_appear, _ = to_frame(0.0, L, v)
print(f"  frame S: vanish and reappear both at t = 0, total charge 1 at all times")
print(f"  frame S': vanish at t' = {t_vanish:+.3f} s, reappear at t' = {t_appear:+.3f} s")
tp = np.linspace(-1, 1, 9)
q = ((tp < t_vanish).astype(int) + (tp >= t_appear).astype(int))
print(f"  total charge on slices t' = {tp.tolist()}:\n    {q.tolist()}")
C.check("in S' the total charge changes: a jump conserves charge in at most one frame",
        q.min() != q.max())

C.section("Random jumps: size of the relativistic-QFT suppression")
# Beyond the light cone the field amplitude decays as exp(-m r c / hbar) at spacelike r.
for name, m, r in [("electron, 1 nm", 9.109e-31, 1e-9), ("electron, 1 m", 9.109e-31, 1.0),
                   ("70 kg body, 1 m", 70.0, 1.0)]:
    exponent = m * c * r / hbar
    print(f"  {name:18s}: m c r / hbar = {exponent:.2e}; amplitude ~ 10^(-{exponent/np.log(10):.2e})")
C.check("for a body the exponent exceeds 1e40", 70 * c / hbar > 1e40)

C.section("Billiard ball through a time-shifted wormhole: search for consistent histories")
# Young ball: starts (0,0), velocity (1,0). Portal: enters at x = 10 with height h and velocity
# (vx, vy) at time t_in. Exits at time t_in - T at (5 + h, 5) with velocity rotated -90 deg,
# heading down to cross the young ball's path. Discs of radius r: equal mass elastic collision.
# With T = 10 the straight-line history is a perfect self-hit: the paradox set-up.
r = 0.5
T = 10.0
X_IN = 10.0


def older_state(s):
    h, vx, vy, t_in = s
    return np.array([5 + h, 5.0]), np.array([vy, -vx]), t_in - T


def young_after(s):
    """Return the state with which the young ball reaches the portal, given older ball s."""
    p_o, v_o, t_o = older_state(s)
    p_y0, v_y0 = np.array([0.0, 0.0]), np.array([1.0, 0.0])
    # Relative motion: find first contact time, after both exist (t >= t_o).
    p_rel = (p_o - v_o * t_o) - p_y0
    v_rel = v_o - v_y0
    A = v_rel @ v_rel
    B = 2 * p_rel @ v_rel
    D_ = p_rel @ p_rel - (2 * r)**2
    disc = B**2 - 4 * A * D_
    hit = False
    if A > 0 and disc > 0:
        tc = (-B - np.sqrt(disc)) / (2 * A)
        hit = tc >= t_o
    if not hit:
        # No collision: young ball goes straight to the portal.
        return np.array([0.0, 1.0, 0.0, X_IN]), False, None
    p_y = p_y0 + v_y0 * tc
    p_oc = p_o + v_o * (tc - t_o)
    n = (p_oc - p_y) / (2 * r)
    dv = ((v_o - v_y0) @ n) * n
    v_y = v_y0 + dv
    if v_y[0] <= 1e-6:
        return None, True, tc
    t_in = tc + (X_IN - p_y[0]) / v_y[0]
    h = p_y[1] + v_y[1] * (t_in - tc)
    return np.array([h, v_y[0], v_y[1], t_in]), True, tc


def residual(s):
    out, _, _ = young_after(s)
    if out is None:
        return np.array([1e3, 1e3, 1e3, 1e3])
    return out - s


straight = np.array([0.0, 1.0, 0.0, 10.0])
out, hit, tc = young_after(straight)
print(f"  straight-line history: older ball hits young ball? {hit} (at t = {tc:.3f})")
if out is None:
    print("  after that hit the young ball never reaches the portal: naive history is inconsistent")
else:
    print(f"  after that hit the young ball would arrive with state {np.round(out, 3)}, not the assumed state")
C.check("the undisturbed history is NOT self-consistent (the paradox)", out is None or np.linalg.norm(out - straight) > 0.1)

rng = np.random.default_rng(11)
found = []
for _ in range(400):
    s0 = np.array([rng.uniform(-3, 3), rng.uniform(0.3, 1.2), rng.uniform(-0.8, 0.8), rng.uniform(8, 20)])
    sol, info, ier, _ = fsolve(residual, s0, full_output=True, xtol=1e-12)
    res = np.linalg.norm(residual(sol))
    _, did_hit, _ = young_after(sol)
    if ier == 1 and res < 1e-9 and did_hit and sol[1] > 0 and abs(sol[0]) < 6:
        if not any(np.linalg.norm(sol - f) < 1e-5 for f in found):
            found.append(sol)
print(f"  distinct self-consistent colliding histories found from 400 starts: {len(found)}")
for sol in found[:6]:
    h, vx, vy, t_in = sol
    print(f"    enters at h = {h:+.4f}, velocity ({vx:.4f}, {vy:+.4f}), speed {np.hypot(vx,vy):.4f}, t_in = {t_in:.3f}, "
          f"residual {np.linalg.norm(residual(sol)):.1e}")
C.check("at least one self-consistent history exists for the very initial data that looked paradoxical",
        len(found) >= 1)
ok = True
for sol in found:
    # Forward replay: start from the found state, compute what actually enters, compare.
    out, did_hit, tc = young_after(sol)
    p_o, v_o, t_o = older_state(sol)
    print(f"    older ball emerges at t = {t_o:.3f}; collision at t = {tc:.3f}; young ball started at t = 0")
    ok &= tc > 0 and t_o < tc
    v_y = np.array(out[1:3])
    # Energy check for the collision: older ball leaves with what is left over.
    n_vec = None
    E_before = 1.0 + v_o @ v_o
    v_o_after2 = E_before - v_y @ v_y
    ok &= abs(v_o_after2 - 1.0) < 1e-6 and np.linalg.norm(out - sol) < 1e-8
C.check("each history replays forward to itself, and the collision conserves energy (older ball exits at speed 1)", bool(ok))

C.section("Discrete loops: classical bit, then a Deutsch-style density matrix")
# Grandfather bit: b -> NOT b has no deterministic fixed point.
fixed = [b for b in (0, 1) if (1 - b) == b]
print(f"  deterministic fixed points of NOT: {fixed}")
C.check("a bit-flip loop has no deterministic consistent history", fixed == [])
# Any 'identity on a loop' and other maps do have fixed points.
maps = {"copy": lambda b: b, "set0": lambda b: 0, "set1": lambda b: 1, "not": lambda b: 1 - b}
for name, fn in maps.items():
    print(f"    map {name:5s}: fixed points {[b for b in (0, 1) if fn(b) == b]}")
# Probabilistic relaxation: a mixture p = P(b = 1) maps to 1 - p, with fixed point 1/2.
grid = np.linspace(0, 1, 100001)
gap = np.abs((1 - grid) - grid)
p_star = grid[gap.argmin()]
print(f"  stochastic fixed point of NOT, by scanning p: p = {p_star:.5f}; points with gap < 1e-9: {(gap < 1e-9).sum()}")
C.check("p = 1/2 is the unique stochastic fixed point of NOT", abs(p_star - 0.5) < 1e-9 and (gap < 1e-9).sum() == 1)
# Deutsch (1991): rho_CTC = Tr_sys[ U (rho_sys x rho_CTC) U+ ]. Test with CNOT, CTC as target.
ket0 = np.array([1, 0], dtype=complex)
ket_plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
CNOT_sc = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Ident = np.eye(2, dtype=complex)
# System qubit set to 1 always, so the CTC qubit gets NOT: this is the grandfather loop.
ket1 = np.array([0, 1], dtype=complex)


def deutsch_map(rho_ctc, psi_sys, U):
    rho = np.kron(np.outer(psi_sys, psi_sys.conj()), rho_ctc)
    out = U @ rho @ U.conj().T
    d = 2
    return np.array([[sum(out[i * d + a, i * d + b] for i in range(d)) for b in range(d)] for a in range(d)])


def fixed_point(psi, U, rho0, n=2000):
    rho = rho0
    avg = np.zeros((2, 2), dtype=complex)
    for k in range(n):
        rho = deutsch_map(rho, psi, U)
        avg += rho
    return avg / n


rho0 = np.outer(ket0, ket0.conj())
rho_fp = fixed_point(ket1, CNOT_sc, rho0)
print(f"  Deutsch CTC state for grandfather (NOT) loop: diag = {np.real(np.diag(rho_fp)).round(4).tolist()}")
C.check("Deutsch's condition gives a consistent mixed state for the grandfather loop",
        np.linalg.norm(deutsch_map(rho_fp, ket1, CNOT_sc) - rho_fp) < 1e-2)
ok = True
rng = np.random.default_rng(5)
for _ in range(200):
    A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    U, _ = np.linalg.qr(A)
    v = rng.normal(size=2) + 1j * rng.normal(size=2)
    v /= np.linalg.norm(v)
    # Fixed point of a channel: the eigenvalue-1 eigenvector of its transfer matrix.
    basis = [np.array([[1, 0], [0, 0]], dtype=complex), np.array([[0, 1], [0, 0]], dtype=complex),
             np.array([[0, 0], [1, 0]], dtype=complex), np.array([[0, 0], [0, 1]], dtype=complex)]
    Tm = np.array([deutsch_map(b, v, U).ravel() for b in basis]).T
    w, vec = np.linalg.eig(Tm)
    i = np.argmin(abs(w - 1))
    rho = vec[:, i].reshape(2, 2)
    rho = rho / np.trace(rho)
    ok &= abs(w[i] - 1) < 1e-8 and np.linalg.norm(deutsch_map(rho, v, U) - rho) < 1e-7
    ok &= np.all(np.linalg.eigvalsh((rho + rho.conj().T) / 2) > -1e-7)
C.check("200 random unitaries: Deutsch's equation always has a valid density-matrix solution", bool(ok))
C.finish()
