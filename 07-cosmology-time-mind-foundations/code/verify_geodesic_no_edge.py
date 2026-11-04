"""C04: a ball in a finite space without an edge, tested on the 3-sphere and a flat torus.

The 3-sphere geodesic is integrated from the Christoffel symbols of the metric alone, so
the code works intrinsically, with no outside space. The embedding in R^4 is used only
afterwards, to check where the ball ended up.
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad, solve_ivp
from _harness import Checks
from common.constants import c, year, ly

C = Checks("Finite space without an edge")
Gly = 1e9 * ly

C.section("Geodesic on a 3-sphere, built from the metric only")
chi, th, ph = sp.symbols("chi theta phi")
X = [chi, th, ph]
R0 = 1.0
g = sp.diag(R0**2, R0**2 * sp.sin(chi)**2, R0**2 * sp.sin(chi)**2 * sp.sin(th)**2)
ginv = g.inv()
Gamma = [[[sum(ginv[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l]))
               for l in range(3)) / 2 for k in range(3)] for j in range(3)] for i in range(3)]
Gfun = sp.lambdify(X, Gamma, "numpy")


def rhs(s, y):
    q, v = y[:3], y[3:]
    G = np.array(Gfun(*q), dtype=float)
    return np.concatenate([v, -np.einsum("ijk,j,k->i", G, v, v)])


q0 = np.array([1.1, 1.3, 0.2])
v0 = np.array([0.3, 0.5, 0.7])
gm = np.array(sp.lambdify(X, g, "numpy")(*q0), dtype=float)
v0 = v0 / np.sqrt(v0 @ gm @ v0)  # unit proper speed
L = 2 * np.pi * R0
sol = solve_ivp(rhs, [0, L], np.concatenate([q0, v0]), method="DOP853", rtol=1e-11, atol=1e-12,
                dense_output=True)


def embed(q):
    a, b, d = q[0], q[1], q[2]
    return R0 * np.array([np.cos(a), np.sin(a) * np.cos(b), np.sin(a) * np.sin(b) * np.cos(d),
                          np.sin(a) * np.sin(b) * np.sin(d)])


s_grid = np.linspace(0, L, 2001)
path = np.array([embed(sol.sol(s)[:3]) for s in s_grid])
start = embed(q0)
sep = np.linalg.norm(path - start, axis=1)
print(f"  path length walked = 2 pi R = {L:.4f} R; distance from start at the end = {sep[-1]:.2e} R")
print(f"  farthest from start along the way = {sep.max():.6f} R (antipode is 2 R away, reached at s = {s_grid[sep.argmax()]:.4f} R)")
C.check("forward-moving ball returns to its starting point after exactly 2 pi R", sep[-1] < 1e-6)
C.check("it goes as far as the antipode, so it really traversed the space", abs(sep.max() - 2) < 1e-5)
steps = np.linalg.norm(np.diff(path, axis=0), axis=1)
C.check("no jump anywhere on the path (step sizes uniform to 0.01%)",
        (steps.max() - steps.min()) / steps.mean() < 1e-4)
# Covariant acceleration: the geodesic equation residual measured by finite differences.
s_m = 1.7
h = 1e-3
qm, qp, q_c = sol.sol(s_m - h)[:3], sol.sol(s_m + h)[:3], sol.sol(s_m)[:3]
acc = (qp - 2 * q_c + qm) / h**2
vel = (qp - qm) / (2 * h)
G = np.array(Gfun(*q_c), dtype=float)
cov = acc + np.einsum("ijk,j,k->i", G, vel, vel)
print(f"  covariant acceleration at mid-path = {np.abs(cov).max():.1e} (zero means no force acted)")
C.check("no force acted on the ball: covariant acceleration vanishes", np.abs(cov).max() < 1e-4)
# Compare with a different direction: every direction closes, so no edge in any direction.
ok = True
rng = np.random.default_rng(1)
for _ in range(3):
    v = rng.normal(size=3)
    v = v / np.sqrt(v @ gm @ v)
    s2 = solve_ivp(rhs, [0, L], np.concatenate([q0, v]), method="DOP853", rtol=1e-8, atol=1e-9)
    ok &= np.linalg.norm(embed(s2.y[:3, -1]) - start) < 1e-4
C.check("three random throw directions all close after 2 pi R", bool(ok))

C.section("Measuring curvature from inside: triangle angle excess")
a = np.array([1.0, 0, 0])
b = np.array([0, 1.0, 0])
d = np.array([0, 0, 1.0])


def angle(p, q, r):
    t1 = q - p * (p @ q)
    t2 = r - p * (p @ r)
    return np.arccos(t1 @ t2 / np.linalg.norm(t1) / np.linalg.norm(t2))


sumang = angle(a, b, d) + angle(b, a, d) + angle(d, a, b)
print(f"  octant triangle on a unit sphere: angle sum = {np.degrees(sumang):.1f} deg; "
      f"excess = {sumang-np.pi:.4f}; area/R^2 = {4*np.pi/8:.4f}")
C.check("angle excess equals area/R^2, detectable without any outside space",
        abs((sumang - np.pi) - 4 * np.pi / 8) < 1e-12)

C.section("Flat torus: zero curvature, still finite and edgeless")
Lx, Ly = 10.0, 7.0
vx, vy = 2.0, 1.4
pos = np.array([1.0, 2.0])
start_t = pos.copy()
dt = 1e-3
T_close = 5.0  # vx*T = Lx, vy*T = Ly
for _ in range(int(round(T_close / dt))):
    pos = (pos + dt * np.array([vx, vy])) % np.array([Lx, Ly])
print(f"  after T = {T_close}: position {pos}, started at {start_t}")
C.check("straight line on a flat torus also returns with no force and no wall",
        np.linalg.norm(pos - start_t) < 1e-6)

C.section("Measured curvature and what it allows")
H0 = 67.4e3 / 3.0856775814913673e22
RH = c / H0
for Ok in [0.005, 0.003]:
    print(f"  |Omega_K| = {Ok}: curvature radius c/(H0 sqrt|Omega_K|) = {RH/np.sqrt(Ok)/Gly:.0f} Gly")
R_lim = RH / np.sqrt(0.003)
print(f"  Planck+BAO Omega_K = 0.001 +/- 0.002: at -2 sigma (closed side) Omega_K = -0.003, "
      f"R > {R_lim/Gly:.0f} Gly, circumference > {2*np.pi*R_lim/Gly:.0f} Gly")

C.section("Can a real ball circle a closed universe? Two tests of light")
# (1) Closed matter-only universes: light's comoving angle from bang to crunch.
for Om in [1.5, 2.0, 5.0]:
    amax = Om / (Om - 1)
    sqrtk = np.sqrt(Om - 1)
    # d(angle) = sqrt(Om-1) H0 dt / a, with da/dt = H0 sqrt(Om/a + 1 - Om); expansion plus recollapse.
    ang = 2 * sqrtk * quad(lambda x: 1 / (x * np.sqrt(Om / x - (Om - 1))), 0, amax, limit=500)[0]
    print(f"  closed matter universe, Omega_m = {Om}: light's angle over bang-to-crunch = {ang/np.pi:.5f} pi")
    C.check(f"Omega_m = {Om}: light goes round exactly once in the universe's whole life",
            abs(ang - 2 * np.pi) < 1e-5)
# (2) LCDM with accelerated expansion: light launched today is limited by the event horizon.
OmL = 0.685
Omm = 0.315
chi_e = quad(lambda u: np.exp(-u) / np.sqrt(Omm * np.exp(-3 * u) + OmL), 0, 40, limit=500)[0] * RH
circ = 2 * np.pi * R_lim
print(f"  light launched today covers at most {chi_e/Gly:.1f} Gly of comoving distance; "
      f"needed for a lap > {circ/Gly:.0f} Gly")
C.check("with Lambda and the Planck curvature bound, not even light launched now can finish a lap",
        chi_e < circ)
# A ball at 0.1 c today: peculiar momentum falls as 1/a.
v_now = 0.1


def ball_rate(u):
    a = np.exp(u)
    p = (v_now / np.sqrt(1 - v_now**2)) / a
    v = p / np.sqrt(1 + p**2)
    return v / (a * np.sqrt(Omm / a**3 + OmL))


chi_ball = quad(ball_rate, 0, 40, limit=500)[0] * RH
print(f"  ball thrown today at 0.1 c: total comoving distance ever = {chi_ball/Gly:.2f} Gly (it slows as space expands)")
C.check("the ball's lifetime travel is finite and smaller than light's", chi_ball < chi_e)
C.finish()
