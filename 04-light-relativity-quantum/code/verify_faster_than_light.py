"""L12: things that look faster than light, and why they don't break special relativity.

- Blob in a jet: simulate pulse arrival times, extract apparent speed.
- Dispersive wave: FFT propagation with omega^2 = omega_p^2 + c^2 k^2; crest speed,
  envelope speed, and a strict causal front.
- Entanglement: Bob's reduced state for random operations by Alice.
- Alcubierre: energy density from the 3+1 constraint with sympy, sign and total.
- Superluminal signals lead to reply-before-send: explicit event calculation.
"""
import numpy as np
import sympy as sp
from scipy.integrate import quad

from checklib import K, check, finish, section

c = K.c
rng = np.random.default_rng(13)

section("Apparent superluminal motion")


def apparent_speed(beta, theta):
    """Pulses emitted from a blob moving at beta, angle theta to our line of sight,
    received at t + (D - z)/c. Returns observed transverse rate in units of c."""
    t = np.linspace(0, 10.0, 400)
    z = beta * np.cos(theta) * t
    xt = beta * np.sin(theta) * t
    t_obs = t - z  # D/c constant removed, c = 1
    return np.polyfit(t_obs, xt, 1)[0]


for beta, deg in [(0.9, 30), (0.99, 8), (0.995, 5)]:
    th = np.radians(deg)
    a = apparent_speed(beta, th)
    f = beta * np.sin(th) / (1 - beta * np.cos(th))
    check(f"beta = {beta}, theta = {deg} deg: simulated rate matches the formula", abs(a / f - 1) < 1e-9, f"apparent {a:.3f} c")
for beta in [0.9, 0.99, 0.999]:
    ths = np.linspace(1e-4, np.pi / 2, 200001)
    vals = beta * np.sin(ths) / (1 - beta * np.cos(ths))
    i = vals.argmax()
    g = 1 / np.sqrt(1 - beta**2)
    check(f"beta = {beta}: maximum is gamma beta at cos(theta) = beta",
          abs(vals[i] / (g * beta) - 1) < 1e-6 and abs(np.cos(ths[i]) - beta) < 1e-3, f"max {vals[i]:.3f} c at {np.degrees(ths[i]):.2f} deg")
bmin = 6 / np.sqrt(37)
print(f"  An apparent 6c (M87 jet) needs gamma beta >= 6: beta >= {bmin:.4f}, gamma >= {1 / np.sqrt(1 - bmin**2):.2f}")

section("Phase, group and front velocity (c = 1)")
N = 2**15
L = 4000.0
dx = L / N
x = (np.arange(N) - N // 2) * dx
k = 2 * np.pi * np.fft.fftfreq(N, dx)
wp = 1.0
w = np.sqrt(wp**2 + k**2)
k0 = 2.0
sg = 25.0
psi0 = np.exp(-(x**2) / (2 * sg**2) + 1j * k0 * x)
F = np.fft.fft(psi0)


def at_time(t):
    return np.fft.ifft(F * np.exp(-1j * w * t))


def peak_pos(psi):
    a = np.abs(psi) ** 2
    i = a.argmax()
    y1, y2, y3 = np.log(a[i - 1]), np.log(a[i]), np.log(a[i + 1])
    return x[i] + 0.5 * dx * (y1 - y3) / (y1 - 2 * y2 + y3)


ts = np.array([100.0, 200.0, 300.0, 400.0])
pk = np.array([peak_pos(at_time(t)) for t in ts])
vg = np.polyfit(ts, pk, 1)[0]
w0 = np.sqrt(wp**2 + k0**2)
check("envelope moves at the group velocity k/omega < c", abs(vg - k0 / w0) < 2e-3, f"measured {vg:.5f} c, theory {k0 / w0:.5f} c")
# Crest speed: track the phase at the envelope peak over a short step
t1 = 200.0
dt = 0.5
p1 = peak_pos(at_time(t1))
i1 = np.argmin(np.abs(x - p1))
psi_a = at_time(t1)
psi_b = at_time(t1 + dt)
dphi = np.angle(psi_b[i1] / psi_a[i1])  # phase change at a fixed point = -omega dt
om_meas = -dphi / dt
vp = om_meas / k0
check("crests move at the phase velocity omega/k > c", vp > 1 and abs(vp - w0 / k0) < 2e-2, f"measured {vp:.4f} c, theory {w0 / k0:.4f} c")
check("v_phase * v_group = c^2", abs(vp * vg - 1) < 0.02, f"{vp * vg:.4f}")
# Front: real Klein-Gordon field with compact support and zero initial velocity.
a = 20.0


def bump(y):
    out = np.zeros_like(y)
    ins = np.abs(y) < a
    out[ins] = np.exp(1 - 1 / (1 - (y[ins] / a) ** 2))
    return out


psi = bump(x)
G = np.fft.fft(psi)
tf = 300.0
field = np.fft.ifft(G * np.cos(w * tf)).real
outside = np.abs(field[np.abs(x) > a + tf + 5]).max()
inside = np.abs(field[np.abs(x) <= a + tf]).max()
check("disturbance from a compact source stays inside x = a + c t", outside < 1e-6 * inside, f"outside {outside:.1e}, inside {inside:.1e}")

section("Cherenkov threshold")
n_w = 1.33
b_th = 1 / n_w
g_th = 1 / np.sqrt(1 - b_th**2)
T_th = (g_th - 1) * K.m_e * c**2 / K.eV / 1e6
print(f"  In water n = {n_w}: light moves at {1 / n_w:.3f} c, electron threshold kinetic energy {T_th:.3f} MeV")
print(f"  Cone angle for beta -> 1: {np.degrees(np.arccos(1 / n_w)):.1f} deg. Particle still below vacuum c")
check("electron threshold in water near 0.26 MeV", 0.25 < T_th < 0.27)

section("Expansion")
H0 = 67.4e3 / 3.0857e22
print(f"  H0 = 67.4 km/s/Mpc: distance receding at c = {c / H0 / K.ly / 1e9:.1f} Gly")

section("Entanglement cannot signal")
phi = np.array([1, 0, 0, 1], complex) / np.sqrt(2)  # Bell state
rho = np.outer(phi, phi.conj())


def bob_reduced(r):
    r = r.reshape(2, 2, 2, 2)
    return np.einsum("ijik->jk", r)


def rand_unitary(d):
    q, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    return q


worst = 0.0
for _ in range(3000):
    U = np.kron(rand_unitary(2), np.eye(2))
    worst = max(worst, np.abs(bob_reduced(U @ rho @ U.conj().T) - np.eye(2) / 2).max())
check("Bob's state unchanged by any unitary on Alice's particle", worst < 1e-12, f"worst {worst:.1e}")
worst = 0.0
for _ in range(3000):
    # Alice measures in a random basis, outcomes unselected
    V = rand_unitary(2)
    out = np.zeros((4, 4), complex)
    for i in range(2):
        P = np.kron(np.outer(V[:, i], V[:, i].conj()), np.eye(2))
        out += P @ rho @ P
    worst = max(worst, np.abs(bob_reduced(out) - np.eye(2) / 2).max())
check("Bob's state unchanged by any measurement by Alice", worst < 1e-12, f"worst {worst:.1e}")


def corr(a, b):
    """Correlation of spin measurements at angles a and b in the Bell state."""
    A = np.cos(a) * np.array([[1, 0], [0, -1]]) + np.sin(a) * np.array([[0, 1], [1, 0]])
    B = np.cos(b) * np.array([[1, 0], [0, -1]]) + np.sin(b) * np.array([[0, 1], [1, 0]])
    return (phi.conj() @ np.kron(A, B) @ phi).real


S = corr(0, np.pi / 4) + corr(0, -np.pi / 4) + corr(np.pi / 2, np.pi / 4) - corr(np.pi / 2, -np.pi / 4)
check("the same state still gives CHSH = 2 sqrt 2", abs(abs(S) - 2 * np.sqrt(2)) < 1e-12, f"S = {S:.4f}")

section("Superluminal signal leads to a reply before the question")
Lp = 1.0  # distance, in light-years, in the base frame
u_sig = 2.0  # signal speed in c units, in any inertial frame (hypothetical)


def lt(t, x, b):
    g = 1 / np.sqrt(1 - b**2)
    return g * (t - b * x), g * (x - b * t)


def reply_time(v_f):
    """A sends at (0, 0) at speed u. B, moving at v_f, receives at (Lp/u, Lp) and replies
    at speed u in B's own frame toward A. Returns the A-frame event where it reaches x = 0."""
    tb1, xb1 = lt(Lp / u_sig, Lp, v_f)
    # In B's frame A sits at x' = -v t'; the reply is x' = xb1 - u (t' - tb1).
    tpr = (xb1 + u_sig * tb1) / (u_sig - v_f)
    xpr = -v_f * tpr
    ga = 1 / np.sqrt(1 - v_f**2)
    return ga * (tpr + v_f * xpr), ga * (xpr + v_f * tpr)


v_crit = 2 * u_sig / (1 + u_sig**2)
print(f"  Signals at u = {u_sig} c: reply comes back early only if B moves faster than 2u/(1+u^2) = {v_crit:.3f} c")
print("  v/c    time the reply reaches A (yr; A asked at t = 0)")
for v_f in [0.5, v_crit, 0.9, 0.99]:
    print(f"  {v_f:.3f}  {reply_time(v_f)[0]:+.3f}")
t_return, x_return = reply_time(0.9)
check("at v = 0.9c the reply reaches A before A asked", t_return < 0 and abs(x_return) < 1e-12, f"t = {t_return:.3f} yr")
check("at v = 0.5c it doesn't, so the crossover is real", reply_time(0.5)[0] > 0)
check("crossover at v = 2u/(1+u^2) found by the event calculation", abs(reply_time(v_crit)[0]) < 1e-12)

section("Alcubierre metric: sign of the energy density")
X, Y, Z, vs, Gs = sp.symbols("x y z v_s G", real=True)
r = sp.sqrt(X**2 + Y**2 + Z**2)
f = sp.Function("f")(r)
beta_x = -vs * f
coords = (X, Y, Z)
beta_vec = [beta_x, 0, 0]
Kij = sp.Matrix(3, 3, lambda i, j: (sp.diff(beta_vec[i], coords[j]) + sp.diff(beta_vec[j], coords[i])) / 2)
trK = Kij.trace()
KK = sum(Kij[i, j] ** 2 for i in range(3) for j in range(3))
rho_e = sp.simplify((trK**2 - KK) / (16 * sp.pi * Gs))
target = -vs**2 * sp.diff(sp.Function("f")(sp.Symbol("q")), sp.Symbol("q")) ** 2
# evaluate for a concrete f to keep the check clear
sigma, Rb = sp.symbols("sigma R", positive=True)
q = sp.symbols("q", positive=True)
fconc = (sp.tanh(sigma * (q + Rb)) - sp.tanh(sigma * (q - Rb))) / (2 * sp.tanh(sigma * Rb))
fprime = sp.diff(fconc, q)
beta_concrete = -vs * fconc.subs(q, r)
Kij2 = sp.Matrix(3, 3, lambda i, j: (sp.diff([beta_concrete, 0, 0][i], coords[j]) + sp.diff([beta_concrete, 0, 0][j], coords[i])) / 2)
rho2 = (Kij2.trace() ** 2 - sum(Kij2[i, j] ** 2 for i in range(3) for j in range(3))) / (16 * sp.pi * Gs)
expected = -(vs**2) * (Y**2 + Z**2) / (32 * sp.pi * Gs * r**2) * fprime.subs(q, r) ** 2
num = {X: 3.0, Y: 1.2, Z: -0.4, vs: 1.0, Gs: 1.0, sigma: 2.0, Rb: 4.0}
d_val = float((rho2 - expected).subs(num))
val = float(rho2.subs(num))
check("rho = -(v^2/32 pi G)(y^2+z^2)/r^2 f'(r)^2 (from K^2 - K_ij K^ij)", abs(d_val) < 1e-12 * abs(val) + 1e-30, f"sample value {val:.4e}")
samples = rng.normal(size=(200, 3)) * 5
s_fun = sp.lambdify((X, Y, Z), rho2.subs({vs: 1.0, Gs: 1.0, sigma: 2.0, Rb: 4.0}), "numpy")
vals = s_fun(samples[:, 0], samples[:, 1], samples[:, 2])
check("energy density is never positive at 200 random points", np.all(vals <= 1e-14), f"most negative {vals.min():.3e}")


def warp_energy(R, sg_, v_over_c=1.0):
    """Total Eulerian energy (J), integrating over angles: <sin^2 theta> = 2/3."""
    def sech2(z):
        e = np.exp(-2 * abs(z))
        return 4 * e / (1 + e) ** 2  # safe form of 1/cosh^2

    def fp(q_):
        return (sg_ / (2 * np.tanh(sg_ * R))) * (sech2(sg_ * (q_ + R)) - sech2(sg_ * (q_ - R)))
    I = quad(lambda q_: q_**2 * fp(q_) ** 2, 0, R + 40 / sg_, points=[R], limit=500)[0]
    return -(c**4 / K.G) * v_over_c**2 / (32 * np.pi) * (4 * np.pi * 2 / 3) * I


E_w = warp_energy(100.0, 1.0)
print(f"  R = 100 m bubble, wall ~ 1 m, v = c: E = {E_w:.2e} J, i.e. {abs(E_w) / c**2 / K.M_sun:.2f} M_sun of negative energy")
check("matches the thin-wall estimate -(c^4/G)(v/c)^2 R^2 sigma/36", abs(E_w / (-(c**4 / K.G) * 100.0**2 * 1.0 / 36) - 1) < 0.02)
E_w2 = warp_energy(100.0, 100.0)
print(f"  With a 1 cm wall: {abs(E_w2) / c**2 / K.M_sun:.1f} M_sun. Thinner walls cost more")

finish()
