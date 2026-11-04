"""F05: centre of mass of a donut, by three independent routes, plus what the point means.

Routes: exact integration (sympy), Monte Carlo over the solid, and a simulated tossed ring
of particles whose centre of mass must follow a parabola though no particle is there.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.spatial import ConvexHull, Delaunay

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
fails = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        fails.append(name)


Rv, rv = 0.04, 0.015  # donut: ring radius R = 4 cm, tube radius r = 1.5 cm

# 1. Exact integration in toroidal coordinates with sympy.
R, r, s, phi, th = sp.symbols("R r s phi theta", positive=True)
rho_cyl = R + s * sp.cos(phi)
jac = s * rho_cyl
x = rho_cyl * sp.cos(th)
y = rho_cyl * sp.sin(th)
z = s * sp.sin(phi)


def integ(f):
    return sp.simplify(sp.integrate(sp.expand(f * jac), (th, 0, 2 * sp.pi), (phi, 0, 2 * sp.pi), (s, 0, r)))


V = integ(1)
print("    exact volume V =", V)
com_exact = [integ(q) for q in (x, y, z)]
print("    exact integrals of x, y, z dV =", com_exact)
check("volume is 2 pi^2 R r^2", sp.simplify(V - 2 * sp.pi ** 2 * R * r ** 2) == 0)
check("first moments vanish: COM at the geometric centre", all(c == 0 for c in com_exact))

# 2. Monte Carlo over a bounding box, knowing nothing of the symmetry.
rng = np.random.default_rng(5)
N = 2_000_000
box = Rv + rv
pts = rng.uniform([-box, -box, -rv], [box, box, rv], size=(N, 3))


def in_torus(p):
    return (np.hypot(p[..., 0], p[..., 1]) - Rv) ** 2 + p[..., 2] ** 2 <= rv ** 2


inside = pts[in_torus(pts)]
vol_mc = (2 * box) ** 2 * (2 * rv) * len(inside) / N
com_mc = inside.mean(axis=0)
err = inside.std(axis=0) / np.sqrt(len(inside))
print(f"    Monte Carlo: {len(inside)} points, V = {vol_mc * 1e6:.3f} cm^3 "
      f"(exact {2 * np.pi ** 2 * Rv * rv ** 2 * 1e6:.3f}), COM = {np.round(com_mc * 1e3, 3)} mm "
      f"+- {np.round(err * 1e3, 3)} mm")
check("Monte Carlo COM within 4 standard errors of the centre", np.all(np.abs(com_mc) < 4 * err))
check("Monte Carlo volume within 0.2%", abs(vol_mc / (2 * np.pi ** 2 * Rv * rv ** 2) - 1) < 2e-3)
check("the COM point itself is not in the material", not in_torus(np.zeros(3)))

# 3. COM lies in the convex hull but need not lie in the body.
sub = inside[rng.choice(len(inside), 4000, replace=False)]
hull = Delaunay(sub[ConvexHull(sub).vertices])
check("COM lies inside the convex hull of the donut", hull.find_simplex(com_mc) >= 0)

# Half-ring (a thin C shape): exact COM is 2R/pi from the centre, on the empty side of the wire.
theta = np.linspace(-np.pi / 2, np.pi / 2, 200001)
xc = np.trapezoid(Rv * np.cos(theta), theta) / np.pi
check("half ring COM at 2R/pi, off the wire", abs(xc - 2 * Rv / np.pi) < 1e-9,
      f"{xc * 100:.3f} cm from the centre; the wire is {Rv * 100:.1f} cm away at that angle")

# 4. Non-uniform donut: icing on top moves the COM up the axis, but it stays in the hole.
icing = 1.5  # top half 1.5 times denser
w = np.where(inside[:, 2] > 0, icing, 1.0)
com_ice = (inside * w[:, None]).sum(axis=0) / w.sum()
print(f"    iced donut COM = {np.round(com_ice * 1e3, 2)} mm")
check("icing shifts COM along the axis only, still not in material",
      com_ice[2] > 1e-3 and np.hypot(*com_ice[:2]) < 1e-4 and not in_torus(com_ice))

# 5. Toss a ring of particles joined by springs, spinning and wobbling, and track the COM.
n = 72
ang = np.arange(n) * 2 * np.pi / n
pos = np.column_stack((Rv * np.cos(ang), Rv * np.sin(ang)))
m = np.ones(n) * (0.06 / n)
m[:12] *= 3.0  # lopsided, so the COM is not the geometric centre
com0 = (m[:, None] * pos).sum(0) / m.sum()
v0 = np.array([1.5, 3.0])
spin = 25.0
rel = pos - com0
vel = v0 + spin * np.column_stack((-rel[:, 1], rel[:, 0]))
vel += rng.normal(0, 0.02, vel.shape)  # internal wobble
vel -= (m[:, None] * vel).sum(0) / m.sum() - v0
k = 400.0
g = np.array([0.0, -9.81])


# Massless cross-links keep the ring's shape; they hold no mass near the centre.
links = (1, 6, 18)
rests = {q: np.linalg.norm(pos - np.roll(pos, -q, axis=0), axis=1) for q in links}


def acc(p):
    a_ = np.zeros_like(p)
    for q in links:
        d = np.roll(p, -q, axis=0) - p
        L = np.linalg.norm(d, axis=1)
        f = (k * (L - rests[q]) / L)[:, None] * d
        a_ += f - np.roll(f, q, axis=0)
    return a_ / m[:, None] + g


dt, T, steps = 4e-5, 0.5, 12500
a = acc(pos)
for _ in range(steps):
    vel += 0.5 * dt * a
    pos += dt * vel
    a = acc(pos)
    vel += 0.5 * dt * a
com_T = (m[:, None] * pos).sum(0) / m.sum()
pred = com0 + v0 * T + 0.5 * g * T ** 2
nearest = np.linalg.norm(pos - com_T, axis=1).min()
check("tossed lopsided ring: COM follows the parabola to under 1 micron",
      np.linalg.norm(com_T - pred) < 1e-6, f"miss {np.linalg.norm(com_T - pred):.1e} m")
check("no particle sits at the COM at the end of the toss", nearest > 0.02,
      f"nearest particle {nearest * 100:.2f} cm away")

# 6. Hang the uniform ring from a nail on the inner rim: the rest position puts the centre below.
beta = np.linspace(-np.pi, np.pi, 3601)  # tilt of the pivot-to-centre line from straight down
ring = np.column_stack((Rv * np.cos(ang), Rv * np.sin(ang))) - [Rv, 0.0]  # pivot at origin


def pot(b):
    c, s_ = np.cos(b), np.sin(b)
    # Rotate so that at b=0 the ring's centre lies straight below the nail, then sum m g y.
    rotx = c * ring[:, 1] + s_ * ring[:, 0]
    roty = c * ring[:, 0] - s_ * ring[:, 1]
    return (m.mean() * 9.81 * roty).sum(), rotx


U = np.array([pot(b)[0] for b in beta])
b_min = beta[np.argmin(U)]
print(f"    ring on a nail, potential minimum at tilt {np.degrees(b_min):.2f} deg from straight down")
check("hanging ring rests with its empty centre directly below the nail", abs(b_min) < 2e-3)

# 7. The COM is not where gravity 'lives': field of the donut at its own centre vs far away.
G = 6.674e-11
Mass = 0.06
ms = Mass / len(inside)
far = np.array([3.0, 0.0, 0.0])
d_far = far - inside
F_far = -(G * ms * d_far / np.linalg.norm(d_far, axis=1)[:, None] ** 3).sum(0)
d_c = -inside
F_c = -(G * ms * d_c / np.linalg.norm(d_c, axis=1)[:, None] ** 3).sum(0)
scale = G * Mass / (Rv ** 2)
print(f"    field at centre / G M R^-2 = {np.linalg.norm(F_c) / scale:.1e}; "
      f"field at 3 m / (G M / 9) = {np.linalg.norm(F_far) / (G * Mass / 9):.5f}")
check("field at the donut's centre is near zero, far field looks like a point at the COM",
      np.linalg.norm(F_c) / scale < 0.01 and abs(np.linalg.norm(F_far) / (G * Mass / 9) - 1) < 1e-3)

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All F05 checks passed.")
