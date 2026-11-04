"""L05-L07: gravitational field, shell theorem, centre of mass and a rattling ball.

Direct sums of Newtonian point masses are compared with the shell theorem and with
the multipole expansion. A rattling ball is simulated by integrating each body
separately with impulses at contact, so the COM motion is an output, not an input.
"""
import numpy as np

from checklib import K, check, finish, section

G = K.G
rng = np.random.default_rng(7)


def fib_sphere(n):
    i = np.arange(n) + 0.5
    z = 1 - 2 * i / n
    phi = np.pi * (1 + 5**0.5) * i
    s = np.sqrt(1 - z**2)
    return np.column_stack([s * np.cos(phi), s * np.sin(phi), z])


def field_at(points, masses, x):
    """Newtonian field at x from a set of point masses."""
    d = points - x
    r = np.linalg.norm(d, axis=1)
    return G * ((masses / r**3)[:, None] * d).sum(axis=0)


section("Shell theorem by direct summation")
M = 1.0
R = 0.1
Ns = 400000
shell = R * fib_sphere(Ns)
ms = np.full(Ns, M / Ns)
for f in [0.0, 0.3, 0.6]:
    gv = field_at(shell, ms, np.array([0, 0, f * R]))
    ref = G * M / R**2
    check(f"inside a shell at r = {f} R: field vanishes", np.linalg.norm(gv) / ref < 1e-4, f"|g|/(GM/R^2) = {np.linalg.norm(gv) / ref:.1e}")
for f in [1.5, 3.0, 10.0]:
    r = f * R
    gv = field_at(shell, ms, np.array([0, 0, r]))
    check(f"outside a shell at r = {f} R: |g| = GM/r^2", abs(np.linalg.norm(gv) / (G * M / r**2) - 1) < 1e-4,
          f"ratio {np.linalg.norm(gv) / (G * M / r**2):.7f}")

section("Solid uniform ball: Gauss's law check")
# Point masses on a cubic lattice; random points are too noisy for a field probe inside.
h = R / 50
axis = (np.arange(-51, 51) + 0.5) * h
gx, gy, gz = np.meshgrid(axis, axis, axis, indexing="ij")
lat = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
lat = lat[np.linalg.norm(lat, axis=1) <= R]
print(f"  {len(lat)} lattice points, spacing R/50")
for f in [0.25, 0.5, 0.8, 2.0]:
    # Probe at a cell centre (integer multiple of h), equidistant from its 8 nearest points,
    # so that the graininess of nearby point masses cancels by symmetry.
    probe = np.round(np.array([1, 2, 2]) / 3 * f * R / h) * h
    r = np.linalg.norm(probe)
    pred = G * M * r / R**3 if f < 1 else G * M / r**2
    gv = field_at(lat, np.full(len(lat), M / len(lat)), probe)
    ratio = -(gv @ probe / r) / pred
    check(f"uniform ball at r = {f} R: inward field matches Gauss's law", abs(ratio - 1) < 5e-3, f"ratio {ratio:.5f}")

section("L06: the field is not the field of a point at the COM")
m_each = 0.5
sep = 1.0  # m, dumbbell with masses at +- sep/2 on the x axis
db = np.array([[-sep / 2, 0, 0], [sep / 2, 0, 0]])
dm = np.array([m_each, m_each])
for point, label in [(np.array([0.0, 0, 0]), "at the COM itself"),
                     (np.array([1.0, 0, 0]), "on the axis at 1 m"),
                     (np.array([0.0, 1.0, 0]), "broadside at 1 m"),
                     (np.array([0.7, 0.7, 0]), "diagonal at 1 m")]:
    gv = field_at(db, dm, point)
    r = np.linalg.norm(point)
    if r == 0:
        check("dumbbell field at its COM is zero, not infinite", np.linalg.norm(gv) < 1e-20, f"|g| = {np.linalg.norm(gv):.1e}")
        continue
    gm = G * 1.0 / r**2
    angle = np.degrees(np.arccos(np.clip(-(gv @ point) / (np.linalg.norm(gv) * r), -1, 1)))
    print(f"  {label}: |g|/(GM/r^2) = {np.linalg.norm(gv) / gm:.4f}, field misses the COM direction by {angle:.2f} deg")
gv_far = field_at(db, dm, np.array([0.0, 100.0, 0]))
check("at 100 separations the point-at-COM picture is good to 1e-4",
      abs(np.linalg.norm(gv_far) / (G * 1.0 / 100**2) - 1) < 1e-4)

section("L07: multipole expansion against direct summation")
M_s, m_i = 1.0, 0.1  # outer shell and inner ball, kg
R_s, r_i = 0.10, 0.02
d_vec = np.array([0.05, 0.03, -0.02])  # inner ball offset from shell centre, m
Ntot = M_s + m_i
com = m_i * d_vec / Ntot


def monopole(x):
    return -G * Ntot * x / np.linalg.norm(x) ** 3


def quad_moment(centres, masses):
    Q = np.zeros((3, 3))
    for p, mi in zip(centres, masses):
        Q += mi * (3 * np.outer(p, p) - (p @ p) * np.eye(3))
    return Q


def quad_field(x, Q):
    """g from the l=2 term of the potential Phi = -G/(2 r^5) x.Q.x"""
    r = np.linalg.norm(x)
    return G * (Q @ x / r**5 - 2.5 * (x @ Q @ x) * x / r**7)


c_s = -com  # shell centre relative to the total COM
c_i = d_vec - com
Q = quad_moment([c_s, c_i], [M_s, m_i])
dip = M_s * c_s + m_i * c_i
check("dipole about the total COM is zero", np.linalg.norm(dip) < 1e-16, f"{np.linalg.norm(dip):.1e} kg m")
check("quadrupole is traceless", abs(np.trace(Q)) < 1e-15)
mu = M_s * m_i / Ntot
eig = np.linalg.eigvalsh(Q)
print(f"  reduced mass {mu:.4f} kg, d = {np.linalg.norm(d_vec) * 100:.2f} cm, Q eigenvalues {np.round(eig * 1e3, 4)} x 1e-3 kg m^2")
check("largest Q eigenvalue = 2 mu d^2 (two-point analytic)", abs(eig[-1] / (2 * mu * d_vec @ d_vec) - 1) < 1e-10)

# Bodies built of many points, so that the direct sum is not a 2-point sum in disguise.
N1, N2 = 200000, 100000
P1 = R_s * fib_sphere(N1) + c_s
P2 = r_i * fib_sphere(N2) * rng.random(N2)[:, None] ** (1 / 3) + c_i
pts = np.vstack([P1, P2])
mass = np.concatenate([np.full(N1, M_s / N1), np.full(N2, m_i / N2)])
print("  r (m)  |g_direct|     residual after monopole   residual after monopole+quadrupole")
rs = np.array([0.4, 0.6, 1.0, 1.6, 2.5])
dirn = np.array([0.6, -0.48, 0.64])
res1, res2, res_d = [], [], []
for r in rs:
    x = r * dirn
    gd = field_at(pts, mass, x)
    e1 = gd - monopole(x)
    e2 = e1 - quad_field(x, Q)
    # Same field, but with the monopole wrongly placed at the shell centre
    e_d = gd - monopole(x - c_s)
    res1.append(np.linalg.norm(e1))
    res2.append(np.linalg.norm(e2))
    res_d.append(np.linalg.norm(e_d))
    print(f"  {r:4.1f}   {np.linalg.norm(gd):.5e}    {np.linalg.norm(e1):.4e}              {np.linalg.norm(e2):.4e}")
slope1 = np.polyfit(np.log(rs), np.log(res1), 1)[0]
slope2 = np.polyfit(np.log(rs), np.log(res2), 1)[0]
slope_d = np.polyfit(np.log(rs), np.log(res_d), 1)[0]
check("after removing the monopole the residual falls as r^-4 (quadrupole)", abs(slope1 + 4) < 0.15, f"fitted power {slope1:.3f}")
check("after removing the quadrupole too it falls as r^-5 (octupole)", abs(slope2 + 5) < 0.3, f"fitted power {slope2:.3f}")
check("quadrupole removes most of the residual at r = 1 m", res2[2] < 0.2 * res1[2], f"{res2[2] / res1[2]:.3f} remaining")
check("expanding about the shell centre instead leaves a dipole, r^-3", abs(slope_d + 3) < 0.15,
      f"fitted power {slope_d:.3f}, size {res_d[2]:.3e} m/s^2 at 1 m, dipole m d = {m_i * np.linalg.norm(d_vec):.2e} kg m")

section("L07: a thrown ball with a rattling ball inside")


def rattle(T=2.0, gz=-9.80665, seed=11, e=1.0):
    """Shell (centre X, vel V) and inner ball (x, v) in uniform gravity, separate
    integration of each, with an impulse at each contact conserving momentum."""
    r = np.random.default_rng(seed)
    X = np.zeros(3)
    V = np.array([8.0, 0.0, 8.0])
    x = np.array([0.03, 0.0, 0.01])
    v = V + r.normal(0, 1.5, 3)
    gv = np.array([0, 0, gz])
    gap = R_s - r_i
    t = 0.0
    log = []
    n_hits = 0
    dt_log = 0.005
    next_log = 0.0
    while t < T:
        rel, vr = x - X, v - V
        A = vr @ vr
        B = 2 * rel @ vr
        C = rel @ rel - gap**2
        s = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
        s = max(s, 1e-12)
        while next_log <= min(t + s, T) and next_log <= T:
            h = next_log - t
            Xn = X + V * h + 0.5 * gv * h**2
            xn = x + v * h + 0.5 * gv * h**2
            log.append((next_log, Xn, xn))
            next_log += dt_log
        if t + s > T:
            break
        X = X + V * s + 0.5 * gv * s**2
        x = x + v * s + 0.5 * gv * s**2
        V = V + gv * s
        v = v + gv * s
        t += s
        rel = x - X
        n = rel / np.linalg.norm(rel)
        vn = (v - V) @ n
        Jmag = (1 + e) * vn * M_s * m_i / (M_s + m_i)
        v = v - Jmag / m_i * n
        V = V + Jmag / M_s * n
        n_hits += 1
    return log, n_hits, np.array([8.0, 0.0, 8.0]), gv


log, hits, V0, gv = rattle()
ts = np.array([q[0] for q in log])
Xs = np.array([q[1] for q in log])
xs = np.array([q[2] for q in log])
coms = (M_s * Xs + m_i * xs) / Ntot
v_com0 = (coms[1] - coms[0] - 0.5 * gv * (ts[1] ** 2 - ts[0] ** 2)) / (ts[1] - ts[0])
para = coms[0] + np.outer(ts, v_com0) + 0.5 * np.outer(ts**2, gv)
dev = np.linalg.norm(coms - para, axis=1).max()
check("COM stays on one parabola through all the collisions", dev < 1e-9, f"{hits} collisions, worst deviation {dev:.1e} m")
print(f"  shell centre wiggles about the COM by up to {np.linalg.norm(Xs - coms, axis=1).max() * 1e3:.2f} mm")
# field at a point 1 m from the COM, direction fixed, as a function of time
probe = np.array([1.0, 0, 0])
g_all = []
for X, x_ in zip(Xs, xs):
    com_now = (M_s * X + m_i * x_) / Ntot
    p = com_now + probe
    g_all.append(field_at(np.array([X, x_]), np.array([M_s, m_i]), p))
g_all = np.array(g_all)
g_mag = np.linalg.norm(g_all, axis=1)
g_mono = G * Ntot
dg = np.linalg.norm(g_all - g_all.mean(axis=0), axis=1).max()
print(f"  At 1 m from the COM: monopole field {g_mono:.4e} m/s^2")
print(f"  Variation from rattling {dg:.3e} m/s^2, fraction {dg / g_mono:.2e}")
gap = R_s - r_i
Qmax = 3 * G * mu * gap**2 / 1.0**4
check("rattling variation is of order the quadrupole estimate 3 G mu d^2 / r^4", dg < 2 * Qmax and dg > 1e-3 * Qmax,
      f"estimate {Qmax:.2e} m/s^2")
check("rattling changes the field by much less than the steady monopole", dg / g_mono < 5e-3)

section("Propagation and radiation")
print(f"  Field change reaches 1 m after {1 / K.c * 1e9:.2f} ns (retardation at c, GR; Newton's formula has none)")
d, om = gap, 2 * np.pi * 10
Qdd = mu * d**2 * om**2
Qddd = Qdd * om
h = 2 * G / K.c**4 * Qdd / 1.0
Pgw = G / (5 * K.c**5) * Qddd**2
print(f"  Quadrupole radiation for a 10 Hz rattle: strain at 1 m ~ {h:.1e}, power ~ {Pgw:.1e} W (order of magnitude only)")
check("gravitational-wave power from the rattle is absurdly small", Pgw < 1e-40)

finish()
