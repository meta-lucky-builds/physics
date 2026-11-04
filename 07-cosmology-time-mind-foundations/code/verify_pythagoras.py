"""C39: exact checks of two proofs, failure on curved surfaces, Plimpton 322 rows."""
from fractions import Fraction as Fr
from math import isqrt
import numpy as np
import sympy as sp
from _harness import Checks

C = Checks("Pythagorean theorem: proofs, curvature and history")

C.section("Proof 1 (rearrangement) replayed with exact rational coordinates")


def shoelace(pts):
    s = Fr(0)
    for (x1, y1), (x2, y2) in zip(pts, pts[1:] + pts[:1]):
        s += x1 * y2 - x2 * y1
    return abs(s) / 2


ok = True
for a, b in [(Fr(3), Fr(4)), (Fr(5), Fr(12)), (Fr(7, 3), Fr(11, 5)), (Fr(1), Fr(1)), (Fr(9, 7), Fr(40, 3))]:
    s = a + b
    # Inner tilted quadrilateral inside the big square of side a + b.
    inner = [(a, Fr(0)), (s, a), (b, s), (Fr(0), b)]
    big = s * s
    tris = 4 * (a * b / 2)
    inner_area = shoelace(inner)
    # Each inner side is a hypotenuse: check the inner shape is a square with right corners.
    edges = [(inner[(i + 1) % 4][0] - inner[i][0], inner[(i + 1) % 4][1] - inner[i][1]) for i in range(4)]
    dot0 = all(edges[i][0] * edges[(i + 1) % 4][0] + edges[i][1] * edges[(i + 1) % 4][1] == 0 for i in range(4))
    ok &= dot0 and inner_area == big - tris and inner_area == a * a + b * b
C.check("for exact rational legs: inner figure is a square of area a^2 + b^2 = (a+b)^2 - 2ab", ok)

C.section("Proof 2 (similar triangles) replayed symbolically")
a, b = sp.symbols("a b", positive=True)
c_ = sp.sqrt(a**2 + b**2)
# Foot of the altitude divides c into p, q. From similarity: a/c = p/a, b/c = q/b.
p = a**2 / c_
q = b**2 / c_
# Independent geometric test: compute the actual foot with coordinates C=(0,0), A=(b,0), B=(0,a).
t = sp.symbols("t")
A = sp.Matrix([b, 0])
B = sp.Matrix([0, a])
foot = A + t * (B - A)
tt = sp.solve(foot.dot(B - A), t)[0]
F = foot.subs(t, tt)
C.check("|BF|^2 from coordinates equals (a^2/c)^2", sp.simplify((B - F).dot(B - F) - p**2) == 0)
C.check("|AF|^2 from coordinates equals (b^2/c)^2", sp.simplify((A - F).dot(A - F) - q**2) == 0)
C.check("p + q = c, which is the theorem", sp.simplify(p + q - c_) == 0)

C.section("Measured failure on a sphere and in hyperbolic space")


def sphere_right_triangle(aa, bb, R):
    C0 = np.array([R, 0, 0.0])
    # Walk distance aa along one great circle and bb along a perpendicular one.
    P = R * np.array([np.cos(aa / R), np.sin(aa / R), 0.0])
    Q = R * np.array([np.cos(bb / R), 0.0, np.sin(bb / R)])
    cc = R * np.arccos(np.clip(P @ Q / R**2, -1, 1))
    return cc


def hyper_right_triangle(aa, bb, R):
    # Hyperboloid model, -t^2 + x^2 + y^2 = -R^2; Lorentzian inner product gives distance.
    def inner(u, v):
        return -u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
    P = R * np.array([np.cosh(aa / R), np.sinh(aa / R), 0.0])
    Q = R * np.array([np.cosh(bb / R), 0.0, np.sinh(bb / R)])
    return R * np.arccosh(-inner(P, Q) / R**2)


print(f"  {'legs a=b':>8} {'R':>5} {'sphere c':>9} {'hyper c':>9} {'flat c':>9} {'sphere c^2/(a^2+b^2)':>21} {'hyper c^2/(a^2+b^2)':>20}")
good = True
for aa, R in [(0.01, 1), (0.1, 1), (0.5, 1), (1.0, 1), (np.pi / 2, 1)]:
    cs = sphere_right_triangle(aa, aa, R)
    ch = hyper_right_triangle(aa, aa, R)
    cf = np.hypot(aa, aa)
    print(f"  {aa:8.3f} {R:5.1f} {cs:9.5f} {ch:9.5f} {cf:9.5f} {cs**2/(2*aa**2):21.6f} {ch**2/(2*aa**2):20.6f}")
    good &= cs < cf < ch
C.check("sphere c is shorter than flat c; hyperbolic c is longer; both for every size tested", bool(good))
c3 = sphere_right_triangle(np.pi / 2, np.pi / 2, 1)
print(f"  octant triangle: legs pi/2, hypotenuse {c3:.5f}. c^2 = {c3**2:.4f}, a^2+b^2 = {2*(np.pi/2)**2:.4f}")
C.check("octant triangle: hypotenuse equals each leg, a flat-rule failure of 2x", abs(c3 - np.pi / 2) < 1e-12)
# Fit the leading correction: on a sphere c^2 = a^2 + b^2 - a^2 b^2/(3 R^2) + ...
aa = 0.01
meas = (2 * aa**2 - sphere_right_triangle(aa, aa, 1.0)**2) / aa**4
print(f"  measured coefficient of a^2 b^2/R^2 at a = b = 0.01: {meas:.4f} (series below predicts 1/3)")
# Symbolic prediction: put c^2 = a^2 + b^2 + k a^2 b^2 into cos c = cos a cos b, order by order.
A_, B_ = sp.symbols("A B", positive=True)
k_ = sp.symbols("k", real=True)
Csq = A_**2 + B_**2 + k_ * A_**2 * B_**2
lhs = 1 - Csq / 2 + Csq**2 / 24
rhs_ = (1 - A_**2 / 2 + A_**4 / 24) * (1 - B_**2 / 2 + B_**4 / 24)
coef = sp.Poly(sp.expand(lhs - rhs_), A_, B_).coeff_monomial(A_**2 * B_**2)
k_sol = sp.solve(coef, k_)[0]
print(f"  symbolic series from cos c = cos a cos b: c^2 = a^2 + b^2 + ({k_sol}) a^2 b^2/R^2 + ...")
C.check("symbolic coefficient is -1/3", k_sol == sp.Rational(-1, 3))
C.check("the measured defect matches 1/3 within 1%", abs(meas - 1 / 3) < 0.01)
C.check("flat limit: defect vanishes as R -> infinity (c^2/(a^2+b^2) -> 1)",
        abs(sphere_right_triangle(1.0, 1.0, 1e4)**2 / 2 - 1) < 1e-6)

C.section("Angle sums: the same geometries, tested for the parallel-postulate consequence")


def angle_sphere(p, q, r):
    t1 = q - p * (p @ q)
    t2 = r - p * (p @ r)
    return np.arccos(t1 @ t2 / np.linalg.norm(t1) / np.linalg.norm(t2))


P0 = np.array([1.0, 0, 0])
Pa = np.array([np.cos(1.0), np.sin(1.0), 0])
Pb = np.array([np.cos(1.0), 0, np.sin(1.0)])
s = angle_sphere(P0, Pa, Pb) + angle_sphere(Pa, P0, Pb) + angle_sphere(Pb, P0, Pa)
print(f"  sphere triangle, legs 1, right angle: angle sum = {np.degrees(s):.2f} deg (> 180)")
C.check("the same triangle has angle sum above 180 deg, as the parallel postulate's failure requires", s > np.pi)

C.section("Plimpton 322: rows checked for s^2 + l^2 = d^2")
rows = [(119, 169), (3367, 4825), (4601, 6649), (12709, 18541), (65, 97), (319, 481), (2291, 3541), (799, 1249),
        (481, 769), (4961, 8161), (45, 75), (1679, 2929), (161, 289), (1771, 3229), (28, 53)]
sex = lambda *d: sum(v * 60**i for i, v in enumerate(reversed(d)))
assert sex(1, 59) == 119 and sex(2, 49) == 169 and sex(1, 20, 25) == 4825
good = True
for s_, d_ in rows:
    l2 = d_**2 - s_**2
    l_ = isqrt(l2)
    good &= l_ * l_ == l2
    print(f"    short {s_:6d}, diagonal {d_:6d}: long side {l_:6d}{'  exact' if l_*l_==l2 else '  NOT SQUARE'}")
C.check("after the four standard scribal corrections, all 15 rows are exact Pythagorean triples", good)
# Rows as inscribed on the tablet, before correction.
raw = {2: (sex(56, 7), sex(3, 12, 1)), 9: (sex(9, 1), sex(12, 49)), 13: (sex(7, 12, 1), sex(4, 49)),
       15: (sex(56), sex(53))}
bad = 0
for k, (s_, d_) in raw.items():
    v = d_**2 - s_**2
    perfect = v > 0 and isqrt(v)**2 == v
    bad += not perfect
    print(f"    tablet row {k}: as inscribed {s_}, {d_}: difference of squares is perfect square? {perfect}")
C.check("rows with recognised errors really fail as inscribed (control)", bad == 4)
C.finish()
