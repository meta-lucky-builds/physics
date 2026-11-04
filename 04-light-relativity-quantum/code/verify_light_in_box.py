"""L01: does light in a box add mass, and what would your hand notice?

Tests, each independent of the formula it checks:
1. Invariant mass of photon sets from raw (E, p) sums, in two frames.
2. A photon bouncing in an accelerated box, with mirror Doppler shifts applied at each
   reflection: the force the box feels gives inertia E/c^2.
3. Monte Carlo ray tracing of light lifetime in a mirrored cube, with and without an
   absorbing "hand", against the mean-chord formula used for the real numbers.
Then the real numbers for a million flashlights.
"""
import numpy as np

from checklib import K, check, finish, section

c = K.c
rng = np.random.default_rng(1)


def invariant_mass(E, p):
    """M from total energy and total momentum vector of a set of photons."""
    Et = E.sum()
    Pt = p.sum(axis=0)
    return np.sqrt(max(Et**2 - c**2 * Pt @ Pt, 0.0)) / c**2


def random_dirs(n):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def boost_photons(E, p, beta):
    """Lorentz boost along z of photon energies and momenta."""
    g = 1 / np.sqrt(1 - beta**2)
    Ep = g * (E - beta * c * p[:, 2])
    pp = p.copy()
    pp[:, 2] = g * (p[:, 2] - beta * E / c)
    return Ep, pp


section("1. Invariant mass of photons (each photon has E = pc, so m = 0 alone)")
E0 = 3e-19  # about a 660 nm photon, J
one = invariant_mass(np.array([E0]), np.array([[0, 0, E0 / c]]))
check("one photon has zero invariant mass", one < 1e-6 * E0 / c**2, f"M = {one:.2e} kg")
same = invariant_mass(np.array([E0, E0]), np.array([[0, 0, E0 / c], [0, 0, E0 / c]]))
check("two parallel photons have zero invariant mass", same < 1e-6 * E0 / c**2, f"M = {same:.2e} kg")
back = invariant_mass(np.array([E0, E0]), np.array([[0, 0, E0 / c], [0, 0, -E0 / c]]))
check("two opposite photons: M = 2E/c^2", abs(back / (2 * E0 / c**2) - 1) < 1e-12, f"M = {back:.4e} kg")
for theta in [0.3, 1.0, 2.0]:
    E1, E2 = 2e-19, 5e-19
    p = np.array([[0, 0, E1 / c], [E2 / c * np.sin(theta), 0, E2 / c * np.cos(theta)]])
    M = invariant_mass(np.array([E1, E2]), p)
    formula = np.sqrt(2 * E1 * E2 * (1 - np.cos(theta))) / c**2
    check(f"two photons at {theta} rad match sqrt(2 E1 E2 (1-cos))/c^2", abs(M / formula - 1) < 1e-9)
N = 200000
E = rng.uniform(1e-19, 5e-19, N)
p = (E / c)[:, None] * random_dirs(N)
M_gas = invariant_mass(E, p)
ratio = M_gas / (E.sum() / c**2)
check("isotropic photon gas: M -> (sum E)/c^2", abs(ratio - 1) < 1e-4, f"M/(sum E/c^2) = {ratio:.8f} for {N} photons")
Eb, pb = boost_photons(E, p, 0.9)
check("M does not change under a boost at 0.9c", abs(invariant_mass(Eb, pb) / M_gas - 1) < 1e-8,
      f"energy changed by factor {Eb.sum() / E.sum():.3f}, M by {invariant_mass(Eb, pb) / M_gas - 1:.1e}")

section("2. A photon's inertia, from mirror kicks in an accelerated box")


def bouncing_photon_force(a, L=1.0, cycles=200, E_start=1.0, c_=1.0):
    """Event-driven photon between a floor at a t^2/2 and a ceiling at L + a t^2/2.

    Each reflection applies the moving-mirror Doppler factor and hands the box the
    momentum change. Returns mean force on box and time-mean photon energy.
    """
    t = 0.0
    x = 0.0
    E = E_start
    impulse = 0.0
    E_t = 0.0
    t0 = 0.0
    for _ in range(cycles):
        # up: x + c(s-t) = L + a s^2/2; smaller positive root
        A, B, C = a / 2, -c_, L + c_ * t - x
        s = (-B - np.sqrt(B**2 - 4 * A * C)) / (2 * A)
        E_t += E * (s - t)
        v = a * s
        Enew = E * (1 - v / c_) / (1 + v / c_)  # receding mirror redshifts
        impulse += (E + Enew) / c_  # ceiling is pushed up
        E, t, x = Enew, s, L + a * s**2 / 2
        # down: x - c(s-t) = a s^2/2
        A, B, C = a / 2, c_, -(x + c_ * t)
        s = (-B + np.sqrt(B**2 - 4 * A * C)) / (2 * A)
        E_t += E * (s - t)
        v = a * s
        Enew = E * (1 + v / c_) / (1 - v / c_)  # approaching mirror blueshifts
        impulse -= (E + Enew) / c_  # floor is pushed down
        E, t, x = Enew, s, a * s**2 / 2
    return impulse / (t - t0), E_t / (t - t0)


for a in [1e-5, 1e-6, 1e-7]:
    F, Ebar = bouncing_photon_force(a)
    r = -F / (Ebar * a)  # c = 1 units: expect F = -(E/c^2) a
    check(f"mean reaction force = -(E/c^2) a at a L/c^2 = {a:g}", abs(r - 1) < 2e-3, f"ratio {r:.5f}")
print("  By the equivalence principle the same force is the photon's weight (E/c^2) g on a scale.")
g = K.g_0
h = 1.0
print(f"  Check in gravity terms: redshift over 1 m is g h/c^2 = {g * h / c**2:.3e}; "
      "floor kick minus ceiling kick over one round trip gives E g/c^2.")

section("3. Ray-traced light lifetime in a mirrored cube")


def trace_cube(R, L=1.0, sphere_radius=0.0, n=40000, max_bounce=200000):
    """Mean path length before absorption. Diffuse walls with reflectivity R; optional
    black sphere at the centre. Returns mean path in units of metres."""
    pos = rng.uniform(0, L, (n, 3))
    if sphere_radius:
        d = np.linalg.norm(pos - L / 2, axis=1)
        pos = pos[d > sphere_radius]
        n = len(pos)
    u = random_dirs(n)
    path = np.zeros(n)
    alive = np.ones(n, bool)
    centre = np.full(3, L / 2)
    for _ in range(max_bounce):
        idx = np.nonzero(alive)[0]
        if idx.size == 0:
            break
        P, U = pos[idx], u[idx]
        with np.errstate(divide="ignore", invalid="ignore"):
            tt = np.where(U > 0, (L - P) / U, np.where(U < 0, -P / U, np.inf))
        axis = tt.argmin(axis=1)
        s = tt[np.arange(idx.size), axis]
        hit_sphere = np.zeros(idx.size, bool)
        if sphere_radius:
            oc = P - centre
            b = (oc * U).sum(axis=1)
            disc = b**2 - ((oc**2).sum(axis=1) - sphere_radius**2)
            ok = disc > 0
            s2 = np.where(ok, -b - np.sqrt(np.where(ok, disc, 0)), np.inf)
            hit_sphere = (s2 > 0) & (s2 < s)
            s = np.where(hit_sphere, s2, s)
        path[idx] += s
        newP = P + U * s[:, None]
        # walls absorb with probability 1-R
        survive = (~hit_sphere) & (rng.random(idx.size) < R)
        alive[idx[~survive]] = False
        keep = np.nonzero(survive)[0]
        # Lambertian reflection off the wall that was hit
        W = random_dirs(keep.size)
        ax = axis[keep]
        P_keep = newP[keep]
        rows = np.arange(keep.size)
        P_keep[rows, ax] = np.where(U[keep, ax] > 0, L, 0.0)
        inward = np.where(U[keep, ax] > 0, -1.0, 1.0)
        # cosine-weighted direction about the inward normal
        r1 = rng.random(keep.size)
        cos_t = np.sqrt(r1)
        phi = 2 * np.pi * rng.random(keep.size)
        sin_t = np.sqrt(1 - r1)
        newdir = np.zeros((keep.size, 3))
        a1 = (ax + 1) % 3
        a2 = (ax + 2) % 3
        newdir[rows, ax] = inward * cos_t
        newdir[rows, a1] = sin_t * np.cos(phi)
        newdir[rows, a2] = sin_t * np.sin(phi)
        pos[idx[keep]] = P_keep
        u[idx[keep]] = newdir
        del W
    return path.mean(), path.std() / np.sqrt(n)


L = 1.0
for R in [0.9, 0.99]:
    mc, err = trace_cube(R, n=20000)
    # 4V/S/(1-R) counts the path correctly when photons start uniform and isotropic
    theory = (2 * L / 3) / (1 - R)
    check(f"bare cube R = {R}: mean path = (4V/S)/(1-R)", abs(mc / theory - 1) < 0.03,
          f"MC {mc:.2f} +- {err:.2f} m, formula {theory:.2f} m")
rs = 0.10
mc, err = trace_cube(0.99, sphere_radius=rs, n=20000)
A_s = 4 * np.pi * rs**2
V = L**3 - 4 / 3 * np.pi * rs**3
theory = 4 * V / (A_s + 0.01 * 6 * L**2)
check("cube + black sphere: mean path = 4V/(A_abs + (1-R)S)", abs(mc / theory - 1) < 0.05,
      f"MC {mc:.2f} +- {err:.2f} m, formula {theory:.2f} m")

section("4. Real numbers")
N_f = 1e6
P_each = 1.0  # W of light per flashlight (assumption)
P = N_f * P_each
print(f"  {N_f:.0e} flashlights at {P_each} W of light each: P = {P:.0e} W of light")

print("  (a) Ideal 1 m cube fed through a small port, no hand in it:")
print("      R            tau_light (s)   stored U (J)    mass of light (kg)")
S_cube, V_cube = 6.0, 1.0
results = {}
for R in [0.95, 0.99, 1 - 1.6e-6]:
    tau = (4 * V_cube / S_cube) / (c * (1 - R))
    U = P * tau
    results[R] = (tau, U, U / c**2)
    print(f"      {R:<12.7f} {tau:<15.3e} {U:<15.3e} {U / c**2:.3e}")
tau_best, U_best, dm_best = results[1 - 1.6e-6]

scale_best = 1e-10  # kg, best comparators, about 0.1 ug (general knowledge)
scale_kitchen = 1e-4  # kg, 0.1 g bench scale (general knowledge)
print(f"      best-case mass {dm_best:.2e} kg vs 0.1 ug comparator {scale_best:.0e} kg "
      f"(short by {scale_best / dm_best:.0f}x) and 0.1 g scale {scale_kitchen:.0e} kg")
check("best-case light mass is below the best balance", dm_best < scale_best)
U_needed = scale_best * c**2
print(f"      stored light needed to register on a 0.1 ug balance: {U_needed:.2e} J")

print("  (b) Now put a hand in (area 0.04 m^2, absorbs 60% of visible, general-knowledge figures):")
A_hand = 0.04
a_hand = 0.6
A_eff_h = a_hand * A_hand
for R in [0.99, 1 - 1.6e-6]:
    A_eff = A_eff_h + (1 - R) * S_cube
    u = 4 * P / (c * A_eff)
    U = u * V_cube
    frac = A_eff_h / A_eff
    P_hand = P * frac
    layer_C = A_hand * 1e-3 * 1000 * 3500  # 1 mm of skin, water-like, J/K
    t_burn = 12 * layer_C / P_hand  # seconds to raise that layer by 12 K (33 C to 45 C)
    press = u / 6  # pressure of isotropic light on a black surface
    F = press * A_hand
    print(f"      R = {R:.7f}: stored {U:.3g} J, mass {U / c**2:.3g} kg, hand absorbs {P_hand:.3g} W "
          f"({frac:.0%}), skin +12 K in {t_burn * 1e3:.2g} ms, light pressure {press:.3g} Pa, squeeze {F:.3g} N")
    if R > 0.999:
        check("with good mirrors the hand, not the mirrors, limits the stored light",
              frac > 0.99 and U < U_best, f"stored energy falls from {U_best:.3g} J to {U:.3g} J")
        U_hand_best = U
        F_max = P / c
        t_burn_best = t_burn
check("skin overheats in well under a second at 1 MW", t_burn_best < 0.1, f"{t_burn_best * 1e3:.2g} ms")
print(f"      maximum possible push if all {P:.0e} W came from one side: P/c = {F_max:.2e} N "
      f"= the weight of {F_max / K.g_0 * 1e3:.2f} g. Heating wins by far.")
check("radiation force is below 5 mN even in the one-sided limit", F_max < 5e-3, f"{F_max * 1e3:.2f} mN")

print("  (c) A million real flashlights sealed inside the box:")
V_each = 1e-4  # m^3
m_each = 0.15  # kg
A_each = 0.03  # m^2 of surface
a_body = 0.5
batt = 10 * 3600.0  # 10 Wh per flashlight, J
P_elec = N_f * 3.0  # 3 W electric per flashlight
V_f = N_f * V_each
Mbox = N_f * m_each
A_abs = N_f * A_each * a_body
u = 4 * P / (c * A_abs)
U = u * V_f
print(f"      flashlights occupy at least {V_f:.0f} m^3 and weigh {Mbox / 1e3:.0f} t; "
      "they themselves absorb the light")
print(f"      light in flight (with {V_f:.0f} m^3 of free space): {U:.2e} J, mass {U / c**2:.2e} kg")
E_tot = N_f * batt
dm_total = E_tot / c**2
print(f"      battery energy {E_tot:.2e} J, E/c^2 = {dm_total * 1e6:.2f} mg of {Mbox / 1e3:.0f} t "
      f"= {dm_total / Mbox:.1e} of the total")
print("      sealed: that energy only moves from chemical to light to heat, so total mass is unchanged")
print(f"      heat raises the box by about {E_tot / (Mbox * 800):.0f} K in {batt / 3.0 / 3600:.1f} h "
      f"({P_elec:.0e} W electric) if nothing leaves; flashlights would fail first")

section("5. Sealed box bookkeeping with four-vectors (units c = 1, exaggerated energies)")
# A source of rest mass 1 emits photons one at a time, recoiling each time. The walls
# (rest mass 1) then absorb those photons. Each body's rest mass comes from its own (E, p).
M0 = 1.0
src_E, src_p = M0, np.zeros(3)
n_ph = 2000
E_ph = rng.uniform(1e-5, 2e-5, n_ph)
dir_ph = random_dirs(n_ph)
ph_E, ph_p = E_ph, E_ph[:, None] * dir_ph
src_E_after = src_E - ph_E.sum()
src_p_after = src_p - ph_p.sum(axis=0)
src_rest = np.sqrt(src_E_after**2 - src_p_after @ src_p_after)
wall_E, wall_p = 1.0, np.zeros(3)
sys_E = src_E_after + ph_E.sum() + wall_E
sys_p = src_p_after + ph_p.sum(axis=0) + wall_p
sys_M_light = np.sqrt(sys_E**2 - sys_p @ sys_p)
wall_E_after = wall_E + ph_E.sum()
wall_p_after = wall_p + ph_p.sum(axis=0)
wall_rest = np.sqrt(wall_E_after**2 - wall_p_after @ wall_p_after)
sys_E2 = src_E_after + wall_E_after
sys_p2 = src_p_after + wall_p_after
sys_M_heat = np.sqrt(sys_E2**2 - sys_p2 @ sys_p2)
print(f"  energy emitted {ph_E.sum():.5f}; source rest mass {M0:.5f} -> {src_rest:.5f}; "
      f"walls 1.00000 -> {wall_rest:.5f} after absorbing")
check("source loses rest mass about equal to the energy it emits", abs((M0 - src_rest) / ph_E.sum() - 1) < 1e-3)
check("walls gain rest mass about equal to the energy they absorb", abs((wall_rest - 1) / ph_E.sum() - 1) < 1e-3)
check("system mass is 2 before, with light in flight, and after absorption",
      abs(sys_M_light - 2) < 1e-12 and abs(sys_M_heat - 2) < 1e-12,
      f"{sys_M_light:.12f} and {sys_M_heat:.12f}")
print("  The parts change; the total doesn't. With the box sealed, no scale can tell on from off.")

finish()
