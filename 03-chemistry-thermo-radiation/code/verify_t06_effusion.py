"""T06: effusion simulated with particles, not with the formula being tested.

Two simulators:

1. Exact ballistic tracing with no time step, in a unit box. Each flight runs from wall to wall. Walls scatter
   diffusely (cosine law, speed kept), standing in for the mixing that collisions and rough walls give a real gas.
   A flight ending in the hole is an escape; the particle is put back at a random place with its own speed and a
   new direction, so each species keeps its number and speed distribution. This is the limit of a hole far
   smaller than the mean free path. A mirror-wall version is also run, because it fails instructively.

2. Time-stepped hard spheres with real collisions, for a hole far wider than the mean free path.

Tests:
  A. fitted exponent of rate against mass. Rival claims: -1/2 and -1.
  B. absolute rate per area and mean escaping energy, against the kinetic-theory derivation.
  C. equal masses, different diameters, sphere must clear the rim.
  D. does the mass separation survive when collisions control the flow?
"""
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

FAILS = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def random_dirs(rng, n):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def fold(u):
    """Position in [0, 1] of a coordinate unfolded across mirror walls."""
    u = np.mod(u, 2.0)
    return np.where(u > 1.0, 2.0 - u, u)


def next_wall(pos, vel):
    """Time and axis of the next wall hit for flights inside the unit box."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t_axis = np.where(vel > 0, (1.0 - pos) / vel, np.where(vel < 0, -pos / vel, np.inf))
    axis = np.argmin(t_axis, axis=1)
    return t_axis[np.arange(len(pos)), axis], axis


def ballistic(masses, counts, hole_radius, t_end, seed, diameters=None, walls="diffuse"):
    """Exact collisionless effusion. Returns escapes and mean kinetic energy of escapers per species."""
    rng = np.random.default_rng(seed)
    species = np.concatenate([np.full(c, k) for k, c in enumerate(counts)])
    N = len(species)
    mass = np.asarray(masses, float)[species]
    limit = np.full(N, float(hole_radius))
    if diameters is not None:
        limit = hole_radius - 0.5 * np.asarray(diameters, float)[species]
    vel = rng.normal(size=(N, 3)) * np.sqrt(1.0 / mass)[:, None]  # kT = 1
    speed = np.linalg.norm(vel, axis=1)
    pos = rng.random((N, 3))
    clock = np.zeros(N)
    escapes = np.zeros(len(masses), int)
    ke_sum = np.zeros(len(masses))
    active = np.arange(N)
    while len(active):
        dt, axis = next_wall(pos[active], vel[active])
        keep = clock[active] + dt < t_end
        active, dt, axis = active[keep], dt[keep], axis[keep]
        if len(active) == 0:
            break
        clock[active] += dt
        pos[active] += vel[active] * dt[:, None]
        at_x0 = (axis == 0) & (vel[active, 0] < 0)
        out = np.zeros(len(active), bool)
        out[at_x0] = np.hypot(pos[active[at_x0], 1] - 0.5, pos[active[at_x0], 2] - 0.5) < limit[active[at_x0]]
        ids = active[out]
        if len(ids):
            np.add.at(escapes, species[ids], 1)
            np.add.at(ke_sum, species[ids], 0.5 * mass[ids] * speed[ids] ** 2)
            pos[ids] = rng.random((len(ids), 3))
            vel[ids] = random_dirs(rng, len(ids)) * speed[ids][:, None]
        stay = active[~out]
        ax = axis[~out]
        sides = pos[stay, ax] > 0.5
        pos[stay, ax] = np.where(sides, 1.0, 0.0)
        if walls == "mirror":
            vel[stay, ax] = -vel[stay, ax]
        else:
            # cosine law about the inward wall normal, same speed
            u1, u2 = rng.random(len(stay)), rng.random(len(stay))
            r = np.sqrt(u1)
            phi = 2 * np.pi * u2
            d = np.stack([r * np.cos(phi), r * np.sin(phi), np.sqrt(1 - u1)], axis=1)  # third entry along the normal
            v_new = np.empty((len(stay), 3))
            for k in range(3):
                mask = ax == k
                others = [m for m in range(3) if m != k]
                v_new[mask, k] = d[mask, 2] * np.where(sides[mask], -1.0, 1.0)
                v_new[mask, others[0]] = d[mask, 0]
                v_new[mask, others[1]] = d[mask, 1]
            vel[stay] = v_new * speed[stay][:, None]
    return escapes, ke_sum / np.maximum(escapes, 1)


def collide(pos, vel, mass, d, pairs):
    """Elastic hard-sphere impulse for approaching overlapping pairs; each particle at most once per step."""
    if len(pairs) == 0:
        return 0
    i, j = pairs[:, 0], pairs[:, 1]
    dr = pos[j] - pos[i]
    dist = np.linalg.norm(dr, axis=1)
    n = dr / dist[:, None]
    vn = np.sum((vel[i] - vel[j]) * n, axis=1)
    keep = vn > 0
    i, j, n, vn = i[keep], j[keep], n[keep], vn[keep]
    # a pair is used only if both its particles have not yet collided this step
    seen = np.zeros(len(pos), dtype=bool)
    use = []
    for k, (a, b) in enumerate(zip(i.tolist(), j.tolist())):
        if not (seen[a] or seen[b]):
            seen[a] = seen[b] = True
            use.append(k)
    if not use:
        return 0
    i, j, n, vn = i[use], j[use], n[use], vn[use]
    imp = (2 * vn * mass[i] * mass[j] / (mass[i] + mass[j]))[:, None] * n
    vel[i] -= imp / mass[i][:, None]
    vel[j] += imp / mass[j][:, None]
    return len(i)


def hard_spheres(masses, counts, d, hole_radius, t_end, dt, seed):
    rng = np.random.default_rng(seed)
    species = np.concatenate([np.full(c, k) for k, c in enumerate(counts)])
    N = len(species)
    mass = np.asarray(masses, float)[species]
    pos = rng.random((N, 3))
    vel = rng.normal(size=(N, 3)) * np.sqrt(1.0 / mass)[:, None]
    E0 = 0.5 * np.sum(mass[:, None] * vel**2)
    escapes = np.zeros(len(masses), int)
    n_coll = 0
    for _ in range(int(round(t_end / dt))):
        pos += vel * dt
        n_coll += collide(pos, vel, mass, d, cKDTree(pos).query_pairs(d, output_type="ndarray"))
        low = pos[:, 0] < 0
        if low.any():
            idx = np.nonzero(low)[0]
            out = np.hypot(pos[idx, 1] - 0.5, pos[idx, 2] - 0.5) < hole_radius
            esc = idx[out]
            np.add.at(escapes, species[esc], 1)
            if len(esc):
                speed = np.linalg.norm(vel[esc], axis=1)
                vel[esc] = random_dirs(rng, len(esc)) * speed[:, None]
                pos[esc] = rng.random((len(esc), 3))
            bounce = idx[~out]
            pos[bounce, 0] = -pos[bounce, 0]
            vel[bounce, 0] = -vel[bounce, 0]
        for ax in range(3):
            hi = pos[:, ax] > 1
            pos[hi, ax] = 2 - pos[hi, ax]
            vel[hi, ax] = -vel[hi, ax]
            if ax > 0:
                lo = pos[:, ax] < 0
                pos[lo, ax] = -pos[lo, ax]
                vel[lo, ax] = -vel[lo, ax]
    E1 = 0.5 * np.sum(mass[:, None] * vel**2)
    return escapes, 2 * n_coll / N, E1 / E0


def main():
    t0 = time.time()

    print("--- A. Collisionless gas, kT = 1, hole radius 0.05 in a unit box, exact tracing, diffuse walls ---")
    masses = [1.0, 2.0, 4.0, 8.0, 16.0, 32.0]
    counts = [1500] * 6
    a = 0.05
    t_end = 800.0
    esc, ke = ballistic(masses, counts, a, t_end, 1)
    print("  mass   escapes   rate      kinetic theory n vbar A/4   ratio")
    area = np.pi * a**2
    ratios = []
    for m, n_esc, c in zip(masses, esc, counts):
        rate = n_esc / t_end
        theory = c * np.sqrt(8 / (np.pi * m)) * area / 4
        ratios.append(rate / theory)
        print(f"  {m:5.0f} {n_esc:8d} {rate:8.3f} {theory:12.3f} {rate / theory:17.3f}")
    x = np.log(masses)
    slope, _ = np.polyfit(x, np.log(esc), 1)
    err = np.sqrt(np.sum((x - x.mean()) ** 2 / esc)) / np.sum((x - x.mean()) ** 2)
    print(f"  fitted exponent of rate against mass = {slope:.4f} +/- {err:.4f} (Poisson)")
    check("exponent consistent with -1/2", abs(slope + 0.5) < 4 * err)
    check("exponent excludes the -1 of 'rate inversely proportional to mass'", abs(slope + 1.0) > 20 * err)
    check("absolute rate agrees with n vbar A / 4 within 4 percent for every mass", all(abs(r - 1) < 0.04 for r in ratios),
          f"worst {max(abs(r - 1) for r in ratios):.3f}")
    print(f"  ratio light/heavy for masses 4 and 32: measured {esc[2] / esc[5]:.2f}, sqrt law {np.sqrt(8):.2f}, 1/M law {8:.2f}")

    print("\n  Same run with perfect mirror walls and no collisions (nothing mixes directions):")
    esc_m, _ = ballistic(masses, counts, a, 500.0, 1, walls="mirror")
    ratios_m = [e / 500.0 / (c * np.sqrt(8 / (np.pi * m)) * area / 4) for e, c, m in zip(esc_m, counts, masses)]
    print("  rate / (n vbar A / 4) by mass: " + ", ".join(f"{r:.2f}" for r in ratios_m))
    slope_m, _ = np.polyfit(x, np.log(esc_m), 1)
    print(f"  exponent {slope_m:.3f}: orbits that can reach the hole are used up, so the formula's assumption of an")
    print("  isotropic gas at the hole has failed. The law needs something to keep the gas mixed.")
    check("mirror box falls below the kinetic-theory rate, showing the isotropy assumption does work", ratios_m[0] < 0.9)

    print("\n--- B. Energy carried by escaping particles ---")
    for m, k in zip(masses, ke):
        print(f"  mass {m:4.0f}: mean KE of escapers = {k:.3f} kT (bulk gas 1.5 kT, flux theory 2 kT)")
    mean_ke = float(np.average(ke, weights=esc))
    check("escapers carry 2 kT, not 1.5 kT", abs(mean_ke - 2.0) < 0.05, f"{mean_ke:.3f} kT")
    check("escaper energy is independent of mass", np.ptp(ke) < 0.25, f"range {np.ptp(ke):.3f}")

    print("\n--- C. Equal mass, different diameter, sphere must clear the rim by d/2 (dilute limit) ---")
    print("  The rim rule is a geometric input here; the run shows what follows from it and how small it is.")
    for d1, d2 in ((0.001, 0.004), (0.001, 0.05)):
        a2 = 0.05
        esc2, _ = ballistic([1.0, 1.0], [1500, 1500], a2, 400.0, 7, diameters=[d1, d2])
        expected = ((a2 - d1 / 2) / (a2 - d2 / 2)) ** 2
        ratio = esc2[0] / esc2[1]
        unc = ratio * np.sqrt(1 / esc2[0] + 1 / esc2[1])
        print(f"  d = {d1} and {d2}, hole radius {a2}: escapes {esc2[0]} and {esc2[1]}, ratio {ratio:.3f} +/- {unc:.3f}, "
              f"geometry alone {expected:.3f}")
        check(f"ratio agrees with aperture geometry for d2 = {d2}", abs(ratio - expected) < 3.5 * unc)
    print("  Real gases: a 0.3 nm molecule at a 1 um hole gives d/(2a) = 1.5e-4, a correction of 3e-4 in the rate.")

    print("\n--- D. Mixture of masses 1 and 16: separation against Knudsen number Kn = mean free path / hole diameter ---")
    results = []
    esc3, _ = ballistic([1.0, 16.0], [1500, 1500], 0.05, 400.0, 11)
    sep = esc3[0] / esc3[1]
    unc = sep * np.sqrt(1 / esc3[0] + 1 / esc3[1])
    results.append((sep, unc))
    print(f"  no collisions (Kn infinite): escapes {esc3[0]} / {esc3[1]}, separation {sep:.2f} +/- {unc:.2f} (sqrt law 4.00, 1/M law 16.00)")
    N, d, a3, tend, dt = 10000, 0.033, 0.20, 2.0, 0.003
    lam = 1 / (np.sqrt(2) * N * np.pi * d**2)
    esc4, coll, er = hard_spheres([1.0, 16.0], [N // 2, N // 2], d, a3, tend, dt, 5)
    sep2 = esc4[0] / esc4[1]
    unc2 = sep2 * np.sqrt(1 / esc4[0] + 1 / esc4[1])
    phi = N * np.pi * d**3 / 6
    print(f"  hard spheres: N = {N}, d = {d}, volume fraction {phi:.2f}, mean free path {lam:.3f}, hole diameter {2 * a3:.2f}, "
          f"Kn = {lam / (2 * a3):.3f}")
    print(f"    {coll:.1f} collisions per particle, energy ratio final/initial {er:.4f}")
    print(f"    escapes {esc4[0]} / {esc4[1]}, separation {sep2:.2f} +/- {unc2:.2f}")
    check("collisionless limit gives the sqrt(M) separation", abs(sep - 4) < 3.5 * unc)
    check("collisions conserve energy within 2 percent (re-insertion also keeps speeds)", abs(er - 1) < 0.02)
    check("a hole wide against the mean free path separates far less", sep2 + 3 * unc2 < sep - 3 * unc, f"{sep2:.2f} vs {sep:.2f}")

    print(f"\nrun time {time.time() - t0:.0f} s, {len(FAILS)} checks failed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
