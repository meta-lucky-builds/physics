"""T10-T17: does radiation come back? Tests, all simulated, none assuming the answer.

  A. Self-absorption of a shape (Monte Carlo rays): convex sphere, hemispherical cup, grey-cup emissivity.
  B. Kirchhoff: photons exchanged between two grey plates. Equal T stays balanced only if absorptivity = emissivity.
  C. Einstein detailed balance: stochastic atoms and one light mode. With stimulated emission the mode reaches
     the Bose-Einstein occupation; without it, the Wien form.
  D. Direction: Fresnel emissivity of water at 10 um against angle, polarisation.
  E. Photon diffusion in absorbing, re-emitting air: chance of returning to a plane or a small body.
  F. Kac ring: deterministic, reversible, yet relaxes; exact recurrence at 2N steps.
  G. Ehrenfest urn: mean time to return to all-left against 2^N. Then the scale for real photons.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.special import gammaln

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import c, h, k_B, sigma_SB

FAILS = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def lambert_dirs(rng, normals):
    """Cosine-weighted directions about the given unit normals."""
    n = len(normals)
    u1, u2 = rng.random(n), rng.random(n)
    r = np.sqrt(u1)
    phi = 2 * np.pi * u2
    local = np.stack([r * np.cos(phi), r * np.sin(phi), np.sqrt(1 - u1)], axis=1)
    a = np.where(np.abs(normals[:, 0:1]) < 0.9, np.array([[1.0, 0, 0]]), np.array([[0, 1.0, 0]]))
    t1 = np.cross(normals, a)
    t1 /= np.linalg.norm(t1, axis=1, keepdims=True)
    t2 = np.cross(normals, t1)
    return local[:, 0:1] * t1 + local[:, 1:2] * t2 + local[:, 2:3] * normals


def part_a(rng):
    print("--- A. Does an emitting body catch its own radiation? ---")
    N = 2_000_000
    # Convex sphere radius 1: launch outward from surface, test for any later hit on the same sphere.
    p = rng.normal(size=(N, 3))
    p /= np.linalg.norm(p, axis=1, keepdims=True)
    d = lambert_dirs(rng, p)
    b = np.sum(p * d, axis=1)
    # second root of |p + t d| = 1 is t = -2 b; b >= 0 outward means no positive hit
    returned = np.mean(-2 * b > 1e-9)
    print(f"  convex sphere: fraction of rays that return to the body = {returned:.6f}")
    check("convex body in empty space never catches its own rays", returned == 0)

    # Hemispherical cup radius 1, z < 0 part of a sphere, opening at z = 0.
    p = rng.normal(size=(N, 3))
    p[:, 2] = -np.abs(p[:, 2])
    p /= np.linalg.norm(p, axis=1, keepdims=True)
    d = lambert_dirs(rng, -p)  # inward normal
    t = -2 * np.sum(p * d, axis=1)
    hit_z = p[:, 2] + t * d[:, 2]
    f_self = np.mean(hit_z < 0)
    print(f"  hemispherical cup, inner surface: self view factor by rays = {f_self:.4f}; reciprocity gives 1 - A_open/A_cup = 0.5")
    check("cup view factor agrees with area reciprocity", abs(f_self - 0.5) < 0.002)

    # Grey cup: follow rays with absorption probability eps at each wall hit, launched from the opening.
    for eps in (0.05, 0.3, 0.9):
        M = 400_000
        r = np.sqrt(rng.random(M))
        ph = 2 * np.pi * rng.random(M)
        o = np.stack([r * np.cos(ph), r * np.sin(ph), np.zeros(M)], axis=1)
        down = np.tile([0.0, 0.0, -1.0], (M, 1))
        d = lambert_dirs(rng, down)
        alive = np.ones(M, bool)
        absorbed = np.zeros(M, bool)
        pos = o.copy()
        for bounce in range(400):
            if not alive.any():
                break
            idx = np.nonzero(alive)[0]
            bq = np.sum(pos[idx] * d[idx], axis=1)
            cq = np.sum(pos[idx] ** 2, axis=1) - 1
            t = -bq + np.sqrt(np.clip(bq**2 - cq, 0, None))
            hit = pos[idx] + t[:, None] * d[idx]
            esc = hit[:, 2] >= 0
            alive[idx[esc]] = False
            keep = idx[~esc]
            hit = hit[~esc]
            hit /= np.linalg.norm(hit, axis=1, keepdims=True)
            ab = rng.random(len(keep)) < eps
            absorbed[keep[ab]] = True
            alive[keep[ab]] = False
            surv = keep[~ab]
            pos[surv] = hit[~ab]
            d[surv] = lambert_dirs(rng, -hit[~ab])
        sim = absorbed.mean()
        F = 0.5
        analytic = eps / (1 - (1 - eps) * F)
        print(f"  grey cup surface emissivity {eps:.2f}: opening looks like emissivity {sim:.4f}; "
              f"multiple-reflection formula eps/(1-(1-eps)F) = {analytic:.4f}")
        check(f"cavity effect for eps = {eps}", abs(sim - analytic) < 0.004)


def part_b(rng):
    print("\n--- B. Kirchhoff by photon bookkeeping: two large parallel plates, same T ---")
    print("  Plate 1 emits per unit area eps1*E and absorbs with probability a1. Plate 2 is grey, eps2 = a2 = 0.6.")

    def net_power(eps1, a1, eps2=0.6, a2=0.6, n=600_000):
        # photons launched in proportion to emissivity; each bounces until absorbed
        rates = np.array([eps1, eps2])
        start = rng.random(n) < rates[1] / rates.sum()  # True = starts on plate 2
        target = np.where(start, 0, 1)   # 0 = plate 1, 1 = plate 2
        alive = np.ones(n, bool)
        absorbed = np.array([0, 0])
        while alive.any():
            idx = np.nonzero(alive)[0]
            p_abs = np.where(target[idx] == 0, a1, a2)
            ab = rng.random(len(idx)) < p_abs
            absorbed[0] += np.sum(target[idx][ab] == 0)
            absorbed[1] += np.sum(target[idx][ab] == 1)
            alive[idx[ab]] = False
            target[idx[~ab]] = 1 - target[idx[~ab]]
        emitted1 = np.sum(~start)
        return (emitted1 - absorbed[0]) / n * rates.sum()  # net loss of plate 1, per unit E

    for eps1, a1 in ((0.2, 0.2), (0.9, 0.9), (0.5, 0.5), (0.8, 0.4), (0.3, 0.7)):
        val = net_power(eps1, a1)
        print(f"  eps1 = {eps1}, a1 = {a1}: net power out of plate 1 = {val:+.4f} E")
        if eps1 == a1:
            check(f"eps = a = {eps1}: no net flow at equal temperature", abs(val) < 0.004)
        else:
            sign_ok = (val > 0.02) if eps1 > a1 else (val < -0.02)
            check(f"eps {eps1} != a {a1}: heat would flow between equal temperatures, a second-law violation", sign_ok)


def part_c(rng):
    print("\n--- C. Einstein detailed balance, stochastic simulation of atoms and one mode ---")

    def simulate(x, stimulated, n_events=300_000):
        """Gillespie: N1, N2 atoms held by a thermal bath at ratio exp(-x); photon number n."""
        ratio = np.exp(-x)
        N_total = 1000.0
        N2 = N_total * ratio / (1 + ratio)
        N1 = N_total - N2
        n = 0
        t_tot = 0.0
        acc = 0.0
        waits = rng.exponential(1.0, n_events).tolist()
        picks = rng.random(n_events).tolist()
        burn = n_events // 10
        for i in range(n_events):
            up = N2 * (n + 1) if stimulated else N2
            down = N1 * n
            total = up + down
            dt = waits[i] / total
            if i > burn:
                acc += n * dt
                t_tot += dt
            n += 1 if picks[i] < up / total else -1
        return acc / t_tot

    print("   h nu / kT   sim (with stimulated)  Bose 1/(e^x-1)   sim (spontaneous only)  Wien e^-x")
    for x in (0.5, 1.0, 2.0, 4.0):
        with_st = simulate(x, True)
        without = simulate(x, False)
        bose = 1 / np.expm1(x)
        wien = np.exp(-x)
        print(f"   {x:6.2f} {with_st:17.4f} {bose:16.4f} {without:20.4f} {wien:14.4f}")
        check(f"x = {x}: with stimulated emission gives Planck/Bose occupation", abs(with_st / bose - 1) < 0.08)
        check(f"x = {x}: without it, Wien form", abs(without / wien - 1) < 0.08)
    # Rate asymmetry: emission rate ~ n+1, absorption ~ n. At n = 0 only emission survives.
    print("   per-atom rates into a mode holding n photons: absorb = A n; emit = A (n + 1).")
    print("   Spontaneous emission is the '+1'. It is the only term with no absorption counterpart.")
    for T, lam in ((307.0, 10e-6), (5778.0, 500e-9), (3000.0, 589e-9)):
        x = h * c / (lam * k_B * T)
        print(f"   T = {T:.0f} K, {lam * 1e9:.0f} nm: x = {x:.2f}, mean photons per mode {1 / np.expm1(x):.3g}")


def part_d():
    print("\n--- D. Direction: water at 10 um (n = 1.218, k = 0.0508) ---")

    def fresnel(theta, n_complex):
        ci = np.cos(theta)
        s2 = np.sin(theta) ** 2
        ct = np.sqrt(n_complex**2 - s2)
        rs = (ci - ct) / (ci + ct)
        rp = (n_complex**2 * ci - ct) / (n_complex**2 * ci + ct)
        return np.abs(rs) ** 2, np.abs(rp) ** 2

    nc = 1.218 + 0.0508j
    print("   angle   emissivity s   emissivity p   mean")
    for ang in (0, 30, 60, 75, 85, 89):
        rs, rp = fresnel(np.deg2rad(ang), nc)
        print(f"   {ang:5d} {1 - rs:12.4f} {1 - rp:14.4f} {1 - 0.5 * (rs + rp):8.4f}")
    th = np.linspace(0, np.pi / 2, 20001)
    rs, rp = fresnel(th, nc)
    hemi = 2 * np.trapezoid((1 - 0.5 * (rs + rp)) * np.cos(th) * np.sin(th), th)
    rs0, rp0 = fresnel(0.0, nc)
    print(f"   hemispherical emissivity {hemi:.4f}, normal {1 - rs0:.4f}")
    check("emissivity depends strongly on direction near grazing", (1 - rs[-2]) < 0.1 and (1 - rs[0]) > 0.97)
    check("s and p agree at normal incidence, split at oblique incidence", abs(rs0 - rp0) < 1e-12 and
          abs(fresnel(np.deg2rad(75), nc)[0] - fresnel(np.deg2rad(75), nc)[1]) > 0.02)


def part_e(rng):
    print("\n--- E. Radiative random walk: a photon absorbed and re-emitted isotropically in gas ---")

    def slab(tau, n=200_000):
        """Start at z=0 going into the medium, cosine-distributed. Absorbing surface at z=0, free escape at z=tau."""
        cos0 = np.sqrt(rng.random(n))
        z = -np.log(1 - rng.random(n)) * cos0
        alive = (z < tau) & (z > 0)
        returned = 0
        transmitted = np.sum(z >= tau)
        z = z[alive]
        while len(z):
            mu = 2 * rng.random(len(z)) - 1
            z = z - np.log(1 - rng.random(len(z))) * mu
            returned += np.sum(z <= 0)
            transmitted += np.sum(z >= tau)
            z = z[(z > 0) & (z < tau)]
        return returned / n, transmitted / n

    print("   optical depth   back to source   through    two-stream 1/(1 + 3 tau/4) for transmission")
    results = []
    for tau in (0.1, 1, 5, 20, 100):
        back, thru = slab(tau, n=120_000 if tau < 50 else 8_000)
        results.append((tau, back, thru))
        est = 1 / (1 + 0.75 * tau)
        print(f"   {tau:9.1f} {back:15.4f} {thru:10.4f} {est:12.4f}")
        check(f"tau = {tau}: every photon ends somewhere", abs(back + thru - 1) < 1e-9)
    check("plane source: return probability tends to 1 as the gas gets thicker", results[-1][1] > 0.95 and results[0][1] < 0.15)

    def sphere(a, n=20_000, r_max=40.0):
        """Small absorbing sphere radius a (mean free paths) in a uniform re-emitting gas."""
        p = rng.normal(size=(n, 3))
        p /= np.linalg.norm(p, axis=1, keepdims=True)
        d = lambert_dirs(rng, p)
        pos = p * a + d * (-np.log(1 - rng.random(n)))[:, None]
        back = 0
        outward = 0
        alive = np.ones(n, bool)
        for _ in range(6000):
            if not alive.any():
                break
            idx = np.nonzero(alive)[0]
            step = rng.normal(size=(len(idx), 3))
            step /= np.linalg.norm(step, axis=1, keepdims=True)
            L = -np.log(1 - rng.random(len(idx)))
            o = pos[idx]
            b = np.sum(o * step, axis=1)
            disc = b**2 - (np.sum(o**2, axis=1) - a**2)
            t = -b - np.sqrt(np.clip(disc, 0, None))
            hit = (disc > 0) & (t > 0) & (t < L)
            back += hit.sum()
            alive[idx[hit]] = False
            rem = idx[~hit]
            pos[rem] = o[~hit] + (L[~hit])[:, None] * step[~hit]
            far = np.linalg.norm(pos[rem], axis=1) > r_max
            outward += far.sum()
            alive[rem[far]] = False
        return back / n, outward / n, alive.sum() / n, r_max

    print("   body radius (mean free paths)   returned by 40 mfp   reached 40 mfp   still walking   corrected to infinity")
    vals = []
    for a in (0.1, 1.0, 10.0):
        ret, out, un, rm = sphere(a)
        # A walker at radius rm later returns with probability about a / rm (harmonic function in 3D).
        total = ret + (out + un) * a / rm
        vals.append((a, total))
        print(f"   {a:10.1f} {ret:26.3f} {out:16.3f} {un:15.4f} {total:18.3f}")
    est = 0.75 * 10 / (1 + 0.75 * 10)
    print(f"   diffusion limit for a >> 1, by analogy with the slab: 0.75 a / (1 + 0.75 a) = {est:.3f} at a = 10")
    check("small body: most photons never come back", vals[0][1] < 0.3)
    check("large body relative to mean free path: return becomes likely", vals[2][1] > 0.7)
    check("a = 10 agrees with the diffusion limit within 0.03", abs(vals[2][1] - est) < 0.03)
    print("   Earth's plane-like surface: back radiation about 333 W/m2 (Trenberth et al. 2009), i.e. this effect at planetary scale.")


def part_f():
    print("\n--- F. Kac ring: deterministic reversible dynamics ---")

    def Hbin(m):
        p = np.clip((1 + m) / 2, 1e-12, 1 - 1e-12)
        return -(p * np.log(p) + (1 - p) * np.log(1 - p))

    for N, seed, flip in ((5000, 101, False), (5000, 101, True)):
        rng = np.random.default_rng(seed)
        marks = rng.random(N) < 0.1
        if flip:  # make the number of marks have the opposite parity
            marks[np.argmin(marks)] = True
        s0 = np.ones(N, dtype=int)
        sign = np.where(marks, -1, 1)

        def step(s):
            return np.roll(s * sign, 1)

        def back(s):
            return np.roll(s, -1) * sign

        s = s0.copy()
        mags = []
        first_return = None
        for t in range(1, 2 * N + 1):
            s = step(s)
            mags.append(s.mean())
            if first_return is None and np.array_equal(s, s0):
                first_return = t
            if t == 150:
                s150 = s.copy()
        mags = np.array(mags)
        mu = marks.mean()
        n_marks = int(marks.sum())
        print(f"   N = {N}, {n_marks} marks ({'even' if n_marks % 2 == 0 else 'odd'}): magnetisation at t = 10, 50, 150, N, 2N: "
              f"{mags[9]:.3f}, {mags[49]:.3f}, {mags[149]:.3f}, {mags[N - 1]:+.3f}, {mags[2 * N - 1]:+.3f}; first exact return t = {first_return}")
        print(f"     molecular-chaos prediction (1 - 2 mu)^t at t = 10, 50: {(1 - 2 * mu) ** 10:.3f}, {(1 - 2 * mu) ** 50:.3f}")
        print(f"     coarse-grained entropy per site at t = 1, 10, 150, first return: {Hbin(mags[0]):.3f}, {Hbin(mags[9]):.3f}, "
              f"{Hbin(mags[149]):.3f}, {Hbin(mags[first_return - 1]):.3f} (max ln 2 = {np.log(2):.3f})")
        expected = N if n_marks % 2 == 0 else 2 * N
        check("exact return at N for an even number of marks, 2N for odd (each site meets every mark once per lap)",
              first_return == expected, f"t = {first_return}")
        check("early relaxation follows the molecular-chaos prediction", abs(mags[9] - (1 - 2 * mu) ** 10) < 0.05)
        check("magnetisation has decayed by t = 150", abs(mags[149]) < 0.1)
        x = s150.copy()
        for _ in range(150):
            x = back(x)
        check("reversing the map from t = 150 restores the start exactly", np.array_equal(x, s0))


def part_g(rng):
    print("\n--- G. Ehrenfest urn: mean time to return to all particles in one box ---")
    print("   N    measured mean return   2^N")
    good = True
    for N in (2, 4, 6, 8, 10, 12):
        n_returns = 2000 if N <= 8 else 400
        lefts = N
        steps = 0
        returns = 0
        # vectorised over batches of random draws
        while returns < n_returns:
            draws = rng.integers(0, N, 100000)
            for d in draws.tolist():
                steps += 1
                if d < lefts:
                    lefts -= 1
                else:
                    lefts += 1
                if lefts == N:
                    returns += 1
                    if returns >= n_returns:
                        break
        mean = steps / returns
        print(f"   {N:2d} {mean:15.1f} {2**N:12d}")
        good &= abs(mean / 2**N - 1) < 0.12
    check("mean return time matches Kac's lemma, 1/P(state) = 2^N", good)

    print("\n--- Scale for real radiation ---")
    V = 30.0  # m^3 room
    lam1, lam2 = 8e-6, 12e-6
    nu = np.linspace(c / lam2, c / lam1, 2000)
    M = np.trapezoid(8 * np.pi * nu**2 / c**3, nu) * V
    print(f"   electromagnetic modes in a {V:.0f} m3 room between 8 and 12 um: {M:.2e}")
    T = 307.0
    A = 1.35
    x = np.linspace(0.01, 60, 200000)
    photon_flux = A * 2 * np.pi * (k_B * T) ** 3 / (h**3 * c**2) * np.trapezoid(x**2 / np.expm1(x), x)
    print(f"   photons emitted per second by {A} m2 of 307 K skin (black): {photon_flux:.2e}")
    log10_one = -np.log10(M)
    log10_all = log10_one * photon_flux
    print(f"   chance a given photon sits in its original mode after mixing: 10^{log10_one:.1f}")
    print(f"   chance all of one second's photons do so at once: 10^({log10_all:.2e})")
    # Entropy route: S = (4/3) U / T, recurrence time ~ tau exp(S/k)
    U = sigma_SB * T**4 * A
    S_over_k = 4 / 3 * U / T / k_B
    print(f"   entropy of one second's emission: S/k = {S_over_k:.2e}; recurrence time ~ exp(S/k): "
          f"10^({S_over_k / np.log(10):.2e}) of any time unit")
    log_age = np.log10(4.35e17 / 1e-14)  # age of universe in optical periods
    print(f"   for comparison, age of the universe is only 10^{log_age:.0f} optical periods")
    check("the two routes agree that the exponent itself is above 1e22", abs(log10_all) > 1e22 and S_over_k / np.log(10) > 1e22)
    print("   Single-photon case: a lossless mirror box returns a photon pulse every round trip, about")
    print(f"   {2 * 3 / c * 1e9:.0f} ns for 3 m. Real wall reflectance near 0.1 at 10 um: survival after n bounces 0.1^n.")
    for n in (1, 3, 10):
        print(f"     n = {n:2d}: {0.1 ** n:.0e}")
    # Ratio of emitted energy per photon
    E_ph = sigma_SB * T**4 * A / photon_flux
    print(f"   mean photon energy {E_ph:.2e} J = {E_ph / (k_B * T):.2f} kT")
    check("mean thermal photon energy equals the known 2.70 kT", abs(E_ph / (k_B * T) - 2.70) < 0.02)


def main():
    rng = np.random.default_rng(2024)
    part_a(rng)
    part_b(rng)
    part_c(rng)
    part_d()
    part_e(rng)
    part_f()
    part_g(rng)
    print(f"\n{len(FAILS)} checks failed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
