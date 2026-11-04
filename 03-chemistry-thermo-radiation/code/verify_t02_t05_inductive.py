"""T02-T05: how far a substituent's effect on acidity reaches, tested against aqueous pKa data.

Two rival descriptions are fitted to the same data:
  (a) empirical through-bond attenuation: dpKa = A * f^(n-1), n = carbons from carboxyl;
  (b) through-space electrostatics: a C-Cl dipole acting on the carboxylate charge,
      averaged over chain conformations, with one fitted effective dielectric constant.

Model (b) has one overall scale, so its decay shape is a prediction, though the conformer
treatment was selected after trying three (all three are printed). Data: Engineering ToolBox pKa table (25 C, water), OpenStax Organic Chemistry 20.4.
"""
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import N_A, R, e as e_charge
from scipy.constants import epsilon_0

FAILS = []


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


PKA = {
    "formic": 3.74, "acetic": 4.76, "propanoic": 4.87, "butanoic": 4.82, "pentanoic": 4.86, "hexanoic": 4.87,
    "2-methylpropanoic": 4.84, "pivalic": 4.78,
    "chloroacetic": 2.87, "dichloroacetic": 1.35, "trichloroacetic": 0.66,
    "2-chloropropanoic": 2.83, "3-chloropropanoic": 3.98,
    "2-chlorobutanoic": 2.86, "3-chlorobutanoic": 4.05, "4-chlorobutanoic": 4.52,
    "acrylic": 4.25, "propiolic": 1.89,
}
# (position n of Cl counted from the carboxyl carbon's neighbour, dpKa relative to its own parent)
SERIES_BUT = [(1, PKA["butanoic"] - PKA["2-chlorobutanoic"]),
              (2, PKA["butanoic"] - PKA["3-chlorobutanoic"]),
              (3, PKA["butanoic"] - PKA["4-chlorobutanoic"])]
ALL_POINTS = SERIES_BUT + [(1, PKA["acetic"] - PKA["chloroacetic"]),
                           (1, PKA["propanoic"] - PKA["2-chloropropanoic"]),
                           (2, PKA["propanoic"] - PKA["3-chloropropanoic"])]


def fit_attenuation(points):
    n = np.array([p[0] for p in points], float)
    y = np.log([p[1] for p in points])
    slope, intercept = np.polyfit(n, y, 1)
    return np.exp(slope), np.exp(intercept + slope)  # f, A (value at n=1)


def rot(axis, angle):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * K @ K


def next_atom(a, b, c, length, bond_angle, dihedral):
    """Place d so that |cd| = length, angle bcd = bond_angle, dihedral abcd = dihedral."""
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    if np.linalg.norm(n) < 1e-9:
        n = np.cross(bc, [0.3, 0.5, 0.8])
    n = n / np.linalg.norm(n)
    # direction making bond_angle with (b - c), then rotated by dihedral about bc
    direction = rot(n, np.pi - bond_angle) @ bc
    direction = rot(bc, dihedral) @ direction
    return c + direction * length


TET = np.deg2rad(111.0)
TRIG = np.deg2rad(122.0)
R_CC = 1.53
R_CCL = 1.79
R_CO2 = 1.53        # C(alpha)-C(carboxyl)
ION_OFFSET = 0.70   # carboxylate charge centre beyond carboxyl C (midpoint of O atoms, about 1.25 cos 56 deg)
MU_CCL = 1.87 * 3.33564e-30  # C-Cl group dipole taken from CH3Cl, C m (see README)


def energy_kJ(rCl, rC, ion, eps):
    """Interaction of unit negative charge at ion with the C(+q)-Cl(-q) dipole, kJ/mol. Negative stabilises."""
    q = MU_CCL / (R_CCL * 1e-10)
    pref = -e_charge * q / (4 * np.pi * epsilon_0 * eps * 1e-10) * N_A / 1000
    r_C = np.linalg.norm(rC - ion)
    r_Cl = np.linalg.norm(rCl - ion)
    # charge -e with +q on C and -q on Cl
    return pref * (1 / r_C - 1 / r_Cl)


def conformers(n_cl, double_bond=False, r23=None):
    """Yield (Cl position, the carbon bearing Cl, ion position) for all staggered conformers.

    Chain placed so the carboxyl lies on the x axis: ion, C1 (carboxyl), C2 ... Cl on carbon C(1+n_cl).
    The ion sits on the C2->C1 axis, so rotation about C1-C2 changes nothing.
    double_bond=True makes C2=C3 trans and planar, the (E)-but-2-enoic case with Cl on C4.
    """
    ion = np.array([R_CO2 + ION_OFFSET, 0.0, 0.0])
    c1 = np.array([R_CO2, 0.0, 0.0])
    c2 = np.array([0.0, 0.0, 0.0])
    ang2 = TRIG if double_bond else TET
    len23 = r23 if r23 is not None else (1.33 if double_bond else R_CC)
    c3 = c2 + len23 * np.array([np.cos(ang2), np.sin(ang2), 0.0])
    dih = [np.pi, np.pi / 3, -np.pi / 3]
    if n_cl == 1:
        # Cl directly on C2; geometry fixed by the bond angle
        rcl = c2 + R_CCL * np.array([np.cos(TET), np.sin(TET), 0.0])
        yield rcl, c2, ion
        return
    if n_cl == 2:
        for d in dih:
            yield next_atom(c1, c2, c3, R_CCL, TET, d), c3, ion
        return
    if n_cl == 3:
        if double_bond:
            c4 = next_atom(c1, c2, c3, 1.50, TRIG, np.pi)  # trans across the double bond
            for d in (np.pi, np.pi / 3, -np.pi / 3):
                yield next_atom(c2, c3, c4, R_CCL, TET, d), c4, ion
            return
        for d1, d2 in itertools.product(dih, dih):
            c4 = next_atom(c1, c2, c3, R_CC, TET, d1)
            yield next_atom(c2, c3, c4, R_CCL, TET, d2), c4, ion
        return
    if n_cl >= 4:
        for ds in itertools.product(dih, repeat=n_cl - 1):
            atoms = [c1, c2, c3]
            for d in ds[:-1]:
                atoms.append(next_atom(atoms[-3], atoms[-2], atoms[-1], R_CC, TET, d))
            yield next_atom(atoms[-3], atoms[-2], atoms[-1], R_CCL, TET, ds[-1]), atoms[-1], ion


K_LN10 = np.log(10) * R * 298.15 / 1000  # kJ/mol per pKa unit


def shifts_unit_eps(n_cl, double_bond=False, r23=None):
    """Interaction energies at eps = 1 for each staggered conformer, kJ/mol. First is all-anti."""
    return np.array([energy_kJ(rcl, rc, ion, 1.0) for rcl, rc, ion in conformers(n_cl, double_bond, r23)])


def dpka_field(n_cl, eps, mode, double_bond=False, r23=None):
    """pKa shift for conformer treatment `mode`:
    equal: every staggered conformer counted equally, populations not altered by the charge;
    anti: extended chain only;
    boltzmann: conformers re-weighted by the same screened energy in the anion.
    """
    U = shifts_unit_eps(n_cl, double_bond, r23) / eps
    if mode == "equal":
        dG = U.mean()
    elif mode == "anti":
        dG = U[0]
    else:
        dG = -(R * 298.15 / 1000) * np.log(np.mean(np.exp(-U / (R * 298.15 / 1000))))
    return -dG / K_LN10


def fit_eps(points, mode):
    grid = np.linspace(1.5, 40, 771)
    errs = [sum((dpka_field(n, e, mode) - y) ** 2 for n, y in points) for e in grid]
    i = int(np.argmin(errs))
    return grid[i], errs[i]


def main():
    print("--- T02: replacing H by alkyl (data, pKa in water at 25 C) ---")
    for name in ("formic", "acetic", "propanoic", "butanoic", "pentanoic", "hexanoic"):
        print(f"  {name:12s} {PKA[name]:.2f}")
    d1 = PKA["acetic"] - PKA["formic"]
    d2 = abs(PKA["propanoic"] - PKA["acetic"])
    d3 = max(PKA[k] for k in ("propanoic", "butanoic", "pentanoic", "hexanoic")) - \
        min(PKA[k] for k in ("propanoic", "butanoic", "pentanoic", "hexanoic"))
    print(f"  H->CH3 on the carboxyl: dpKa = {d1:+.2f}")
    print(f"  CH3->C2H5:              dpKa = {d2:+.2f}")
    print(f"  spread of C3-C6 acids:  {d3:.2f}")
    check("first alkyl matters, later ones do not (factor over 5)", d1 > 5 * d2 and d3 <= 0.06)
    print(f"  electronegativity (Pauling): H 2.20, C 2.55 -- not equal")
    print("  hybridisation of the attached carbon (propanoic sp3, acrylic sp2, propiolic sp):")
    print(f"    {PKA['propanoic']:.2f} -> {PKA['acrylic']:.2f} -> {PKA['propiolic']:.2f}")
    check("acidity rises with s character", PKA["propanoic"] > PKA["acrylic"] > PKA["propiolic"])

    print("\n--- T04/T05: attenuation fit, Cl position along the chain ---")
    for n, y in SERIES_BUT:
        print(f"  Cl {n} carbon(s) from COOH, butanoic series: dpKa = {y:.2f}")
    f_but, A_but = fit_attenuation(SERIES_BUT)
    f_all, A_all = fit_attenuation(ALL_POINTS)
    print(f"  butanoic only:  dpKa(n) = {A_but:.2f} x {f_but:.3f}^(n-1)")
    print(f"  all six points: dpKa(n) = {A_all:.2f} x {f_all:.3f}^(n-1)")
    r12 = SERIES_BUT[1][1] / SERIES_BUT[0][1]
    r23 = SERIES_BUT[2][1] / SERIES_BUT[1][1]
    print(f"  successive ratios: {r12:.3f}, {r23:.3f}")
    check("decay is close to geometric (ratios within 0.03)", abs(r12 - r23) < 0.03)

    rng = np.random.default_rng(1)
    boots = []
    for _ in range(4000):
        noisy = [(n, y + rng.normal(0, 0.02 * np.sqrt(2))) for n, y in SERIES_BUT]
        boots.append(fit_attenuation(noisy)[0])
    print(f"  f with pKa errors of +/-0.02: {np.mean(boots):.3f} +/- {np.std(boots):.3f}")
    check("attenuation factor differs from 1/2 given +/-0.02 pKa errors",
          abs(np.mean(boots) - 0.5) > 3 * np.std(boots))

    print("  extrapolated shifts (no cut-off is built in):")
    for n in range(1, 8):
        v = A_but * f_but ** (n - 1)
        print(f"    n = {n}: dpKa = {v:.3f}  Ka ratio = {10**v:.3g}")
    v4 = A_but * f_but**3
    v5 = A_but * f_but**4
    check("n=4 is about the size of table scatter, n=5 is below +/-0.05",
          0.05 < v4 < 0.2 and v5 < 0.05, f"n=4 {v4:.3f}, n=5 {v5:.3f}")

    print("\n--- T03: same distance, different skeleton ---")
    print(f"  Cl on C2: acetic {ALL_POINTS[3][1]:.2f}, propanoic {ALL_POINTS[4][1]:.2f}, butanoic {SERIES_BUT[0][1]:.2f}")
    print(f"  Cl on C3: propanoic {ALL_POINTS[5][1]:.2f}, butanoic (Cl-bearing carbon also carries CH3) {SERIES_BUT[1][1]:.2f}")
    alpha = [ALL_POINTS[3][1], ALL_POINTS[4][1], SERIES_BUT[0][1]]
    check("Cl shift at fixed distance depends weakly on the rest of the skeleton (spread < 0.2)",
          max(alpha) - min(alpha) < 0.2 and abs(ALL_POINTS[5][1] - SERIES_BUT[1][1]) < 0.2)

    print("  additivity of branches, chloroacetic acids:")
    base = PKA["acetic"]
    s = [base - PKA["chloroacetic"], base - PKA["dichloroacetic"], base - PKA["trichloroacetic"]]
    for k, v in enumerate(s, 1):
        print(f"    {k} Cl: dpKa = {v:.2f}, additive prediction {k * s[0]:.2f}, ratio {v / (k * s[0]):.2f}")
    check("second Cl adds about 80 percent as much as the first; third about 35 percent",
          0.7 < (s[1] - s[0]) / s[0] < 0.9 and 0.3 < (s[2] - s[1]) / s[0] < 0.5)

    print("  through-bond path counting for the bicyclo[2.2.2]octane cage (three 4-bond paths, X on C4):")
    per = A_all * f_all**3
    print(f"    single chain at that distance dpKa = {per:.2f}; three paths in parallel = {3 * per:.2f} (Cl-equivalent)")

    print("\n--- Through-space test: C-Cl dipole acting on carboxylate charge ---")
    print("  One fitted constant (effective dielectric constant) for each of three conformer treatments.")
    results = {}
    for mode in ("boltzmann", "anti", "equal"):
        eps, err = fit_eps(SERIES_BUT, mode)
        pred = [dpka_field(n, eps, mode) for n in range(1, 6)]
        results[mode] = (eps, err, pred)
        print(f"  {mode:9s}: eps = {eps:5.2f}, shifts n=1..5 = " + ", ".join(f"{x:.2f}" for x in pred)
              + f", worst error at n=1..3 = {max(abs(pred[i] - SERIES_BUT[i][1]) for i in range(3)):.2f}")
    print("  data     :               shifts n=1..3 = " + ", ".join(f"{y:.2f}" for _, y in SERIES_BUT))
    eq = results["equal"]
    check("equal-weight electrostatics matches all three positions within 0.1 with one constant",
          all(abs(eq[2][i] - SERIES_BUT[i][1]) < 0.1 for i in range(3)))
    check("result is sensitive to conformer treatment: Boltzmann and anti-only both miss by more than 0.25",
          max(abs(results["boltzmann"][2][i] - SERIES_BUT[i][1]) for i in range(3)) > 0.25
          and max(abs(results["anti"][2][i] - SERIES_BUT[i][1]) for i in range(3)) > 0.25)
    check("every fitted dielectric constant lies far below bulk water (78)", all(r[0] < 15 for r in results.values()))
    rms_a = np.sqrt(np.mean([(A_but * f_but ** (n - 1) - y) ** 2 for n, y in SERIES_BUT]))
    rms_b = np.sqrt(eq[1] / 3)
    print(f"  rms misfit: 2-parameter exponential {rms_a:.3f}, 1-parameter electrostatic {rms_b:.3f}.")
    print("  The conformer treatment was chosen after seeing the misfits, so this shows electrostatics is")
    print("  sufficient in size and decay; it does not prove the through-bond picture wrong.")

    print("\n--- T05: does a C=C between Cl and COOH extend the reach? (equal-weight model, same eps) ---")
    eps = eq[0]
    base = dpka_field(3, eps, "equal")
    short = dpka_field(3, eps, "equal", False, 1.33)
    uns = dpka_field(3, eps, "equal", True)
    print(f"  Cl-CH2-CH2-CH2-COOH, all staggered conformers           : dpKa = {base:.3f}")
    print(f"  same, C2-C3 shortened from 1.53 to 1.33 A, nothing else  : dpKa = {short:.3f}")
    print(f"  Cl-CH2-CH=CH-COOH (E): shorter, planar, no folding       : dpKa = {uns:.3f}")
    for label, flag, r in (("saturated", False, None), ("shortened", False, 1.33), ("E-alkene", True, None)):
        rs = [np.linalg.norm(rcl - ion) for rcl, _, ion in conformers(3, flag, r)]
        print(f"  {label:9s}: Cl to charge centre min/mean/max = {min(rs):.2f} / {np.mean(rs):.2f} / {max(rs):.2f} A")
    print(f"  gain from bond length alone {short - base:+.3f}; total gain with the rigid alkene geometry {uns - base:+.3f}")
    check("bond shortening alone cannot supply a whole extra carbon of reach (needs factor 1/f = 2.6)",
          short / base < 1 / f_but, f"factor {short / base:.2f} vs {1 / f_but:.2f}")
    check("alkene geometry transmits at least as well as the saturated chain in this model", uns >= base)

    print("\n--- Resonance: why a distance rule alone fails (para-benzoic acids, OpenStax Table 20.4) ---")
    benz = {"H": 4.19, "NO2": 3.41, "CN": 3.55, "CHO": 3.75, "Br": 3.96, "Cl": 4.00, "CH3": 4.34, "OCH3": 4.46, "OH": 4.48}
    for k, v in benz.items():
        print(f"    p-{k:5s} {v:.2f}  dpKa = {benz['H'] - v:+.2f}")
    check("p-OCH3 and p-OH weaken the acid despite electronegative O (resonance donation beats induction)",
          benz["OCH3"] > benz["H"] and benz["OH"] > benz["H"])
    meta_cl, para_cl = 3.84, 4.00
    print(f"    Cl: meta {4.20 - meta_cl:+.2f}, para {4.20 - para_cl:+.2f} (Engineering ToolBox, benzoic 4.20)")
    check("Cl stronger at meta than para", meta_cl < para_cl)

    print(f"\n{len(FAILS)} checks failed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
