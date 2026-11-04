"""T09: colour of a heated body from Planck's law, and what a coloured lamp does to melanopsin.

Colour: Planck spectrum -> CIE 1931 XYZ using the multi-lobe Gaussian fit of Wyman, Sloan and Shirley
(JCGT 2013) -> chromaticity, dominant wavelength, linear sRGB hue. Anchor check: CIE Illuminant A is defined
as a Planckian radiator at about 2856 K with published chromaticity (0.44758, 0.40745).

Melanopsin: Govardovskii et al. (2000) A1 alpha-band template with peak 480 nm, no lens or age correction.
Used only for ratios between lamps; this is not the CIE S 026 calculation.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from common.constants import c, h, k_B

FAILS = []
OUT = Path(__file__).resolve().parents[1] / "figures"


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


WL = np.arange(360.0, 831.0, 1.0)


def lobe(w, mu, s1, s2):
    t = (w - mu) / np.where(w < mu, s1, s2)
    return np.exp(-0.5 * t**2)


def cmf(w):
    x = 1.056 * lobe(w, 599.8, 37.9, 31.0) + 0.362 * lobe(w, 442.0, 16.0, 26.7) - 0.065 * lobe(w, 501.1, 20.4, 26.2)
    y = 0.821 * lobe(w, 568.8, 46.9, 40.5) + 0.286 * lobe(w, 530.9, 16.3, 31.1)
    z = 1.217 * lobe(w, 437.0, 11.8, 36.0) + 0.681 * lobe(w, 459.0, 26.0, 13.8)
    return x, y, z


XB, YB, ZB = cmf(WL)


def planck(w_nm, T):
    w = w_nm * 1e-9
    return 2 * h * c**2 / w**5 / np.expm1(h * c / (w * k_B * T))


def rayleigh_jeans(w_nm):
    return (w_nm * 1e-9) ** -4.0


def to_xyz(spd):
    return np.trapezoid(spd * XB, WL), np.trapezoid(spd * YB, WL), np.trapezoid(spd * ZB, WL)


def xy(spd):
    X, Y, Z = to_xyz(spd)
    return X / (X + Y + Z), Y / (X + Y + Z)


M = np.array([[3.2406, -1.5372, -0.4986], [-0.9689, 1.8758, 0.0415], [0.0557, -0.2040, 1.0570]])


def linear_rgb(spd):
    X, Y, Z = to_xyz(spd)
    return M @ np.array([X, Y, Z]) / Y


def hue_deg(rgb):
    r, g, b = np.clip(rgb, 0, None)
    mx, mn = max(r, g, b), min(r, g, b)
    if mx - mn < 1e-12:
        return float("nan")
    if mx == r:
        hue = ((g - b) / (mx - mn)) % 6
    elif mx == g:
        hue = (b - r) / (mx - mn) + 2
    else:
        hue = (r - g) / (mx - mn) + 4
    return 60 * hue


LOCUS = np.array([(XB[i] / (XB[i] + YB[i] + ZB[i]), YB[i] / (XB[i] + YB[i] + ZB[i])) for i in range(len(WL))])
WHITE = np.array([0.3127, 0.3290])


def dominant(pt):
    """Dominant wavelength: the spectrum colour on the ray from white through pt."""
    d = np.asarray(pt) - WHITE
    ang = np.arctan2(d[1], d[0])
    # The Gaussian fit is unreliable in the far tails, where the true CMFs are tiny and the fitted
    # chromaticity wanders off the real locus, so the search stays within 430-680 nm.
    use = (WL >= 430) & (WL <= 680)
    v = LOCUS[use] - WHITE
    a = np.arctan2(v[:, 1], v[:, 0])
    diff = np.abs(np.angle(np.exp(1j * (a - ang))))
    i = int(np.argmin(diff))
    return WL[use][i], float(diff[i]), np.linalg.norm(d) / np.linalg.norm(v[i])


def describe(T):
    spd = planck(WL, T)
    pt = xy(spd)
    rgb = linear_rgb(spd)
    w, _, purity = dominant(pt)
    return pt, rgb, w, purity, hue_deg(rgb)


def srgb_gamma(u):
    u = np.clip(u, 0, 1)
    return np.where(u <= 0.0031308, 12.92 * u, 1.055 * u ** (1 / 2.4) - 0.055)


def melanopsin(w):
    lam = 480.0
    x = lam / w
    a = 0.8795 + 0.0459 * np.exp(-((lam - 300) ** 2) / 11940)
    return 1 / (np.exp(69.7 * (a - x)) + np.exp(28 * (0.922 - x)) + np.exp(-14.9 * (1.104 - x)) + 0.674)


def main():
    print("--- Anchor checks of the colour pipeline ---")
    e = np.ones_like(WL)
    px, py = xy(e)
    print(f"  equal-energy spectrum: xy = ({px:.4f}, {py:.4f}) (exact CIE value 1/3, 1/3)")
    check("equal-energy white within 0.005", abs(px - 1 / 3) < 0.005 and abs(py - 1 / 3) < 0.005)
    ax, ay = xy(planck(WL, 2856.0))
    print(f"  Planck at 2856 K: xy = ({ax:.4f}, {ay:.4f}); CIE Illuminant A (0.44758, 0.40745)")
    check("Illuminant A within 0.005", abs(ax - 0.44758) < 0.005 and abs(ay - 0.40745) < 0.005)
    grid = np.linspace(380, 780, 4001)
    num = np.trapezoid(planck(grid, 5778.0), grid) * 1e-9
    frac = num * np.pi / (5.670374e-8 * 5778.0**4)
    print(f"  fraction of 5778 K blackbody power in 380-780 nm: {frac:.3f}")
    check("visible fraction near the textbook 0.45", 0.38 < frac < 0.52)

    print("\n--- Colour of a heated body ---")
    print("   T/K    x      y     dominant nm  purity  hue deg   linear R:G:B      peak nm")
    table = [800, 1000, 1500, 1900, 2700, 3500, 4500, 5500, 6500, 8000, 10000, 15000, 30000, 100000]
    for T in table:
        pt, rgb, w, pur, hue = describe(T)
        norm = rgb / np.max(rgb)
        peak = 2.897771955e-3 / T * 1e9
        print(f"{T:7d} {pt[0]:.3f} {pt[1]:.3f} {w:9.0f} {pur:9.2f} {hue:7.0f}   "
              f"{norm[0]:5.2f}:{norm[1]:5.2f}:{norm[2]:5.2f} {peak:9.0f}")
    lim = xy(rayleigh_jeans(WL))
    print(f"  T -> infinity (Rayleigh-Jeans) xy = ({lim[0]:.3f}, {lim[1]:.3f}), dominant {dominant(lim)[0]:.0f} nm")

    Ts = np.geomspace(800, 1e6, 500)
    doms, hues, gfrac, chrom, purs = [], [], [], [], []
    for T in Ts:
        pt, rgb, w, pur, hue = describe(T)
        doms.append(w)
        hues.append(hue)
        purs.append(pur)
        clipped = np.clip(rgb, 0, None)
        gfrac.append(clipped[1] / clipped.max())
        chrom.append(pt)
    doms, hues, gfrac, chrom, purs = map(np.array, (doms, hues, gfrac, chrom, purs))
    # Near the white crossing the direction from white is ill-defined, so state the claim quantitatively:
    # how saturated does a blackbody ever get while its direction lies in the green or violet sector?
    green_sector = (doms >= 490) & (doms <= 565)
    violet_sector = (doms < 450) | (np.isfinite(hues) & (hues > 250) & (hues < 330))
    hue_green = np.isfinite(hues) & (hues > 75) & (hues < 170)
    print(f"  samples in the green sector: {green_sector.sum() + hue_green.sum()} of {len(Ts)}")
    print(f"  samples in the violet sector: {violet_sector.sum()} of {len(Ts)}, at T = "
          f"{Ts[violet_sector].min():.0f} to {Ts[violet_sector].max():.0f} K, highest purity there {purs[violet_sector].max():.3f}")
    tinted = purs > 0.05
    warm = hues[tinted & (hues < 180)]
    cool = hues[tinted & (hues >= 180)]
    print(f"  bodies with purity above 0.05: hue {warm.min():.0f}-{warm.max():.0f} deg (red-orange) or {cool.min():.0f}-{cool.max():.0f} deg (blue)")
    check("no temperature falls in the green sector, at any purity", green_sector.sum() + hue_green.sum() == 0)
    check("the only violet-side samples are near-white, purity below 0.05 (locus passes D65 on the magenta side)",
          purs[violet_sector].max() < 0.05)
    check("linear sRGB green never exceeds both red and blue", np.all(gfrac < 0.999999), f"max G / max(R,G,B) = {gfrac.max():.4f}")
    d_green = np.min(np.linalg.norm(chrom - LOCUS[WL.tolist().index(530.0)], axis=1))
    d_violet = np.min(np.linalg.norm(chrom - LOCUS[WL.tolist().index(400.0)], axis=1))
    print(f"  closest chromaticity distance to pure 530 nm: {d_green:.3f}; to pure 400 nm: {d_violet:.3f}")
    # Monotonic path: blue/red ratio always increases.
    br = np.array([np.trapezoid(planck(WL, T) * ZB, WL) / np.trapezoid(planck(WL, T) * XB, WL) for T in Ts])
    check("Z/X (blue against red) rises monotonically with T, so colour order cannot reverse", np.all(np.diff(br) > 0))
    # The student's "peak in green" trap: at the T whose peak is 530 nm the body looks white.
    T_green = 2.897771955e-3 / 530e-9
    pt, rgb, w, pur, hue = describe(T_green)
    print(f"  T whose Planck peak is at 530 nm (green) = {T_green:.0f} K: purity {pur:.2f}, xy = ({pt[0]:.3f}, {pt[1]:.3f})")
    check("a blackbody peaking in the green looks nearly white (purity below 0.2)", pur < 0.2)

    # White crossing: smallest distance to D65 white.
    dist = np.linalg.norm(chrom - WHITE, axis=1)
    T_white = Ts[int(np.argmin(dist))]
    print(f"  closest approach to D65 white: {dist.min():.4f} at {T_white:.0f} K")

    print("\n--- Why a thermal red lamp is useless: luminous efficacy of a blackbody ---")
    wide = np.geomspace(100.0, 2e6, 200000)
    eff = {}
    for T in (800, 1000, 1500, 2700, 3500, 6500):
        total = np.trapezoid(planck(wide, float(T)), wide)
        vis = np.trapezoid(planck(WL, float(T)) * YB, WL)
        eff[T] = 683 * vis / total
        print(f"  {T:5d} K: {eff[T]:9.3g} lumens per watt radiated")
    check("efficacy rises steeply with temperature in this range", all(a < b for a, b in zip(list(eff.values())[:-1], list(eff.values())[1:])))
    check("2700 K efficacy is near the 10-20 lm/W of a filament lamp", 8 < eff[2700] < 25, f"{eff[2700]:.1f}")

    print("\n--- Coloured lamps are not thermal: chromaticity and melanopsin per lumen ---")
    mel = melanopsin(WL)
    print(f"  template peak {WL[np.argmax(mel)]:.0f} nm, half-maximum at "
          f"{WL[np.argmax(mel > 0.5)]:.0f} and {WL[len(WL) - 1 - np.argmax(mel[::-1] > 0.5)]:.0f} nm")

    def gaussian_led(peak, fwhm):
        return np.exp(-0.5 * ((WL - peak) / (fwhm / 2.355)) ** 2)

    def filt(spd, cut):
        return spd * (WL >= cut)

    sources = {
        "candle-like 1900 K": planck(WL, 1900.0),
        "incandescent 2700 K": planck(WL, 2700.0),
        "daylight-like 6500 K": planck(WL, 6500.0),
        "2700 K with 600 nm long-pass filter": filt(planck(WL, 2700.0), 600.0),
        "amber LED 590 nm (FWHM 15)": gaussian_led(590, 15),
        "red LED 630 nm (FWHM 20)": gaussian_led(630, 20),
        "blue LED 460 nm (FWHM 20)": gaussian_led(460, 20),
    }
    ref = None
    ratios = {}
    for name, spd in sources.items():
        Y = np.trapezoid(spd * YB, WL)
        M_ = np.trapezoid(spd * mel, WL)
        ratios[name] = M_ / Y
    ref = ratios["daylight-like 6500 K"]
    print("  (dominant wavelength is poorly resolved above about 620 nm, where the spectrum locus is nearly a straight line)")
    print("  source                                xy              dominant  purity  melanopsin per lumen (6500 K = 1)")
    for name, spd in sources.items():
        pt = xy(spd)
        w, _, pur = dominant(pt)
        print(f"  {name:36s} ({pt[0]:.3f},{pt[1]:.3f}) {w:7.0f} {pur:7.2f} {ratios[name] / ref:10.4f}")
    check("red LED gives under 1 percent of the melanopsin drive per lumen of 6500 K light",
          ratios["red LED 630 nm (FWHM 20)"] / ref < 0.01)
    check("melanopsin per lumen falls with lower colour temperature",
          ratios["candle-like 1900 K"] < ratios["incandescent 2700 K"] < ratios["daylight-like 6500 K"])
    check("a real red lamp exists at purity above 0.9, outside any blackbody",
          dominant(xy(sources["red LED 630 nm (FWHM 20)"]))[2] > 0.9)
    bb_max_pur_red = max(describe(T)[3] for T in (800, 1000, 1500))
    print(f"  blackbody at 800-1500 K purity up to {bb_max_pur_red:.2f}: deep but dim; and the 2700 K filament has "
          f"{np.trapezoid(sources['2700 K with 600 nm long-pass filter'] * YB, WL) / np.trapezoid(planck(WL, 2700.0) * YB, WL):.2f} "
          "of its luminous output left after the filter")

    print("\n--- Does the bulb's colour change room temperature? ---")
    V = 40.0
    ach = 0.5
    ua_walls = 100.0  # W/K, assumed: about 100 m2 of ordinary room surface at 1 W/m2/K
    vent = 1.2 * 1005 * V * ach / 3600
    G_loss = vent + ua_walls
    for P, label in ((5.0, "5 W LED"), (60.0, "60 W filament lamp"), (80.0, "sleeping person")):
        print(f"  {label:20s}: room warms by {P / G_loss:.2f} K at steady state ({vent:.1f} + {ua_walls:.0f} W/K loss)")
    check("5 W lamp changes room temperature by under 0.1 K", 5.0 / G_loss < 0.1)

    # Figure: colour strip with spectra.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 1, figsize=(8, 4.3), gridspec_kw={"height_ratios": [1, 2]})
    Tstrip = np.geomspace(800, 40000, 300)
    colours = []
    for T in Tstrip:
        rgb = np.clip(linear_rgb(planck(WL, T)), 0, None)
        colours.append(srgb_gamma(rgb / rgb.max()))
    ax[0].imshow(np.array(colours)[None, :, :], aspect="auto", extent=[np.log10(800), np.log10(40000), 0, 1])
    ticks = [800, 1500, 3000, 6000, 10000, 20000, 40000]
    ax[0].set_xticks(np.log10(ticks))
    ax[0].set_xticklabels([str(t) for t in ticks])
    ax[0].set_yticks([])
    ax[0].set_xlabel("temperature, K (brightness normalised at each T)")
    for T, col in ((1500, "#b5442f"), (3000, "#c98a2b"), (6000, "#777777"), (15000, "#2f5da8")):
        spd = planck(WL, T)
        ax[1].plot(WL, spd / spd.max(), color=col, label=f"{T} K")
    ax[1].plot(WL, mel, "k--", lw=1, label="melanopsin template")
    ax[1].set_xlabel("wavelength, nm")
    ax[1].set_ylabel("relative power per nm")
    ax[1].legend(fontsize=8, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.18), frameon=False)
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    path = OUT / "t09_blackbody_colour.png"
    fig.savefig(path, dpi=110)
    size = path.stat().st_size
    print(f"\nfigure {path.name}: {size / 1024:.0f} KB")
    check("figure under 300 KB", size < 300 * 1024)

    print(f"\n{len(FAILS)} checks failed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
