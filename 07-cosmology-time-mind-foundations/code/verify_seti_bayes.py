"""C41: glass versus ocean, and a Bayesian posterior for a rate after zero detections."""
import numpy as np
from scipy.integrate import quad
from _harness import Checks

C = Checks("SETI: search fraction and zero-detection Bayes")

C.section("Volumes")
ocean_L = 1.335e21  # Eakins & Sharman 2010, as used by Wright et al. 2018
glass_L = 6 * 0.0295735  # 6 US fluid ounces in litres
f_wright = 6.0e-18  # Wright, Kanodia & Lubar 2018, Table 1 total
print(f"  ocean = {ocean_L:.3e} L, 6 oz glass = {glass_L:.3f} L, ratio = {glass_L/ocean_L:.1e}")
print(f"  Wright et al. fraction {f_wright:.1e} of ocean = {f_wright*ocean_L:.0f} L (paper says about 8000 L, hot tub to small pool)")
C.check("fraction x ocean reproduces the paper's roughly 8000 L", 7000 < f_wright * ocean_L < 9000)
print(f"  on Wright's numbers the search so far is {f_wright*ocean_L/glass_L:.0f} of your glasses, a bigger sample than your glass, still tiny")

C.section("Ocean analogy with fish actually simulated")
rng = np.random.default_rng(1)
f = 1e-3  # Use a searchable scale: each trial, each of N fish is in our sample with chance f
print("  fish N, chance that a sample of fraction f = 1e-3 catches none; analytic exp(-N f); simulated")
good = True
for Nf in [10, 100, 1000, 5000]:
    catches = rng.binomial(Nf, f, size=200000)
    sim = np.mean(catches == 0)
    ana = (1 - f)**Nf
    print(f"    N = {Nf:5d}: analytic {ana:.4f}; simulated {sim:.4f}")
    good &= abs(sim - ana) < 5e-3
C.check("simulated empty-sample rate agrees with (1-f)^N", bool(good))

C.section("Bayesian update on zero detections")
# Rate model: N transmitters in the haystack, Poisson mean N f detections if we search fraction f.
# Likelihood of zero: exp(-N f). Prior for N: log-uniform on [1, N_max] (wide ignorance).


def posterior_summary(f, Nmax=1e22, p_none=0.0):
    u = np.linspace(0, np.log(Nmax), 200001)
    like = np.exp(-np.exp(u) * f)
    post = like / np.trapezoid(like, u)
    cdf = np.cumsum(post) * (u[1] - u[0])
    return np.exp(u[np.searchsorted(cdf, 0.5)]), np.exp(u[np.searchsorted(cdf, 0.95)]), np.trapezoid(like, u) / np.log(Nmax)


print(f"  {'search fraction f':>18} {'N median':>10} {'N 95% upper':>12} {'Bayes factor (0 detected) vs prior':>36}")
for f in [1.3e-22, 6.0e-18, 1e-6, 1e-2, 0.5]:
    med, up, bf = posterior_summary(f)
    print(f"  {f:18.1e} {med:10.1e} {up:12.1e} {bf:36.3f}")
med0, up0, bf0 = posterior_summary(6e-18)
C.check("at the Wright fraction zero detections barely shift a wide prior", bf0 > 0.75)
# Conjugate test: flat prior on N >= 0 gives posterior exp(-N f) f, 95% limit -ln(0.05)/f = 3/f
limit = -np.log(0.05) / 6e-18
num = quad(lambda x: np.exp(-x), 0, -np.log(0.05))[0]
print(f"  flat prior: 95% upper limit N < 3/f = {limit:.1e} transmitters of the haystack's kind; integral check {num:.3f}")
C.check("flat-prior 95% credible region integrates to 0.95", abs(num - 0.95) < 1e-9)
print(f"  compare stars in the Milky Way (order 1e11): upper limit is {limit/1e11:.0e} per star. Not a useful bound")
C.check("limit far exceeds the number of stars, so no claim of absence follows", limit > 1e11)

C.section("Hypothesis 'none exist' versus 'exist', both fairly")
# Prior 50/50. H_some: N log-uniform in [1, Nmax]; H_none: N = 0.
for f in [6.0e-18, 1e-6, 1e-2, 0.5]:
    bf_some = posterior_summary(f)[2]  # P(zero | some)
    post_none = 0.5 * 1 / (0.5 * 1 + 0.5 * bf_some)
    print(f"  f = {f:.1e}: P(none | zero detected) = {post_none:.3f} from a 0.500 prior")
bf_stars = posterior_summary(6.0e-18, Nmax=1e11)[2]
post_stars = 0.5 / (0.5 + 0.5 * bf_stars)
print(f"  prior capped at one transmitter per star (N <= 1e11): P(none | zero detected) = {post_stars:.6f}")
print("  so the small rise above depends on allowing N far above the number of stars")
C.check("with N capped near the star count, the Wright-sized search moves the odds by under 0.01%",
        abs(post_stars - 0.5) < 1e-4)
post_small = 0.5 / (0.5 + 0.5 * posterior_summary(6e-18)[2])
C.check("posterior for 'none' rises only slightly at the real search fraction", post_small < 0.6)
C.check("absence IS evidence, just weak: posterior for 'none' rises at any f > 0", post_small > 0.5)

C.section("Where zero detections do bite: civilisations that are easy to see")
n_gal = 1e5
upper = -np.log(0.05) / n_gal
print(f"  G-hat: about {n_gal:.0e} galaxies, none reprocessing > 85% of starlight (Griffith et al. 2015)")
print(f"  if detection is reliable: 95% limit on that kind of civilisation = {upper:.1e} per galaxy")
C.check("a big fish that fills its ocean is constrained even by a modest search", upper < 1e-4)
C.finish()
