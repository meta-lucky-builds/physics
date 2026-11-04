# 7. Cosmology, time, foundations and mind

This section covers the questions furthest from the lab bench: why the universe expands, what a
time machine would need, whether physics is "real", and what a thought is. The answers vary a lot
in how firm they are. Some are derivations with code behind them. Some are conditional
thought experiments. Some cannot be settled by any experiment we know how to do, and the status
line says so. Read the status before you trust the answer.

**Scripts** (run from the repo root, each exits non-zero on failure):

| Script | Questions |
|---|---|
| [code/verify_cosmology.py](code/verify_cosmology.py) | C01, C02, C03, C27 |
| [code/verify_geodesic_no_edge.py](code/verify_geodesic_no_edge.py) | C04 |
| [code/verify_time_dilation.py](code/verify_time_dilation.py) | C06, C07, C08 |
| [code/verify_ctc_solutions.py](code/verify_ctc_solutions.py) | C09 to C17 |
| [code/verify_block_time.py](code/verify_block_time.py) | C18 to C25 |
| [code/verify_boltzmann_brains.py](code/verify_boltzmann_brains.py) | C31 to C34 |
| [code/verify_effective_theories.py](code/verify_effective_theories.py) | C35, C36, C37 |
| [code/verify_pythagoras.py](code/verify_pythagoras.py) | C39 |
| [code/verify_brain_numbers.py](code/verify_brain_numbers.py) | C40 |
| [code/verify_seti_bayes.py](code/verify_seti_bayes.py) | C41 |
| [code/verify_dna_bayes.py](code/verify_dna_bayes.py) | C42 to C46 |

Sources are listed in [References](#references), which also says which ones could not be confirmed.

---

### C01: Why the universe expands, into what, and how we know

**Status:** Derived + Code + Established + Open

**Short answer:** No. The universe does not expand in order to raise entropy. Expansion follows
from initial conditions plus the field equations, and nothing in those equations aims at anything.
It doesn't need to expand "into" anything, because the equations describe distances inside space
and need no outside. The Big Bang model has no particle at a location: it says all of space was
hotter and denser, the part we can now see was once very small, and the infinite-density moment is
an extrapolation at which the theory stops working.

**Starting from:** a uniform distribution of matter on large scales (the cosmological principle);
Newton's gravity with the shell theorem for the dust derivation; general relativity for the
pressure term and for the meaning of curvature.

**Derivation:**

*Friedmann equations.* Take a uniform ball of density $\rho$ and a shell at radius $R(t)=a(t)\chi$,
with $\chi$ fixed. The shell theorem says only the mass inside acts, $M=\tfrac{4}{3}\pi R^3\rho$, and
$M$ stays constant as dust moves with the shell. So

$$\ddot R=-\frac{GM}{R^2}\quad\Rightarrow\quad\frac{\ddot a}{a}=-\frac{4\pi G}{3}\rho.$$

Multiply by $\dot R$ and integrate. The integration constant is the shell's energy per unit mass,
which we write $-\tfrac12 kc^2\chi^2$:

$$\left(\frac{\dot a}{a}\right)^2=\frac{8\pi G}{3}\rho-\frac{kc^2}{a^2}.$$

General relativity gives the same equations, adds $\Lambda$, replaces $\rho$ with
$\rho+3p/c^2$ in the acceleration equation, and makes $k$ the spatial curvature:

$$H^2=\frac{8\pi G}{3}\rho-\frac{kc^2}{a^2}+\frac{\Lambda c^2}{3},\qquad
\frac{\ddot a}{a}=-\frac{4\pi G}{3}\left(\rho+\frac{3p}{c^2}\right)+\frac{\Lambda c^2}{3},\qquad
\dot\rho=-3H\left(\rho+\frac{p}{c^2}\right).$$

*Entropy.* No entropy term appears in these equations. Today's $\dot a>0$ is carried from an earlier
$\dot a>0$, which is an initial condition. Here is a sharper test. Take a comoving volume $V=a^3$ of a
perfect fluid, so its energy is $E=\rho c^2V$. Then

$$T\,dS=dE+p\,dV=c^2V\,d\rho+(\rho c^2+p)\,dV=0,$$

where the last step uses the continuity equation above with $dV/V=3\,da/a$. Smooth expansion of a
perfect fluid leaves its entropy unchanged. So expansion isn't in itself an entropy-making process.
Entropy grows through stars shining, gas clumping under gravity and black holes forming, and why the
early universe had so little entropy available is a separate [open question](#c27-what-exists-before-the-origin-of-time).

You're onto a real connection, though. Expansion cools the universe faster than some reactions can
keep up (that's why nucleosynthesis stopped with helium mostly, not iron), and this keeps matter out
of equilibrium. Expansion helps entropy production happen. It's not there for that purpose.

*Into what?* Everything measurable is in the metric $ds^2=-c^2dt^2+a(t)^2d\Sigma^2$. It says the
distance between two galaxies at fixed coordinates grows. That's a statement about lengths inside
space, checkable from inside, and no surrounding space is used. A balloon's surface is a useful
picture only if you ignore the air around it. Nothing forbids a larger space, but the equations don't
need one, and no measurement has required one.

*The "particle".* This is the correction that matters most. The model says density was higher
everywhere: $\rho_m\propto a^{-3}$, $\rho_r\propto a^{-4}$. If space is infinite now, it was
infinite at every earlier time with $a>0$. What was small was the patch that grew into our observable
universe. Scaling by $a\approx T_0/T$, today's 46 Gly radius was roughly 100 light-years when nuclei
formed at $T\sim10^9$ K (Estimate). As $a\to0$ the equations give $\rho\to\infty$ at a finite
past time. Most physicists read this as a sign that the theory has left its domain, not as a physical
infinite-density object: general relativity ignores quantum effects, and these should matter long
before that.

*What is the universe?* All of spacetime and everything in it. The observable universe is the ball
around us from which light has had time to arrive. We don't know how far the whole extends.

*How we know it expands.* Three independent lines:

1. **Redshift against distance.** Hubble (1929) reported a linear relation using 24 galaxies, with
   velocities largely from Slipher, and Lemaître (1927) had already derived expansion and suggested a
   similar relation. Davis & Lineweaver (2004) analysed supernova data and report treating redshift as
   a plain special-relativistic Doppler shift ruled out at 23 sigma.
2. **Cosmic microwave background.** Found by Penzias & Wilson (1965). Its measured temperature is
   $2.72548\pm0.00057$ K (Fixsen 2009). A hot, dense past predicts such a relic.
3. **Light-element abundances.** Fusion in the first minutes predicts helium and deuterium abundances.
   Fields et al. (2020) compare these with Planck data. The abstract reports good precision and also a
   live issue: one theoretical rate change it examined moves D/H out of agreement.

*Hubble tension.* Planck gives $H_0=67.4\pm0.5$ km/s/Mpc, assuming ΛCDM. SH0ES gives
$73.04\pm1.04$ km/s/Mpc from Cepheids and supernovae (Riess et al. 2022), which those authors call a
5-sigma difference with unknown cause. This is about the rate, not whether expansion happens.

**Check:** [code/verify_cosmology.py](code/verify_cosmology.py) integrates the shell's acceleration
with no use of the first Friedmann equation and finds that equation holds as a constant of motion to
$10^{-9}$. It says nothing about the pressure term, which comes from general relativity. The output is
under [C02](#c02-93-billion-light-years-in-138-billion-years).

**Limits:** the Newtonian route works for dust, with $k$ only an integration constant. The entropy
result assumes a perfect fluid with no dissipation. Times before roughly nucleosynthesis are inferred
with less direct evidence.

### C02: 93 billion light-years in 13.8 billion years

**Status:** Derived + Code + Established

**Short answer:** Space kept stretching while the light was on its way, so the place that emitted it
is now much farther than light-speed times travel time. The age is about 13.8 billion years (13.6 is
an older or rounded figure). The computed radius today is about 46 Gly, so the diameter is about
92 Gly. Nothing moved faster than light through the space around it. And that's the observable part,
not necessarily the whole universe.

**Starting from:** the Friedmann equation above with flat ΛCDM; Planck 2018 values $H_0=67.4$ km/s/Mpc
and $\Omega_m=0.315$ (Planck base-ΛCDM age: $13.797\pm0.023$ Gyr); radiation from $T_0=2.72548$ K and
$N_\text{eff}=3.046$.

**Derivation:** with $E(a)=\sqrt{\Omega_r a^{-4}+\Omega_m a^{-3}+\Omega_\Lambda}$,

$$t_0=\int_0^1\frac{da}{aH_0E(a)}.$$

Light obeys $ds=0$, so $c\,dt=a\,d\chi$ and the comoving distance covered is

$$\chi=\int\frac{c\,dt}{a}=\int_0^1\frac{c\,da}{a^2H_0E(a)}.$$

Because $a<1$ in the past, $1/a>1$, so every early stretch of the path counts for more today than
$c\,dt$. That's why $\chi>ct_0$.

Now the speed limit. A photon heading to us has proper distance $d=a\chi$, and differentiating gives

$$\dot d=Hd-c.$$

The $-c$ is the photon moving at $c$ through local space, and the $Hd$ is the stretching of everything
between. The script confirms this by finite differences along a real ray. When $Hd>c$, the
photon's distance from us grows, even as it heads our way. The light of the CMB was emitted just 41 Mly
from where we would be, and it spent ages being carried away before it could make headway. Relativity's
limit applies between objects at the same place. Recession rate isn't a local speed, and for distant
objects it can exceed $c$ without any conflict (Davis & Lineweaver 2004).

**Check:** [code/verify_cosmology.py](code/verify_cosmology.py) gets age 13.79 Gyr, within 0.05% of
Planck's value, and matches the analytic no-radiation formula. The horizon comes out at 46.13 Gly by two
integration routes, so the diameter is 92.3 Gly. The often-quoted 93 Gly differs from this only by
choice of parameters and rounding. Actual output:

<details><summary>verify_cosmology.py output</summary>

```text
Cosmology: Friedmann, age, horizons, reachability
=================================================

-- Newtonian shell versus Friedmann equation (matter + Lambda)
  a after 40 Gyr = 10.464; max |k c^2 / H0^2| over the run = 1.47e-09
  [PASS] a'^2 - (8 pi G rho/3) a^2 - Lambda a^2/3 stays constant (0 for flat)

-- Age and horizons (Planck 2018 parameters)
  Omega_r = 9.209e-05, Omega_Lambda = 0.6849, 1/H0 = 14.51 Gyr
  age t0 = 13.791 Gyr (Planck base-LCDM table: 13.797 +/- 0.023)
  [PASS] age agrees with Planck to 0.3%
  [PASS] numeric age without radiation matches the analytic asinh formula (13.7962 Gyr)
  particle horizon = 46.13 Gly; diameter = 92.3 Gly; c*t0 = 13.79 Gly
  [PASS] particle horizon: adaptive quadrature and fixed-grid march agree to 0.01%
  [PASS] diameter lies in 90-95 Gly and exceeds 2 c t0
  CMB light: source now at 45.2 Gly; proper distance at emission 41.4 Mly
  its recession rate: 66.4 c at emission, 3.12 c now (coordinate rates, not local speeds)
  [PASS] that source recedes faster than c, yet its light reached us
  at a = 0.5: d(proper distance)/dt = -0.3154 c; H d - c = -0.3154 c
  [PASS] proper-distance rate = H d - c: the photon moves at c through local space

-- Hubble tension sensitivity
  c/H0 = 14.51 Gly (Planck), 13.39 Gly (SH0ES)
  age with H0 = 73.04, Om = 0.315 held fixed = 12.73 Gyr (illustrative, not a fit)
  tension = 5.64 km/s/Mpc = 4.9 sigma from the quoted errors
  [PASS] the two H0 values differ by close to the 5 sigma Riess et al. quote

-- Lookback, Hubble radius, event horizon, the 94% figure
  z =     1: lookback   7.95 Gyr, comoving distance  11.1 Gly
  z =     3: lookback  11.65 Gyr, comoving distance  21.2 Gly
  z =    10: lookback  13.32 Gyr, comoving distance  31.4 Gly
  z = 14.44: lookback  13.51 Gyr, comoving distance  33.8 Gly
  Hubble radius c/H0 = 14.51 Gly; comoving event horizon = 16.68 Gly
  unreachable fraction of the observable VOLUME = 1-(chi_e/chi_p)^3 = 95.3%
  counting only out to the most distant confirmed galaxy, z = 14.44: 88.0%
  [PASS] student's 94% lies between the two counting assumptions
  a photon launched now and followed to a = 1e8 covers 16.68 Gly of comoving distance
  [PASS] marched photon saturates at the computed event horizon (agree to 0.1%)

-- Future reachability of today's observable volume
   cosmic age (Gyr)      a  event horizon (Gly)  reachable % of volume
               13.8      1                16.68                 4.7261
               24.9      2                 8.70                 0.6715
               37.0      4                 4.38                 0.0855
               49.1      8                 2.19                 0.0107
               61.2     16                 1.10                 0.0013
               85.5     64                 0.27                 0.0000
  99% unreachable at cosmic age 23 Gyr; 99.99% unreachable at age 49 Gyr
  late-time e-folding of the scale factor = 17.5 Gyr
  [PASS] 99.99% unreachable arrives in tens of Gyr, not trillions of years

-- Andromeda: light delay, expansion and binding
  light delay = 2.5e+06 yr; Hubble flow at 2.5 Mly would be 52 km/s recession
  measured radial velocity = -109.3 km/s (van der Marel et al. 2012), the opposite sign
  mass that makes -109 km/s bound at this separation = 1.06e+12 Msun; measured M_LG = 3.17e12 Msun
  [PASS] Local Group mass exceeds the binding mass, so the Hubble flow does not apply here
  fractional stretch of M31's light if expansion did apply, H0 D/c = 1.7e-04
  fractional change of the CMB temperature over 2000 yr = 1.4e-07

-- Backward extrapolation: scale factor reaches 0 at finite time
  age at a = 1e-12 is 7.56e-13 yr
  [PASS] the integral of da/(aH) converges as a -> 0: finite age within this model

13/13 checks passed
```

</details>

**Limits:** assumes exact flat ΛCDM through all epochs. The numbers shift with $H_0$: the Hubble tension
alone changes $c/H_0$ by about 8%.

### C03: Are we seeing the past, does expansion count, and is 94% out of reach?

**Status:** Derived + Code + Established

**Short answer:** Yes, everything we see is the past, and a star could have exploded with the news
still on its way. For Andromeda it's a 2.5-million-year delay, and expansion doesn't count, because
Andromeda is approaching and gravitationally bound to us. Your 94% is close: I get 95.3% by volume of
the observable universe, or 88% if you count only out to the farthest confirmed galaxy. 99.99%
unreachable comes after about 35 more billion years, so between your "few billion" and "trillion".

**Starting from:** the same model, plus $c$ as the limit on signals. Andromeda: 2.5 Mly (as you say),
radial velocity $-109.3\pm4.4$ km/s and Local Group mass $(3.17\pm0.57)\times10^{12}\,M_\odot$ (van der
Marel et al. 2012).

**Derivation:**

*The past.* Light from distance $D$ left $D/c$ ago. Proxima Centauri is about 4 years ago, Andromeda
2.5 million. An observer sees events on their past light cone only, so a star could have died
without us knowing, and until the news arrives "has it exploded?" has no observational answer (see
[C25](#c25-why-the-block-universe-was-proposed) for why "now" over there is frame-dependent anyway).

*What did people in the past see?* The same sky, from an earlier light cone. Someone 2000 years ago
saw Andromeda as it was 2000 years before the image we get. Their ancestors saw earlier still. The
sky barely changes over such times: Andromeda moved under one light-year in that span at 109 km/s, and
the CMB was warmer by one part in $10^7$. The chain of "the past of the past" ends at the CMB,
because the universe before it was opaque to light.

*Does expansion count for Andromeda?* Hubble's law would give 52 km/s of recession at this distance.
The measured value is $-109$ km/s, the opposite sign. It takes only about $10^{12}\,M_\odot$ to bind
at this speed and separation, and the measured mass is three times that. Expansion applies to
scales where matter is not bound. Bound systems don't take part.

*Hubble radius and event horizon.* These differ, and they're often confused.

- **Hubble radius** $c/H_0=14.5$ Gly: where the recession rate equals $c$.
- **Event horizon** $\chi_e=\int_1^\infty c\,da/(a^2H)=16.7$ Gly today: any signal we send now can reach
  no farther. The script confirms it by following a photon launched today, which saturates at this value.
- **Particle horizon** 46.1 Gly: the farthest light that has reached us by now.

*The 94% check.* If galaxies have uniform comoving density, the reachable share is the volume ratio:

$$1-\left(\frac{16.68}{46.13}\right)^3=95.3\%.$$

This needs care, because galaxies didn't exist all the way out to the particle horizon. If you stop at
the most distant confirmed galaxy, $z=14.44$ (Naidu et al. 2025), then 33.8 Gly gives 88%. So the
number depends on what you count, and your 94% lies inside that range. The assumption is that ΛCDM
continues forever with constant $\Lambda$. It's also for light; ships are worse.

"Out of reach" doesn't mean out of sight. We keep seeing those galaxies, with their light redder and
older, but what we send now will never arrive.

*The future.* As $a$ grows, the comoving event horizon shrinks like $1/a$, so reachable volume falls
quickly:

| Cosmic age (Gyr) | Scale factor | Event horizon (Gly) | Reachable % of today's observable volume |
|---|---|---|---|
| 13.8 | 1 | 16.68 | 4.73 |
| 24.9 | 2 | 8.70 | 0.67 |
| 37.0 | 4 | 4.38 | 0.086 |
| 49.1 | 8 | 2.19 | 0.011 |
| 61.2 | 16 | 1.10 | 0.0013 |

99.99% unreachable arrives at age 49 Gyr, about 35 Gyr from now. Only bound systems such as the Local
Group stay with us. Krauss & Scherrer (2007) argue that far-future observers in such an island may lose
the evidence for expansion altogether, which is a sobering thought given how sure we feel now.

**Check:** [code/verify_cosmology.py](code/verify_cosmology.py), with output under
[C02](#c02-93-billion-light-years-in-138-billion-years). The table above is copied from it.

**Limits:** assumes $\Lambda$ never changes, which isn't known. If dark energy evolves, the far-future
numbers change. The binding estimate treats the pair as two point masses on a radial orbit.

### C04: What happens to a ball thrown at the edge of the universe?

**Status:** Derived + Code + Established + Open + Conceptual

**Short answer:** Your dilemma is a good and very old argument, and it's right against an edge with
a wall. But it has a hidden third case: the ball can come back with nothing pushing it, from the
opposite direction, because space closes on itself. That's what "no edge" means. It isn't known
whether our universe is shaped that way. The measured curvature is consistent with flat, and what lies
beyond the horizon cannot be observed.

**Starting from:** Riemannian geometry. A force-free object follows a geodesic, the straightest
possible path. Measurements: $\Omega_K=0.001\pm0.002$ from Planck plus BAO (Planck 2018).

**Derivation:**

*Your argument has a pedigree.* Simplicius reports, from Eudemus, that Archytas of Tarentum (4th century
BCE) asked whether a person at the edge of the heavens could stretch out a hand or a staff beyond it.
Either something blocks it or there is more space, and you repeat at the new limit. Lucretius later made
it a spear. Your ball is the same argument, and you reached it on your own. You also asked the
correct follow-up: what does "beyond" mean?

*Testing each branch fairly.*

| What happens | Your reading | Status |
|---|---|---|
| Ball is stopped | Something's there, so that wasn't the end of everything | Sound. A wall is more stuff. |
| Ball carries on | More space, so that wasn't the edge | Sound. |
| Ball comes back | "It bounced off something" | **The hidden assumption.** Return needs no bounce. |
| Ball's path just ends | Not on your list | Mathematically possible (a space with boundary, needing nothing outside it). No evidence, and physics gives no rule for it. |

*How return works with no force.* Picture a very large sphere and walk straight ahead. There is no
edge at any point, every direction looks the same, and you come back from behind. The surface has a
finite area and no boundary. Now do the same one dimension up. A 3-sphere of radius $R$ has metric

$$dl^2=R^2\left[d\chi^2+\sin^2\chi\,(d\theta^2+\sin^2\theta\,d\phi^2)\right].$$

A throw is a geodesic of this metric. It goes through the antipode at $\pi R$ and gets back at $2\pi R$.
Its covariant acceleration stays zero the entire way, so no force ever acts. Your intuition that
"coming back means it was pushed" is true in flat, infinite space, and this geometry is where it fails.

So when you say "space-time warps itself and there's no edge, then what happens when you throw a
ball?" the answer is: it keeps going straight, and you might one day see it arrive from behind you.

*Nothing outside is needed.* Curvature is measurable from inside. A triangle on a sphere has an angle
sum over 180°, by exactly area$/R^2$. The script checks this on the octant triangle, with 270° and
excess $\pi/2$. No outside point enters the calculation, and so "what's outside the sphere?" isn't a
question this geometry asks. In this picture "beyond the universe" is like "north of the North Pole".

*Curved isn't required for finite.* A flat torus has zero curvature, finite volume and no edge. It's
like a video-game screen that wraps. So flatness wouldn't prove infinity.

*What we measure.* $\Omega_K=0.001\pm0.002$ is consistent with flat. If closed at about $-2\sigma$
($\Omega_K=-0.003$), the curvature radius is at least 265 Gly and the circumference at least about
1660 Gly. Searches for repeated patterns in the CMB set $R_i>0.97\,\chi_\text{rec}$ for a cubic torus
(Planck 2015). Anything beyond the particle horizon can't be measured, even in principle, today.

*Would a real ball make it round?* Here's a check against the tidy picture, and the answer is no. In
closed matter-only models, light goes exactly once round over the universe's entire life from bang to
crunch, whatever the density (the script gets $2.00000\pi$ for three densities). With measured
$\Lambda$, light sent today covers 16.7 Gly at most, and a lap needs over 1660 Gly. A ball at $0.1c$
slows as momentum falls like $1/a$, and covers just 0.81 Gly forever. So the return works in the
geometry, but not in our actual expanding universe.

**Check:** [code/verify_geodesic_no_edge.py](code/verify_geodesic_no_edge.py) computes the geodesic from
Christoffel symbols of the metric alone, then uses the embedding only to see where the ball is. It shows
how such a geometry behaves. It can't tell which geometry we live in.

<details><summary>verify_geodesic_no_edge.py output</summary>

```text
Finite space without an edge
============================

-- Geodesic on a 3-sphere, built from the metric only
  path length walked = 2 pi R = 6.2832 R; distance from start at the end = 3.01e-12 R
  farthest from start along the way = 2.000000 R (antipode is 2 R away, reached at s = 3.1416 R)
  [PASS] forward-moving ball returns to its starting point after exactly 2 pi R
  [PASS] it goes as far as the antipode, so it really traversed the space
  [PASS] no jump anywhere on the path (step sizes uniform to 0.01%)
  covariant acceleration at mid-path = 4.0e-08 (zero means no force acted)
  [PASS] no force acted on the ball: covariant acceleration vanishes
  [PASS] three random throw directions all close after 2 pi R

-- Measuring curvature from inside: triangle angle excess
  octant triangle on a unit sphere: angle sum = 270.0 deg; excess = 1.5708; area/R^2 = 1.5708
  [PASS] angle excess equals area/R^2, detectable without any outside space

-- Flat torus: zero curvature, still finite and edgeless
  after T = 5.0: position [1. 2.], started at [1. 2.]
  [PASS] straight line on a flat torus also returns with no force and no wall

-- Measured curvature and what it allows
  |Omega_K| = 0.005: curvature radius c/(H0 sqrt|Omega_K|) = 205 Gly
  |Omega_K| = 0.003: curvature radius c/(H0 sqrt|Omega_K|) = 265 Gly
  Planck+BAO Omega_K = 0.001 +/- 0.002: at -2 sigma (closed side) Omega_K = -0.003, R > 265 Gly, circumference > 1664 Gly

-- Can a real ball circle a closed universe? Two tests of light
  closed matter universe, Omega_m = 1.5: light's angle over bang-to-crunch = 2.00000 pi
  [PASS] Omega_m = 1.5: light goes round exactly once in the universe's whole life
  closed matter universe, Omega_m = 2.0: light's angle over bang-to-crunch = 2.00000 pi
  [PASS] Omega_m = 2.0: light goes round exactly once in the universe's whole life
  closed matter universe, Omega_m = 5.0: light's angle over bang-to-crunch = 2.00000 pi
  [PASS] Omega_m = 5.0: light goes round exactly once in the universe's whole life
  light launched today covers at most 16.7 Gly of comoving distance; needed for a lap > 1664 Gly
  [PASS] with Lambda and the Planck curvature bound, not even light launched now can finish a lap
  ball thrown today at 0.1 c: total comoving distance ever = 0.81 Gly (it slows as space expands)
  [PASS] the ball's lifetime travel is finite and smaller than light's

12/12 checks passed
```

</details>

**Limits:** whether the whole universe is finite is Open and may stay so. The ±2σ curvature figures assume
ΛCDM with a Gaussian error. Whether "edge" could be meaningful at all in some deeper theory is Conceptual.

### C05: If the arrow of time exists, why consider past time travel at all?

**Status:** Conceptual + Established

**Short answer:** Because the arrow of time is a statistical fact about entropy, and it isn't a law
that forbids any path through spacetime. The field equations of general relativity have solutions
with closed timelike curves, so physicists have to ask whether those solutions can be real. Studying
them is largely how we test where the theory breaks.

**Starting from:** the second law as a statistical statement; the time-symmetric form of the known
microscopic laws apart from small weak-interaction effects (general knowledge); general relativity.

**Derivation:** a traveller on a closed timelike curve moves forward in their own light cone all the
way (shown explicitly in [C10](#c10-could-rotation-destroy-the-distinction-between-past-and-future)).
Their own clock runs forward and their own entropy can rise, so the local arrow is untouched. The real
tension comes in at the loop itself: a system that returns exactly to its starting state must have the
same entropy at the end as at the start. So a global arrow can't hold on a closed loop for an isolated
system. That's a constraint on what can circulate (see
[C20](#c20-if-one-lands-in-1984-do-people-of-1984-exist-normally)), and it isn't a proof that loops
can't exist. The arrow's existence tells us about the boundary condition of our universe. Whether that
condition rules out loops everywhere is what the chronology protection conjecture is about
([C17](#c17-could-wormholes-with-time-shifted-mouths-enable-past-travel)).

**Limits:** no experiment has seen any travel to the past. This answers why it's worth asking, not
whether it can happen.

### C06: Time passed at 99.95% of light speed

**Status:** Derived + Code + Established

**Short answer:** At $v=0.9995c$, $\gamma=31.6$. Each year on the traveller's clock is 31.6 years on
Earth. A 100 light-year trip takes 100.05 years by Earth clocks and 3.16 years on the traveller's.

**Starting from:** the speed of light is the same for all inertial observers, and the laws of physics
take the same form in all inertial frames.

**Derivation:** make a clock from a light pulse bouncing between mirrors a distance $L$ apart. At rest it
ticks every $2L/c$. Moving sideways at $v$, the pulse must follow a diagonal path, still at $c$, so in
half a tick $(c\,\Delta t/2)^2=L^2+(v\,\Delta t/2)^2$. Solving:

$$\Delta t=\frac{2L/c}{\sqrt{1-v^2/c^2}}=\gamma\,\frac{2L}{c},\qquad\gamma(0.9995)=31.63.$$

| Distance | Earth time | Traveller's time |
|---|---|---|
| 4.37 ly (Alpha Centauri) | 4.37 yr | 0.14 yr |
| 100 ly | 100.05 yr | 3.16 yr |
| 2.5 Mly (Andromeda) | 2.5 million yr | 79,000 yr |

A real trip has to speed up, and at a comfortable 1 g that takes 4 ship-years to reach this speed. A
round trip with four such burns and two one-year coasts, integrated step by step, ages the traveller
18.07 years and Earth 185.7 years, with a turnaround at 90.9 ly.

The muon experiment at CERN gives direct evidence at almost the same $\gamma$: lifetimes of stored muons
were measured at $\gamma=29.33$ and agreed with special relativity (Bailey et al. 1977). At the other
extreme, optical clocks detect the effect below 10 m/s (Chou et al. 2010).

**Check:** [code/verify_time_dilation.py](code/verify_time_dilation.py) counts light-clock bounces event by
event and integrates rapidity along the whole round trip, comparing with the closed form.

<details><summary>verify_time_dilation.py output</summary>

```text
Time dilation at 0.9995 c and the real hazards
==============================================
  gamma(0.9995) = 31.627

-- Light clock simulated event by event, no gamma in the code
  200 round trips: Earth-frame time 4.21982e-05 s; clock's own count 1.33426e-06 s
  measured ratio = 31.627
  [PASS] light-clock ratio matches 1/sqrt(1-beta^2) to 0.1%

-- Lorentz transform checks
  10 ship-years correspond to 316.3 Earth years and 316.1 ly
  [PASS] interval is invariant
  [PASS] Earth time = gamma * ship time

-- Coasting trips at 0.9995 c
  distance       4.37 ly: Earth         4.37 yr, traveller       0.14 yr
  distance        100 ly: Earth       100.05 yr, traveller       3.16 yr
  distance  2500000.0 ly: Earth   2501250.63 yr, traveller   79086.60 yr

-- A round trip with real acceleration: 1 g, integrated numerically
  peak gamma = 31.63, peak speed = 0.99950 c
  traveller aged 18.07 yr; Earth aged 185.74 yr; turned around at 90.9 ly
  back at start? final x = 2.49e-14 ly; acceleration phase alone: 4.02 yr proper time
  [PASS] numerical trip returns to its start
  [PASS] Earth ages far more than the traveller, as in the twin effect
  [PASS] integrated Earth time matches closed form (185.74 yr) to 0.1%

-- Biology check: nothing local changes (a chemical rate model, boosted)
  speed 0.0    : Earth sees half-life     0.693; proper time between events 0.693
  [PASS] proper half-life unchanged at beta = 0.0
  speed 0.9    : Earth sees half-life     1.590; proper time between events 0.693
  [PASS] proper half-life unchanged at beta = 0.9
  speed 0.9995 : Earth sees half-life    21.922; proper time between events 0.693
  [PASS] proper half-life unchanged at beta = 0.9995

-- Hazards at gamma = 31.6
  n = 0.1 H/cm^3: proton energy 28.7 GeV, flux 9.48e+10 /cm^2/s, beam power 4.4 MW/m^2, ionisation dose 30 Gy/s
  [PASS] unshielded dose at n = 0.1 is above 4 Gy (roughly lethal) within one second
  n = 1.0 H/cm^3: proton energy 28.7 GeV, flux 9.48e+11 /cm^2/s, beam power 43.6 MW/m^2, ionisation dose 304 Gy/s
  [PASS] unshielded dose at n = 1.0 is above 4 Gy (roughly lethal) within one second
  1 microgram dust grain: kinetic energy 2.75e+09 J = 0.66 tonnes TNT
  ideal photon rocket mass ratio per speed change = 63.2; round trip with stops = 1.60e+07
  kinetic energy per kg of ship = 2.75e+18 J/kg
  CMB head-on temperature in the ship frame = 172 K (a minor effect)
  Poisson count on 1 cm^2 per ship microsecond (ISM frame sweep) = 948020; expected 947671
  [PASS] frame-sweep count agrees with gamma*n*v flux

12/12 checks passed
```

</details>

**Limits:** ignores gravity and everything in the way. Those appear in [C08](#c08-does-physics-dominate-life-at-extreme-speed).

### C07: Would biology be overridden by physics in such travel?

**Status:** Derived + Code + Established

**Short answer:** No. Nothing is overridden, because the traveller's body is at rest in its own frame,
and everything local runs at the normal rate. Your heart beats at the usual rate by your watch. It's the
comparison with Earth's clocks that differs. Biology was never separate from physics, so physics has
nothing to override.

**Starting from:** the principle of relativity: no experiment inside a uniformly moving ship can detect
the motion.

**Derivation:** chemical reactions, nerve pulses and cell division are physics, and all of it obeys laws
that take the same form in every inertial frame. Therefore a body in the ship has the same chemistry as a
body at rest. The proper time between two events on the body's own path,
$\tau=\sqrt{t^2-x^2/c^2}$, is invariant. Take a reaction half-life and Lorentz-transform its start and end
events: Earth sees it stretched by $\gamma$ (a half-life of 0.693 becomes 21.92), and
the proper time recovered is still 0.693. You wouldn't feel the speed, since speed has no local meaning. You
can feel acceleration, which is why the trip above stays at 1 g.

**Check:** [code/verify_time_dilation.py](code/verify_time_dilation.py), output under
[C06](#c06-time-passed-at-9995-of-light-speed). The reaction check there confirms invariance of proper time,
and it isn't a biological experiment.

**Limits:** no organism has been tested at these speeds. The muon results show the effect for a decay
process, and the claim for biology rests on the principle of relativity holding for every interaction.

### C08: Does physics dominate life at extreme speed?

**Status:** Derived + Code + Estimate

**Short answer:** Not through speed itself. What kills the traveller is what they hit. At $\gamma=31.6$ the
thin gas of space arrives as a beam of 28.7 GeV protons, giving an ionising dose of order 30 to 300 gray per
second without shielding, which is lethal within a second. Then there's fuel.

**Starting from:** interstellar hydrogen density 0.1 to 0.2 per cm³ in the local cloud (arXiv:1910.01243), with 1 per cm³ as a round higher figure. Water stopping power about 2 MeV cm²/g and 4 Gy as a rough
lethal whole-body dose (both from general knowledge).

**Derivation:** in the ship's frame each proton has energy $(\gamma-1)m_pc^2=28.7$ GeV. The density is
Lorentz-contracted to $\gamma n$, so the flux is $\gamma n\beta c$.

| Quantity | $n=0.1$/cm³ | $n=1$/cm³ |
|---|---|---|
| Proton flux | $9.5\times10^{10}$ /cm²/s | $9.5\times10^{11}$ /cm²/s |
| Beam power on the bow | 4.4 MW/m² | 44 MW/m² |
| Ionising dose | 30 Gy/s | 304 Gy/s |

The dose uses ionisation only. Nuclear cascades add more, and shielding at this energy is a heavy problem of
its own (Semyonov 2009 treats this as the main radiation hazard).

Other costs: a 1 microgram dust grain carries $2.75\times10^9$ J, or about 0.66 tonnes of TNT. Kinetic energy
is $2.75\times10^{18}$ J per kg of ship. An ideal photon rocket needs a mass ratio
$\sqrt{(1+\beta)/(1-\beta)}=63$ for one speed change, and $1.6\times10^7$ for a round trip with stops. The CMB,
blueshifted to about 172 K, is a small effect. The acceleration needs are mild (1 g for 4 ship-years). The
body can take that; the engineering can't yet.

So your question has a real answer, in a different place: physics limits the trip, through radiation and
energy, and it does this by damaging biology as it would at any speed.

**Check:** [code/verify_time_dilation.py](code/verify_time_dilation.py), output under
[C06](#c06-time-passed-at-9995-of-light-speed). The flux is cross-checked by counting from the gas's frame as
well as the ship's.

**Limits:** order-of-magnitude. Density varies along any route. No shielding design is assessed.

### C09: What if the entire universe is rotating?

**Status:** Derived + Code + Established + Speculative

**Short answer:** It would be measurable, since rotation isn't relative the way velocity is. A gyroscope
would turn against distant galaxies. Measurements show no such rotation. Gödel (1949) found an exact solution
of Einstein's equations with a rotating universe, but it has features ours doesn't: no expansion and a
negative cosmological constant.

**Starting from:** general relativity. The Gödel metric, signature $(-,+,+,+)$, $c=G=1$:

$$ds^2=4a^2\left[-dt^2+dr^2+dz^2+(\sinh^2r-\sinh^4r)\,d\phi^2+2\sqrt2\,\sinh^2r\,d\phi\,dt\right].$$

**Derivation:** the script computes the full Ricci tensor symbolically, then finds all 16 components of
$G_{\mu\nu}+\Lambda g_{\mu\nu}=8\pi\rho\,u_\mu u_\nu$ vanish with dust $u=(1/2a)\partial_t$,

$$\Lambda=-\frac{1}{2a^2},\qquad 8\pi\rho=\frac{1}{a^2},\qquad\omega^2=\frac{1}{2a^2}=4\pi\rho=-\Lambda.$$

Gödel's universe is homogeneous, so there's no center, and every point is on an axis of rotation.

If it had our density (critical density $8.5\times10^{-27}$ kg/m³), then $\omega=2.7\times10^{-18}$ rad/s, a
period of 74 Gyr, and CTCs would open at a proper radius of 29.5 Gly. For scale, the rate is $\omega/H_0=1.2$.

Measured limits:

- Planck 2015: $(\omega/H)_0<7.6\times10^{-10}$ at 95% in a Bianchi VII$_h$ fit. That is $1.6\times10^{9}$ times
  smaller than the hypothetical above.
- Saadeh et al. 2016: vector-mode anisotropic expansion $(\sigma_V/H)_0<4.7\times10^{-11}$ at 95%, and overall
  odds of 121,000:1 against anisotropic expansion. This limits a shear mode, not rotation alone.
- Su & Chu 2009: under $10^{-9}$ rad/yr at last scattering, from the second-order Sachs-Wolfe effect.

**Check:** [code/verify_ctc_solutions.py](code/verify_ctc_solutions.py). Symbolic and numeric results:

<details><summary>verify_ctc_solutions.py output</summary>

```text
Solutions of GR with closed timelike curves
===========================================

-- Goedel metric: symbolic field-equation check (signature -+++, units c = G = 1)
  u^mu u_mu = -1, Ricci scalar R = -1/a**2
  dust plus Lambda solves the field equations with Lambda = -1/(2*a**2), 8 pi G rho = a**(-2)
  [PASS] all 16 field-equation components vanish
  [PASS] Lambda is negative: this is not the positive Lambda we measure
  vorticity^2 = 1/(2*a**2); 4 pi G rho = 1/(2*a**2)
  [PASS] Goedel's relation omega^2 = 4 pi G rho = -Lambda holds

-- Goedel: sign change of g_phiphi and the tilting light cone
  g_phiphi = 4 a^2 sinh^2 r (1 - sinh^2 r): zero at sinh r = 1, r_c = ln(1+sqrt 2) = 0.881374
  [PASS] numerical root of g_phiphi equals ln(1+sqrt 2)
    r/r_c  sinh r  g_phiphi/4a^2         dt/dphi cone edges  phi circle
     0.30   0.268         0.0664      -0.1757,       0.3781  spacelike
     0.60   0.554         0.2126      -0.1993,       1.0668  spacelike
     0.90   0.879         0.1756      -0.0776,       2.2633  spacelike
     1.00   1.000         0.0000       0.0000,       2.8284  null
     1.10   1.129        -0.3490       0.0996,       3.5037  TIMELIKE
     1.50   1.742        -6.1789       0.7929,       7.7929  TIMELIKE
     2.00   2.828       -56.0000       2.8284,      19.7990  TIMELIKE
  [PASS] lower cone edge crosses dt/dphi = 0 exactly at r_c: the cone has tipped past the phi circle

-- Goedel: an explicit closed curve, tested for timelike and future-directed
  r = 1.5 r_c, curve t = const, phi from 0 to -2 pi: g(T,T) = -24.7157 (negative = timelike), g(T,u) = -8.5858 (negative = future-directed)
  [PASS] the closed curve is timelike everywhere
  [PASS] it is future-directed: the traveller always moves forward in their own light cone
  [PASS] it closes: same t, r, z and phi mod 2 pi
  [PASS] the same loop at r = 0.8 r_c is spacelike, so no such loop there

-- Goedel: numbers if our universe's density were used
  density = critical = 8.53e-27 kg/m^3 -> omega = 2.68e-18 rad/s = 8.44e-11 rad/yr
  rotation period 2 pi/omega = 74 Gyr
  Goedel scale a = 16.8 Gly; proper CTC radius 2 a r_c = 29.5 Gly
  omega/H0 for this hypothetical = 1.22; Planck 2015 limit (omega/H)_0 < 7.6e-10; observed ratio is below hypothetical by 1.6e+09
  that CMB limit in today's rate: omega < 5.2e-20 rad/yr (Bianchi VII_h fit, model dependent)
  [PASS] a Goedel rotation at the real density is excluded by CMB limits by a factor > 1e8
  cheapest constant-r t=const loop: sinh^2 r = 1.366, total integrated acceleration = 27.67 (units of c)
  ideal photon-rocket mass ratio = exp(TA) = 1.04e+12 (Natario 2012 gives a lower-cost non-circular candidate)
  the traveller's own clock runs 149 Gyr per loop at that r
  [PASS] even an ideal rocket faces a large mass ratio

-- Kerr (M = 1, spin a): horizons and where g_phiphi < 0 on the equatorial plane
  a = 0.5: r+ = 1.866, r- = 0.134; g_phiphi < 0 for -0.689 < r < 0
  [PASS] a = 0.5: CTC region lies at r < 0, inside r- and r+, behind the horizons
  a = 0.9: r+ = 1.436, r- = 0.564; g_phiphi < 0 for -0.948 < r < 0
  [PASS] a = 0.9: CTC region lies at r < 0, inside r- and r+, behind the horizons
  a = 0.998: r+ = 1.063, r- = 0.937; g_phiphi < 0 for -0.999 < r < 0
  [PASS] a = 0.998: CTC region lies at r < 0, inside r- and r+, behind the horizons
  Over-spun a > M: 1 - a^2 < 0, no real horizon, so the same region would be exposed.

-- Gott strings: holonomy test in 2+1 dimensions
  G mu = 0.1: holonomy turns boost-like at beta = 0.30901699; cos(4 pi G mu) = 0.30901699
  [PASS] G mu = 0.1: numerical threshold reproduces Gott's beta > cos(4 pi G mu)
  G mu = 0.01: holonomy turns boost-like at beta = 0.99211470; cos(4 pi G mu) = 0.99211470
  [PASS] G mu = 0.01: numerical threshold reproduces Gott's beta > cos(4 pi G mu)
  G mu = 0.001: holonomy turns boost-like at beta = 0.99992104; cos(4 pi G mu) = 0.99992104
  [PASS] G mu = 0.001: numerical threshold reproduces Gott's beta > cos(4 pi G mu)
  G mu = 1e-07: need gamma > 1/sin(4 pi G mu) = 7.96e+05, 1-beta < 7.9e-13
  each string's energy must exceed its rest energy by 8.0e+05 (Planck limit G mu < 1.3e-7 used)

-- Morris-Thorne wormhole: symbolic energy condition at the throat
  8 pi G (rho + p_radial) = G_tt + G_ll = -2*Derivative(R(l), (l, 2))/R(l)
  [PASS] rho + p_radial = -(2 R''/R)/(8 pi G): flare-out R'' > 0 forces NEC violation
  Ellis throat R = sqrt(b0^2 + l^2): at l = 0, 8 pi G (rho + p) = -2/b0**2
  [PASS] the sign is negative for a concrete smooth throat
  for a 1 m throat in SI: rho + p/c^2 = -c^2/(4 pi G b0^2) = -1.07e+26 kg/m^3 (compare water, 1e3 kg/m^3)

-- Making a time machine from a wormhole: moving one mouth
  mouth round trip 10 yr outside at 0.9995 c: mouth clock 0.32 yr; time shift 9.68 yr
  mouth round trip 100 yr outside at 0.9995 c: mouth clock 3.16 yr; time shift 96.84 yr
  CTCs appear once shift exceeds the external separation D/c: shift 9.68 yr and D < 9.68 ly.

19/19 checks passed
```

</details>

**Limits:** these limits apply to the classes of model fitted. A small rotation has not been excluded at
zero. Exact solution here means exactly solving the equations, and it doesn't imply it describes anything real.

### C10: Could rotation destroy the distinction between past and future?

**Status:** Derived + Code + Conceptual

**Short answer:** Locally, no. Every traveller in Gödel's universe still has a clear past and future in
their own light cone. Globally, yes, in a precise way: there's no universal time, so no consistent labelling
of events in which every future-directed path has an increasing label.

**Starting from:** the [Gödel metric](#c09-what-if-the-entire-universe-is-rotating). Time-orientation
is given by the dust flow $u$.

**Derivation:** the script takes the loop $t=$ const, $r=1.5r_c$, $\phi$ going from $0$ to $-2\pi$, whose tangent
$T$ satisfies $g(T,T)<0$ (timelike) and $g(T,u)<0$ (future-directed, same side of the cone as the dust). So the
traveller is going forward the entire way and arrives at their starting event.

Now a one-line theorem. Suppose a function $f$ on spacetime increases along every future-directed timelike curve.
Follow a closed one through $p$. Then $f(p)>f(p)$, a contradiction. So **any spacetime with a CTC has no global
time function**. The past/future distinction survives as a local fact, and it fails as a global ordering. Gödel
drew philosophical conclusions from this about time's objectivity (see Savitt's review), which is a separate
argument from the calculation.

**Check:** [code/verify_ctc_solutions.py](code/verify_ctc_solutions.py), output under
[C09](#c09-what-if-the-entire-universe-is-rotating).

**Limits:** this applies within the Gödel solution, which doesn't describe our universe.

### C11: Could tilted light cones allow forward-time travel to earlier events?

**Status:** Derived + Code

**Short answer:** Yes, in those solutions, and "forward-time travel to earlier events" is a fair
description. In Gödel's universe the cones tip over with distance until the circle around the axis lies inside
the cone. Walking round then stays inside your local future and still arrives back.

**Starting from:** a null path satisfies $ds^2=0$. In the Gödel metric at fixed $r,z$:
$-dt^2+2\sqrt2 s^2\,dt\,d\phi+(s^2-s^4)\,d\phi^2=0$, where $s=\sinh r$.

**Derivation:** solving gives the two cone edges

$$\frac{dt}{d\phi}=\sqrt2 s^2\pm\sqrt{2s^4+s^2-s^4}=\sqrt2 s^2\pm s\sqrt{1+s^2}.$$

The lower edge is negative when $s<1$ and passes through zero at $s=1$, that is at $r_c=\ln(1+\sqrt2)=0.8814$.

| $r/r_c$ | $\sinh r$ | $g_{\phi\phi}/4a^2$ | Lower and upper cone edge $dt/d\phi$ | $\phi$ circle |
|---|---|---|---|---|
| 0.60 | 0.554 | +0.2126 | −0.199, 1.067 | spacelike |
| 1.00 | 1.000 | 0 | 0, 2.828 | null |
| 1.50 | 1.742 | −6.179 | 0.793, 7.793 | **timelike** |

Coordinate pictures of tilted cones can mislead. A tilt can be just a bad choice of coordinates, as with
Schwarzschild coordinates at a horizon. What makes this one real is that $\partial_\phi$ has closed orbits and
its own length $g_{\phi\phi}$ is a coordinate-independent scalar. When it's negative, the closed circle really
is timelike.

The price: the simple constant-$r$ loop needs integrated acceleration of 27.7 in units of $c$, or a photon-rocket
mass ratio of $e^{27.7}\approx10^{12}$. At our density it takes 149 Gyr on the traveller's own clock per lap. Natário
(2012) gives a cheaper candidate curve than Malament's conjecture, so this figure is only an upper bound on the
minimum.

**Check:** [code/verify_ctc_solutions.py](code/verify_ctc_solutions.py), output under
[C09](#c09-what-if-the-entire-universe-is-rotating). The root of $g_{\phi\phi}$ and the table come from it.

**Limits:** only this solution is tested. A tilt strong enough for CTCs isn't a generic result of rotation.

### C12: Could an infinite rotating cylinder create time loops?

**Status:** Established + Speculative

**Short answer:** On paper, yes. Tipler (1974) showed that the field of a rapidly rotating infinite
cylinder (the van Stockum family) violates causality, with closed timelike lines able to link any two events.
A real one doesn't exist since nothing is infinite, and a finite one with ordinary matter runs into Hawking's
theorem ([C17](#c17-could-wormholes-with-time-shifted-mouths-enable-past-travel)).

**Starting from:** Tipler's published result. It's cited here, not re-derived.

**Derivation:** the mechanism is frame dragging: the rotating mass tilts nearby light cones around the cylinder.
It's the same effect as Gödel's tilt in a different source. For hardware: a rim moving at speed $v$ has hoop stress
of order $\rho v^2$, so a rim speed near $c$ needs a stress near $\rho c^2$, which is the extreme limit for matter.
Ordinary materials manage $v$ of order a kilometre per second at most (Estimate). Then there's the infinite length,
which the exact solution actually uses.

**Limits:** Tipler's specific parameter values weren't reproduced here. A long but finite cylinder is the physical
question, and it's limited by the theorem cited, whose conclusions depend on its assumptions.

### C13: Do rotating black holes contain regions where causality breaks down?

**Status:** Derived + Code + Established + Open

**Short answer:** In the exact mathematical Kerr solution, yes. Carter (1968) showed that it has closed timelike
curves deep inside, near and past the ring singularity. Whether a real black hole has such a region is doubtful,
because the interior of real collapse is expected to differ from the exact solution.

**Starting from:** the Kerr metric in Boyer-Lindquist coordinates, mass $M=1$, spin $a$. On the equatorial plane,

$$g_{\phi\phi}=r^2+a^2+\frac{2a^2}{r}.$$

**Derivation:** the same test as Gödel: $\phi$ circles are closed, so where $g_{\phi\phi}<0$ they're timelike. For
$r>0$ all terms are positive. For small negative $r$, the last term dominates and is negative. In the analytically
extended solution, $r<0$ is reached through the ring. The script finds the boundary:

| Spin $a$ | $r_+$ | $r_-$ | $g_{\phi\phi}<0$ on the equator for |
|---|---|---|---|
| 0.5 | 1.866 | 0.134 | $-0.689<r<0$ |
| 0.9 | 1.436 | 0.564 | $-0.948<r<0$ |
| 0.998 | 1.063 | 0.937 | $-0.999<r<0$ |

So "causality breaks down" is a fair summary, in a region of the extended math.

**Check:** [code/verify_ctc_solutions.py](code/verify_ctc_solutions.py), output under
[C09](#c09-what-if-the-entire-universe-is-rotating). It checks only this formula's sign; the Kerr metric's validity
as a vacuum solution is taken from the literature.

**Limits:** the interior is likely replaced in real black holes. See [C14](#c14-are-such-regions-sealed-off-by-horizons).

### C14: Are such regions sealed off by horizons?

**Status:** Derived + Code + Established + Open

**Short answer:** For any spin below the maximum, yes. The CTC region lies inside both the outer and inner
horizons, so no signal from it reaches us. It might not exist at all in real black holes, because the inner
horizon seems to be unstable. Whether nature always keeps such things hidden is the cosmic censorship conjecture,
which is unproven.

**Starting from:** the table in [C13](#c13-do-rotating-black-holes-contain-regions-where-causality-breaks-down):
$r_\pm=1\pm\sqrt{1-a^2}$.

**Derivation:** $r<0$ is less than $r_-$ and $r_+$ for every spin $a<1$. If $a>M$ the root becomes complex, with no
horizon, and the region would be exposed. That's one reason censorship matters. Going in, the inner horizon is the
weak point. Poisson & Israel (1989) found mass inflation there. Dafermos & Luk show that, from suitable data inside a perturbed Kerr black hole, the
metric extends continuously across a piece of Cauchy horizon, and they argue that such horizons are generically weak null
singularities. So the smooth continuation needed to reach the CTC region is in doubt, and a revised form of strong censorship
may survive.

**Check:** [code/verify_ctc_solutions.py](code/verify_ctc_solutions.py) tests the location against the horizons.
Output under [C09](#c09-what-if-the-entire-universe-is-rotating).

**Limits:** "sealed" means sealed from outside. An infalling observer is a different question, and the answer depends
on an interior no one can observe. Cosmic censorship is Open.

### C15: Could spacetime topology alone create time loops?

**Status:** Derived

**Short answer:** Yes, as mathematics. Take flat empty space and declare $t$ and $t+T$ to be the same. You get
loops with zero curvature and no matter at all. This is why local curvature alone doesn't prevent them. It says
nothing about whether nature has such a shape, or about how to get there from a normal start.

**Starting from:** Minkowski metric, and the definition of chronological future $I^+(p)$: events reachable from
$p$ by future-directed timelike curves. These sets are open.

**Derivation:** *Example.* Identify $(t,\mathbf x)\sim(t+T,\mathbf x)$. The field equations are local, and the local
geometry is untouched, so this solves them in vacuum. A clock at rest traces a closed timelike curve of length
$cT$. Unwrap the identification and the loop disappears, so the loop is due to global shape only.

*Theorem: every compact spacetime contains a CTC.* Every point lies in some $I^+(p)$ (go slightly back along a timelike
curve), so the open sets $I^+(p)$ cover the space. By compactness finitely many suffice. Take a smallest such list $I^+(p_1),\dots,I^+(p_n)$.
Then $p_1$ lies in some $I^+(p_i)$. If $i\ne1$, then $I^+(p_1)\subset I^+(p_i)$, so $I^+(p_1)$ was redundant, contradicting
the smallest list. So $i=1$: $p_1$ lies in its own future, which is a CTC.

Geroch (1967) studied how topology can change in general relativity, a closely related question. The standard statement
is that such change needs causality violation or singular behaviour (general knowledge; only the citation was confirmed).

**Limits:** this says nothing about reachability. To start normal and end with loops is what
[C17](#c17-could-wormholes-with-time-shifted-mouths-enable-past-travel) addresses. No observation suggests periodic
time.

### C16: Could cosmic strings create closed timelike curves?

**Status:** Derived + Code + Established + Speculative

**Short answer:** Gott (1991) found exact solutions where two infinite strings passing fast enough allow CTCs. The
condition is $\gamma\sin(4\pi G\mu/c^2)>1$. For strings at the observational limit, each one needs $\gamma\sim10^6$.
Follow-up work shows it can't be set up in an open universe with physically allowed total momentum.

**Starting from:** in 2+1 dimensions, a string has a deficit angle $8\pi G\mu/c^2$, so going round it is a rotation by
that angle, which one can compute in flat space.

**Derivation:** going round both strings multiplies two rotations, each boosted to its string's frame. In
$\mathrm{SO}(2,1)$ a matrix is a rotation for trace $<3$ and a boost for trace $>3$. A boost-type total means a
spacelike total momentum, the condition for Gott's CTCs. The script doesn't use Gott's formula. It solves numerically for
when the trace of the product hits 3:

| $G\mu$ | Numerical threshold $\beta$ | $\cos(4\pi G\mu)$ |
|---|---|---|
| 0.1 | 0.30901699 | 0.30901699 |
| 0.01 | 0.99211470 | 0.99211470 |
| 0.001 | 0.99992104 | 0.99992104 |

So $\beta>\cos(4\pi G\mu)$, which is equivalent to $\gamma\sin(4\pi G\mu)>1$. Using the Planck bound
$G\mu/c^2<1.3\times10^{-7}$ for Nambu-Goto strings at about that scale, $G\mu=10^{-7}$ needs $\gamma>8\times10^5$. No
cosmic string has been observed.

Objections are serious. Boost-type total momentum is what Deser, Jackiw & 't Hooft (1992) flag: it's tachyonic.
Carroll, Farhi & Guth (1992) show that Gott's machine can't exist in an open (2+1)-dimensional universe with timelike
total momentum, and that decays don't supply enough mass.

**Check:** [code/verify_ctc_solutions.py](code/verify_ctc_solutions.py), output under
[C09](#c09-what-if-the-entire-universe-is-rotating).

**Limits:** only the holonomy condition is tested. It's a necessary algebraic condition for this construction, not a
full construction of the curves. Closed-universe cases are not ruled out by the same argument.

### C17: Could wormholes with time-shifted mouths enable past travel?

**Status:** Derived + Code + Established + Open + Speculative

**Short answer:** If a traversable wormhole exists, moving one mouth at high speed builds a time shift, and with a
large enough shift it becomes a time machine (Morris, Thorne & Yurtsever 1988). It needs matter that violates the null
energy condition, which is unlike anything we've seen on macroscopic scales. Hawking's theorem generalizes this:
any finite time machine with compact generation needs such a violation. His chronology protection conjecture, that
physics forbids it entirely, isn't proved.

**Starting from:** the wormhole metric $ds^2=-dt^2+dl^2+R(l)^2d\Omega^2$, with a throat where $R$ has a minimum.

**Derivation:**

*Energy condition.* The Einstein tensor computed symbolically gives

$$8\pi G(\rho+p_\text{radial})=-\frac{2R''}{R}.$$

A throat requires $R''>0$, so $\rho+p_\text{radial}<0$: the null energy condition fails at the throat. For a
1-metre Ellis throat this is $\rho+p/c^2=-1.07\times10^{26}$ kg/m³. Visser, Kar & Dadhich (2003) show the total amount
of violating matter can be made arbitrarily small, while Ford & Roman (1996) use quantum inequalities to argue that
macroscopic versions are very unlikely. Quantum fields can do this in small amounts, but nobody knows of a way to
sustain a large throat.

*The time shift.* Move a mouth out and back at $0.9995c$ for 10 years external. Its clock reads 0.32 yr, giving a
shift of 9.68 yr through the throat. Once the shift exceeds the outside travel time $D/c$, a closed path exists. It's
the twin effect from [C06](#c06-time-passed-at-9995-of-light-speed), put to a new use.

*Hawking's theorem.* Hawking (1992) treats causality violation in a finite region without curvature singularities.
There the boundary of the predictable region is a compactly generated Cauchy horizon. If the time machine develops
from a non-compact initial surface, the averaged weak energy condition must fail on the horizon. Sketch, not a full
proof: the horizon's null generators, followed backward, stay trapped in a compact set, and so must focus or wind up
on closed null geodesics. The Raychaudhuri equation
$\dot\theta=-\theta^2/2-\sigma^2-R_{\mu\nu}k^\mu k^\nu$ forbids the area behaviour this needs when
$T_{\mu\nu}k^\mu k^\nu\ge0$ along the generators. So the null contraction of stress-energy has to be negative
somewhere, which means the null energy condition fails. Whatever the mechanism, a finite machine of this kind needs exotic matter.

*Chronology protection.* Hawking also argued that quantum effects grow on the horizon. Kay, Radzikowski & Wald (1997) prove
that for a linear scalar field the stress tensor is ill-defined at base points of a compactly generated horizon, which supports
the idea but also reveals the breakdown of the approximation. Visser (2002) reviews the conjecture. Ori (2005) offers a model
with a compact vacuum core and matter that satisfies the standard energy conditions, regular and topologically trivial, with
stability open. Whether it's truly an exception depends on how its horizon is generated, and I haven't checked this
in detail. Hawking himself notes the final answer awaits quantum gravity.

*Summary of what's real.*

| System | Solves GR? | Known obstacle | Status as real object |
|---|---|---|---|
| Gödel | Yes | Needs negative $\Lambda$, no expansion, conflicts with rotation limits | Not our universe |
| Tipler | Yes | Infinite length, extreme rim speeds | Not buildable |
| Kerr interior | Yes | Behind horizons; inner horizon unstable | Doubtful |
| Periodic time | Yes | No origin from a normal start | No evidence |
| Gott | Yes in 2+1 | Tachyonic total momentum | Unobtainable in open universe |
| Morris-Thorne | Yes with exotic source | NEC violation | Unknown |

**Check:** [code/verify_ctc_solutions.py](code/verify_ctc_solutions.py) gets the throat formula and evaluates it,
output under [C09](#c09-what-if-the-entire-universe-is-rotating). The theorem isn't checked by code.

**Limits:** semiclassical arguments can't be trusted at the Planck scale. So physics hasn't proved time machines
impossible. It hasn't shown any are possible, either.

### C18: Are your ideas connected to block universe theory?

**Status:** Conceptual + Derived

**Short answer:** Yes. Asking "where is 1984 now, so that I can land in it?" treats past events as still being
somewhere, which is the core of the block universe: spacetime as a four-dimensional whole in which all events have
their place. C19 to C24 are easiest to think through in this picture, and that's a reason to use it, not proof that
it's correct.

**Starting from:** Minkowski spacetime with worldlines, and the relativity of simultaneity
([C25](#c25-why-the-block-universe-was-proposed)).

**Derivation:** a person is a worldline, a thin four-dimensional track. In GR, a time machine isn't a device for moving
"backward along" time. It's a track that's arranged so that later points on it (by its own clock) lie at earlier coordinates.
That only makes sense if earlier coordinates aren't erased. A strict presentist, for whom only the present exists, has
no 1984 to land in, unless they say the time-loop world differs in some deep way. Hence the link.

**Limits:** the converse fails. Block-universe views don't need or predict time machines. See
[C30](#c30-can-block-universe-ever-be-experimentally-proven) for what's unsettled.

### C19: Does time travel to the past mean total reversal of reality?

**Status:** Derived + Code + Conceptual

**Short answer:** No. Nothing runs backward. The traveller's clock runs forward and their body ages normally, and the
world around them runs forward too. They just arrive at an earlier date. Physics of this kind would be a loop in
spacetime, not a rewind.

**Starting from:** *assume* a closed timelike curve exists and that local physics holds everywhere along it.

**Derivation:** the explicit Gödel loop in [C10](#c10-could-rotation-destroy-the-distinction-between-past-and-future)
stays future-directed everywhere: $g(T,u)<0$ along it. A film of the whole world playing in reverse would have broken
eggs reassembling, and nothing in the solution does that. Consistency is the real constraint. If there's only one
history, the events of 1984 already include the traveller. Novikov and collaborators' self-consistency approach (Friedman et
al. 1990) says only self-consistent histories occur.

The classic test: a ball enters a wormhole mouth, emerges earlier, and knocks its younger self off course, so it never
enters. Paradox? I set this up in code with the undisturbed path being a perfect self-hit, and the undisturbed history does
fail. But searching for fixed points of the loop map finds a consistent one: the older ball emerges with a slightly different
speed and direction, gives its younger self a glancing hit, and that hit is exactly what sends the younger ball into the mouth
with that state. For these initial conditions:

- enters at height $-2.086$ with velocity $(1.184,-0.320)$, speed 1.2265
- residual $10^{-14}$ when replayed
- the older ball leaves at speed 1, which conserves energy in the collision

So what looks like a paradox has a solution for the same initial conditions, in line with the Echeverria, Klinkhammer &
Thorne work on billiard balls (only its existence was confirmed for this write-up).

Not every system has this. A bit that flips itself has no deterministic fixed point. It has a mixed one, $p=1/2$, and
Deutsch's (1991) quantum model always has a solution (tested on 200 random unitaries, and with the same answer $p=1/2$ for
the grandfather loop). Lloyd et al. (2011) give a physically different post-selection model, so the theory of CTCs isn't
unique.

**Check:** [code/verify_block_time.py](code/verify_block_time.py). It finds at least one solution; uniqueness isn't
tested, and it's a toy model of hard discs with a prescribed portal rule.

<details><summary>verify_block_time.py output</summary>

```text
Block universe inputs, jumps, and time-loop consistency
=======================================================

-- Andromeda calculation (Rietdijk-Putnam-Penrose)
  walker at +1.4 m/s: 'now' on Andromeda is t = +4.264 days in Earth's frame
  [PASS] root agrees with v D/c^2 for v = +1.4
  walker at -1.4 m/s: 'now' on Andromeda is t = -4.264 days in Earth's frame
  [PASS] root agrees with v D/c^2 for v = -1.4
  two walkers passing each other disagree by 8.53 days about 'now' at Andromeda
  [PASS] the disputed event is spacelike separated from the walkers
  Andromeda event at t = 0: t' = -4.264 days for one walker, +4.264 days for the other
  [PASS] the order of spacelike events flips between frames
  [PASS] 10000 random timelike pairs: no boost reverses their order, so causes stay first

-- Why a global jump breaks conservation in another frame (Feynman's argument)
  frame S: vanish and reappear both at t = 0, total charge 1 at all times
  frame S': vanish at t' = +0.000 s, reappear at t' = -0.577 s
  total charge on slices t' = [-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]:
    [1, 1, 2, 2, 1, 1, 1, 1, 1]
  [PASS] in S' the total charge changes: a jump conserves charge in at most one frame

-- Random jumps: size of the relativistic-QFT suppression
  electron, 1 nm    : m c r / hbar = 2.59e+03; amplitude ~ 10^(-1.12e+03)
  electron, 1 m     : m c r / hbar = 2.59e+12; amplitude ~ 10^(-1.12e+12)
  70 kg body, 1 m   : m c r / hbar = 1.99e+44; amplitude ~ 10^(-8.64e+43)
  [PASS] for a body the exponent exceeds 1e40

-- Billiard ball through a time-shifted wormhole: search for consistent histories
  straight-line history: older ball hits young ball? True (at t = 4.293)
  after that hit the young ball never reaches the portal: naive history is inconsistent
  [PASS] the undisturbed history is NOT self-consistent (the paradox)
  distinct self-consistent colliding histories found from 400 starts: 1
    enters at h = -2.0858, velocity (1.1839, -0.3204), speed 1.2265, t_in = 8.803, residual 1.9e-14
  [PASS] at least one self-consistent history exists for the very initial data that looked paradoxical
    older ball emerges at t = -1.197; collision at t = 2.293; young ball started at t = 0
  [PASS] each history replays forward to itself, and the collision conserves energy (older ball exits at speed 1)

-- Discrete loops: classical bit, then a Deutsch-style density matrix
  deterministic fixed points of NOT: []
  [PASS] a bit-flip loop has no deterministic consistent history
    map copy : fixed points [0, 1]
    map set0 : fixed points [0]
    map set1 : fixed points [1]
    map not  : fixed points []
  stochastic fixed point of NOT, by scanning p: p = 0.50000; points with gap < 1e-9: 1
  [PASS] p = 1/2 is the unique stochastic fixed point of NOT
  Deutsch CTC state for grandfather (NOT) loop: diag = [0.5, 0.5]
  [PASS] Deutsch's condition gives a consistent mixed state for the grandfather loop
  [PASS] 200 random unitaries: Deutsch's equation always has a valid density-matrix solution

14/14 checks passed
```

</details>

**Limits:** conditional on CTCs existing. What doesn't follow: that any of it is possible, or that a traveller could
change what happened.

### C20: If one lands in 1984, do people of 1984 exist normally?

**Status:** Conceptual + Derived

**Short answer:** Within the assumption, yes. They're ordinary people living ordinary lives forward, unaware unless you
show up. They aren't ghosts or recordings. If you met them, that meeting was always part of 1984 in the single-history
version.

**Starting from:** *assume* one self-consistent history (Novikov type), and local physics unchanged.

**Derivation:** events in that region obey the same laws as any others. A person in 1984 has an ordinary worldline, and
your arrival is a new input to their region, just as a stranger walking in the door is. In the single-history version, the
1984 you remember always contained your visit, whether or not you knew. There's a bookkeeping point: in the billiard run an
extra ball is present between $t=-1.2$ and $t=8.8$, and collisions conserve energy locally. In a wormhole version, the mouths
account for this, by mass changes (general knowledge; not verified).

In a Deutsch-type branching model, you could land in a 1984 which diverges from the one you remember. Those people would
still live normally, and you couldn't count on your memories of it.

**Limits:** purely conditional. There's no evidence that anyone has done this, and the choice between models isn't settled.

### C21: If you weren't born yet, do you exist there?

**Status:** Conceptual + Derived

**Short answer:** Within the assumption, yes, you'd exist as a visitor. The baby you doesn't exist yet, and the adult you
standing there does. No contradiction, because "you exist at date $t$" gets answered by whether your worldline passes through
that date, and yours does, by the route of the time machine.

**Starting from:** *assume* a CTC-style trip from after your birth to before it, with one history.

**Derivation:** your age at an event is measured along your worldline, not by calendar date. So you can be 30 by your own
clock at a date earlier than your birth. Your atoms are a nice detail: atoms in your body in 1984 are also elsewhere in 1984,
as air, soil and food, on their way to you. In that case, the same atoms are present twice on that time slice. Because
identical atoms are indistinguishable and the conservation law is local, this is allowed, with the time machine accounting for
the books (see [C20](#c20-if-one-lands-in-1984-do-people-of-1984-exist-normally)).

Doing something that would stop your birth is the grandfather paradox. Per
[C19](#c19-does-time-travel-to-the-past-mean-total-reversal-of-reality), in the single-history version it doesn't happen.
Your attempts fail for ordinary reasons. Nothing special stops you; it's just that this history isn't one of the consistent
ones.

**Limits:** conditional. What doesn't follow: that it's possible, or that free attempts are "prevented" by some extra force.

### C22: What happens to the 2022 version of you?

**Status:** Conceptual

**Short answer:** Nothing happens to it. In the single-history version the 2022 you who steps into the machine is you, at
an earlier stage on your worldline. It isn't replaced, frozen or deleted. If you go from 2022 to 1984 and then live on, you
reach 2022 again as an old person, and you can watch your younger self leave. For a while, both exist.

**Starting from:** *assume* a trip from 2022 to 1984, and that you survive afterward. One self-consistent history.

**Derivation:** the worldline has two stretches that overlap in date. Birth to 2022, then 1984 onward. Take a person who's 20
at departure, born in 2002. The traveller is 20 at arrival in 1984, then 38 when the baby is born, and 58 in 2022 when the
20-year-old leaves. Both are present in 2002 to 2022 as two bodies. After the 20-year-old departs, only the older one
continues, and nothing else changes. A branching (Deutsch-type) version is different: in the branch you arrive in, the
younger you may never be born, while the branch you left carries on without you.

**Limits:** conditional. Both versions can't be chosen between without knowing whether CTCs exist and what theory applies.

### C23: Does future travel imply existing in a time where you are already dead?

**Status:** Derived + Conceptual

**Short answer:** No, unless you took a shortcut around your own death. Forward travel by time dilation means you skipped
ahead while staying alive. You simply haven't died yet: your death, if it comes, is later on your worldline. The date on a
calendar doesn't decide when you die. Your own clock does.

**Starting from:** the trip in [C06](#c06-time-passed-at-9995-of-light-speed): proper time differs between paths.

**Derivation:** the 1 g round trip gets you back after 18 years on your clock, while 186 Earth years have passed. In the
Earth records you were never recorded as having died at age 80, because that event didn't occur. Your friends are gone, and
people may presume you dead. That's the real cost. A genuinely strange case needs a CTC: if you first lived to old age and
died, and a younger stage of you later travelled forward, then you could see your grave. In the block picture that's not a
contradiction, because your younger stage and your death are different points on one track, but it needs a machine beyond
dilation, and it is purely conditional.

**Limits:** forward travel is tested physics, on clocks and particles. No human has gone more than a small fraction of a second
ahead this way (general knowledge).

### C24: Why don't random spacetime jumps happen?

**Status:** Derived + Code + Estimate

**Short answer:** Because conserved things can only move continuously. A jump conserves in at most one frame: in another
frame, the object has vanished before it reappears, or reappears before it vanishes. For a large body, quantum effects
are suppressed by an exponent of about $10^{44}$. Worldlines being continuous is built into all the physics we've tested.

**Starting from:** local conservation $\partial_\mu J^\mu=0$; Lorentz transformations; the known form of quantum propagation.

**Derivation:** suppose a charge vanishes at $x=0$ and appears at $x=L$ at the same time in frame $S$. The total is conserved in
$S$. In a frame moving at $v$ those two events are at different times: with $L=1$ light-second and $v=0.5c$, reappearance is
0.577 s before vanishing. On slices in between, the charge is 2. So a jump breaks conservation in every frame other than one, and
no frame is special. This is Feynman's argument for why conservation has to be local (Feynman Lectures, general knowledge). The
only way to have the charge disappear from here and appear there is by flowing in between, which is continuity.

Quantum mechanics doesn't rescue random jumps. For separation $r$ outside the light cone, the free field amplitude falls as about
$e^{-mcr/\hbar}$. For an electron at 1 nm the exponent is $2.6\times10^3$, and for a 70 kg body at 1 m it's $2\times10^{44}$, so the
amplitude is about $10^{-8.6\times10^{43}}$. Treating a body as one particle is crude, though it's generous.

Also, jumping to another *time* would need the track to go somewhere outside the light cone, which brings its own causal problems.

**Check:** [code/verify_block_time.py](code/verify_block_time.py) counts charge on slices in a boosted frame and evaluates the
exponents. Output under [C19](#c19-does-time-travel-to-the-past-mean-total-reversal-of-reality).

**Limits:** this describes tested physics. A deeper theory could allow new things, but it would have to reproduce conservation where
it's measured.

### C25: Why the block universe was proposed

**Status:** Derived + Code + Established + Conceptual

**Short answer:** It grew from relativity of simultaneity. Minkowski (1908) rewrote relativity as four-dimensional geometry. Then
Rietdijk (1966) and Putnam (1967) argued: observers in relative motion disagree about which distant events are "now", so if "now" is
what's real, then future events must already be real. Penrose (1989) made the striking Andromeda version. This is an argument with
philosophical premises, not a consequence of the equations alone.

**Starting from:** the Lorentz transformation $t'=\gamma(t-vx/c^2)$.

**Derivation:** set $t'=0$ at a distant place $x=D$ and solve: $t=vD/c^2$. Different $v$ means different "now".

Take two people walking past each other at $\pm1.4$ m/s, with Andromeda at $D=2.5$ Mly. Then

$$\frac{vD}{c^2}=4.26\text{ days},$$

so they disagree by 8.53 days about what's happening "now" in Andromeda. One's "now" could include a fleet launch that for the other
is still days ahead. Neither can see or affect it, since the events are spacelike separated, and that means nothing observable is in
conflict. The script solves the transform by root-finding, and shows order flips between frames for such events. It also tests 10,000
random timelike pairs: no boost reverses cause and effect.

The further step: (a) what's simultaneous with me is real for me, (b) "real for" is transitive, (c) no observer is privileged. Then every
event is real. Weyl's line sums up the attitude: "The objective world simply is, it does not happen" (Weyl 1949, quoted from secondary
sources). Stein (1968) challenged the premises: he defines a relation of "having become" using the past light cone. It is transitive, yet it
doesn't hold between all pairs of events, so it doesn't make everything real. Savitt's review reports this as formally undercutting the argument,
and notes objections to it. He also says eternalism is "optional" in a related argument. So relativity removes a global, frame-free "now". Whether that yields a block is a further step.

**Check:** [code/verify_block_time.py](code/verify_block_time.py), output under
[C19](#c19-does-time-travel-to-the-past-mean-total-reversal-of-reality). It checks the arithmetic of simultaneity, not the philosophy.

**Limits:** Minkowski's text is cited from general knowledge. General relativity has no global Lorentz frames, so the argument has to be
restated there.

### C26: If block universe is true, what do choices mean?

**Status:** Conceptual

**Short answer:** Physics can't settle this. A common answer is that choices remain real events, part of what causes what comes after.
In a block, your deliberation is in the block, and tomorrow's outcome depends on it. Whether that amounts to the freedom people care
about is a debate in philosophy, and there's no result that ends it.

**Starting from:** a block universe, taken as an assumption. Also a distinction: the block is about what *exists*, and determinism is
about what laws *fix*.

**Derivation:** the strongest argument is by counterfactual dependence. If you'd chosen differently, later events would differ, by the same
laws that link cause to effect in any picture. So the fatalist line "it'll happen whatever I do" doesn't follow. It'll happen *because* of what you
do. People who hold this view are called compatibilists, and others argue that it's not enough, since real choice needs open alternatives.
Block also isn't the same as determinism: a block could contain outcomes that no law fixed in advance. Quantum measurements, in some interpretations, are like this. So
"the future exists" doesn't mean "the future can be computed". This is an argument, and not a proof.

**Limits:** no experiment distinguishes these views. Neuroscience of decisions addresses timing, not this conceptual point.

### C27: What exists "before" the origin of time?

**Status:** Open + Speculative + Derived + Code

**Short answer:** Nobody knows, and maybe the question lacks meaning. In the simple model, time runs out at a finite past, so there's no
"before" within the model. But the model isn't reliable there. Serious proposals include time having no boundary, a bounce, or a larger
pre-existing spacetime.

**Starting from:** the Friedmann equation running backward.

**Derivation:** what's solid: $t=\int da/(aH)$ converges as $a\to0$ for matter and radiation, so the model has finite age (script: at
$a=10^{-12}$, age $7.6\times10^{-13}$ yr). Inflation doesn't by itself remove the problem: Borde, Guth & Vilenkin (2003) show a space that
expands quickly enough on average must be incomplete to the past, and physics beyond inflation is required for that boundary. Their argument
is kinematic, needing no energy condition.

Options, all unverified by data:

- **No boundary** (Hartle & Hawking 1983): time becomes space-like near the start, so asking what's before is like asking what's north of the
  North Pole.
- **Bounce or cycles:** an earlier contracting phase. It has to deal with the entropy issue and the BGV condition.
- **Larger spacetime:** our region started from something older, with the Big Bang local.

There's also a low-entropy puzzle: whatever came "first" needs to explain why entropy started low, which [C31](#c31-what-if-we-are-boltzmann-brains)
shows is a serious constraint. That's really the question hidden in [C01](#c01-why-the-universe-expands-into-what-and-how-we-know).

**Check:** [code/verify_cosmology.py](code/verify_cosmology.py) confirms the finite integral only. Output under
[C02](#c02-93-billion-light-years-in-138-billion-years).

**Limits:** none of these is established. A full theory of quantum gravity may change the question itself.

### C28: Does block universe imply parallel universes?

**Status:** Conceptual

**Short answer:** No. A block universe is one spacetime with all its events. Parallel universes are additional spacetimes or branches.
You can accept either without the other, or both, or neither.

**Starting from:** the definitions. Block: one four-dimensional manifold. Many-worlds: a branching quantum state. Multiverse: other regions of
a bigger spacetime.

**Derivation:** the logic of it is independent. Relativity's argument for the block never mentions other worlds. The Everett view is a response to
quantum measurement, not simultaneity. Eternal inflation's many regions depend on inflation's dynamics. If all were true they'd fit together: an
Everett universe could be a block in a much larger sense. But none implies another. People confuse them, I think, because "other times exist" sounds like
"other worlds exist".

**Limits:** interpretation questions, with no experiment that distinguishes them currently accepted.

### C29: Does block universe validate religious ideas of a written destiny?

**Status:** Conceptual

**Short answer:** Physics can't validate or refute this. The block view says all events have a place. It doesn't say anyone wrote them, read them or
planned them, and a religious idea of destiny usually involves such things. Traditions differ, and many believers and non-believers have reasoned
carefully about this. It isn't for physics to settle.

**Starting from:** a block universe assumed true, and no assumption about religion either way.

**Derivation:** there's a similarity worth noting honestly: both say the future isn't nothing. There are differences too. The block has no author,
purpose or message in its equations. It's not readable from within either: chaos, quantum uncertainty and light cones prevent anyone from looking
ahead, and nothing allows you to look up your future. Also the block alone doesn't imply outcomes are fixed by law
([C26](#c26-if-block-universe-is-true-what-do-choices-mean)). So a person of faith can see in it a picture that fits their view, and a skeptic can
see a mere geometry. Neither is shown right by physics. And the block view itself isn't established
([C30](#c30-can-block-universe-ever-be-experimentally-proven)), so it couldn't serve as proof anyway.

**Limits:** this addresses what physics implies, and it doesn't evaluate any faith's claims.

### C30: Can block universe ever be experimentally proven?

**Status:** Conceptual + Open

**Short answer:** Probably not directly. Eternalism and presentism are widely regarded as hard to separate empirically, because each tends to be
stated so as to predict the same observations. Experiment can constrain particular versions, for example ones with a detectable preferred frame.

**Starting from:** the distinction between a physical theory and its metaphysical reading.

**Derivation:** relativity's successes confirm the geometry of spacetime, and they don't tell you if other times "exist". A presentist can keep a
hidden preferred frame, and a neo-Lorentzian version reproduces special relativity's predictions. Balashov & Janssen (2003) argue against this approach, on grounds
of explanation, not new data. Dorato (2006) argues the debate is ill-posed and there's no clear
way to state a disagreement between the two. Savitt's review reports no consensus, with options such as local becoming.

What could count: a found preferred frame would hurt a pure block motivation, though it wouldn't refute it logically. A confirmed CTC would strengthen it.
A theory of quantum gravity might make time emergent, reshaping the question. None of this is in hand.

**Limits:** "can't be tested" is a claim about current formulations, and it might not be permanent.

### C31: What if we are Boltzmann brains?

**Status:** Estimate + Code + Conceptual + Open

**Short answer:** It can't be refuted by experiment, since a fluctuation could produce a brain with your exact memories of every experiment. It's rejected
for a reason about reasoning: if you were one, the evidence you'd use for believing it, including all of physics, would be unreliable, so the belief
undermines itself. That is a good argument, but it isn't a proof that you aren't one.

**Starting from:** Boltzmann statistics $P\propto e^{\Delta S/k}$, measured $\Lambda$ with $\Omega_\Lambda=0.685$, and a de Sitter horizon with temperature
$T_\text{dS}=\hbar H/2\pi k$ and entropy $S_\text{dS}/k=\pi L^2/\ell_P^2$. Those last two are cited general results. Brain: 1.4 kg treated as water, 70 J/mol/K
(general knowledge).

**Derivation:**

| Quantity | Value |
|---|---|
| Horizon radius | $1.66\times10^{26}$ m |
| $T_\text{dS}$ | $2.2\times10^{-30}$ K |
| $S_\text{dS}/k$ | $3.3\times10^{122}$ |
| Brain's own thermal entropy | $3.9\times10^{26}\,k$ |
| Entropy cost of making that mass-energy from horizon, $E/kT_\text{dS}$ | $4.2\times10^{69}$ |
| Probability per trial | $\exp(-4\times10^{69})$, so $\log_{10}P\approx-1.8\times10^{69}$ |
| Recurrence time | $\exp(3.3\times10^{122})$, any units |

The numbers are absurd, but a ratio is what counts. Exponent of the recurrence exceeds the brain's by a factor of $8\times10^{52}$. So if the phase lasts anywhere near
recurrence, then freak brains vastly outnumber ordinary observers. A lower "information-only" cost, such as $10^{16}$ bits or $\exp(-7\times10^{15})$, changes the exponent but
doesn't change that conclusion. This figure is a rough scale.

The logic that smaller fluctuations are more common is testable in a toy. The Ehrenfest urn measured below matches $\exp(\Delta S)$ over three decades. A dip to 12 is 63 times more
common than a dip to 8. Applied to us: a lone brain costs much less than an entire ordered universe. So if observers arise by fluctuation, you should expect to be a lone brain with
a false memory. This goes back to Boltzmann's ideas (general knowledge), and in modern form to Albrecht & Sorbo (2004) and Dyson, Kleban & Susskind (2002).

*Why people reject it.* Carroll (2017) calls such theories "cognitively unstable": they can't be both true and justifiably believed. If you were a Boltzmann brain, your reasons
for thinking so would be random noise. Equal likelihood means no observation can shift the odds, as the check shows. So models predicting many freak observers get rejected as
unsuited to be believed, not as contradicted by data. That is a judgement about reasoning, and honest people could reject the premise.

**Check:** [code/verify_boltzmann_brains.py](code/verify_boltzmann_brains.py). It tests statistics in a toy, and evaluates the formulas. It doesn't test whether de Sitter space does this.

<details><summary>verify_boltzmann_brains.py output</summary>

```text
Boltzmann brain estimates
=========================

-- Late de Sitter phase
  H_inf = 1.808e-18 /s; horizon radius = 1.658e+26 m; T_dS = 2.20e-30 K
  de Sitter entropy S_dS/k = pi L^2/l_P^2 = 3.31e+122
  [PASS] S_dS is about 3e122 for the measured Lambda

-- Brain estimates (inputs stated, order of magnitude only)
  thermal entropy of a 1.4 kg brain as water: 3.94e+26 k
  energy 1.26e+17 J at T_dS removes 4.15e+69 k from the horizon
  P(one brain in the bath per trial) ~ exp(-4.1e+69): log10 P = -1.80e+69
  [PASS] horizon cost dwarfs the brain's own entropy
  atoms in brain ~ 1.4e+26
  low-end 'information only' cost, 1e16 bits: exp(-6.9e+15) (rough, shows the range)

-- Comparison with recurrence time
  recurrence time ~ exp(S_dS): log10 t_rec ~ 1.44e+122
  brain waiting time per horizon volume: log10 t_BB ~ 1.80e+69
  recurrence/brain exponent ratio = 8.0e+52
  brains per recurrence: exp(3.31e+122 - 4.15e+69) = exp(3.31e+122), which is huge
  [PASS] if the phase lasts a recurrence time, freak brains vastly outnumber ordinary observers
  brains avoided only if the vacuum decays in time < exp(4.1e+69) (any unit), far longer than 13.8 Gyr

-- Direct test of the statistics: Ehrenfest urn
  N = 40, steps = 3000000
   n up     S/k   measured P   exp(S-N ln2)
     20   25.65    1.258e-01      1.254e-01
     15   24.42    3.658e-02      3.658e-02
     12   22.44    5.066e-03      5.081e-03
     10   20.56    7.788e-04      7.709e-04
      8   18.16    8.083e-05      6.994e-05
  [PASS] measured macrostate frequency matches exp(Delta S) over 3+ decades (within 25%)
  smaller dip (n=12) is 63 times as common as larger dip (n=8); predicted 73
  [PASS] smaller fluctuations are far more common than larger ones

-- The undefeatable part, as a Bayes check
  priors [0.001, 0.5, 0.999] -> posteriors [0.001, 0.5, 0.999] when likelihoods are equal
  [PASS] equal likelihood means no experiment can shift the odds

6/6 checks passed
```

</details>

**Limits:** order-of-magnitude. Whether $\Lambda$ is constant, and whether de Sitter vacuum fluctuates dynamically, are unknown.

### C32: Does block universe increase this possibility?

**Status:** Conceptual

**Short answer:** No. How many Boltzmann brains there are depends on the dynamics and how long the universe lasts, which block-or-not doesn't change. A block view only says that
all of them, past and future, exist in the same sense. The count doesn't change.

**Starting from:** [C31](#c31-what-if-we-are-boltzmann-brains), [C25](#c25-why-the-block-universe-was-proposed).

**Derivation:** there's one way they do interact. "How many observers like me exist?" is a tenseless question, and a block view makes it natural to count across all time. A
presentist might say future brains don't exist yet. But then the reasoning "which observer am I likely to be?" gets awkward in that view too, and most use time-summed counts
regardless. So the block changes the framing and not the physics.

**Limits:** the reasoning from counts to "what should I expect?" is contested in itself.

### C33: Could Boltzmann brains arise from spacetime or QFT fluctuations?

**Status:** Open + Speculative

**Short answer:** It depends on an unresolved question. Thermal fluctuations of a de Sitter horizon are the usual route. Boddy, Carroll & Pollack argue that in many cases, de Sitter
vacuum is quiet: stationary states have no dynamical fluctuations producing recordable events. If so, brains don't form. If the patch has a finite-dimensional Hilbert space, they say
that recurrences and such fluctuations do occur.

**Starting from:** quantum field theory in de Sitter space.

**Derivation:** the "vacuum fluctuations" of QFT aren't things happening over time in a stationary state. They're correlations, shown by the Casimir effect. For an event to occur
there has to be something that evolves and records. Boddy et al. argue that's the missing part, and that a field relaxes to the vacuum by the cosmic no-hair theorem. De Simone et al.
(2010) estimate nucleation rates under various assumptions, treating brains as real. So experts disagree about this. The Estimate in C31 assumes the fluctuating view.

**Limits:** decided by quantum gravity, which is not known.

### C34: Could parallel universes generate such observers?

**Status:** Speculative + Open

**Short answer:** In eternal-inflation multiverse models, yes, and it's a known problem for them. The ratio of Boltzmann brains to normal observers depends on the rule used to count in
infinite spacetimes, and this "measure problem" isn't solved. In fact this ratio is used as a test of proposed rules.

**Starting from:** eternal inflation as an assumption, with its infinite volume.

**Derivation:** an infinite multiverse has infinitely many of every kind of observer, so ratios aren't defined until a cutoff is chosen. De Simone et al. (2010) use a
scale-factor cutoff, get a finite ratio, express it through nucleation and vacuum-decay rates, and propose that a viable measure must not be dominated by freak observers. Note
how this works: the absence of a Boltzmann-brain majority is an assumption made to select theories, on the reasoning in C31. It's not observed. The Everett view is separate, and it
doesn't add new fluctuating brains by branching alone (Boddy et al. argue that perturbations in a vacuum don't decohere).

**Limits:** every step depends on inflationary and landscape physics that hasn't been tested.

### C35: If a foundation goes wrong, does 90% of physics collapse?

**Status:** Derived + Code + Established + Conceptual

**Short answer:** Not in practice. Your logic is right for a chain of pure deductions, but physics results are held up by measurements as well as by their premises. When Newton's
foundations failed (absolute time, instant gravity), almost every number he gave stayed right in the range where it had been tested. The picture changed, and the measured
numbers stayed.

**Starting from:** a theory's tested domain. The relevant expansion parameters are $v^2/c^2$ and $GM/rc^2$.

**Derivation:**

*The old theory comes out as a limit.* Expanding:

$$E=\frac{mc^2}{\sqrt{1-v^2/c^2}}=mc^2+\frac12mv^2+\frac38\frac{mv^4}{c^2}+\dots$$

and the Lorentz transform becomes the Galilean one as $c\to\infty$. So Einstein has to contain Newton. This is a requirement any replacement has to meet.

*Size of the corrections.*

| Case | Relativistic correction |
|---|---|
| Car at 30 m/s | $7.5\times10^{-15}$ |
| Earth's orbit | $9.9\times10^{-9}$ |
| Sun's surface potential | $2.1\times10^{-6}$ |
| Hydrogen electron | $5.3\times10^{-5}$ |
| GPS clocks | $4.5\times10^{-10}$, which still matters because they demand it |
| LHC proton | $7\times10^{3}$, Newton fails totally |

This table is hand-picked to show how it works. It isn't a survey, and so it doesn't give a figure to set against your 90%.

*Worked example: GPS.* With $a=26561.75$ km and the geoid potential at the equator:

$$\frac{\Delta f}{f}=\frac{\Phi_\text{sat}-\Phi_\text{geoid}}{c^2}-\frac{v^2}{2c^2}=+5.30\times10^{-10}-0.83\times10^{-10}=+4.465\times10^{-10},$$

which is $+45.8-7.2=+38.6$ µs/day. Physics Today's Ashby article gives 4.4645 parts in $10^{10}$. Left uncorrected, that's about 11.6 km per day of ranging error. Simulating an eccentric
orbit with $e=0.01$ and accumulating proper time gives the same mean, and a periodic term of 22.9 ns matching the formula.

*Worked example: Mercury.* Integrating the relativistic orbit equation

$$u''+u=\frac{GM}{L^2}+\frac{3GM}{c^2}u^2$$

gives $5.02\times10^{-7}$ rad per orbit, or 42.98 arcsec/century, which matches Will (2014). With the term off, it gives no shift. The measured total is $575.31\pm0.0015$ arcsec/century (Park
et al. 2017). So the Newtonian part supplies 532.33, or 92.5%, by subtraction. Einstein didn't overturn this: he fixed the remaining 7.5%.

*What survived Newton → Einstein.* Numbers for the orbits, tides, engineering, and thermodynamics. Conservation of energy and momentum, restated. Maxwell's equations, which were already
relativistic. *What didn't:* the concepts of absolute time, gravity as force, and additive velocity at high speed. So a dependency graph has two kinds of nodes: concepts and numbers. A fallen
concept drags the ones above it, but the numbers that were measured stay measured.

Where your worry does apply: in deductions with no independent test. A long chain of theory resting on a single untested assumption can fall all at once.

**Check:** [code/verify_effective_theories.py](code/verify_effective_theories.py) integrates both systems instead of just reusing formulas.

<details><summary>verify_effective_theories.py output</summary>

```text
Effective theories: Newton inside Einstein, and why patterns matter
===================================================================

-- Newtonian limit by Taylor expansion
  E = m*(8*c**4*(2*c**2 + v**2) + 6*c**2*v**4 + 5*v**6)/(16*c**4)
  [PASS] E = m c^2 + m v^2/2 + (3/8) m v^4/c^2 + ...
  [PASS] Lorentz transform becomes the Galilean transform as c -> infinity

-- GPS: relativistic rate offset, from the orbit itself
  period check: a = 26561.8 km gives T = 11.967 h
  gravity: +5.2996e-10 = +45.79 us/day; speed: -8.3485e-11 = -7.21 us/day
  net rate offset +4.4647e-10 = +38.58 us/day (Physics Today, Ashby 2002: 4.4645e-10)
  ignoring this: 11.6 km of ranging error per day
  [PASS] computed net offset agrees with the published 4.4645e-10 within 0.2%
  simulated orbit (e = 0.01): mean clock rate offset = +4.4647e-10
  periodic term amplitude: simulated 22.9 ns; 2 sqrt(GM a) e/c^2 = 22.9 ns
  [PASS] simulated mean rate agrees with the closed form to 1%
  [PASS] simulated eccentricity term matches formula to 2%

-- Mercury's perihelion: integrate the relativistic orbit, not the formula
  Newton (no 3GMu^2/c^2 term): shift per orbit = 8.03e-14 rad
  relativistic orbit integrated: 5.0186e-07 rad/orbit = 42.98 arcsec/century
  formula 6 pi G M/(c^2 a (1-e^2)) = 5.0186e-07 rad/orbit = 42.98 arcsec/century (Will 2014: 42.98)
  [PASS] Newtonian control gives no shift
  [PASS] integrated value reproduces 42.98 arcsec/century within 0.1%
  measured total 575.31 arcsec/century (Park et al. 2017): Newtonian perturbations keep 532.33, or 92.5%

-- What the change from Newton to Einstein actually touched
  quantity                                    relativistic correction
  Car at 30 m/s, kinetic energy                              7.51e-15
  Airliner at 250 m/s, kinetic energy                        5.22e-13
  Earth orbit, v^2/c^2                                       9.87e-09
  Sun's surface potential GM/(R c^2)                         2.12e-06
  Earth's surface potential GM/(R c^2)                       6.96e-10
  Hydrogen electron, (alpha)^2                               5.33e-05
  GPS clocks, net offset                                     4.46e-10
  Proton at the LHC, 6.8 TeV, gamma - 1                      7.25e+03
  7/8 of these hand-picked rows change by under 0.1%. This sample shows the mechanism;
  it is not a census of physics and gives no percentage for 'all of physics'.
  [PASS] the same replacement is negligible for slow, weak-field systems and huge at the LHC

-- A law versus a lookup table: predicting data held back
  fit on the 4 inner planets only: exponent = 1.4993 (Kepler's 3/2)
    Jupiter : predicted      4331 d, observed      4331 d, error 0.01%
    Saturn  : predicted     10799 d, observed     10747 d, error 0.48%
    Uranus  : predicted     30577 d, observed     30589 d, error 0.04%
    Neptune : predicted     60409 d, observed     59800 d, error 1.02%
  'every planet has its own law' has no prediction; using the inner mean gives errors 94-100%
  [PASS] a law fit on inner planets predicts the outer ones to under 2%
  [PASS] the fitted exponent is 3/2 within 1%
  control, periods shuffled among planets: median worst error 738%, fraction beating the real law 0.0000
  [PASS] control: shuffled data almost never predicts this well

11/11 checks passed
```

</details>

**Limits:** no percentage of "all physics" exists to compute. The GPS semi-major axis is from memory of Ashby's paper, with the period check shown in the output.

### C36: Physics is not the reality, it's what we make up to understand reality

**Status:** Conceptual

**Short answer:** You're half right, and it's the half many people forget. Equations and concepts are human-made descriptions, like a map. But "made up" can't mean arbitrary, because
reality can reject them. Mercury's 43 arcseconds was not up to us. Physicists and philosophers disagree on how far a good theory captures what's really there.

**Starting from:** a difference between our description (a theory) and what it's a description of.

**Derivation:** evidence for the human-made side: we've replaced concepts while keeping numbers. Gravity was a force, then curvature, and equivalent theories can use very different
pictures. Laudan's historical argument (general knowledge; not confirmed for this write-up) says past successful theories had central terms that we now reject, so perhaps today's do too.

Evidence for the "not arbitrary" side: the [Kepler test](#c37-what-proved-means-and-whether-laws-are-only-local). A rule fit to four planets predicts the outer four within about 1%, and the same procedure
on shuffled data fails with median error over 700%. So the pattern wasn't invented by the fit. Theories also predict things nobody put in. Realists take that as a sign that theories
get something right, Instrumentalists as a sign that they're useful tools. This argument doesn't end it.

**Limits:** this is a philosophical stance, and it's not a physics result. A reasonable working stance is to treat a theory's tested numbers as firm, and its picture as provisional.

### C37: What "proved" means, and whether laws are only local

**Status:** Conceptual + Code + Established

**Short answer:** You're right about fallibilism: no physical law is proved in the way a theorem is, and tested theories have been replaced before. Quantum mechanics and special
relativity aren't proved. They've survived very severe tests. But "maybe there are no laws at all, only special cases" has to explain what we observe, and it doesn't: the
same small rule keeps predicting data it was never fitted to. A future theory will contain today's as a limit, rather than make it false in its tested range.

**Starting from:** a mathematical proof derives from axioms with certainty relative to those axioms. A physical claim is tested by observation and can only be supported, or refuted.

**Derivation:**

*Credit where it's due.* Hume's problem of induction is real: no finite set of observations proves the next one will fit. "Here we know" is exactly the right caution. We already know
that general relativity and quantum theory don't fit together at Planck scales, so at least one will change. Your point stands.

*History needs one correction.* Relativity didn't come from nowhere with everyone laughing. There were measured problems: Mercury's anomaly, known since the 1800s, and the failure to detect
motion through an ether. A new idea wins by matching old data and adding more. "Newton's wrong" alone doesn't do that. Asimov (1989) has the right idea: flat Earth and spherical Earth are both
"wrong" but not equally wrong.

*Why "no laws" isn't the default.* Compare a lookup table with a law. A table storing each case has no predictive power, and needs as many entries as cases. A law is a compression. If the universe had
a special law per case, compression would fail on new data. Test it: fit period against distance for the four inner planets using NASA data. The exponent is 1.4993, and the outer planets come out
within 1.02%. Shuffled data: median worst error 738%, and in 500 tries, never that good. So your hypothesis does make a prediction, and it's one that has failed so far in every domain where we've looked.
It could still be true beyond what we've looked at. That's logically open, but it needs its own evidence. And the "undergrad till your funeral" worry is itself a measure of compression:
the reason you can learn so much in four years is that it compresses so well.

*What a future discovery will do.* If a law fails somewhere, it'll be in a new regime, and it'll keep its accuracy where it's already been tested. See [C35](#c35-if-a-foundation-goes-wrong-does-90-of-physics-collapse).
And the universe laughing? I'd guess it's more that we're slowly getting less wrong.

**Check:** [code/verify_effective_theories.py](code/verify_effective_theories.py), output under [C35](#c35-if-a-foundation-goes-wrong-does-90-of-physics-collapse). This shows compression in one data set, and can't show laws hold everywhere.

**Limits:** induction can't be proved. Eight planets and one shuffle control illustrate a principle, and they don't establish it for all of physics.

### C38: What do blind people see in dreams?

**Status:** Established + Open

**Short answer:** It depends mainly on when sight was lost. People blind from birth or very early generally report dreams made of sound, touch, smell, taste and movement, with emotions and stories
like anyone else's. People who lost sight later often keep visual dreams, which fade with years of blindness. Whether people blind from birth have anything truly visual in their dreams is disputed.

**Starting from:** dream reports, which are the main evidence, along with some EEG. These are small studies.

**Derivation:**

- **Hurovitz, Dunn, Domhoff & Fiss (1999)** studied 15 blind adults and 372 dream reports (a corpus summary lists 381). Those blind from birth or very early had no visual imagery and more taste, smell and touch. The paper notes that blindness before about age 5 seldom leaves visual dream imagery.
- **Meaidi, Jennum, Ptito & Kupers (2014)** collected reports over 4 weeks: 11 congenitally blind, 14 late blind, 25 sighted. All blind participants had fewer visual impressions than sighted controls. In the late blind, duration of blindness correlated negatively with visual duration, clarity and color. The congenitally blind reported more sound, touch, taste and smell, and the late blind more touch only. Emotions and themes didn't differ. The congenitally blind reported more nightmares, perhaps tied to threatening experiences in daily life. Secondary accounts give figures for the congenitally blind of 93% hearing, 67% touch, 40% smell and 26% taste in at least one dream; I couldn't confirm those in the paper, so I'm not relying on them.
- **The disagreement.** Bértolo et al. (2003) studied 10 congenitally blind and 9 sighted people, and reported visual content in dream reports and drawings, and links with EEG alpha power. Kerr (2004) objects that "visual" there may be spatial imagery that doesn't depend on sight, and that visual-cortex activity isn't restricted to visual processing. Lopes da Silva (2003) also commented (only its existence is confirmed here).

So the premise "dreams are visions" is about sighted dreamers. A dream is built from what the brain has. Without visual experience, it uses other senses and spatial understanding. As for "what happens in their mind", the honest limit is that we only get reports, and the meaning of "see" is what's in dispute.

**Limits:** small samples and self-report. Numbers from one study may not generalize. I read abstracts and excerpts, not all full papers.

### C39: How the Pythagorean theorem was found, and rigorous proofs

**Status:** Derived + Code + Established

**Short answer:** It started as practical knowledge of certain triangles and became a rigorous theorem, with real proofs. It's as certain as anything in mathematics, relative to Euclid's axioms, which include the parallel postulate. It fails on curved surfaces, where the parallel postulate fails. Pythagoras himself probably didn't prove it, and earlier cultures knew the rule.

**Starting from:** Euclidean plane geometry, including the parallel postulate. Its consequence that a triangle's angles sum to 180° and that similar triangles exist.

**Derivation:**

*Proof 1, rearrangement.* Take a right triangle with legs $a,b$ and hypotenuse $c$. Place four copies inside a square of side $a+b$, with one on each side, each leg lying along the outer edges. The inner figure has four sides of length $c$. At each of its corners, the two acute angles of the triangles and the inner angle make a straight line. Since the acute angles of a right triangle add to 90° (this is where the angle sum, and so the parallel postulate, gets used), the inner angle is 90°, so the inner figure is a square. Equating areas:

$$(a+b)^2=c^2+4\cdot\tfrac12ab\quad\Rightarrow\quad a^2+2ab+b^2=c^2+2ab\quad\Rightarrow\quad a^2+b^2=c^2.\qquad\blacksquare$$

*Proof 2, similar triangles.* Drop the altitude from the right angle to the hypotenuse, splitting $c$ into $p+q$, with $p$ next to $a$. Each small triangle shares an acute angle with the big one and has a right angle. Since angle sums are fixed, the third angles agree, so each is similar to the big one. Therefore $a/c=p/a$ and $b/c=q/b$, giving $a^2=cp$ and $b^2=cq$. Adding: $a^2+b^2=c(p+q)=c^2$. $\blacksquare$

*Which axiom matters.* The script measures it: on a unit sphere, legs $a=b=1$ give $c=1.2746$ against flat $1.4142$, and hyperbolic space gives $1.5134$. The octant triangle has all sides $\pi/2$, so $c^2$ is half of $a^2+b^2$. Substituting $c^2=a^2+b^2+ka^2b^2/R^2$ into the spherical law $\cos c=\cos a\cos b$ gives $k=-1/3$ exactly, and the measured defect at small size agrees. As $R\to\infty$ the flat rule returns. In a suitable axiom system (Millman and Parker's, according to a source I saw only as an excerpt), the theorem is equivalent to the parallel postulate. The equivalence depends on the axioms: a non-Archimedean counterexample exists, the same excerpt says. So mathematics proves the theorem in Euclidean geometry, and whether real space is Euclidean is a measurement. It's very close: with $R>265$ Gly from [C04](#c04-what-happens-to-a-ball-thrown-at-the-edge-of-the-universe), corrections on any humanly accessible scale are immeasurable.

*History.*

- **Plimpton 322**, a Babylonian tablet from about 1800 BCE, lists numbers that fit the rule. After correcting the four rows with recognized errors, all 15 rows give exact triples, and the rows as inscribed fail, as a control (script). Robson (2001) argues against reading it as a trigonometric table, and connects it to reciprocal pairs and scribal school practice. So it's knowledge of the relation, but it's not a proof.
- **Baudhayana Sulba Sutra**, conventionally about 800 BCE, with uncertain dating. MacTutor gives its statement for the diagonal of a square (doubling the area), with the rectangle version attributed to Katyayana. Other sources read Baudhayana 1.48 as already stating the rectangle case. The Sulbasutras give rules for altars, without proofs.
- **Euclid** I.47, about 300 BCE, is the oldest surviving axiomatic proof.
- **Pythagoras.** The Stanford Encyclopedia entry says "there is not a jot of evidence for a proof by Pythagoras", and the attribution comes from late sources such as Proclus.

So your two options, "rigorous proof or random observation", were both true at different stages.

**Check:** [code/verify_pythagoras.py](code/verify_pythagoras.py) replays the steps with exact rationals and symbolic coordinates. A computed check isn't the proof; the proof is above.

<details><summary>verify_pythagoras.py output</summary>

```text
Pythagorean theorem: proofs, curvature and history
==================================================

-- Proof 1 (rearrangement) replayed with exact rational coordinates
  [PASS] for exact rational legs: inner figure is a square of area a^2 + b^2 = (a+b)^2 - 2ab

-- Proof 2 (similar triangles) replayed symbolically
  [PASS] |BF|^2 from coordinates equals (a^2/c)^2
  [PASS] |AF|^2 from coordinates equals (b^2/c)^2
  [PASS] p + q = c, which is the theorem

-- Measured failure on a sphere and in hyperbolic space
  legs a=b     R  sphere c   hyper c    flat c  sphere c^2/(a^2+b^2)  hyper c^2/(a^2+b^2)
     0.010   1.0   0.01414   0.01414   0.01414              0.999983             1.000017
     0.100   1.0   0.14130   0.14154   0.14142              0.998331             1.001664
     0.500   1.0   0.69172   0.72121   0.70711              0.956948             1.040281
     1.000   1.0   1.27456   1.51337   1.41421              0.812246             1.145150
     1.571   1.0   1.57080   2.52669   2.22144              0.500000             1.293702
  [PASS] sphere c is shorter than flat c; hyperbolic c is longer; both for every size tested
  octant triangle: legs pi/2, hypotenuse 1.57080. c^2 = 2.4674, a^2+b^2 = 4.9348
  [PASS] octant triangle: hypotenuse equals each leg, a flat-rule failure of 2x
  measured coefficient of a^2 b^2/R^2 at a = b = 0.01: 0.3333 (series below predicts 1/3)
  symbolic series from cos c = cos a cos b: c^2 = a^2 + b^2 + (-1/3) a^2 b^2/R^2 + ...
  [PASS] symbolic coefficient is -1/3
  [PASS] the measured defect matches 1/3 within 1%
  [PASS] flat limit: defect vanishes as R -> infinity (c^2/(a^2+b^2) -> 1)

-- Angle sums: the same geometries, tested for the parallel-postulate consequence
  sphere triangle, legs 1, right angle: angle sum = 213.24 deg (> 180)
  [PASS] the same triangle has angle sum above 180 deg, as the parallel postulate's failure requires

-- Plimpton 322: rows checked for s^2 + l^2 = d^2
    short    119, diagonal    169: long side    120  exact
    short   3367, diagonal   4825: long side   3456  exact
    short   4601, diagonal   6649: long side   4800  exact
    short  12709, diagonal  18541: long side  13500  exact
    short     65, diagonal     97: long side     72  exact
    short    319, diagonal    481: long side    360  exact
    short   2291, diagonal   3541: long side   2700  exact
    short    799, diagonal   1249: long side    960  exact
    short    481, diagonal    769: long side    600  exact
    short   4961, diagonal   8161: long side   6480  exact
    short     45, diagonal     75: long side     60  exact
    short   1679, diagonal   2929: long side   2400  exact
    short    161, diagonal    289: long side    240  exact
    short   1771, diagonal   3229: long side   2700  exact
    short     28, diagonal     53: long side     45  exact
  [PASS] after the four standard scribal corrections, all 15 rows are exact Pythagorean triples
    tablet row 2: as inscribed 3367, 11521: difference of squares is perfect square? False
    tablet row 9: as inscribed 541, 769: difference of squares is perfect square? False
    tablet row 13: as inscribed 25921, 289: difference of squares is perfect square? False
    tablet row 15: as inscribed 56, 53: difference of squares is perfect square? False
  [PASS] rows with recognised errors really fail as inscribed (control)

12/12 checks passed
```

</details>

**Limits:** the converse direction of the equivalence (the theorem implies the parallel postulate) is only checked here for constant-curvature geometries. Plimpton row values come from a published listing, not the tablet itself.

### C40: How can thoughts come from atoms and not be atoms?

**Status:** Conceptual + Open + Established + Code

**Short answer:** A thought seems best understood as a process that atoms carry out, not as a thing made of additional stuff. A wave, a song and a running program aren't extra atoms either. So you aren't making something non-physical out of physical things. The part still unexplained is why some of these processes are felt from the inside. That's called the hard problem, and it's unsolved. Nobody has solved it here.

**Starting from:** the observation that thoughts have measurable physical costs and correlates. Herculano-Houzel (2012): 86 billion neurons, 516 kcal/day. Lennie (2003): $2.4\times10^9$ ATP per spike and about 0.16 spikes/s in cortex. Azevedo et al. (2009): 19% of neurons in cortex. A round 50 kJ/mol per ATP (general knowledge).

**Derivation:**

*Thinking has a bill.* The brain runs at 25 W, or about 0.29 nW per neuron, and turns over about $3\times10^{20}$ ATP per second. One spike costs $2\times10^{-10}$ J, which is $7\times10^{10}$ times the Landauer limit at body temperature. Cortical spikes at Lennie's rate need about 0.5 W, which matches his point that activity must be sparse, since firing every neuron at 5 per second at this cost would use 86 W. Attwell & Laughlin (2001) gave an earlier signaling budget, and later revisions changed the share for action potentials, so such budgets carry large uncertainty.

*A process is not an extra substance.* The script simulates 2000 masses joined by springs under Newton's laws. A pulse moves 998 sites while no mass moves more than 0.01 of a spacing. The number of masses doesn't change. Is the pulse "made of atoms"? It's something atoms are doing. That's the sense in which your premise "not made of atoms" is correct, and it isn't mysterious: "running" isn't made of legs either. And "we can't touch it" applies to any pattern. You could also say thoughts *can* be measured, through energy use and brain imaging, which suggests that they aren't non-physical in the sense of having no physical footprint.

*What's left.* Chalmers (1995) separates easy problems (how the brain discriminates, integrates, reports: hard work, but we know how to go about it) from the hard problem (why this is accompanied by experience at all). The pulse says nothing about that. Physicalists expect an explanation in physical terms, or say the question will dissolve. Others think something's missing. You've got a real question, and the strongest honest statement is that the process view handles the "made of" puzzle and doesn't handle the "what it's like" one.

**Check:** [code/verify_brain_numbers.py](code/verify_brain_numbers.py) checks arithmetic and simulates a pulse. It's no model of a neuron or a thought.

<details><summary>verify_brain_numbers.py output</summary>

```text
Brain numbers and signal-as-process
===================================

-- Energy budget
  brain power = 516.0 kcal/day = 25.0 W; per neuron = 0.29 nW
  [PASS] brain power is in the familiar 20-25 W range
  ATP turnover at 50 kJ/mol = 3.0e+20 ATP/s; mass of ATP cycled per day = 22 kg
  one spike at 2.4e+09 ATP = 2.0e-10 J
  cortex neurons (19% of 86 billion) = 1.6e+10
  cortical spikes/s at 0.16/s: 2.6e+09; power for those spikes = 0.52 W
  [PASS] spikes alone use a modest part of the total, consistent with budgets where other costs matter
    if all 86e9 neurons fired 0.16/s at that cost: 2.7 W
    if all 86e9 neurons fired 1/s at that cost: 17.1 W
    if all 86e9 neurons fired 5/s at that cost: 85.7 W
  Landauer minimum at 310 K = 2.97e-21 J/bit; a spike costs 6.7e+10 times that
  [PASS] a real spike is far above the thermodynamic floor
  atoms in a 1.4 kg water-like brain = 1.4e+26; a thought adds none of its own

-- A pulse travels; the atoms do not: mass-spring chain, Newton's laws only
  pulse moved 998 lattice sites in time 1000; speed = 0.998 (continuum prediction 1)
  largest displacement of any mass in the whole run = 0.0100 lattice spacings
  ratio of distance signal went to distance any mass went = 1.0e+05
  [PASS] pulse speed matches the wave speed to 1%
  [PASS] signal crosses about 1000 sites while no mass moves one site
  [PASS] number of masses is unchanged: the signal is a pattern, not a new object

6/6 checks passed
```

</details>

**Limits:** the numbers are estimates with wide ranges, and rest on modeling. Correlation between brain states and thoughts doesn't explain experience.

### C41: A glass of water from the ocean, and "there's nobody else"

**Status:** Derived + Code + Established + Open

**Short answer:** Your logic is right, and it's a well-known analogy. Wright, Kanodia & Lubar (2018) credit Tarter et al. (2010) with the comparison of a drinking glass of seawater, and update it: about 8000 litres, a large hot tub or small pool. Zero detections in such a sample barely changes the odds. And nobody serious can claim that there aren't any. But zero detections do constrain some things: a civilization that's very easy to see, like a fish filling the sea.

**Starting from:** Wright et al.'s search fraction $6.0\times10^{-18}$ of an 8-dimensional haystack, and ocean volume $1.335\times10^{21}$ L. A Poisson detection model.

**Derivation:**

*Your glass.* A 6 oz glass is 0.177 L, so the ratio is $1.3\times10^{-22}$. Wright's number is 45,000 of those glasses, which is 8010 L by the script, and still tiny.

*Bayes.* Let $N$ be transmitters in the haystack, and $f$ the searched fraction. $P(\text{zero}\mid N)=e^{-Nf}$. With a flat prior, the 95% limit is $N<3/f=5\times10^{17}$. Compare with about $10^{11}$ stars in our galaxy: that's no useful bound at all.

With a log-uniform prior on $N$ from 1 to $N_\text{max}$:

| Searched fraction $f$ | Bayes factor for zero detections | $P(\text{none})$ from a 50% prior ($N_\text{max}=10^{22}$) |
|---|---|---|
| $6\times10^{-18}$ | 0.771 | 0.565 |
| $10^{-6}$ | 0.261 | 0.793 |
| $10^{-2}$ | 0.080 | 0.926 |
| $0.5$ | 0.011 | 0.989 |

That first row depends on allowing $N$ far beyond the star count. Cap $N$ at $10^{11}$, and the shift is below 0.01%. So absence of evidence is evidence of absence, strictly speaking, with a strength set by how much you've searched. Right now it's near zero.

*Being fair to the other side.* Two things cut against "the ocean must be full". First, Sandberg, Drexler & Ord (2018) show that once you use real uncertainty in the origin-of-life steps, there's a substantial prior probability that we're alone in the observable universe, so silence isn't that surprising. Second, some searches are strong for their target: G-hat looked at about $10^5$ galaxies and found none with over 85% of starlight reprocessed (Griffith et al. 2015). That implies below $3\times10^{-5}$ of that type per galaxy, if detection is reliable. So a big fish is limited even though small ones aren't.

So: "we've found none" is true, "there are none" isn't known, and "there must be many" isn't known either.

**Check:** [code/verify_seti_bayes.py](code/verify_seti_bayes.py) simulates fish in a sample and computes posteriors numerically.

<details><summary>verify_seti_bayes.py output</summary>

```text
SETI: search fraction and zero-detection Bayes
==============================================

-- Volumes
  ocean = 1.335e+21 L, 6 oz glass = 0.177 L, ratio = 1.3e-22
  Wright et al. fraction 6.0e-18 of ocean = 8010 L (paper says about 8000 L, hot tub to small pool)
  [PASS] fraction x ocean reproduces the paper's roughly 8000 L
  on Wright's numbers the search so far is 45142 of your glasses, a bigger sample than your glass, still tiny

-- Ocean analogy with fish actually simulated
  fish N, chance that a sample of fraction f = 1e-3 catches none; analytic exp(-N f); simulated
    N =    10: analytic 0.9900; simulated 0.9902
    N =   100: analytic 0.9048; simulated 0.9051
    N =  1000: analytic 0.3677; simulated 0.3667
    N =  5000: analytic 0.0067; simulated 0.0069
  [PASS] simulated empty-sample rate agrees with (1-f)^N

-- Bayesian update on zero detections
   search fraction f   N median  N 95% upper   Bayes factor (0 detected) vs prior
             1.3e-22    6.1e+10      3.3e+20                                0.981
             6.0e-18    3.1e+08      1.4e+16                                0.771
             1.0e-06    7.5e+02      4.3e+05                                0.261
             1.0e-02    8.1e+00      1.0e+02                                0.080
             5.0e-01    1.7e+00      4.8e+00                                0.011
  [PASS] at the Wright fraction zero detections barely shift a wide prior
  flat prior: 95% upper limit N < 3/f = 5.0e+17 transmitters of the haystack's kind; integral check 0.950
  [PASS] flat-prior 95% credible region integrates to 0.95
  compare stars in the Milky Way (order 1e11): upper limit is 5e+06 per star. Not a useful bound
  [PASS] limit far exceeds the number of stars, so no claim of absence follows

-- Hypothesis 'none exist' versus 'exist', both fairly
  f = 6.0e-18: P(none | zero detected) = 0.565 from a 0.500 prior
  f = 1.0e-06: P(none | zero detected) = 0.793 from a 0.500 prior
  f = 1.0e-02: P(none | zero detected) = 0.926 from a 0.500 prior
  f = 5.0e-01: P(none | zero detected) = 0.989 from a 0.500 prior
  prior capped at one transmitter per star (N <= 1e11): P(none | zero detected) = 0.500000
  so the small rise above depends on allowing N far above the number of stars
  [PASS] with N capped near the star count, the Wright-sized search moves the odds by under 0.01%
  [PASS] posterior for 'none' rises only slightly at the real search fraction
  [PASS] absence IS evidence, just weak: posterior for 'none' rises at any f > 0

-- Where zero detections do bite: civilisations that are easy to see
  G-hat: about 1e+05 galaxies, none reprocessing > 85% of starlight (Griffith et al. 2015)
  if detection is reliable: 95% limit on that kind of civilisation = 3.0e-05 per galaxy
  [PASS] a big fish that fills its ocean is constrained even by a modest search

9/9 checks passed
```

</details>

**Limits:** Wright et al. say the numbers are very sensitive to the haystack definition (volume goes as $d_\text{max}^5$), which covers particular radio signals only. It doesn't cover life without technology. Targeted surveys are more complete for nearby stars.

### C42: 98% DNA, and a civilization 2% above us

**Status:** Established + Code + Conceptual

**Short answer:** Your sense of how large the gap is, between chimps and us, is fair. The arithmetic behind "2% more" isn't: the percent is a count of letters, and it doesn't measure ability. One letter can destroy a protein, and 2% of letters can change nothing. Chimps aren't only "smashing boxes" either. In some tests they match human toddlers or beat adults. Nobody can rule out what you're imagining at the end, though there's no evidence for it ([C45](#c45-could-such-manipulation-already-have-happened)).

**Starting from:** Chimpanzee Sequencing and Analysis Consortium (2005), plus the standard genetic code.

**Derivation:**

*The number.* Mean substitution rate 1.23% between copies of the two genomes, with 1.06% or less fixed between species. About 35 million substitutions. Indels: about 90 Mb in total in about 5 million events, or about 3% more. So by events, substitutions lead 7 to 1. By bases, indels lead 3 to 1. About 29% of orthologous proteins are identical, with the typical one differing by two amino acids. Complete genomes (Yoo et al. 2025) show that 12.5-27.3% of an ape genome fails to align or isn't consistent with simple one-to-one alignment, covering all apes. So the "difference" is around 1.2%, 4% or above 12%, depending on the definition. These aren't conflicting. They count different things. Some people misuse the larger number to say that apes aren't closely related; gap measures count structural changes where single events can span many thousands of bases.

*Percent doesn't equal function.* In a random 500-codon gene:

| Change | Bases altered | Protein identical |
|---|---|---|
| Synonymous only | 2.0% | 100% |
| Random positions | 2.0% | 61% |
| 1-base frameshift | 0.07% | 2.2% |
| 1-base stop | 0.07% | 0.4% |

King & Wilson (1975) argued long ago that much of the difference may lie in regulation. And $86/28$: humans have about 3.1 times the predicted neuron count of chimps (Herculano-Houzel 2012). So a tiny sequence change can be tied to big differences in development.

*Intelligence has no 2% scale.* IQ is an interval scale built by norming, so adding a constant changes any ratio (script). There's no meaningful sense of "2% smarter", and comparing Einstein to a chimp doesn't reduce to a number.

*What chimps can do.* Herrmann et al. (2007): human 2-year-olds and apes scored about the same on physical-world tasks, with toddlers far ahead on social tasks. Inoue & Matsuzawa (2007): young chimps beat human adults on one rapid numeral memory task. Krupenye et al. (2016): apes look ahead to where an agent falsely believes something to be, implicit false-belief understanding, though interpretations are debated. Tennie, Call & Tomasello (2009): cumulative culture, where each generation builds on the last, seems to be the key human distinction. This is what your "exploring space" is really about. No individual human could do it. Newton didn't invent algebra. So the gap between us and them may be largely about accumulation, and that changes the analogy: we might accumulate past a threshold that doesn't have a next level in the same sense. That part's speculation.

**Check:** [code/verify_dna_bayes.py](code/verify_dna_bayes.py) uses the real genetic code on a random sequence. It's illustration, not a model of human evolution.

<details><summary>verify_dna_bayes.py output</summary>

```text
Sequence percent and undetectable manipulation
==============================================

-- Genome arithmetic (Chimpanzee Sequencing and Analysis Consortium 2005)
  aligned sequence 2.4e+09 bp, substitutions at 1.23%: 3.0e+07 (paper: about 35 million)
  [PASS] 1.23% of the aligned sequence reproduces about 30-35 million
  indel sequence total about 90 Mb in about 5 million events = 18 bp per event on average
  paper's total: 1.23% by substitutions, about 3% more by indels. Events: 7 times more substitutions, bases: 3.0 times more in indels
  [PASS] event count and base count rank the two kinds of change oppositely: the percent depends on the unit
  Yoo et al. 2025 complete genomes: 12.5-27.3% of an ape genome fails simple one-to-one alignment (all apes)
  Same pair, three legitimate numbers: about 1.2%, about 4%, above 12%. Which is 'the' difference? None is.

-- Same percent, opposite function: real genetic code
  synonymous only, 2% of bases        : base differences  30/1503, protein identity 100.00%
  random, 2% of bases                 : base differences  30/1503, protein identity  61.48%
  1-base frameshift, 0.07% of bases   : base differences   1/1503, protein identity   2.20%
  1-base stop gain, 0.07% of bases    : base differences   1/1503, protein identity   0.40%
  [PASS] 2% synonymous change leaves the protein 100% identical
  [PASS] one base (0.07%) can destroy over 90% of a protein

-- Brains are not 2% apart either
  neurons: human 86 billion, chimpanzee predicted 28 billion (Herculano-Houzel 2012): ratio 3.1
  So 1.23% sequence divergence can coexist with a threefold neuron count difference.
  [PASS] the neuron ratio is not 1.02
  Scale test: IQ is an interval scale built by norming, not a ratio. 2% above something depends on the zero:
    add    0 to every score: 100 -> 100, 130 -> 130, ratio 1.300
    add  100 to every score: 100 -> 200, 130 -> 230, ratio 1.150
    add  -50 to every score: 100 -> 50, 130 -> 80, ratio 1.600
  [PASS] ratios on a shiftable scale are not meaningful

-- Bayes: hypotheses and their predictions
  H_mimic: prior 0.0100; after 10000 observations posterior 0.0100
  [PASS] the mimicking hypothesis can't be moved in either direction
  H_any: spreads over all outcomes; log10 Bayes factor after 10000 observations = -1694.8
  [PASS] the vague 'could do anything' version loses to a sharper model
  H_trace: manipulators leave an anomaly with chance q = 0.001 per independent observation
        100 clean observations: Bayes factor 9.05e-01; posterior from 0.01 prior 9.06e-03
       1000 clean observations: Bayes factor 3.68e-01; posterior from 0.01 prior 3.70e-03
      10000 clean observations: Bayes factor 4.52e-05; posterior from 0.01 prior 4.56e-07
     100000 clean observations: Bayes factor 3.54e-44; posterior from 0.01 prior 3.57e-46
  [PASS] any version with detectable traces gets testable: clean data lowers it

9/9 checks passed
```

</details>

**Limits:** the chimp neuron number is a prediction, not a count. The cognition studies have specific samples and tasks.

### C43: Could a species only 2% more intelligent dominate us without us knowing?

**Status:** Conceptual + Speculative

**Short answer:** "2% more intelligent" has no defined meaning, so the question needs restating. Could a far more capable species secretly influence us? It can't be ruled out, and there's no evidence for it. Whether it could stay undetected depends on how much it does and how clean its traces are.

**Starting from:** [C42](#c42-98-dna-and-a-civilization-2-above-us), for why percent doesn't work.

**Derivation:** a capability gap could matter a lot. History has examples of groups being dominated by others with more technology, and those being dominated generally did know they were. Secret influence is harder than open power: acting on a whole civilization means many interventions, and each risks a trace. Whether that's small depends on capability, and that's unknown. See [C46](#c46-would-we-even-be-capable-of-recognizing-it) for how detection scales.

**Limits:** no data. This is a discussion of what's conceivable.

### C44: Are humans to such beings what chimps are to humans?

**Status:** Conceptual + Speculative + Open

**Short answer:** Possibly in capability, but the analogy has a weak point. Chimps have no science. We do, and that changes what can stay hidden. Whether there are ideas beyond human reach is an open question, and I can't answer it.

**Starting from:** the cumulative culture point in [C42](#c42-98-dna-and-a-civilization-2-above-us).

**Derivation:** where the analogy works: a chimp doesn't understand why it's in a research center. So understanding can have limits for a given mind, and we can't assume ours is the exception. Where it fails: a chimp does notice a cage, and it can't build instruments, write down data, or pass on a question for 400 years. Humans outsource thinking to records, math and machines, and that may extend what we can reach. Whether this hits a limit is unknown. Also, a chimp does know many things in its world, and the "knows nothing" picture is too low.

**Limits:** no way to test this without meeting such beings.

### C45: Could such manipulation already have happened?

**Status:** Conceptual + Code

**Short answer:** It can't be ruled out, and there's no evidence that it happened. Don't call it impossible. But the version that predicts every observation can't gain support from any observation, which means there's no reason to believe it, and no way to refute it. Credit to you: you saw that "we'd notice" isn't guaranteed. That's a good epistemic point.

**Starting from:** Bayes' theorem, posterior odds = prior odds × likelihood ratio.

**Derivation:** split the claim into three, since they behave differently.

| Version | What it predicts | What happens with data |
|---|---|---|
| **Mimic**: they produce exactly what an unmanipulated world would | $P(D\mid H)=P(D\mid H_0)$ for every $D$ | Posterior stays at prior forever. From 0.01, after 10,000 observations: 0.0100. |
| **Anything goes**: they could cause any outcome | Likelihood spread thin over outcomes | Loses to any sharper model. $\log_{10}$ Bayes factor $=-1695$ over 10,000 observations in the toy. |
| **Leaves traces** with chance $q$ per independent observation | Some anomalies | Testable. For $q=0.001$, 1000 clean observations give a factor 0.37, and $10^5$ give $3.5\times10^{-44}$. |

The Bayes argument for the first is one line. If likelihoods are equal for all $D$, the ratio is 1, so the posterior equals the prior. The prior is all you have, and that's no evidence. It works just the same for rival stories that the same mimic could tell, such as a different species or different motives, and those conflict with each other, so none gets a boost.

This has a known cousin, Ball's (1973) zoo hypothesis, that others may be leaving us alone, on purpose, as in a reserve.

Saying "unfalsifiable" isn't saying "false". It says no evidence can help. Best stance: can't be ruled out, no evidence, and not a basis for action.

**Check:** [code/verify_dna_bayes.py](code/verify_dna_bayes.py), output under [C42](#c42-98-dna-and-a-civilization-2-above-us). It verifies the structure of the update. It can't assign a real prior.

**Limits:** Occam-type reasons to give such hypotheses low priors are judgement calls. The toy distribution is arbitrary, and only the general behaviour matters.

### C46: Would we even be capable of recognizing it?

**Status:** Conceptual + Code + Open

**Short answer:** It depends on the manipulation. Anything that leaves measurable anomalies could be found, given enough independent checks, even if we didn't understand why. A perfect mimic couldn't be. And a mind can be limited in ways it can't see, so we can't be sure we'd know what to look for.

**Starting from:** the third row of the table in [C45](#c45-could-such-manipulation-already-have-happened), with $q$ the chance of a trace per independent observation.

**Derivation:** after $n$ clean observations the Bayes factor against a trace-leaving version is $(1-q)^n$. So to reach a factor of 0.05 you need about $3/q$ independent observations. Science has many such measurements: conservation laws, clocks, genetics, cross-checks between labs. So even subtle interference has to be very careful. That's a practical point in favour of recognizing it, unlike a chimp. A chimp recognizes a cage, and it can't tell why. We could likewise find an anomaly and not know what it means.

The catches:

1. Data can be selected or fed: a clever manipulator could make the "independent" observations less independent, so the count $n$ would mean less.
2. Recognizing a pattern depends on the concepts you already have.
3. False positives. Most anomalies have ordinary explanations, and a mind that wants to see patterns finds them. Any serious test has to be stated in advance.

**Check:** [code/verify_dna_bayes.py](code/verify_dna_bayes.py), output under [C42](#c42-98-dna-and-a-civilization-2-above-us). It tests only the arithmetic of the updating.

**Limits:** independent observations and $q$ can't be estimated for an unknown agent. So this gives the form of the answer, not a number.

---

## References

"Confirmed" means I saw the source's record, abstract or text while writing. Details in brackets say how far that went.

**Cosmology**

- Planck Collaboration (2020). Planck 2018 results VI. A&A 641, A6. [arXiv:1807.06209](https://arxiv.org/abs/1807.06209). Abstract read for $H_0$, $\Omega_m$, $\Omega_K$; age from [the A&A table](https://www.aanda.org/articles/aa/full_html/2020/09/aa33910-18/T1.html) via a search excerpt.
- Planck Collaboration (2016). Planck 2015 results XVIII. A&A 594, A18. [arXiv:1502.01593](https://arxiv.org/abs/1502.01593).
- Planck Collaboration (2014). Planck 2013 results XXV. [arXiv:1303.5085](https://arxiv.org/abs/1303.5085).
- Riess et al. (2022). ApJL. [arXiv:2112.04510](https://arxiv.org/abs/2112.04510).
- Davis & Lineweaver (2004). Expanding Confusion. PASA. [arXiv:astro-ph/0310808](https://arxiv.org/abs/astro-ph/0310808).
- Fields, Olive, Yeh & Young (2020). JCAP 03, 010. [arXiv:1912.01132](https://arxiv.org/abs/1912.01132). Abstract only.
- Fixsen (2009). ApJ 707, 916. [arXiv:0911.1955](https://arxiv.org/abs/0911.1955).
- Hubble (1929). PNAS 15, 168. [Reprint](https://apod.nasa.gov/diamond_jubilee/1996/hub_1929.html). On Slipher and Lemaître, through [a historical review](https://pmc.ncbi.nlm.nih.gov/articles/PMC314128).
- van der Marel et al. (2012). ApJ 753, 8. [arXiv:1205.6864](https://arxiv.org/abs/1205.6864).
- Krauss & Scherrer (2007). GRG 39, 1545. [arXiv:0704.0221](https://arxiv.org/abs/0704.0221).
- Naidu et al. (2025). [arXiv:2505.11263](https://arxiv.org/abs/2505.11263).
- Saadeh et al. (2016). PRL 117, 131302. [arXiv:1605.07178](https://arxiv.org/abs/1605.07178).
- Su & Chu (2009). ApJ 703, 354. [arXiv:0902.4575](https://arxiv.org/abs/0902.4575).
- Borde, Guth & Vilenkin (2003). PRL 90, 151301. [arXiv:gr-qc/0110012](https://arxiv.org/abs/gr-qc/0110012).
- Hartle & Hawking (1983). PRD 28, 2960. [INSPIRE](https://inspirehep.net/literature/192909). Citation only.
- Archytas: [Javelin argument](https://en.wikipedia.org/wiki/Javelin_argument) and Schofield & Huffman's [chapter](https://cambridge.org/core/books/history-of-pythagoreanism/archytas/07B9D035C82DBA6C37CEC45B8B09C03F). Search excerpts, not the Greek source.
- Redfield et al. (2019) on local density. [arXiv:1910.01243](https://arxiv.org/abs/1910.01243). Excerpt.

**Relativity and time machines**

- Bailey et al. (1977). Nature 268, 301. [Record](https://www.nature.com/articles/268301a0.pdf). Excerpt.
- Chou, Hume, Rosenband & Wineland (2010). Science 329, 1630. [Record](https://vivo.colorado.edu/display/pubid_113209).
- Semyonov. Radiation hazard of relativistic interstellar flight. [arXiv:physics/0610030](https://arxiv.org/abs/physics/0610030).
- Gödel (1949). RMP 21, 447. Citation confirmed through a reference list.
- Natário (2012). GRG 44, 855. [arXiv:1105.6197](https://arxiv.org/abs/1105.6197).
- Tipler (1974). PRD 9, 2203. [Record](https://scilit.com/publications/e36364e9afb2cb00927787d5ff9c4f6b). Abstract excerpt.
- Carter (1968). Phys. Rev. 174, 1559. [PDF](https://swx.sinp.msu.ru/panov/Lib/Papers/BHI/PhysRev.174.1559.pdf). Excerpt.
- Poisson & Israel (1989). PRL 63, 1663. [PubMed](https://pubmed.ncbi.nlm.nih.gov/10040638/). Citation only.
- Dafermos & Luk. [arXiv:1710.01722](https://arxiv.org/abs/1710.01722).
- Geroch (1967). J. Math. Phys. 8, 782. Citation only.
- Gott (1991). PRL 66, 1126. [NASA record](https://ntrs.nasa.gov/citations/19910042584). Volume and page, and the $\gamma\sin$ form, are from general knowledge. The threshold is reproduced by code.
- Carroll, Farhi & Guth (1992). PRL 68, 263. [APS](https://link.aps.org/doi/10.1103/PhysRevLett.68.263). Abstract excerpt.
- Deser, Jackiw & 't Hooft (1992). PRL 68, 267. Volume and page confirmed through secondary sources; title not confirmed.
- Morris, Thorne & Yurtsever (1988). PRL 61, 1446. [Caltech](https://authors.library.caltech.edu/9262/).
- Hawking (1992). PRD 46, 603. [APS](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.46.603). Abstract; the sketch in C17 is my summary.
- Kay, Radzikowski & Wald (1997). CMP 183, 533. [arXiv:gr-qc/9603012](https://arxiv.org/abs/gr-qc/9603012).
- Visser (2002). [arXiv:gr-qc/0204022](https://arxiv.org/abs/gr-qc/0204022).
- Visser, Kar & Dadhich (2003). PRL 90, 201102. [arXiv:gr-qc/0301003](https://arxiv.org/abs/gr-qc/0301003).
- Ford & Roman (1996). PRD 53, 5496. [arXiv:gr-qc/9510071](https://arxiv.org/abs/gr-qc/9510071).
- Ori (2005). PRL 95, 021101. [arXiv:gr-qc/0503077](https://arxiv.org/abs/gr-qc/0503077).
- Friedman et al. (1990). PRD 42, 1915. Citation only.
- Echeverria, Klinkhammer & Thorne: billiard balls in wormhole spacetimes. [Archive record](https://public.ucrlib.aspace.cdlib.org/repositories/3/archival_objects/15583). Manuscript record only.
- Deutsch (1991). PRD 44, 3197. [PDF](https://gwern.net/doc/science/physics/1991-deutsch.pdf). Citation.
- Lloyd et al. (2011). PRL 106, 040403. [arXiv:1005.2219](https://arxiv.org/abs/1005.2219).

**Block universe**

- Savitt. [Being and Becoming in Modern Physics](https://plato.stanford.edu/entries/spacetime-bebecome/). Stanford Encyclopedia of Philosophy.
- Rietdijk (1966). Phil. Sci. 33, 341. Putnam (1967). J. Phil. 64, 240. Stein (1968). J. Phil. 65, 5. Citations confirmed.
- Penrose (1989). The Emperor's New Mind. Through Savitt's bibliography.
- Balashov & Janssen (2003). BJPS 54, 327. Dorato (2006). [PhilSci](https://philsci-archive.pitt.edu/2214). Abstracts.
- Weyl (1949). Philosophy of Mathematics and Natural Science. The quote is confirmed only in secondary sources, such as [this chapter](https://springerlink.bibliotecabuap.elogim.com/content/pdf/10.1007/978-3-540-68001-7_8.pdf).

**Boltzmann brains**

- Carroll (2017). [arXiv:1702.00850](https://arxiv.org/abs/1702.00850).
- Albrecht & Sorbo (2004). PRD 70, 063528. [arXiv:hep-th/0405270](https://arxiv.org/abs/hep-th/0405270).
- Dyson, Kleban & Susskind (2002). JHEP 0210, 011. [arXiv:hep-th/0208013](https://arxiv.org/abs/hep-th/0208013).
- De Simone et al. (2010). PRD 82, 063520. [arXiv:0808.3778](https://arxiv.org/abs/0808.3778).
- Boddy, Carroll & Pollack. [arXiv:1405.0298](https://arxiv.org/abs/1405.0298).

**Foundations**

- Ashby (2002). [Relativity and the Global Positioning System](https://physicstoday.aip.org/features/relativity-and-the-global-positioning-system). Physics Today. Ashby (2003). Living Rev. Relativity 6, 1. [Record](https://www.maths.tcd.ie/EMIS/journals/LRG/Articles/lrr-2003-1/). The full text wasn't readable, so the 4.4645e-10 comes from the Physics Today article.
- Will (2014). [arXiv:1403.7377](https://arxiv.org/abs/1403.7377). 42.98 arcsec/century.
- Park et al. (2017). AJ 153, 121. [MIT](https://dspace.mit.edu/handle/1721.1/109312). Abstract through a search excerpt.
- NASA [Planetary Fact Sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/).
- Asimov (1989). The Relativity of Wrong. Skeptical Inquirer 14(1), 35. [Text](https://web.williams.edu/Mathematics/sjmiller/public_html/238/handouts/asimovrelativityofwrong.htm).

**Dreams, brain, mathematics**

- Hurovitz, Dunn, Domhoff & Fiss (1999). Dreaming 9, 183. [Library copy](https://dreams.ucsc.edu/Library/hurovitz_1999a.html). Search excerpts only.
- Meaidi, Jennum, Ptito & Kupers (2014). Sleep Med. 15, 586. [Abstract](https://research.regionh.dk/en/publications/the-sensory-construction-of-dreams-and-nightmare-frequency-in-con/).
- Bértolo et al. (2003). Cogn. Brain Res. 15, 277. Search excerpts. Lopes da Silva (2003). Trends Cogn. Sci. 7, 328. Existence only. Kerr (2004). [Text excerpt](https://dreams.ucsc.edu/Library/kerr_2004.html).
- Robson (2001). Historia Math. 28, 167. [Record](https://keibi.uni-tuebingen.de/Record/KEI00000991). Abstract.
- [Plimpton 322](https://en.wikipedia.org/wiki/Plimpton_322) for the row values.
- MacTutor. [The Indian Sulbasutras](https://mathshistory.st-andrews.ac.uk/HistTopics/Indian_sulbasutras/).
- Stanford Encyclopedia. [Pythagoras](https://plato.stanford.edu/entries/pythagoras/).
- Parallel-postulate equivalence: [Cut-the-Knot](https://www.cut-the-knot.org/triangle/pythpar/PTimpliesPP.shtml), and a Journal of Automated Reasoning paper, DOI 10.1007/s10817-017-9422-8, seen only as an excerpt (authors not confirmed).
- Herculano-Houzel (2012). PNAS 109 (Suppl 1), 10661. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC3386878).
- Azevedo et al. (2009). J. Comp. Neurol. 513, 532. [PubMed](https://pubmed.ncbi.nlm.nih.gov/19226510/).
- Lennie (2003). Curr. Biol. 13, 493. [PDF](https://www2.bcs.rochester.edu/sites/plennie/pdfs/Lennie03a.pdf). Excerpt.
- Attwell & Laughlin (2001). JCBFM 21, 1133. [PubMed](https://pubmed.ncbi.nlm.nih.gov/11598490/). Abstract excerpt.
- Chalmers (1995). JCS 2(3), 200.

**SETI, genomes, cognition**

- Wright, Kanodia & Lubar (2018). AJ 156, 260. [arXiv:1809.07252](https://arxiv.org/abs/1809.07252). Paper text read for the numbers.
- Griffith et al. (2015). ApJS 217, 25. [arXiv:1504.03418](https://arxiv.org/abs/1504.03418).
- Sandberg, Drexler & Ord (2018). [arXiv:1806.02404](https://arxiv.org/abs/1806.02404).
- Chimpanzee Sequencing and Analysis Consortium (2005). Nature 437, 69. [PDF](https://www.genome.gov/Pages/Research/DIR/Chimp_Analysis.pdf). Text read.
- Yoo et al. (2025). Nature. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12058530/). Main text excerpt for 12.5-27.3%; that human-chimp gap figures of about 12 to 13% are from secondary discussion.
- King & Wilson (1975). Science 188, 107. Idea confirmed through a textbook; volume and page from general knowledge.
- Herrmann et al. (2007). Science 317, 1360. Search excerpts; volume and page from general knowledge.
- Inoue & Matsuzawa (2007). Curr. Biol. 17, R1004.
- Krupenye et al. (2016). Science 354, 110.
- Tennie, Call & Tomasello (2009). Phil. Trans. R. Soc. B 364, 2405. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC2865079).
- Ball (1973). Icarus 19, 347. [ADS](https://ui.adsabs.harvard.edu/abs/1973Icar...19..347B/abstract).

**Not confirmed.** The search budget ran out before I could check these, so they stand as general knowledge: Penzias & Wilson (1965, ApJ 142, 419); Minkowski (1908, "Raum und Zeit"); Laudan (1981); Feynman Lectures on local conservation; Boltzmann's fluctuation idea; the Novikov credit beyond his co-authorship of Friedman et al.; the 4 Gy lethal dose, 2 MeV cm²/g stopping power, and water's molar entropy.
