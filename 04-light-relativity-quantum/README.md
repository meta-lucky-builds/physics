# Light, relativity, quantum mechanics and measurement

Sixteen questions, L01 to L16. They share a few threads, and the answers point to each
other where they meet:

- A box of light shows that mass is not additive (L01). The same box reappears when
  Einstein and Bohr argued about uncertainty (L10).
- Only comparisons can be measured. That is why the one-way speed of light is a convention
  (L02) and why a world with one kind of clock can't see its own drift (L14).
- "The fastest thing" is more careful than it sounds: crests, shadows and images can beat
  $c$, while energy and information can't (L12).
- "Classical" fails in surprising places: thermal expansion itself freezes out for
  quantum reasons (L15).

Scripts are in [`code/`](code/), and each prints its checks and exits non-zero on a failure.
Run them from the repo root, for example `python 04-light-relativity-quantum/code/verify_eclipse.py`.
The output blocks below are pasted from real runs and cut to the relevant part.
Facts I could not confirm online are marked as general knowledge and are listed at the end.

---

### L01: A million flashlights in a box

**Status:** Derived + Code + Estimate

**Short answer:** Your hand would burn within milliseconds, while the mass difference would be far too small for any scale. Light in a box really does add mass $m = E/c^2$ even though each photon has zero mass, but a megawatt of light stored even between record mirrors weighs about $10^{-14}$ kg. If the flashlights and batteries are sealed inside the box, its total mass doesn't change at all when they switch on.

**Starting from:** Special relativity's energy-momentum relation $E^2 = (pc)^2 + (mc^2)^2$, with $E = pc$ for a photon. Conservation of energy and momentum. The inputs below are assumptions I chose: 1 W of light and 3 W of electricity per flashlight, 10 Wh battery, 100 cm³ and 150 g per flashlight, a hand of area 0.04 m² absorbing 60% of visible light.

**Derivation:**

*Why zero-mass photons can make mass.* Invariant mass comes from the *totals* of a system, not from adding up each part's mass:

$$
M^2c^4 = \Big(\sum E_i\Big)^2 - \Big|\sum \vec p_i\Big|^2c^2 .
$$

For two photons with energies $E_1, E_2$ at angle $\theta$, using $|\vec p_i| = E_i/c$:

$$
M^2c^4 = E_1^2 + E_2^2 + 2E_1E_2 - (E_1^2 + E_2^2 + 2E_1E_2\cos\theta) = 2E_1E_2(1-\cos\theta).
$$

Parallel photons ($\theta = 0$) have $M = 0$, however bright the beam. Opposite photons ($\theta = \pi$) with equal energy have $M = 2E/c^2$. In a box, photons head in all directions, so the total momentum averages to zero and $M = \sum E_i / c^2$.

*Why a scale would notice.* Suppose one photon bounces vertically between floor and ceiling a height $h$ apart. Each bounce kicks the mirror by $2E/c$. Climbing redshifts the photon, $E_{\rm top} = E(1 - gh/c^2)$, so the floor's downward kick beats the ceiling's upward kick. The net downward impulse per round trip is $2Egh/c^3$, and a round trip takes $2h/c$:

$$
F = \frac{2Egh/c^3}{2h/c} = \frac{E}{c^2}\,g .
$$

That is exactly the weight of a mass $E/c^2$. The script tests the equivalent version, a photon in an accelerated box with Doppler shifts at each mirror, without using $m = E/c^2$ anywhere.

*How long the light lasts.* In a box of volume $V$ and wall area $S$, the average straight path between hits is $4V/S$ (a standard result, tested below by ray tracing). With reflectivity $R$:

$$
\tau = \frac{4V/S}{c(1-R)}, \qquad U = P\tau, \qquad \Delta m = \frac{P\tau}{c^2}.
$$

Add any other absorber of effective area $A_{\rm abs}$ and $1-R$ is replaced by $(A_{\rm abs} + (1-R)S)/S$.

*The three cases.*

1. **Light fed into an empty 1 m mirrored cube at 1 MW.** With $R = 0.99$ the light survives 0.22 µs and weighs $2.5\times10^{-18}$ kg. With the best mirrors ever measured, $R = 1 - 1.6\times10^{-6}$ (Rempe et al. 1992), it survives 1.4 ms, stores 1.4 kJ and weighs $1.5\times10^{-14}$ kg. Those mirrors work only near one wavelength at near-normal incidence, so this is a generous upper limit for white light from a flashlight. A 0.1 µg comparator balance (general knowledge) would still need 6,500 times more.

2. **Now put your hand in.** The hand becomes the main loss, so it collapses the stored light by a factor of about 2,500 and takes in essentially the whole megawatt. Heating 1 mm of skin by 12 K (to about the point of pain) takes about 2 ms. Radiation pressure is real but tiny: at most $P/c = 3.3$ mN, the weight of 0.34 g, if all the light came from one side. In a box it comes from all sides, so there's no net push. You would feel a burn, not a push.

3. **A million real flashlights sealed in the box.** They need at least 100 m³ and weigh 150 t, so a "box" of them is bigger than a shipping container. Their bodies absorb light so well that only about $10^{-4}$ J of light is in flight at any moment. And the batteries settle the mass question: their $3.6\times10^{10}$ J is already part of the box's mass ($E/c^2 = 0.4$ mg, or $3\times10^{-12}$ of the total). Switching on moves energy from chemical to light to heat, all inside. Total mass stays constant, and only drops when heat escapes the box. Without any leak, that heat would raise the box by about 300 K.

So your question had a real answer hiding in it. The two boxes differ in mass only if energy came from outside.

**Check:** [`code/verify_light_in_box.py`](code/verify_light_in_box.py) computes invariant mass from raw totals in two frames, runs the moving-mirror photon, ray traces a cube and tracks every four-vector in a sealed box.

```text
--- 1. Invariant mass of photons (each photon has E = pc, so m = 0 alone)
[PASS] one photon has zero invariant mass: M = 0.00e+00 kg
[PASS] two parallel photons have zero invariant mass: M = 0.00e+00 kg
[PASS] two opposite photons: M = 2E/c^2: M = 6.6759e-36 kg
[PASS] two photons at 0.3 rad match sqrt(2 E1 E2 (1-cos))/c^2
[PASS] two photons at 1.0 rad match sqrt(2 E1 E2 (1-cos))/c^2
[PASS] two photons at 2.0 rad match sqrt(2 E1 E2 (1-cos))/c^2
[PASS] isotropic photon gas: M -> (sum E)/c^2: M/(sum E/c^2) = 0.99999899 for 200000 photons
[PASS] M does not change under a boost at 0.9c: energy changed by factor 2.291, M by -6.1e-14

--- 2. A photon's inertia, from mirror kicks in an accelerated box
[PASS] mean reaction force = -(E/c^2) a at a L/c^2 = 1e-05: ratio 1.00000
[PASS] mean reaction force = -(E/c^2) a at a L/c^2 = 1e-06: ratio 1.00000
[PASS] mean reaction force = -(E/c^2) a at a L/c^2 = 1e-07: ratio 1.00000
  By the equivalence principle the same force is the photon's weight (E/c^2) g on a scale.
  Check in gravity terms: redshift over 1 m is g h/c^2 = 1.091e-16; floor kick minus ceiling kick over one round trip gives E g/c^2.
```

```text
--- 3. Ray-traced light lifetime in a mirrored cube
[PASS] bare cube R = 0.9: mean path = (4V/S)/(1-R): MC 6.51 +- 0.05 m, formula 6.67 m
[PASS] bare cube R = 0.99: mean path = (4V/S)/(1-R): MC 66.15 +- 0.46 m, formula 66.67 m
[PASS] cube + black sphere: mean path = 4V/(A_abs + (1-R)S): MC 21.02 +- 0.15 m, formula 21.45 m
```

```text
--- 4. Real numbers
  1e+06 flashlights at 1.0 W of light each: P = 1e+06 W of light
  (a) Ideal 1 m cube fed through a small port, no hand in it:
      R            tau_light (s)   stored U (J)    mass of light (kg)
      0.9500000    4.448e-08       4.448e-02       4.949e-19
      0.9900000    2.224e-07       2.224e-01       2.474e-18
      0.9999984    1.390e-03       1.390e+03       1.546e-14
      best-case mass 1.55e-14 kg vs 0.1 ug comparator 1e-10 kg (short by 6467x) and 0.1 g scale 1e-04 kg
[PASS] best-case light mass is below the best balance
      stored light needed to register on a 0.1 ug balance: 8.99e+06 J
  (b) Now put a hand in (area 0.04 m^2, absorbs 60% of visible, general-knowledge figures):
      R = 0.9900000: stored 0.159 J, mass 1.77e-18 kg, hand absorbs 2.86e+05 W (29%), skin +12 K in 5.9 ms, light pressure 0.0265 Pa, squeeze 0.00106 N
      R = 0.9999984: stored 0.556 J, mass 6.18e-18 kg, hand absorbs 1e+06 W (100%), skin +12 K in 1.7 ms, light pressure 0.0926 Pa, squeeze 0.0037 N
[PASS] with good mirrors the hand, not the mirrors, limits the stored light: stored energy falls from 1.39e+03 J to 0.556 J
[PASS] skin overheats in well under a second at 1 MW: 1.7 ms
      maximum possible push if all 1e+06 W came from one side: P/c = 3.34e-03 N = the weight of 0.34 g. Heating wins by far.
[PASS] radiation force is below 5 mN even in the one-sided limit: 3.34 mN
  (c) A million real flashlights sealed inside the box:
      flashlights occupy at least 100 m^3 and weigh 150 t; they themselves absorb the light
      light in flight (with 100 m^3 of free space): 8.90e-05 J, mass 9.90e-22 kg
      battery energy 3.60e+10 J, E/c^2 = 0.40 mg of 150 t = 2.7e-12 of the total
      sealed: that energy only moves from chemical to light to heat, so total mass is unchanged
      heat raises the box by about 300 K in 3.3 h (3e+06 W electric) if nothing leaves; flashlights would fail first
```

```text
--- 5. Sealed box bookkeeping with four-vectors (units c = 1, exaggerated energies)
  energy emitted 0.03007; source rest mass 1.00000 -> 0.96992; walls 1.00000 -> 1.03007 after absorbing
[PASS] source loses rest mass about equal to the energy it emits
[PASS] walls gain rest mass about equal to the energy they absorb
[PASS] system mass is 2 before, with light in flight, and after absorption: 2.000000000000 and 2.000000000000
  The parts change; the total doesn't. With the box sealed, no scale can tell on from off.
```

The cube test shows the mean-path formula is good to about 3%, enough for these estimates. The script's moving-mirror test treats the box as rigid, which is valid only to first order in $ah/c^2$.

**Limits:** Flashlight, hand and scale numbers are round-number assumptions, so read the results as orders of magnitude. The mirror result is a best-case limit. Heating ignores blood flow and conduction, which don't matter over 2 ms. $m = E/c^2$ itself has been tested directly to $4.4\times10^{-7}$ (Rainville et al. 2005).

**Sources:** [Rempe et al., Opt. Lett. 17, 363 (1992)](https://authors.library.caltech.edu/records/t3g05-4xw31); [Rainville et al., Nature 438, 1096 (2005)](https://ideas.repec.org/a/nat/nature/v438y2005i7071d10.1038_4381096a.html).

### L02: How did we measure the speed of light?

**Status:** Derived + Code + Established + Conceptual

**Short answer:** You don't need to walk millions of kilometres. You can use distances the solar system already provides, or chop or deflect light over a few kilometres, or measure frequency and wavelength in a lab. Today $c$ isn't measured at all, because the metre is defined from it. Nearly all of these methods measure the *round trip*, and the one-way speed turns out to depend on how you choose to synchronise distant clocks.

**Starting from:** speed = distance/time, plus each experiment's geometry. For the one-way section: the laws of special relativity and a free choice of coordinates.

**Derivation:**

| Who | Method | Key numbers | Result (km/s) |
|---|---|---|---|
| Rømer, 1676 | Io's eclipses run late when Earth is farther from Jupiter | 22 min across Earth's orbit | about 227,000 with the modern AU (he didn't compute this figure himself) |
| Bradley, 1729 | Stellar aberration: rain seems to slant when you run | $\tan\alpha = v/c$, $\alpha = 20.5''$ | 299,690 with modern $v$ |
| Fizeau, 1849 | Toothed wheel | 720 teeth, 8,633 m each way, dark at 12.6 rev/s | 313,274 |
| Foucault, 1862 | Rotating mirror | Mirror turns while light is away | 298,000 |
| Rosa & Dorsey, 1907 | Electrical units, $c = 1/\sqrt{\mu_0\epsilon_0}$ (see L03) | No light used | $299{,}710 \pm 22$ |
| Michelson, 1926 | Rotating octagon, Mt Wilson to Mt San Antonio | About 35 km baseline | $299{,}796 \pm 4$ |
| Essen & Gordon-Smith, 1940s | Microwave cavity resonance, $c = f\lambda$ with $\lambda$ set by the cavity | | $299{,}792 \pm 9$ |
| Evenson et al., 1972 | Laser frequency × wavelength | 88.376181627 THz × 3.392231376 µm | $299{,}792.4562 \pm 0.0011$ |

*Fizeau's wheel.* A pulse leaves through a gap and returns after $2d/c$. If the wheel has turned half a tooth period by then, a tooth blocks it. With $N$ teeth at $f$ rev/s that takes $1/(2Nf)$, so

$$
\frac{2d}{c} = \frac{1}{2Nf} \quad\Rightarrow\quad c = 4dNf = 4(8633)(720)(12.6)\ \text{m/s} = 313{,}274\ \text{km/s}.
$$

That is 4.5% high. The true dark point is 12.06 rev/s, so an error of about half a revolution per second explains all of it. Finding the darkest point by eye is hard.

*Why Rømer could see it.* Each Io orbit shifts by at most $(v/c)P = 15$ s, too small to see once. But the shifts add up over months to about 16.6 min across the full orbit. This is his trick: he used a clock (Io) and a ruler (Earth's orbit) that were already there.

*Definition.* The 17th CGPM (1983) defined the metre as the path of light in $1/299\,792\,458$ s. The wording was revised in 2019 but $c$ stays fixed. So 299,792,458 m/s is exact by definition, and an experiment like Evenson's now calibrates *metres*, not $c$.

*The one-way speed is a convention.* To time light from A to B you need synchronised clocks at A and B, and synchronising them needs a signal whose one-way speed you already know. Reichenbach's way to make this precise is to use new coordinates $t' = t + \kappa x/c$ with any $-1 < \kappa < 1$. Then light moving in $+x$ and $-x$ has

$$
c_+ = \frac{c}{1+\kappa}, \qquad c_- = \frac{c}{1-\kappa},
$$

while the round trip is unchanged:

$$
\frac{L}{c_+} + \frac{L}{c_-} = \frac{L(1+\kappa) + L(1-\kappa)}{c} = \frac{2L}{c}.
$$

This is the same physics in different labels, so it can't predict anything different. What round-trip measurements leave open is the value of $\kappa$. Setting $\kappa = 0$ (Einstein synchronisation) is a choice, a very convenient one because it makes the laws look the same in all directions. "Light has speed $c$ both ways" isn't a measured fact in the same sense that "the round trip takes $2L/c$" is.

Doesn't Rømer measure one way? It looks like it, because the light only goes from Jupiter to Earth. But the observer's clock moves around the orbit, and in $\kappa$ coordinates that clock's reading is tied to $\kappa$ too. The script runs this: an observer moving on a set path against their own clock, light travelling at the anisotropic speed, for several $\kappa$. The arrival times on the observer's clock come out identical.

**Check:** [`code/verify_speed_of_light.py`](code/verify_speed_of_light.py) fits synthetic eclipse timings for $c$, simulates the chopped beam and locates the dark point, and solves for arrival times in $\kappa$ coordinates, keeping all the $\kappa$ terms explicit.

```text
--- Fizeau: toothed wheel
[PASS] simulated wheel goes dark where c = 4 d N f gives back c_true = 2.9979e+08: first dark at 12.05 rev/s, min throughput 0.001
[PASS] simulated wheel goes dark where c = 4 d N f gives back c_true = 3.1330e+08: first dark at 12.60 rev/s, min throughput 0.000
  Fizeau's own 12.6 rev/s: c = 313,274 km/s (+4.5% from modern)
  Sensitivity: true c would be dark at 12.06 rev/s, so a 0.5 rev/s reading error matches the 4.5% offset
[PASS] Fizeau's number reproduces about 313,000 km/s
```

```text
--- Evenson 1972: laser frequency times wavelength
[PASS] frequency x wavelength reproduces the reported value to within its 1.1 m/s error: 8.837618163e+13 Hz x 3.392231376e-06 m = 299,792,456.2 m/s
  Difference from today's exact value: 1.8 m/s

--- Reichenbach convention: nothing observable changes
[PASS] kappa = +0.0: one-way speeds 1.000c and 1.000c, round trip unchanged
[PASS] kappa = +0.5: one-way speeds 0.667c and 2.000c, round trip unchanged
[PASS] kappa = -0.7: one-way speeds 3.333c and 0.588c, round trip unchanged
[PASS] kappa = +0.9: one-way speeds 0.526c and 10.000c, round trip unchanged
[PASS] Roemer-style fit recovers c in standard coordinates: fit 2.99795e+08 m/s
[PASS] Roemer arrival times at kappa = +0.5 match kappa = 0: worst difference 3.7e-09 s
[PASS] Roemer arrival times at kappa = -0.7 match kappa = 0: worst difference 3.7e-09 s
[PASS] Roemer arrival times at kappa = +0.9 match kappa = 0: worst difference 3.7e-09 s
  So a Roemer-style observation also yields the two-way-equivalent c. The moving observer's
  clock term cancels the anisotropy, the same way slow clock transport does.
```

The Rømer fit recovers an injected $c$ to within its 1% test band with 30 s of timing noise. That tests the method, not Rømer's actual data. The $\kappa$ test shows this model is self-consistent, and I can't test it against nature, because by construction no experiment can tell these coordinates apart.

**Limits:** Rømer's own numbers and how the 22 min was reached are simplified here. The 35 km Michelson baseline is general knowledge. The convention argument assumes the usual relativity with its coordinates relabelled. A test that really found direction-dependent physics, such as a difference in round-trip times with orientation, would be real new physics and is a separate question.

**Sources:** [Fizeau's measurement (Wikipedia)](https://en.wikipedia.org/wiki/Fizeau%27s_measurement_of_the_speed_of_light_in_air); [Foucault's measurements (Wikipedia)](https://en.wikipedia.org/wiki/Foucault%27s_measurements_of_the_speed_of_light); [Speed of light (Wikipedia)](https://en.wikipedia.org/wiki/Speed_of_light) for Rømer, Bradley, Rosa & Dorsey and Essen; [Michelson's Mt Wilson measurements](https://otherhand.org/home-page/physics/historical-speed-of-light-measurements-in-southern-california/the-mount-wilson-mount-san-antonio-measurements-1922-1926/); [Evenson et al., PRL 29, 1346 (1972)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.29.1346); [CGPM 1983 Resolution 1](https://www.bipm.org/en/committees/cg/cgpm/17-1983/resolution-1); [Anderson, Vetharaniam & Stedman, Phys. Rep. 295, 93 (1998)](https://doi.org/10.1016/S0370-1573(97)00051-3).

### L03: What is light?

**Status:** Derived + Code + Established

**Short answer:** Light is a travelling disturbance of the electromagnetic field. Changing electric fields make magnetic fields and vice versa, and Maxwell's equations show that this pair sustains itself as a wave at $c = 1/\sqrt{\mu_0\epsilon_0}$. Quantum theory adds that the field's energy comes in discrete packets called photons. Radio, infrared, visible light and gamma rays are the same thing at different frequencies.

**Starting from:** Maxwell's equations in vacuum, which summarise experiments on charges and magnets:

$$
\nabla\cdot\vec E = 0,\quad \nabla\cdot\vec B = 0,\quad \nabla\times\vec E = -\frac{\partial\vec B}{\partial t},\quad \nabla\times\vec B = \mu_0\epsilon_0\frac{\partial\vec E}{\partial t}.
$$

**Derivation:** Take the curl of the third equation and substitute the fourth:

$$
\nabla\times(\nabla\times\vec E) = -\frac{\partial}{\partial t}(\nabla\times\vec B) = -\mu_0\epsilon_0\frac{\partial^2\vec E}{\partial t^2}.
$$

The identity $\nabla\times(\nabla\times\vec E) = \nabla(\nabla\cdot\vec E) - \nabla^2\vec E$ with $\nabla\cdot\vec E = 0$ gives

$$
\nabla^2\vec E = \mu_0\epsilon_0\frac{\partial^2\vec E}{\partial t^2}.
$$

This is a wave equation with speed $v = 1/\sqrt{\mu_0\epsilon_0}$. Nothing in it was about light: $\mu_0$ and $\epsilon_0$ come from forces between currents and charges. The historical surprise was that electrical measurements (Weber and Kohlrausch, 1850s, about $3.1\times10^8$ m/s) matched the measured speed of light (L02), which led Maxwell to identify light with this wave.

A plane wave $\vec E = E_0\cos(kx-\omega t)\,\hat y$ with $\vec B = B_0\cos(kx-\omega t)\,\hat z$ solves all four equations only if $\omega/k = c$ and $B_0 = E_0/c$. So light is transverse, and $\vec E$, $\vec B$ and the direction of travel are mutually perpendicular.

*Wave or particle?* Neither word quite fits. Interference needs the field description. Detectors register single clicks, and the energy delivered per click is $E = hf$. Quantum electrodynamics joins these: the field is quantised, and a photon is one quantum of excitation in a mode of it. L04 and L10 continue this.

**Check:** [`code/verify_maxwell_photon.py`](code/verify_maxwell_photon.py) does the algebra in sympy for arbitrary fields, then integrates the first-order Maxwell equations numerically (no wave equation is used) and clocks the pulse.

```text
--- 1. Wave equation from Maxwell (symbolic, arbitrary fields)
[PASS] curl curl E = grad div E - laplacian E
[PASS] plane wave solves Maxwell only if omega/k = 1/sqrt(mu eps): omega/k = 1/(sqrt(epsilon)*sqrt(mu))
[PASS] and B0 = E0 * sqrt(mu eps) (B0 = E0/c)
[PASS] div E = div B = 0 for the transverse wave
  mu0 = 1.256637061e-06 H/m, eps0 = 8.854187819e-12 F/m gives 1/sqrt(mu0 eps0) = 299,792,458 m/s
  Since 1983 this agreement is built into the values of the units, so it is not a test today.
  The historical test used separate electrical measurements: Weber and Kohlrausch, about 3.1e8 m/s.

--- 2. FDTD integration of the first-order equations
[PASS] pulse in vacuum travels at 1/sqrt(mu0 eps0): measured 2.997612e+08 m/s, 0.99990 c
[PASS] electric and magnetic energy in the wave are equal: uE/uB = 1.0000
[PASS] pulse with eps_r = 4 travels at c/2: measured 0.49993 c
```

The FDTD speed is slightly below $c$ (0.9999 $c$) because of the grid's known numerical dispersion. As the output states, the agreement of $1/\sqrt{\mu_0\epsilon_0}$ with 299,792,458 m/s isn't an experimental test today, since the SI values are tied to $c$ by construction.

**Limits:** This is classical vacuum electromagnetism. Effects where the field's quantum nature matters, such as single-photon detection, need QED, which I'm citing and not deriving.

**Sources:** [Rudolf Kohlrausch (Wikipedia)](https://en.wikipedia.org/wiki/Rudolf_Kohlrausch) for the 1850s electrical measurement.

### L04: Does light have mass?

**Status:** Derived + Code + Established

**Short answer:** A photon has zero rest mass, and experiments limit it to below $10^{-18}$ eV/$c^2$. But it has energy and momentum, $p = E/c$, which is how it pushes on things. Your first two sentences need correcting: light doesn't *emit* photons, light *is* photons (quanta of the field), and photons don't have rest mass.

**Starting from:** Special relativity's four-momentum $p^\mu = m\,dx^\mu/d\tau$, and its invariant length.

**Derivation:** For a body of rest mass $m$ at speed $v$, $E = \gamma mc^2$ and $\vec p = \gamma m\vec v$ with $\gamma = (1 - v^2/c^2)^{-1/2}$. Then

$$
E^2 - (pc)^2 = \gamma^2m^2c^4\left(1 - \frac{v^2}{c^2}\right) = (mc^2)^2,
$$

so $E^2 = (pc)^2 + (mc^2)^2$, and dividing the definitions gives $v = pc^2/E$. Light travels at $c$, so $p = E/c$, and then the relation forces $m = 0$. The relation is what survives when $m \to 0$ and $\gamma \to \infty$ together, leaving finite $E$ and $p$.

*Photon momentum has independent evidence.* You don't need the photon idea at all. Classically, shine a wave on a damped charge. The electric field shakes it along $y$ at velocity $v_y$ and delivers power $qE_yv_y$. The wave's magnetic field gives a force $qv_yB_z = qv_yE_y/c$ along $x$. Force per power is $1/c$, so absorbed momentum is $E/c$. Then Compton scattering uses it per photon: treating a photon as a particle with $p = h/\lambda$ gives a wavelength shift $\Delta\lambda = (h/m_ec)(1-\cos\theta)$, as measured.

*Rest mass and "relativistic mass".* Older books call $E/c^2$ the relativistic mass. It's frame dependent, and the same number as energy. For a 2 eV photon it is $3.6\times10^{-36}$ kg. Most physicists now use "mass" only for the invariant $m$, following arguments such as Okun's (1989). Either way, it's a naming choice about the same physics. In that language: zero rest mass and non-zero energy.

*Experimental limits.* The Particle Data Group's listed limit is $m_\gamma < 1\times10^{-18}$ eV from solar-wind magnetohydrodynamics (Ryutov 2007), a model-dependent limit. Fast radio burst timing gives a more direct but weaker limit, $m_\gamma < 2.1\times10^{-15}$ eV at 68% confidence (Wang et al. 2023 in the PDG listing). No experiment can prove exactly zero, and only upper limits exist. That zero is also a theoretical expectation from the gauge symmetry of electromagnetism.

*Then why does gravity bend light, and why does a box of light weigh something?* In general relativity, gravity couples to energy and momentum, not rest mass (L16). And a set of photons going in different directions has invariant mass even if each has none (L01). A single beam of exactly parallel photons doesn't.

**Check:** [`code/verify_maxwell_photon.py`](code/verify_maxwell_photon.py). The momentum ratio comes from integrating the equation of motion. The Compton shift comes from numerically solving the three conservation equations with no shift formula put in.

```text
--- 3. A wave carries momentum p = E/c, from classical absorption
[PASS] momentum absorbed / (energy absorbed / c) = 1: ratio = 1.00000

--- 4. E^2 = (pc)^2 + (mc^2)^2
[PASS] (gamma m c^2)^2 - (gamma m v c)^2 = (m c^2)^2
[PASS] v = c exactly when m = 0, v < c for m > 0

--- 5. Compton scattering from conservation laws only
  Compton wavelength h/(m_e c) = 2.4263 pm
[PASS] shift at 30 deg equals lambda_C (1 - cos): 0.3251 pm
[PASS] shift at 90 deg equals lambda_C (1 - cos): 2.4263 pm
[PASS] shift at 135 deg equals lambda_C (1 - cos): 4.1420 pm
[PASS] shift at 180 deg equals lambda_C (1 - cos): 4.8526 pm

--- 6. Photon mass limit in context
  Limit m < 1e-18 eV/c^2 = 1.78e-54 kg = 2.0e-24 electron masses
  If m were at that limit, a 2.0 eV photon would lag light by 1 - v/c = 1.3e-37
  The quantity E/c^2 for that photon ('relativistic mass') is 3.57e-36 kg; it is not the invariant mass. Momentum E/c = 1.07e-27 kg m/s
```

**Limits:** The classical absorption test uses a low-velocity charge. The mass limits are the quoted values, and the PDG listing itself warns that several older limits have been critiqued. I have not re-derived any of those limits.

**Sources:** [PDG photon listing 2024, Navas et al., Phys. Rev. D 110, 030001](https://pdg.lbl.gov/2024/listings/rpp2024-list-photon.pdf); [Okun, Physics Today 42(6), 31 (1989)](https://physicstoday.aip.org/features/the-concept-of-mass); Compton (1923) is cited from general knowledge.

### L05: What is a gravitational field?

**Status:** Derived + Code

**Short answer:** It's the gravitational force per unit mass at each point in space: $\vec g(\vec x)$ tells any mass $m$ placed there that it feels a force $m\vec g$. It's a way to describe gravity at every place without saying what the test object is. In general relativity the same role is played by the geometry of spacetime.

**Starting from:** Newton's law of gravitation, $\vec F = -Gm_1m_2\hat r/r^2$, and the fact that forces add.

**Derivation:** Divide by the test mass to get the field of a point mass $M$:

$$
\vec g = -\frac{GM}{r^2}\hat r, \qquad \vec g(\vec x) = \sum_i -\frac{Gm_i(\vec x-\vec x_i)}{|\vec x-\vec x_i|^3}.
$$

*Gauss's law.* For a point mass inside a closed surface, the flux of $\vec g$ through a small piece of surface is $-GM\,d\Omega$, where $d\Omega$ is the solid angle of that piece, because the area grows as $r^2$ while the field falls as $1/r^2$. Adding over all the surface gives $-4\pi GM$. A mass outside contributes zero, since each cone of solid angle crosses the surface in and out equally. So

$$
\oint\vec g\cdot d\vec A = -4\pi GM_{\rm enclosed}.
$$

*Shell theorem.* By symmetry, the field of a spherical shell is radial and depends only on $r$. Take a sphere of radius $r$ as the surface. Outside, $4\pi r^2g = -4\pi GM$ gives $g = -GM/r^2$, as if all the mass were at the centre. Inside, $M_{\rm enclosed} = 0$ so $g = 0$. A solid uniform ball, being shells, has $g = -GMr/R^3$ inside.

*Is the field real or bookkeeping?* In Newton's theory you could regard it as bookkeeping, since the force changes instantly and the field adds nothing you can test. In general relativity the changes spread at $c$ and carry energy away as gravitational waves, so the field has its own existence. Locally, you can cancel the field by falling freely (the equivalence principle). You can't cancel the tidal differences between nearby points, and those are the geometry.

**Check:** [`code/verify_gravity_field.py`](code/verify_gravity_field.py) adds up Newton's law over point masses, using no shell theorem or Gauss's law.

```text
--- Shell theorem by direct summation
[PASS] inside a shell at r = 0.0 R: field vanishes: |g|/(GM/R^2) = 2.9e-09
[PASS] inside a shell at r = 0.3 R: field vanishes: |g|/(GM/R^2) = 5.7e-09
[PASS] inside a shell at r = 0.6 R: field vanishes: |g|/(GM/R^2) = 2.9e-08
[PASS] outside a shell at r = 1.5 R: |g| = GM/r^2: ratio 1.0000000
[PASS] outside a shell at r = 3.0 R: |g| = GM/r^2: ratio 1.0000000
[PASS] outside a shell at r = 10.0 R: |g| = GM/r^2: ratio 1.0000000

--- Solid uniform ball: Gauss's law check
  523984 lattice points, spacing R/50
[PASS] uniform ball at r = 0.25 R: inward field matches Gauss's law: ratio 0.99927
[PASS] uniform ball at r = 0.5 R: inward field matches Gauss's law: ratio 0.99930
[PASS] uniform ball at r = 0.8 R: inward field matches Gauss's law: ratio 0.99930
[PASS] uniform ball at r = 2.0 R: inward field matches Gauss's law: ratio 0.99999
```

**Limits:** Newtonian: weak fields and speeds far below $c$. The lattice ball has a stepped surface, so agreement is good to the stated tolerance, not exact.

### L06: Is it related to the centre of mass?

**Status:** Derived + Code

**Short answer:** Related, but not the same. The field depends on where all the mass is, not only on the COM. For a sphere it looks exactly like a point mass at the centre, and far from any object it approaches that. Near a non-spherical object it can be quite different, and it need not even point at the COM.

**Starting from:** The superposition sum from L05.

**Derivation:** Put the origin anywhere and expand each term for $r \gg |\vec x_i|$. The potential is

$$
\Phi(\vec x) = -G\sum_i\frac{m_i}{|\vec x-\vec x_i|} = -G\left[\frac{M}{r} + \frac{\vec D\cdot\hat r}{r^2} + \frac{\hat r\cdot Q\cdot\hat r}{2r^3} + \dots\right],
$$

where $M = \sum m_i$, the dipole is $\vec D = \sum m_i\vec x_i$, and the quadrupole is $Q_{jk} = \sum m_i(3x_{ij}x_{ik} - |\vec x_i|^2\delta_{jk})$.

Here is the exact connection: $\vec D = M\vec x_{\rm COM}$, so **the dipole vanishes if the origin is at the COM**. The COM is therefore the unique best single point to stand in for the body. Treating the mass as sitting there makes an error falling as $1/r^4$ in the field. Any other point leaves an error of $1/r^3$.

So the field's *leading term* cares about the COM, and everything after that cares about shape.

A dumbbell shows this: two 0.5 kg masses 1 m apart. At its own COM, the field is exactly zero, not infinite as a point mass there would be. On the axis at 1 m it is 2.22 times $GM/r^2$. Broadside at 1 m it is 0.72 times. On a diagonal, the field misses the COM by 20°. The same idea applies to the donut in F05 in section 01: its COM is in the hole, and the field there is zero.

The COM has another role that's easy to confuse with this one. In a *uniform* external field, the net force on a body acts through its COM. That's about how a body *responds* to an outside field, while the multipole expansion is about the field the body *makes*.

**Check:** [`code/verify_gravity_field.py`](code/verify_gravity_field.py).

```text
--- L06: the field is not the field of a point at the COM
[PASS] dumbbell field at its COM is zero, not infinite: |g| = 0.0e+00
  on the axis at 1 m: |g|/(GM/r^2) = 2.2222, field misses the COM direction by 0.00 deg
  broadside at 1 m: |g|/(GM/r^2) = 0.7155, field misses the COM direction by 0.00 deg
  diagonal at 1 m: |g|/(GM/r^2) = 1.1216, field misses the COM direction by 20.04 deg
[PASS] at 100 separations the point-at-COM picture is good to 1e-4
```

**Limits:** Newtonian. The expansion fails inside the smallest sphere that contains the mass.

### L07: A ball rattling inside a thrown ball

**Status:** Derived + Code + Estimate

**Short answer:** Yes, the field outside changes at every instant, but only slightly, and not because of the COM. The COM keeps moving on the same parabola no matter how the inner ball rattles, and total mass stays fixed. What changes is the quadrupole and higher moments, roughly 0.1% of the field at 1 m in my example, and each change reaches a distant point after a delay of $r/c$.

**Starting from:** The multipole expansion from L06. Newton's third law for internal forces. The fact from general relativity that gravitational changes propagate at $c$. Example: shell 1 kg, radius 10 cm; inner ball 0.1 kg, radius 2 cm.

**Derivation:**

*COM motion can't change.* Internal forces come in equal and opposite pairs, so they cancel in the total: $M\ddot{\vec x}_{\rm COM} = \vec F_{\rm ext} = M\vec g$. Rattling makes the *shell* wiggle but the COM doesn't notice. In a throw, the COM follows one parabola.

*Term by term, about the COM.*

- **Monopole** $M$: constant. So the leading field is unchanged. (A tiny caveat from L01: if the collisions turn kinetic energy into heat, the mass doesn't change until that heat leaves.)
- **Dipole** $\vec D$: zero about the COM at all times. Measured from any fixed point, $\vec D = M\vec x_{\rm COM}$, so it changes only with the COM's parabola. Rattling adds nothing. This is why gravity has no dipole radiation: a changing mass dipole would mean momentum not conserved.
- **Quadrupole**: this one does change. With inner ball offset $\vec d$ from the shell centre, two spherical bodies have $Q_{jk} = \mu(3d_jd_k - d^2\delta_{jk})$ about the COM, where $\mu = Mm/(M+m)$ is the reduced mass. Its field scales as $3G\mu d^2/r^4$, which for $d = 8$ cm at 1 m is about $10^{-13}$ m/s².

Note that if you choose the shell's centre as the origin rather than the COM, a "dipole" $m\vec d$ appears and varies. That's the same physics in a worse coordinate system: it is the shift of the COM away from that centre. The script confirms this: the residual falls as $r^{-3}$ about the shell's centre and as $r^{-4}$ about the COM.

*Delay and radiation.* Newton's formula responds instantly, which is wrong. In general relativity, a change at the source arrives $r/c$ later: 3.3 ns at 1 m. A changing quadrupole also radiates, with strain $h \approx 2G\ddot Q/(c^4r)$ and power $P \approx (G/5c^5)\dddot Q^2$. For a 10 Hz rattle that's $h \sim 10^{-44}$ at 1 m and $P \sim 10^{-49}$ W, as an order-of-magnitude estimate.

So both parts of your question are right, at different orders: the COM and monopole give nearly all of the field, and that part doesn't change. The rattling really does change the field, but only in the smaller terms that depend on shape.

**Check:** [`code/verify_gravity_field.py`](code/verify_gravity_field.py) sums hundreds of thousands of point masses and compares with the multipole form. The rattle integrates each body separately through contact impulses, so the COM's parabola is a result, not something assumed.

```text
--- L07: multipole expansion against direct summation
[PASS] dipole about the total COM is zero: 4.3e-19 kg m
[PASS] quadrupole is traceless
  reduced mass 0.0909 kg, d = 6.16 cm, Q eigenvalues [-0.3455 -0.3455  0.6909] x 1e-3 kg m^2
[PASS] largest Q eigenvalue = 2 mu d^2 (two-point analytic)
  r (m)  |g_direct|     residual after monopole   residual after monopole+quadrupole
   0.4   4.57516e-10    1.3431e-12              1.6897e-13
   0.6   2.03671e-10    2.6629e-13              2.2254e-14
   1.0   7.33828e-11    3.4517e-14              1.7089e-15
   1.6   2.86734e-11    5.2584e-15              1.5775e-16
   2.5   1.17459e-11    8.8005e-16              1.6028e-17
[PASS] after removing the monopole the residual falls as r^-4 (quadrupole): fitted power -4.000
[PASS] after removing the quadrupole too it falls as r^-5 (octupole): fitted power -5.053
[PASS] quadrupole removes most of the residual at r = 1 m: 0.050 remaining
[PASS] expanding about the shell centre instead leaves a dipole, r^-3: fitted power -2.996, size 4.122e-13 m/s^2 at 1 m, dipole m d = 6.16e-03 kg m

--- L07: a thrown ball with a rattling ball inside
[PASS] COM stays on one parabola through all the collisions: 37 collisions, worst deviation 9.0e-15 m
  shell centre wiggles about the COM by up to 7.26 mm
  At 1 m from the COM: monopole field 7.3417e-11 m/s^2
  Variation from rattling 9.528e-14 m/s^2, fraction 1.30e-03
[PASS] rattling variation is of order the quadrupole estimate 3 G mu d^2 / r^4: estimate 1.16e-13 m/s^2
[PASS] rattling changes the field by much less than the steady monopole

--- Propagation and radiation
  Field change reaches 1 m after 3.34 ns (retardation at c, GR; Newton's formula has none)
  Quadrupole radiation for a 10 Hz rattle: strain at 1 m ~ 3.8e-44, power ~ 1.1e-49 W (order of magnitude only)
[PASS] gravitational-wave power from the rattle is absurdly small
```

The COM deviation of about $10^{-14}$ m is rounding error. The direct sums use sampled balls, so residuals are only trustworthy to the fitted precision shown.

**Limits:** Elastic collisions, no spin, and uniform external gravity. In a non-uniform field a tidal coupling does exist between the rattling and the COM, which is minute here. The radiation estimates ignore angular factors. The field numbers are far below what any instrument could measure for such small masses.

### L08: Why still can't we explain the Heisenberg uncertainty principle?

**Status:** Derived + Code + Established

**Short answer:** We can. It's a theorem of quantum mechanics, proved in a few lines from the Cauchy-Schwarz inequality, and it holds even for classical waves as a property of Fourier transforms. It says a state can't have sharply defined position and momentum at once, in the sense that $\sigma_x\sigma_p \ge \hbar/2$. What people still debate is what the state *is*, which is a separate question (L09).

**Starting from:** States are vectors $|\psi\rangle$ in a complex vector space, observables are Hermitian operators, and expectation values are $\langle A\rangle = \langle\psi|A|\psi\rangle$. For position and momentum, $[x,p] = i\hbar$, which follows from $p = -i\hbar\,\partial_x$. These are the postulates, and the theorem doesn't explain why they hold.

**Derivation (Robertson, 1929):** Define shifted operators $\hat a = A - \langle A\rangle$, $\hat b = B - \langle B\rangle$ and vectors $|f\rangle = \hat a|\psi\rangle$, $|g\rangle = \hat b|\psi\rangle$. Then $\sigma_A^2 = \langle f|f\rangle$ and $\sigma_B^2 = \langle g|g\rangle$.

1. Cauchy-Schwarz gives $\sigma_A^2\sigma_B^2 = \langle f|f\rangle\langle g|g\rangle \ge |\langle f|g\rangle|^2$.
2. Split the product: $\langle f|g\rangle = \langle\hat a\hat b\rangle = \tfrac12\langle[\hat a,\hat b]\rangle + \tfrac12\langle\{\hat a,\hat b\}\rangle$.
3. A commutator of Hermitian operators is anti-Hermitian, so its expectation is purely imaginary. An anticommutator is Hermitian, so its expectation is real. Therefore $|\langle f|g\rangle|^2 = \tfrac14|\langle[\hat a,\hat b]\rangle|^2 + \tfrac14\langle\{\hat a,\hat b\}\rangle^2$.
4. Drop the second, non-negative term, and note $[\hat a,\hat b] = [A,B]$:

$$
\sigma_A\sigma_B \ge \tfrac12\big|\langle[A,B]\rangle\big|.
$$

For position and momentum, $\sigma_x\sigma_p \ge \hbar/2$. Equality needs $|f\rangle \propto |g\rangle$ with a purely imaginary constant, and that condition gives Gaussian wave packets.

*The picture behind it.* A wave packet is built from plane waves $e^{ikx}$, and $p = \hbar k$. Making it narrow takes many values of $k$, and using few values of $k$ makes it long. A short click has no definite pitch for the same reason. The trade-off is plain Fourier mathematics that holds for sound. The quantum content is the single extra fact $p = \hbar k$.

*Spread is not disturbance.* The popular version says "measuring position kicks the momentum". That is a different statement. $\sigma_x$ and $\sigma_p$ belong to the *state*, and can be measured on separate copies, with no copy measured twice. Heisenberg's original 1927 argument did concern disturbance, and the naive form of that error-disturbance relation is actually violated, as Ozawa showed theoretically and as neutron (Erhart et al. 2012) and photon (Rozema et al. 2012) experiments confirmed. Robertson's inequality wasn't affected by those experiments, since it's about something else.

**Check:** [`code/verify_uncertainty.py`](code/verify_uncertainty.py) computes spreads from FFTs of sampled wavefunctions, and tests the general inequality on random matrices and states. The last section samples separate copies for $x$ and for $p$.

```text
--- Commutator
[PASS] [x, p] f = i hbar f: result I*hbar*f(x)

--- Fourier spreads (hbar = 1)
[PASS] Gaussian sigma = 0.3: product = 1/2: sx sk = 0.5000000000
[PASS] Gaussian sigma = 1.0: product = 1/2: sx sk = 0.5000000000
[PASS] Gaussian sigma = 4.0: product = 1/2: sx sk = 0.5000000000
  Chirped Gaussian: sx sk = 1.1180 (same |psi|, but phase structure widens momentum)
[PASS] chirp raises the product above 1/2
[PASS] infinite well ground state matches sqrt(pi^2/3 - 2)/2: FFT 0.56818, exact 0.56786
  Two separated lumps: sx sk = 3.571
  sech packet: sx sk = 0.5236
[PASS] sech packet respects the bound
[PASS] 1000 random packets: smallest product is at least 1/2: smallest 0.50000000

--- Robertson bound with random operators and states
[PASS] sigma_A sigma_B >= |<[A,B]>|/2 in 5000 random trials: smallest margin 1.57e-07
[PASS] Schroedinger's stronger form also holds: smallest margin -1.09e-14
[PASS] commuting operators allow both spreads to be zero

--- State spread with no disturbance: separate copies
[PASS] independent ensembles reproduce the state's spread product: sampled 0.5166, exact 0.5167
[PASS] that product is above 1/2 though no copy was measured twice
```

**Limits:** Numerical tests can't replace the proof, and finite samples can't cover all states. For operators in finite dimensions the bound is state dependent and can be zero. For an angle variable or a particle on a ring, the plain form needs care because $x$ isn't a good operator there. The energy-time relation is different, since $t$ isn't an operator in this theory.

**Sources:** [Robertson, Phys. Rev. 34, 163 (1929)](https://journals.aps.org/pr/abstract/10.1103/PhysRev.34.163); [Rozema et al., PRL 109, 100404 (2012)](https://arxiv.org/abs/1208.0034); [Erhart et al., Nature Phys. 8, 185 (2012), via TU Wien](https://repositum.tuwien.at/handle/20.500.12708/161368).

### L09: Why can't we explain quantum mechanics? Understanding against accuracy

**Status:** Derived + Code + Conceptual + Open

**Short answer:** Your understanding-versus-accuracy trade-off is a nice metaphor, but it isn't physics. The uncertainty principle is about spreads of a physical state, and says nothing about how well anyone knows the theory. Quantum mechanics is extremely accurate *and* its rules are clearly understood. What remains unsettled is what those rules say reality is, and physicists do disagree about that.

**Starting from:** The theorem in L08, plus what the word "explain" could mean. What counts as understanding isn't something physics can derive, so that part is philosophy, and I'm flagging it as such.

**Derivation:**

*Why the analogy fails as physics.* The inequality $\sigma_A\sigma_B \ge \tfrac12|\langle[A,B]\rangle|$ needs two operators on a state's vector space, with a non-zero commutator. "Understanding" and "accuracy of a theory" aren't operators acting on a quantum state. They have no commutator, so the theorem says nothing about them. Even for real observables the bound isn't universal: commuting operators allow both spreads to be zero, as the script shows in L08.

There's also direct evidence against it. Better understanding often *raises* accuracy. The electron's magnetic moment has been measured as $g/2 = 1.00115965218059(13)$, to 0.13 parts per trillion (Fan et al. 2023). Theory matches at that level only because the quantum field theory behind it became better understood and could be computed to higher orders.

*Credit where it's due.* You're in good company. Bohr himself loved extending complementarity beyond physics, and a line "truth and clarity are complementary" is attributed to him, though I could only find it on a quotation site with no primary source, so treat it as unverified. As a remark about the limits of plain language, it may be valuable. It's just not a consequence of the uncertainty principle.

*Three meanings of "explain".*

1. **Can we state the rules and calculate?** Yes. This is settled mathematics, and its predictions work extremely well.
2. **Can we picture it using everyday objects?** Not fully, and there's a reason. Our intuition was built on objects that have definite position and momentum at once, and L08 shows quantum states don't. Bell-type correlations also rule out a large family of simple pictures. The script in L12 computes the CHSH value of $2\sqrt2$ that exceeds the classical limit of 2.
3. **What is the wavefunction, physically, and what happens in a measurement?** **Open.** Copenhagen-type, many-worlds, pilot-wave, epistemic and objective collapse views all exist. Most of them give the same predictions for experiments so far. A 2025 *Nature* survey of over 1,100 researchers found no consensus. Press reports of it put Copenhagen first at about 36%, and only about a quarter of respondents were confident in their preferred view.

So your feeling that something isn't understood is reasonable. But the gap is in meaning, not in accuracy, and a better account wouldn't have to lose any accuracy.

**Check:** [`code/verify_uncertainty.py`](code/verify_uncertainty.py) covers the mathematical part (L08's output). No script can test the philosophical claims.

**Limits:** The argument that understanding isn't an operator shows only that the theorem doesn't apply to it. It doesn't prove that no trade-off of any kind could ever exist between simple pictures and precise theories, which is a question about people and models and not a theorem. The survey numbers come from press summaries of the Nature article, and a survey doesn't settle physics.

**Sources:** [Fan et al., PRL 130, 071801 (2023)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.130.071801); [Report on the Nature 2025 survey](https://thequantuminsider.com/2025/08/02/a-century-into-quantum-mechanics-physicists-still-cant-agree-what-it-means-nature-survey-shows/) and [CERN Courier](https://cerncourier.com/a/the-measurement-problem-measured/).

### L10: Could light and black holes unify QM and relativity?

**Status:** Derived + Code + Estimate + Open + Speculative

**Short answer:** Half of this is exactly where serious work is happening. Black holes are the best known theoretical testing ground for quantum gravity. The other half doesn't work as stated: light being "a wave and a particle" is already handled by quantum field theory, and isn't a gap between QM and GR. No unified theory exists yet, so any route is unproven. Not amateur at all, though, and the idea has a useful sharper version.

**Starting from:** Quantum field theory, general relativity and thermodynamics, as separately tested theories. The combination below is a heuristic and extends beyond anything tested.

**Derivation:**

*Where your idea is strong: black holes.* A black hole is where you can't ignore either theory: big curvature *and* quantum fields. Hawking (1975) found that quantum fields near a horizon radiate as a hot body.

Heuristic for the temperature, flagged as such. Let's see what the constants allow:

- *Units.* Asking for $\hbar^ac^bG^dM^e$ with units of energy leaves one free exponent (the script solves this). Requiring the energy to be first order in $\hbar$, a quantum effect on a classical background, fixes $k_BT \sim \hbar c^3/(GM)$. This needs that extra assumption, and units can't give the number in front.
- *Uncertainty picture.* A quantum near the hole is confined to about $\Delta x \sim r_s = 2GM/c^2$. L08 then gives energy $\sim \hbar c/r_s$, the same form. This is only a picture, and the "pair at the horizon" story isn't the real calculation.
- *The number in front.* Remove the conical singularity at the horizon in imaginary time. Near $r_s$ the metric looks like flat polar coordinates only if imaginary time has period $\beta = 4\pi r_s/c$. Identifying period with inverse temperature, as in thermal field theory, gives

$$
T_H = \frac{\hbar c}{4\pi r_sk_B} = \frac{\hbar c^3}{8\pi GMk_B} = 6.17\times10^{-8}\ \text{K}\left(\frac{M_\odot}{M}\right).
$$

This agrees with Hawking's result, which comes from a much fuller field theory calculation. The periodicity argument assumes that identification.

Bigger holes are colder. Stellar-mass holes are far colder than the 2.7 K sky, so none of this has been observed from an astrophysical hole. Hawking radiation has wavelength about $8\pi^2r_s$, some 80 times the hole's size.

*The information paradox.* If evaporation is exactly thermal, a pure state turns into a mixed one, which contradicts unitary quantum mechanics. Recent calculations (Penington 2020; Almheiri et al. 2019) reproduce a Page curve, suggesting information does come out, in particular settings such as AdS/CFT. How this works for real black holes is **Open**. Section B covers this in depth, from B25 in section 06 on.

*Where the idea needs correction: light.* Wave-particle duality isn't a QM-GR gap. It's the old-quantum-theory puzzle that QFT resolved: one quantum field, whose excitations show interference and also deposit energy in lumps. That is compatible with special relativity, and tested very precisely (L09). Superposition is broader than wave-particle duality: it means states add, including for spin or for positions. Light is a good example, but not a special clue.

The actual gap is different: general relativity makes spacetime itself dynamical, while standard QFT assumes a fixed background. And gravity's coupling $G$ has units that make perturbation theory break down near the Planck energy, $\sqrt{\hbar c^5/G} = 1.2\times10^{19}$ GeV.

*The sharper version of your idea: light as a probe.* These are real uses of what you suggested:

- **Photon timing.** Some quantum gravity ideas predict that photon speed depends slightly on energy. A linear effect at the Planck scale would delay a 31 GeV photon from GRB 090510 by about 1 s. Fermi saw no such delay, placing that specific scale above the Planck mass (Fermi collaborations 2009). Other models aren't excluded.
- **A box of light.** The box in L01 was central to the 1930 Einstein-Bohr debate. Einstein weighed a box to get a photon's energy. Bohr replied using gravitational redshift of the box's clock (L16's physics) to restore the energy-time uncertainty. So light, gravity and quantum uncertainty have been linked for almost a century, though later authors question whether Bohr's use of relativity was essential.
- **Entanglement through gravity.** Bose et al. (2017) proposed testing whether gravity can entangle two masses, which would show gravity has quantum features. It's proposed, not yet performed to my knowledge. **Speculative** until then.

**Check:** [`code/verify_hawking.py`](code/verify_hawking.py) solves the unit equations and numerically finds the imaginary-time period that makes the origin smooth.

```text
--- Dimensional analysis
  hbar^a c^b G^d M^e with units of energy: a = 1/2 - e/2, b = 5/2 - e/2, d = e/2 - 1/2, e free
  With energy linear in hbar: {a: 1, b: 3, d: -1, e: -1}
[PASS] kT ~ hbar c^3 / (G M)
  Note: units alone leave one free exponent, because (hbar c^5/G)^(1/2) is itself an energy, the
  Planck energy. The step 'linear in hbar' is an extra assumption, and units cannot give the number in front.

--- Euclidean regularity at the horizon (units c = 1, r_s = 1)
  r = r_s + 0.01: Euclidean period required for a smooth origin = 12.650063 r_s/c
  r = r_s + 0.0001: Euclidean period required for a smooth origin = 12.567208 r_s/c
  r = r_s + 1e-06: Euclidean period required for a smooth origin = 12.566380 r_s/c
[PASS] the period is 4 pi r_s / c: found 4.00000 pi
  T = hbar/(k beta) with beta = 4 pi r_s/c and r_s = 2GM/c^2 gives T = hbar c^3/(8 pi G M k)

--- Numbers
  Solar mass: T = 6.170e-08 K
[PASS] solar-mass hole is near 6e-8 K
  Hole as warm as the CMB (2.725 K): M = 4.50e+22 kg, about 0.61 Moon masses
  So every stellar-mass hole now absorbs more CMB than it radiates: T_sun/T_CMB = 2.3e-08
  10 M_sun: T = 6.17e-09 K
  Sgr A* (4.3e6 M_sun): T = 1.43e-14 K
  1e12 kg: T = 1.23e+11 K
[PASS] thermal wavelength hc/(kT) = 8 pi^2 r_s, far larger than the hole: 79.0 r_s
  Planck length 1.616e-35 m, Planck energy 1.221e+19 GeV

--- Light as a probe: a Planck-suppressed delay from GRB 090510
  Linear Lorentz violation at the Planck energy would delay a 31.0 GeV photon by 1.13 s
  That delay is what makes the burst's sub-second structure a Planck-scale test
[PASS] delay at the Planck scale is order 1 s, so the test really reaches that scale
```

The script can't check Hawking's actual calculation, or the idea that a smooth Euclidean geometry sets the physical temperature. The timing estimate uses simple $\Lambda$CDM values and is only an order of magnitude.

**Limits:** Hawking radiation has never been observed from an astrophysical black hole. The heuristics don't derive it. Nothing here shows that black holes *will* give unification, and that part stays speculative.

**Sources:** [Hawking, Commun. Math. Phys. 43, 199 (1975)](https://link.springer.com/doi/10.1007/BF01608497); [Penington, JHEP 09, 002 (2020)](https://arxiv.org/abs/1905.08255); [Almheiri et al., JHEP 12, 063 (2019)](https://arxiv.org/abs/1905.08762); [Fermi GBM/LAT, Nature 462, 331 (2009)](https://arxiv.org/abs/0908.1832); [Bose et al., PRL 119, 240401 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.240401); [The 1930 photon box revisited, arXiv:1203.1139](https://arxiv.org/abs/1203.1139).

### L11: Would I be invisible in a mirror at light speed?

**Status:** Derived + Code + Established

**Short answer:** No. You can't reach light speed, so the setup as stated has no answer in the theory. At any speed below $c$, if the mirror moves with you, you see yourself perfectly normally: the light leaves your face at $c$ in your own frame and returns in $2d/c$, whatever your speed. You had the classic question, the one Einstein reported wondering about as a teenager, and your reasoning is exactly what Newtonian velocity addition predicts. Nature doesn't work that way.

**Starting from:** Einstein's two postulates: the laws of physics are the same in every inertial frame, and light travels at $c$ in every inertial frame, regardless of its source. The second has direct evidence: $\pi^0$ mesons moving at about $0.99975c$ emit photons at $c$ (Alväger et al. 1964).

**Derivation:**

*In your own frame.* You are at rest in your own frame. The postulates say nothing physical can tell you that you're "moving", so the mirror works as in your bedroom.

*Velocity addition.* From the Lorentz transformation, a speed $u$ in one frame becomes in a frame moving at $v$

$$
u' = \frac{u - v}{1 - uv/c^2}.
$$

Put in $u = c$:

$$
u' = \frac{c - v}{1 - v/c} = c .
$$

So light from your face moves away from you at $c$ for every $v$, not at $c - v$. That's what's wrong in the Newtonian argument that gives closing speed zero.

*What an outside observer sees.* They see your mirror ahead of you, moving at $v$. Light chases it at closing speed $c - v$ and reaches it after a long time. Then the light rushes back to you at closing speed $c + v$. With the contracted distance $d/\gamma$:

$$
t_{\rm out} = \frac{d/\gamma}{c - v} = \frac{\gamma d(1+\beta)}{c}, \qquad t_{\rm back} = \frac{\gamma d(1-\beta)}{c}, \qquad t_{\rm total} = \frac{2\gamma d}{c}.
$$

Your clock runs slow by $\gamma$, so it reads $2d/c$. Everybody agrees that you see your image, and they differ only on the time. At $\beta = 0.99$ with $d = 1$ m: 47 ns out and 0.24 ns back for them, 6.67 ns on your clock. The outsider also sees the forward light blue-shifted by $\sqrt{(1+\beta)/(1-\beta)}$ and the returning light red-shifted by the same factor, while you see no shift at all.

*A mirror fixed to the ground.* Now you run at it. Light still leads you and comes back. But the reflected light arrives Doppler-shifted by $(1+\beta)/(1-\beta)$. At $\beta = 0.99$ that's a factor of 199, pushing 550 nm light to 2.8 nm, soft X-rays. So you wouldn't be invisible; your image would be in the wrong colour for your eyes.

*Why $v = c$ isn't an allowed setup.* The energy of a massive object is $(\gamma - 1)mc^2$, which diverges as $v \to c$, so no finite energy gets you there. A constant force gives $p = Ft$, and $v = p/\sqrt{m^2 + p^2/c^2}$ stays below $c$ for all time. Stacked boosts also stay below $c$. Also, $\gamma$ is infinite at $v = c$, so there's no transformation to such a frame. And that frame would contradict the second postulate, since light in it would be at rest. So "what does the person riding at $c$ see" has no consistent answer. It's like asking for the north of the North Pole.

**Check:** [`code/verify_relativity_motion.py`](code/verify_relativity_motion.py) finds meeting events numerically in the outside frame, then Lorentz transforms them to your frame, and integrates a constant force. The stacked boosts use exact fractions.

```text
--- Mirror carried with you: events in the outside frame, transformed to yours
[PASS] beta = 0.0: you see your image after 2d/c, mirror at d: your frame 6.6713 ns; outside frame: out 3.34 ns, back 3.34 ns, total 6.67 ns
[PASS] beta = 0.5: you see your image after 2d/c, mirror at d: your frame 6.6713 ns; outside frame: out 5.78 ns, back 1.93 ns, total 7.7 ns
[PASS] beta = 0.9: you see your image after 2d/c, mirror at d: your frame 6.6713 ns; outside frame: out 14.5 ns, back 0.765 ns, total 15.3 ns
[PASS] beta = 0.99: you see your image after 2d/c, mirror at d: your frame 6.6713 ns; outside frame: out 47.1 ns, back 0.236 ns, total 47.3 ns
[PASS] beta = 0.9999: you see your image after 2d/c, mirror at d: your frame 6.6713 ns; outside frame: out 472 ns, back 0.0236 ns, total 472 ns
[PASS] beta = 0.999999: you see your image after 2d/c, mirror at d: your frame 6.6713 ns; outside frame: out 4.72e+03 ns, back 0.00236 ns, total 4.72e+03 ns

--- Velocity addition from the Lorentz transformation
[PASS] u' = (u - v)/(1 - uv/c^2)
[PASS] u = c gives u' = c for any v
[PASS] 200 successive boosts of c/2 stay below c (exact arithmetic): 1 - w = 7.53e-96
```

```text
--- Force cannot bring a mass to c
  t (yr)   Newton v/c   true v/c
    0.5       0.516   0.4586559750
    1.0       1.032   0.7182527572
    2.0       2.065   0.8999874888
    5.0       5.161   0.9817441863
   10.0      10.323   0.9953407223
   20.0      20.646   0.9988290488
[PASS] v stays below c for all force and time tried
  Energy for a 70 kg person, (gamma-1) m c^2:
    1 - v/c = 0.1: gamma = 2.294e+00, energy = 8.142e+18 J = 1.36e-02 years of world energy use
    1 - v/c = 0.01: gamma = 7.089e+00, energy = 3.831e+19 J = 6.38e-02 years of world energy use
    1 - v/c = 1e-06: gamma = 7.071e+02, energy = 4.442e+21 J = 7.40e+00 years of world energy use
    1 - v/c = 1e-12: gamma = 7.071e+05, energy = 4.449e+24 J = 7.41e+03 years of world energy use
    1 - v/c = 1e-24: gamma = 7.071e+11, energy = 4.449e+30 J = 7.41e+09 years of world energy use
  Bertozzi-type check, electrons in MeV kinetic energy:
      0.5 MeV: v/c = 0.8629; Newton would say 1.40
      1.5 MeV: v/c = 0.9672; Newton would say 2.42
      4.5 MeV: v/c = 0.9948; Newton would say 4.20
     15.0 MeV: v/c = 0.9995; Newton would say 7.66
```

```text
--- What an outside observer sees (photon four-momentum)
[PASS] beta = 0.5: forward photon blue, returning photon red by the Doppler factor: outside energies 1.732 and 0.5774 of your 1.000
[PASS] beta = 0.99: forward photon blue, returning photon red by the Doppler factor: outside energies 14.107 and 0.0709 of your 1.000

--- Mirror fixed to the ground, you running at it
[PASS] beta = 0.5: photon still leads you and comes back: you see the image at 4.448 us (lab) when you are 333.33 m from the mirror; light blue-shifted x3.0
[PASS] beta = 0.99: photon still leads you and comes back: you see the image at 3.352 us (lab) when you are 5.03 m from the mirror; light blue-shifted x199.0
[PASS] beta = 0.9999: photon still leads you and comes back: you see the image at 3.336 us (lab) when you are 0.05 m from the mirror; light blue-shifted x19999.0
  550 nm reflected at beta = 0.99 arrives at 2.76 nm, in the soft X-ray band

--- No rest frame at v = c
  gamma = 1/sqrt(1 - v^2/c^2) grows as: 7.1e+00, 7.1e+02, 7.1e+04
  At v = c the transformation is undefined. Worse, postulate 2 says any inertial frame measures
  that light at c, so a frame where it stands still contradicts the theory that defines frames.
```

The mirror test checks that the formulas hang together. It's built from the Lorentz transformation, so it can't test the postulates. The experiments cited do that.

**Limits:** Flat spacetime with uniform motion. Accelerating to your speed would matter for you physically but not for a mirror check in steady motion. World energy use in the output is an order-of-magnitude figure.

**Sources:** [Alväger et al., Phys. Lett. 12, 260 (1964)](https://docente.unife.it/guido.zavattini/allegati/alvaeger-pl1964.pdf) (the 0.99975 figure is from general knowledge since the scan's digits were unclear); [Bertozzi, Am. J. Phys. 32, 551 (1964)](https://en.wikipedia.org/wiki/Tests_of_relativistic_energy_and_momentum), via Wikipedia.

### L12: Can we create conditions where we cross light speed?

**Status:** Derived + Code + Established + Speculative

**Short answer:** No known conditions let matter, energy or information pass light locally. Many things *appear* faster than light, such as expanding space, wave crests, jet images and entangled correlations, and all of them are consistent with relativity. A warp drive is allowed in the equations only with negative energy we don't know how to build. You're right that SR has a domain of validity, but where it stops applying (large regions of curved spacetime), general relativity keeps the same local limit.

**Starting from:** Special relativity holds in local inertial frames. Quantum theory is used for the entanglement part and general relativity's 3+1 constraint for the warp part.

**Derivation:**

*What SR's domain is.* SR describes physics in a region small enough that spacetime curvature doesn't matter. In general relativity, every small freely falling lab has SR inside it. So the rule is local: nothing outruns light *at the place it is*. Comparing far-apart speeds isn't even well defined in curved spacetime.

*Why the limit is stronger than a speed record.* If some signal moved at $u > c$ in all frames, then frames moving at $v$ with $uv > c^2$ would see it arrive before it left: $\Delta t' = \gamma\Delta t(1 - uv/c^2) < 0$. Two such signals allow a reply to arrive before the question, once $v > 2u/(1 + u^2/c^2)$. So faster-than-light signalling plus relativity gives causality paradoxes, and that makes it a much more serious change than breaking a record. The script runs this at $u = 2c$.

| Case | What goes faster than $c$ | Why relativity survives |
|---|---|---|
| Expansion of space | Distances between far galaxies grow faster than $c$ beyond about 14.5 Gly | Nothing passes light locally; see C02 in section 07 |
| Phase velocity | Crests in plasma or a waveguide, $v_p = \omega/k > c$ | Crests carry no signal. Envelope moves at $v_g = c^2/v_p < c$ |
| Anomalous dispersion | A pulse peak, with negative group velocity (Wang et al. 2000) | The peak is reshaped from the leading edge. The front still goes at $c$ |
| Cherenkov radiation | A particle beats *light in a medium*, $v > c/n$ | Particle stays below $c$. In water, for electrons above 0.264 MeV |
| Jet images | Blobs seem to cross the sky at up to about $6c$ in M87 | Light-travel effect, derived below |
| Entanglement | Correlations appear immediately | Bob's statistics can't depend on Alice's choice |
| Alcubierre warp | Bubble's effective speed | Requires negative energy density, derived below |
| Tachyons | Hypothetical particles with imaginary mass | None found. Causality problems above apply if they carry signals |

*Superluminal jets.* A blob moves at speed $\beta c$ at angle $\theta$ to our line of sight. In time $\Delta t$ it moves $\beta c\Delta t\sin\theta$ across the sky and $\beta c\Delta t\cos\theta$ toward us. Later light has less far to come, so the arrival interval shrinks to $\Delta t(1 - \beta\cos\theta)$:

$$
\beta_{\rm app} = \frac{\beta\sin\theta}{1 - \beta\cos\theta}.
$$

Setting the derivative to zero gives the maximum at $\cos\theta = \beta$, with $\beta_{\rm app} = \gamma\beta$. This can be any size while $\beta < 1$. An apparent $6c$ only needs $\beta \ge 0.986$. Rees predicted this in 1966, before it was observed.

*Entanglement.* Alice acts on her half of a shared state. Bob's statistics are given by his reduced state, $\rho_B = {\rm Tr}_A[(U\otimes I)\rho(U^\dagger\otimes I)] = {\rm Tr}_A\,\rho$, since the trace undoes $U$ by cyclicity. So no choice by Alice changes anything Bob can see. Only on comparing notes, by slower-than-light communication, do the correlations show. This no-signalling result is general (Eberhard & Ross 1989).

*Alcubierre metric.* It uses flat spatial slices and a moving shift $\beta^x = -v_sf(r)$. For flat slices, the Hamiltonian constraint of general relativity gives an Eulerian energy density $\rho = (K^2 - K_{ij}K^{ij})/16\pi G$ with $K_{ij}$ built from the shift. Carrying that through (in $c = 1$ units):

$$
\rho = -\frac{v_s^2}{32\pi G}\,\frac{y^2 + z^2}{r^2}\left(\frac{df}{dr}\right)^2 \le 0 .
$$

That's negative everywhere for any shape $f$. For a 100 m bubble with a 1 m wall at $v_s = c$, it totals about $-3\times10^{46}$ J, or 0.19 solar masses of negative energy, and thinner walls cost more. Later work finds that generic warp drives violate the null energy condition (Santiago, Schuster & Visser). Earlier analysis with quantum inequalities (Pfenning & Ford 1997, as described in the review cited below) made it worse. Nobody has a way to build it, so **Speculative**.

*Can we test the limit further?* Yes, and it's healthy to try. In 2011 OPERA reported neutrinos arriving early, and traced it to a faulty fibre connection and a timing drift, and later measurements agreed with $c$. Quantum gravity could conceivably alter the rule at extreme energy, but tests have found nothing (L10). So far, the honest statement is "never seen, and very hard to fit into known physics", and not "proved impossible for all time".

**Check:** [`code/verify_faster_than_light.py`](code/verify_faster_than_light.py) simulates each case.

```text
--- Apparent superluminal motion
[PASS] beta = 0.9, theta = 30 deg: simulated rate matches the formula: apparent 2.040 c
[PASS] beta = 0.99, theta = 8 deg: simulated rate matches the formula: apparent 7.017 c
[PASS] beta = 0.995, theta = 5 deg: simulated rate matches the formula: apparent 9.870 c
[PASS] beta = 0.9: maximum is gamma beta at cos(theta) = beta: max 2.065 c at 25.84 deg
[PASS] beta = 0.99: maximum is gamma beta at cos(theta) = beta: max 7.018 c at 8.11 deg
[PASS] beta = 0.999: maximum is gamma beta at cos(theta) = beta: max 22.344 c at 2.56 deg
  An apparent 6c (M87 jet) needs gamma beta >= 6: beta >= 0.9864, gamma >= 6.08

--- Phase, group and front velocity (c = 1)
[PASS] envelope moves at the group velocity k/omega < c: measured 0.89434 c, theory 0.89443 c
[PASS] crests move at the phase velocity omega/k > c: measured 1.1181 c, theory 1.1180 c
[PASS] v_phase * v_group = c^2: 0.9999
[PASS] disturbance from a compact source stays inside x = a + c t: outside 3.4e-12, inside 3.4e-01

--- Cherenkov threshold
  In water n = 1.33: light moves at 0.752 c, electron threshold kinetic energy 0.264 MeV
  Cone angle for beta -> 1: 41.2 deg. Particle still below vacuum c
[PASS] electron threshold in water near 0.26 MeV

--- Expansion
  H0 = 67.4 km/s/Mpc: distance receding at c = 14.5 Gly

--- Entanglement cannot signal
[PASS] Bob's state unchanged by any unitary on Alice's particle: worst 6.7e-16
[PASS] Bob's state unchanged by any measurement by Alice: worst 1.1e-15
[PASS] the same state still gives CHSH = 2 sqrt 2: S = 2.8284

--- Superluminal signal leads to a reply before the question
  Signals at u = 2.0 c: reply comes back early only if B moves faster than 2u/(1+u^2) = 0.800 c
  v/c    time the reply reaches A (yr; A asked at t = 0)
  0.500  +0.500
  0.800  +0.000
  0.900  -0.227
  0.990  -0.470
[PASS] at v = 0.9c the reply reaches A before A asked: t = -0.227 yr
[PASS] at v = 0.5c it doesn't, so the crossover is real
[PASS] crossover at v = 2u/(1+u^2) found by the event calculation

--- Alcubierre metric: sign of the energy density
[PASS] rho = -(v^2/32 pi G)(y^2+z^2)/r^2 f'(r)^2 (from K^2 - K_ij K^ij): sample value -5.1117e-05
[PASS] energy density is never positive at 200 random points: most negative -9.765e-03
  R = 100 m bubble, wall ~ 1 m, v = c: E = -3.36e+46 J, i.e. 0.19 M_sun of negative energy
[PASS] matches the thin-wall estimate -(c^4/G)(v/c)^2 R^2 sigma/36
  With a 1 cm wall: 9.4 M_sun. Thinner walls cost more
```

The apparent-speed simulation uses the same light-delay assumption as the derivation, and tests only the algebra and the maximum. The front test is real: a compact source under the Klein-Gordon equation stays inside its light cone, even though crests move faster than $c$. The entanglement test only covers a qubit pair.

**Limits:** Everything assumes known physics. Hubble distance uses $H_0 = 67.4$ km/s/Mpc. Tachyon and warp statements describe hypothetical physics.

**Sources:** [Rees, Nature 211, 468 (1966); Mirabel & Rodríguez, Nature 371, 46 (1994)](https://en.wikipedia.org/wiki/Superluminal_motion); [Biretta, Sparks & Macchetto, ApJ 520, 621 (1999)](https://www.stsci.edu/science/preprints/prep1335/prep1335.html); [Wang, Kuzmich & Dogariu, Nature 406, 277 (2000)](https://ideas.repec.org/a/nat/nature/v406y2000i6793d10.1038_35018520.html); [Eberhard & Ross, Found. Phys. Lett. 2, 127 (1989)](https://link.springer.com/article/10.1007/BF00696109); [Alcubierre, CQG 11, L73 (1994)](https://arxiv.org/abs/gr-qc/0009013); [Santiago, Schuster & Visser, arXiv:2105.03079](https://arxiv.org/abs/2105.03079); [Lobo & Visser, arXiv:gr-qc/0406083](https://arxiv.org/abs/gr-qc/0406083); [OPERA timing errors, Physics World](https://physicsworld.com/a/doubts-grow-over-superluminal-neutrino-result/) and [arXiv:1203.2585 added note](https://arxiv.org/abs/1203.2585).

### L13: A mug in the rain on a trip around the Earth at light speed

**Status:** Derived + Code + Estimate

**Short answer:** Your setup needs a fix, because a mug can't reach $c$ (L11). So take $v$ just below $c$. The answer then depends on the mug's direction. Facing forward, it collects about 86 litres on one lap (more than it can hold), however close to $c$ it goes, even if its own clock shows under a nanosecond. Facing up, it catches less than one raindrop, tending to zero as $v \to c$. At those speeds the air and rain would hit with nuclear-scale energy, so this is only a thought experiment.

**Starting from:** Special relativity. Rain with volume fraction $\phi = R/u$ where the rainfall rate is $R = 10$ mm/h and drop speed $u = 6.5$ m/s (2 mm drops, general knowledge of Gunn-Kinzer's data), mouth area $A = 50$ cm², lap length $L = 40{,}075$ km, and straight-line kinematics. All are assumed, including rain and sea-level air around the entire lap.

**Derivation:**

Counting drops is counting events, and events don't depend on frame. So any correct frame must give the same answer.

*Earth's frame.* The trip takes $T = L/v \approx 0.134$ s. That's $1.3\times10^8$ ns, not "no time", and rain falls only 0.87 m during it. So the drops are almost frozen in place, and the mug acts like a scoop.

- **Mouth facing forward.** The mouth sweeps a tube of volume $AL = 2\times10^5$ m³, whatever $v$ is. Transverse sizes don't contract. Water caught:

$$
V_w = \phi AL = (4.3\times10^{-7})(2.0\times10^5\ \text{m}^3) = 86\ \text{L} .
$$

- **Mouth facing up.** Rain must come down through the opening. In Earth's frame the opening is contracted to width $D/\gamma$ along the direction of motion:

$$
V_w = R\frac{A}{\gamma}T = \frac{1.9\ \mu\text{L}}{\gamma}.
$$

That's less than half of one 4 µL drop even at low $\gamma$.

*Mug's frame.* The mug is at rest and its clock runs for only $T' = T/\gamma$. Rain arrives with three changes:

- Lengths along the motion are contracted, so the drops are denser: $n' = \gamma n$.
- Transverse velocity is reduced by time dilation: $u_y' = u/\gamma$.
- The rain comes almost head-on (aberration), at angle $\tan\alpha' = u/(\gamma v)$ from the direction of travel. For $\gamma = 22.4$, that's $10^{-9}$ rad.

Forward count: $n' \times A \times (\text{path } L/\gamma) = nAL$. The two factors of $\gamma$ cancel.

Upward count: flux through the top is $n'u_y' = nu$, for time $T/\gamma$, through full area $A$, giving $nuAT/\gamma$. That matches.

These are one result: the four-current $J^\mu$ through a surface element is invariant.

*Your "no nanosecond" idea.* Your intuition is that no time passes for the traveller, so nothing could happen. For the mug's clock to read under 1 ns, $\gamma > 1.3\times10^8$. Even then the forward catch is still 86 L. Short time on your clock doesn't mean few events, because the path is compressed and the rain is dense in your frame. This is a good thing to remember about relativity: elapsed proper time doesn't limit what you encounter.

*Now the catastrophe.* The forward mug meets $\rho_{\rm air}AL = 2.4\times10^5$ kg of air. Kinetic energy $(\gamma-1)mc^2$ per drop at $\gamma = 22.4$ is $8\times10^{12}$ J, about 2 kilotons of TNT per drop, and there'd be $2\times10^7$ drops. In the mug's frame, each metre of air brings $1.2\times10^{16}$ J at that speed, and vaporising the mug takes only about $3\times10^6$ J (rough estimate). So even if only a tiny fraction is deposited, the mug is gone within the first moments. If you ignore that and let a 0.3 kg lump keep what it hits, it slows to $c/2$ after sweeping $\gamma m$ of air, about 1.1 km at $\gamma = 22.4$, and it couldn't complete the lap unless $\gamma \gtrsim 8\times10^5$. The whole column would carry $5\times10^{23}$ J in the mug's frame, but only a mug somehow held at speed would meet all of it. Following Earth's curve at this speed needs $v^2/R \approx 10^9\ g$ in Earth's frame, and $\gamma^2$ times that in the mug's.

So the water count is a clean calculation about geometry, not a prediction of what would happen.

**Check:** [`code/verify_rain_mug.py`](code/verify_rain_mug.py) generates random drops and counts catches. The Earth-frame count uses density $n$ with a contracted mouth. The ship-frame count separately draws fresh drops with density $\gamma n$ and boosted velocities, with a mouth at rest.

```text
--- Monte Carlo counts in both frames (test parameters, c = 1)
[PASS] bucket facing forward, gamma = 1.2: count = n A L in both frames: Earth 7692, ship 7835, n A L = 7854
[PASS] bucket facing forward, gamma = 7.1: count = n A L in both frames: Earth 7940, ship 7954, n A L = 7854
[PASS] bucket facing forward, gamma = 70.7: count = n A L in both frames: Earth 7969, ship 7811, n A L = 7854
[PASS] mug facing up, gamma = 1.25: count = n u A T / gamma in both frames: Earth 12620, ship 12640, formula 12566
[PASS] mug facing up, gamma = 7.09: count = n u A T / gamma in both frames: Earth 1325, ship 1360, formula 1343

--- Real rain, real mug
  Rain: 10 mm/h, volume fraction of water in air 4.27e-07, 102 drops/m^3
  Mug mouth 50.3 cm^2, trip length 40075 km
  Trip time in Earth frame (v -> c): 0.1337 s = 1.337e+08 ns
[PASS] trip time is about 0.13 s
  On board, 'less than a nanosecond' needs gamma > 1.34e+08
  Bucket facing forward: swept volume 2.014e+05 m^3, water 86.1 L, 2.06e+07 drops
  Mug facing up, no contraction (gamma = 1 formula): 1.87 uL = 0.45 drops
  gamma     onboard time      facing forward    facing up
  1.2       1.114e-01 s       86.1 L          1.555e+00 uL
  22.4      5.968e-03 s       86.1 L          8.332e-02 uL
  1e+03     1.337e-04 s       86.1 L          1.866e-03 uL
  1e+06     1.337e-07 s       86.1 L          1.866e-06 uL
  1e+09     1.337e-10 s       86.1 L          1.866e-09 uL
  Standing still in the same rain for 0.134 s catches 1.87 uL, the gamma = 1 value of the facing-up formula.
[PASS] forward catch exceeds the facing-up catch by more than 1e7: ratio 4.6e+07

--- Why this can't actually happen: energy and stopping
  Column along the path: air 2.42e+05 kg, water 86 kg, mug 0.3 kg
  This assumes sea-level air and rain along the entire path; real rain covers far less of it.
  gamma    one drop KE         air column KE
  22.4     8.06e+12 J = 1.9 kt  4.65e+23 J
  1000     3.76e+14 J = 90 kt  2.17e+25 J
[PASS] one drop at gamma 22 carries about a kiloton of energy: 1.9 kt TNT
  In the mug's frame each metre of air brings 1.2e+16 J at gamma 22.4; vaporising it needs about 3e6 J
  If the mug simply keeps the matter it hits, v/c ~ gamma m/(gamma m + swept mass):
  gamma 22.4    : speed halves after 1.11 km of air; trip not completed
  gamma 1e+03   : speed halves after 49.7 km of air; trip not completed
  gamma 1e+06   : speed halves after 4.97e+04 km of air; mug could finish the distance in principle
  Needed gamma to still be above c/2 after the trip: 8.1e+05
  Following the curved path takes Earth-frame acceleration v^2/R = 1.41e+10 m/s^2 = 1.4e+09 g,
  and proper acceleration gamma^2 times that, e.g. 7.2e+11 g at gamma 22.4
```

The counts agree within Poisson noise (5σ band). The real-rain lines are direct evaluations of the verified formulas. Test parameters differ from the real values so that counts are large enough to compare.

**Limits:** It treats drops as points, the path as straight for the sweep (the catch depends only on path length), and no interaction with air. Real rain doesn't cover a full lap. Mug capacity, and whether liquid water could stay in it, are ignored, and the 86 L assumes everything sticks.

### L14: A world with only cesium pendulums

**Status:** Derived + Code + Conceptual + Established

**Short answer:** That world could keep time, but its clocks would run 4.2 s/day slower per kelvin of warming. If every clock is the same type at the same temperature, they'd agree perfectly and no comparison among them could detect the drift. They'd find it only by comparing with other processes, such as the day or the planets. Then everything else would seem to speed up on hot days by the same fraction, and that pattern would point at the clock. This is how we actually choose time standards.

**Starting from:** The pendulum equation $\ddot\theta = -(g/L)\sin\theta$. Linear thermal expansion $L = L_0(1 + \alpha\Delta T)$, with $\alpha_{\rm Cs} = 97\times10^{-6}$ /K at 25 °C (some tables give 92.6 at 20 °C). Time scales are tested only by comparing clocks.

**Derivation:**

*Rate.* Small swings give $T = 2\pi\sqrt{L/g}$, so

$$
\frac{\Delta T}{T} = \frac12\alpha\Delta T \quad\Rightarrow\quad 48.5\ \text{ppm per K} = 4.19\ \text{s/day per K}.
$$

The same scaling holds for any rigid pendulum with all its lengths expanding together. Compare: steel 0.52, invar 0.05 s/day/K. A 10 K swing costs 42 s/day. Temperature isn't even the biggest issue: a 5° swing adds $\theta_0^2/16$, or 41 s/day, and moving from equator to pole changes $g$ for 229 s/day, equal to 55 K of cesium.

*Practical issue.* Cesium melts at 28.4 °C and reacts violently with air and water (general knowledge). So your rod melts on a warm day. It's an excellent choice for showing the effect and a terrible one for a clock.

*Why it expands so much.* The Grüneisen relation (L15) gives $\alpha = \gamma_GC/(3BV)$. Cesium's atoms are large and it is very soft: 1.6 GPa bulk modulus and 69 cm³/mol. That makes $3BV$ nine times smaller than copper's. With the measured $\alpha$, the implied $\gamma_G \approx 1.3$ is an ordinary value. So softness explains it, and the atoms aren't unusually "anharmonic".

*The common-mode theorem.* Suppose every clock obeys $d\tau_i = r(T)\,dt$, with the same $r$ for all and the same $T$. Then $\tau_i - \tau_j$ never changes. Clocks are only measured against each other, so a shared drift has no observable effect within that family. "The clock is drifting" has no meaning unless there's something else to compare it with.

*What they could still find.* Different temperatures in different rooms would make pendulums disagree, so they'd discover the temperature effect as soon as they compare a warm and a cold clock. Then they might define time at a reference temperature, as our clockmakers did with compensated pendulums.

But the *natural* clocks expose it too. In pendulum time, the Earth's spin, the Moon's orbit and Mercury's orbit would all appear to speed up on hot days and by exactly the same fraction. They would have two descriptions:

1. Time is what pendulums say, and the Sun, Moon and Mercury all get a mysterious temperature-dependent force, in exact proportion.
2. Pendulums drift with temperature, and the laws stay simple.

Neither is logically forced. Poincaré made this kind of point in 1898: we choose a definition of time that keeps the laws of nature simple. Option 2 wins by far, but preferring it is a judgement and not a proof, which is why this part is Conceptual.

*How we did it.* Exactly this has happened to us. People used Earth's rotation as the clock, until comparisons with the Moon and planets showed that the Earth was the irregular one (de Sitter 1927; Spencer Jones 1939), and ephemeris time was adopted. Then atomic clocks beat that (Essen & Parry 1955, a cesium resonator). And today, comparing Al⁺, Yb and Sr optical clocks to $6$–$8\times10^{-18}$ (BACON 2021) is the same strategy: different physical processes, compared with each other, because any one kind could drift. Since different atoms respond differently to a change in fundamental constants, those comparisons also limit whether constants change.

There's an irony too. Cesium the pendulum rod is about the worst choice for time, and cesium the atom defines the second, since its hyperfine frequency has nothing to do with a metal's expansion.

**Check:** [`code/verify_pendulum_world.py`](code/verify_pendulum_world.py) integrates the full pendulum, then runs a toy world.

```text
--- Pendulum rate
  Length for a 2 s period: L = 0.9936 m
[PASS] small-swing period = 2 pi sqrt(L/g): 2.00000013 s
  Metal rod: temperature coefficient, seconds/day per K
    cesium  alpha =  97.0 ppm/K -> 4.190 s/day/K
    brass   alpha =  19.0 ppm/K -> 0.821 s/day/K
    steel   alpha =  12.0 ppm/K -> 0.518 s/day/K
    invar   alpha =   1.2 ppm/K -> 0.052 s/day/K
[PASS] simulated dT/T per kelvin = alpha/2: 48.50 ppm/K against 48.50 ppm/K
  So a cesium pendulum changes rate by 4.19 s/day per kelvin; a 10 K swing is 42 s/day
  A compound (uniform rod) pendulum has T = 2 pi sqrt(2L/(3g)): same scale law, same alpha/2
[PASS] amplitude matters: T grows by theta0^2/16: 41.1 s/day at 5 degrees
  Moving between equator and pole: 229 s/day, equal to 55 K of cesium

--- Toy world: only ratios can be measured
[PASS] identical clocks at the same temperature never disagree, despite drift: yet each drifts 5393 s against model time
[PASS] a room 8 K warmer shows up immediately: -33.52 s/day
  Apparent rate in pendulum time against daily temperature, 3 ppm measurement noise on each process:
    Earth spin     fractional slope  48.501 ppm/K
    Moon orbit     fractional slope  48.509 ppm/K
    Mercury orbit  fractional slope  48.492 ppm/K
[PASS] all processes appear to speed up on hot days by the same fraction, alpha/2: alpha/2 = 48.50 ppm/K
  Same fractional slope in unrelated processes is the signature: blame the clock, not the Sun, Moon and Mercury
```

The rate tests are real integrations. The "identical clocks never disagree" test is true by construction in the toy world, and illustrates the theorem without testing it. "Model time" exists in the simulation only: inhabitants never have it. The three-process regression recovers the injected effect from noisy data, which shows that the signature is detectable, not that any real world behaves this way.

**Limits:** Simple model: only temperature and noise. Real pendulums also feel air pressure and support flexing.

**Sources:** [Caesium properties (Wikipedia)](https://en.wikipedia.org/wiki/Caesium); [Poincaré 1898, as discussed in arXiv:gr-qc/0610005](https://arxiv.org/abs/gr-qc/0610005); [Ephemeris time (Wikipedia)](https://en.wikipedia.org/wiki/Ephemeris_time), citing Spencer Jones, MNRAS 99, 541 (1939); [Louis Essen (Wikipedia)](https://en.wikipedia.org/wiki/Louis_Essen) for Essen & Parry, Nature 176, 280 (1955); [BACON collaboration, Nature 591, 564 (2021)](https://arxiv.org/abs/2005.14694).

### L15: The most temperature-sensitive pendulum in the coldest gas ever made

**Status:** Derived + Code + Estimate

**Short answer:** Two surprises. Near absolute zero, thermal expansion disappears, so even "the most temperature-sensitive metal" stops being sensitive: the coefficient goes to zero as $T \to 0$, as the third law requires. And an ultracold gas can't cool a metal at all: the whole cloud holds less thermal energy than a single room-temperature thermal photon. Your pendulum would stay at the temperature of its surroundings and would destroy the gas instantly.

**Starting from:** The third law of thermodynamics, and the Debye model for heat capacity. Grüneisen's relation $\alpha_L = \gamma_GC_V/(3BV_m)$. Values for 100 g of copper as a stand-in for a metal with known low-temperature data. A BEC of $10^6$ Rb atoms at $n = 10^{14}$ cm⁻³ and 100 nK, typical order-of-magnitude values.

**Derivation:**

*Third law: why expansion must vanish.* A Maxwell relation says

$$
\left(\frac{\partial V}{\partial T}\right)_P = -\left(\frac{\partial S}{\partial P}\right)_T.
$$

As $T \to 0$, the third law says entropy approaches a constant independent of pressure, so the right side vanishes, and hence so does thermal expansion. This holds for all materials. The most sensitive metal at room temperature isn't exempt, though softer metals have larger coefficients on the way down.

*How it fades.* Expansion follows heat capacity: $\alpha_L = \gamma_GC_V/(3BV_m)$. With a Debye phonon term ($\propto T^3$) and an electron term ($\propto T$):

| $T$ | $\alpha_L$ (1/K) | clock effect |
|---|---|---|
| 293 K | $1.5\times10^{-5}$ | 0.67 s/day/K |
| 4 K | $3.9\times10^{-9}$ | $1.7\times10^{-4}$ s/day/K |
| 1 mK | $4.6\times10^{-13}$ | $2\times10^{-8}$ s/day/K |
| 1 nK | $4.6\times10^{-19}$ | $2\times10^{-14}$ s/day/K |

As a check on the model: integrating it from 293 K to 4 K gives 0.304% contraction, compared with a published 0.324% for copper. Over 99% of the contraction happens above 30 K. Between 1 nK and absolute zero, the length left to change is about $2\times10^{-28}$ of the length, which has no meaning at all. Cooling a copper pendulum from room temperature to absolute zero would speed it up by 131 s/day, nearly all of it above 30 K.

So at the coldest temperatures, the most sensitive pendulum in the world becomes one of the *least* temperature-sensitive objects. That's useful: low temperature is how precision instruments avoid thermal drift. Cryogenic silicon cavities (general knowledge) use this principle.

*About "not going into quantum."* You can't keep this classical. Classically, $C = 3R$ all the way down, so $\alpha$ would stay at $1.6\times10^{-5}$ /K forever and the metal would contract smoothly to zero temperature. The freeze-out **is** the quantum effect, because lattice vibrations are quantised and run out of energy. The classical picture only survives for the pendulum's swing. Its energy quantum is $\hbar\omega$, with $\hbar\omega/k_B = 2.4\times10^{-11}$ K for a 1 m pendulum. The swing would hold about 41 quanta if it could be held at 1 nK, and about 1 at 38 pK. So the swing mode itself approaches its ground state at these temperatures, though no real pendulum could get close.

*Why the gas can't cool it.* A BEC isn't a cold bath you can dip things in.

- **Mass and density.** $10^6$ atoms have a mass of $1.4\times10^{-19}$ kg, and even the density of $1.4\times10^{-5}$ kg/m³ is a hundred-thousandth that of air. A 100 g bob outweighs the cloud by $10^{18}$.
- **Heat to remove.** The bob holds about 7.3 kJ of thermal energy above 0 K, with a heat capacity of 37 J/K at room temperature. That capacity collapses when cold, to about $10^{-12}$ J/K at 1 nK from the electrons. But that is still about $5\times10^4$ times the cloud's $2\times10^{-17}$ J/K, so the cloud is a tiny reservoir even at the very end.
- **What the cloud could absorb.** If every atom warmed to 293 K it would take $6\times10^{-15}$ J, which would cool the bob by $1.6\times10^{-16}$ K. The ratio is $10^{18}$.
- **Worse: sticking.** An atom landing on metal releases about 0.5 eV (order of magnitude), heating the metal *more* than the atoms could remove. So the gas likely warms it.
- **Thermal energy.** The whole cloud at 100 nK holds $2\times10^{-24}$ J, while one 293 K thermal photon has $4\times10^{-21}$ J. A single photon from the room carries about 2000 times the cloud's thermal energy.
- **Speed.** At that density the bob would catch all of the atoms in a few nanoseconds.

*Steady state.* A real BEC sits in a room-temperature vacuum chamber, kept cold only because almost nothing touches it. Feed in atoms at an optimistic $10^{10}$ per second and they can remove at most $6\times10^{-11}$ W. Meanwhile the walls radiate $2\times10^{-2}$ W into even a polished bob, and a thin wire adds another $10^{-4}$ W. That's a shortfall of over $10^8$. The bob stays at room temperature.

*What actually can cool metal.* A dilution refrigerator holds the large CUORE detector near 10 mK. That's about $10^5$ times warmer than a BEC, but it's a real thermal bath. To go lower in bulk you need nuclear demagnetisation, and that still doesn't bring the lattice and electrons of a whole metal object to nanokelvin (general knowledge).

The 38 pK record (Deppner et al. 2021) is an effective kinetic temperature of a freely expanding cloud in a 110 m drop tower, not a material you could hold anything in. See A04 in section 05 for that idea in another setting.

**Check:** [`code/verify_pendulum_world.py`](code/verify_pendulum_world.py) integrates the Debye model numerically, tests the expansion model against the published contraction and computes the energy budget.

```text
--- Debye model and Grueneisen relation
[PASS] Debye C at 1000 K approaches Dulong-Petit 3R
[PASS] Debye C at 5 K follows the T^3 law
[PASS] Grueneisen integral reproduces copper's 293 K to 4 K contraction within 10%: model 0.304%, published 0.324%
  alpha(293 K) model = 15.5 ppm/K (measured about 16.5 ppm/K, general knowledge)
  Why cesium expands so much: alpha = gamma C/(3 B V). Bulk modulus is the main factor
    Cu: 3BV = 2.98e+06 J/mol; Cs: 3BV = 3.31e+05 J/mol, 9.0x softer
    Cs with alpha = 97 ppm/K and C = 3R implies a Grueneisen parameter of 1.29, an ordinary value
[PASS] cesium's implied Grueneisen parameter is ordinary (0.5 to 3)

--- L15: cold pendulum
  T (K)    alpha_L (1/K)   rate change per K (s/day/K)   contraction left between T and 0
  293      1.547e-05       6.685e-01                 3.042e-03
  77       7.227e-06       3.122e-01                 1.895e-04
  4        3.858e-09       1.667e-04                 5.687e-09
  1        4.890e-10       2.113e-05                 2.287e-10
  0.001    4.573e-13       1.976e-08                 2.287e-16
  1e-06    4.573e-16       1.976e-11                 2.287e-22
  1e-09    4.573e-19       1.976e-14                 2.287e-28
[PASS] alpha falls toward zero (third law) from 293 K to 1 nK: alpha(1 nK) = 4.6e-19 /K
  Total rate gain for a copper pendulum 293 K -> 0 K: 131 s/day; almost all of it occurs above 30 K
[PASS] almost all contraction occurs above 30 K: 99.78%
  Classical (Dulong-Petit down to 0) would keep alpha = 1.64e-05 /K at all T,
  which would give a contraction of 0.48% and no freeze-out: the freeze-out IS the quantum effect

--- L15: energy budget for a 100 g copper bob
  Bob radius 1.39 cm, area 24 cm^2, heat capacity 37.0 J/K, heat to remove 7.27e+03 J
  Ultracold cloud: 1e+06 Rb atoms, total mass 1.4e-19 kg, density 1.4e-05 kg/m^3 (air 1.2 kg/m^3)
  Thermal energy of the whole cloud at 100 nK: 2.1e-24 J; one room-temperature thermal photon: 4.0e-21 J
[PASS] one 293 K thermal photon carries more energy than the entire cold cloud
  If every atom warms to 293 K it absorbs 6.1e-15 J, cooling the bob by 1.6e-16 K
[PASS] cloud's maximum heat sink is 1e17+ times too small: ratio 1.2e+18
  Atoms that stick release roughly 0.5 eV each (order of magnitude): 8.0e-14 J, more than they could remove
[PASS] adsorption heat exceeds the most the atoms could carry off
  Bob heat capacity at 0.001 K: 1.09e-06 J/K; cloud 2.1e-17 J/K at most
  Bob heat capacity at 1e-06 K: 1.09e-09 J/K; cloud 2.1e-17 J/K at most
  Bob heat capacity at 1e-09 K: 1.09e-12 J/K; cloud 2.1e-17 J/K at most
[PASS] even at 1 nK the bob's heat capacity exceeds the cloud's
  Collision rate on the bob at that density: 3.0e+14 atoms/s, cloud lost in 3.4 ns if it enveloped the bob
  Sustained cooling with fresh atoms at 1e+10/s: at most 6.1e-11 W
  Heat coming in: radiation from 293 K walls 2.02e-02 W, 0.2 mm steel wire to room temperature 1.38e-04 W
[PASS] heat load beats cold-atom cooling by more than 1e8: ratio 3.3e+08

--- L15: when does the swing itself stop being classical?
  T = 293      K: mean quanta in the swing mode 1.22e+13, rms position 6.42e-11 m
  T = 0.001    K: mean quanta in the swing mode 4.18e+07, rms position 1.19e-13 m
  T = 1e-09    K: mean quanta in the swing mode 41.3, rms position 1.19e-16 m
  T = 3.8e-11  K: mean quanta in the swing mode 1.14, rms position 2.35e-17 m
  hbar omega / k_B = 2.39e-11 K. Below this, the swing mode is in its quantum ground state
  Zero-point rms = 1.30e-17 m
[PASS] the swing mode nears its ground state only near 20 pK
```

The test against copper's published contraction is the one real test of the model here. The rest are estimates that use round-number inputs. Margins of $10^8$ or $10^{18}$ mean those inputs could be off by a lot without changing the conclusion.

**Limits:** Copper's Debye temperature (343 K), electronic coefficient and bulk modulus are general-knowledge values. The Grüneisen parameter is treated as constant and the same for phonons and electrons, which is a rough assumption at low temperature: the sub-kelvin values are order-of-magnitude only. Thermal contact (Kapitza) resistance would make cooling harder still. I didn't compute cesium's low-temperature curve, and only the trend to zero is general.

**Sources:** [Copper contraction 0.324%, arXiv:2307.07496](https://arxiv.org/abs/2307.07496); [Grüneisen parameter of copper 1.96](https://edusj.uomosul.edu.iq/index.php/edusj/article/view/41813); [Anderson et al., Science 269, 198 (1995)](https://scilit.com/publications/01bb6dc362f9e8af027ad006e196eec7); [Deppner et al., PRL 127, 100401 (2021)](https://physics.aps.org/articles/v14/119); [CUORE, Nature 2022](https://arxiv.org/abs/2104.06906).

### L16: Do eclipses matter to any experiment?

**Status:** Derived + Code + Established + Open

**Short answer:** Yes. The 1919 eclipse tested Einstein's prediction that the Sun bends starlight by 1.75 arcseconds, twice the Newtonian value. Eclipses also revealed helium and the million-degree corona, and they act as a natural switch-off of sunlight for studying the ionosphere. Claims that eclipses change pendulums or gravity itself (the Allais effect) have not been reliably reproduced, and a plain tidal effect is far larger.

**Starting from:** Newton's gravity for a particle moving at $c$, and the Schwarzschild solution of general relativity. An eclipse's use is mostly practical: the Moon blocks the Sun so that dim things near it can be seen.

**Derivation:**

*Light bending.* Stars near the Sun can be photographed only during totality. So you compare their positions with a night-time plate of the same field.

Use the orbit equation with $u = 1/r$ and $\phi$ the polar angle, and $h = r^2\dot\phi$. Newton gives $u'' + u = GM/h^2$. Take a particle passing the Sun at speed $c$ with impact parameter $b$, so $h = cb$. The first-order solution is $u = \sin\phi/b + GM/(c^2b^2)$. It reaches $u = 0$ at $\phi = -GM/(c^2b)$ and at $\pi + GM/(c^2b)$, so the total bend is the Newtonian result

$$
\delta_N = \frac{2GM}{c^2b} = 0.875''.
$$

(Free-fall acceleration doesn't depend on the mass, so the answer doesn't need light to have any. In 1919 Dyson labelled Einstein's 1911 calculation of this size the "Newtonian" value, as Kennefick notes.)

For null rays, general relativity gives $u'' + u = 3GMu^2/c^2$. Start from the straight line $u_0 = \sin\phi/b$ and solve to first order:

$$
u = \frac{\sin\phi}{b} + \frac{GM}{c^2b^2}(1 + \cos^2\phi).
$$

At large distance, $u \to 0$ at small $\phi = -\epsilon$ where $\epsilon/b = 2GM/c^2b^2$, so the total bend is $2\epsilon$:

$$
\delta_{GR} = \frac{4GM}{c^2b} = 1.75''.
$$

*Where the doubling comes from.* In the weak field, $ds^2 = -(1-2U)c^2dt^2 + (1+2\gamma U)d\vec x^2$ with $U = GM/rc^2$. Light behaves as if in a medium of index $n = 1 + (1+\gamma)U$. Time curvature gives the "1" and space curvature gives $\gamma$, with $\gamma = 1$ in general relativity. With $\gamma = 0$ you get the Newtonian half. So the eclipse measured space curvature, which Newton's theory doesn't have.

*The 1919 data and reassessment.*

| Instrument | 1919 reduction | 1979 Greenwich re-reduction |
|---|---|---|
| Sobral 4-inch | $1.98'' \pm 0.18''$ | $1.90'' \pm 0.11''$ |
| Sobral astrographic | $0.93''$ (no error stated) | $1.55'' \pm 0.34''$ |
| Principe | $1.61'' \pm 0.30''$ | not in that table |

Critics later claimed Eddington favoured Einstein by discarding the astrographic plates that favoured the Newtonian value. Kennefick's reassessment goes the other way. Dyson made that decision, he wasn't biased toward Einstein, and the problem with those plates was a separate scale change. The 1979 reanalysis (Harvey, Observatory 99, 195) moved the discarded plates to 1.55″, close to the alternative 1.52″ that Dyson and Crommelin had reported. So the historical choice is defensible, even though the data alone in 1919 weren't as strong as newspapers suggested. Modern radio measurements have fixed it: $\gamma = 0.99983 \pm 0.00045$ from VLBI (Shapiro et al. 2004) and $\gamma - 1 = (2.1 \pm 2.3)\times10^{-5}$ from Cassini (Bertotti et al. 2003, by time delay rather than bending).

*The Sun.*

- **Helium.** Janssen, at the 18 August 1868 eclipse in Guntur, saw a yellow line that matched no known element. Lockyer also identified it, and named the element helium. On Earth it was isolated in 1895.
- **Corona.** A green line seen in the 1869 eclipse was blamed on a new element, "coronium". Grotrian (1939) and Edlén (1943) showed it comes from highly ionised iron (Fe XIV), which requires a million-degree corona, hotter than the surface below. Why the corona is that hot is still an active research question.

*The ionosphere.* Sunlight maintains the ionised upper atmosphere. An eclipse switches it off along a known track. In 2017, GPS receivers saw total electron content fall by more than 60% in some regions, exceeding the initial model prediction of 35% (Coster et al. 2017).

*Pendulum and gravity anomalies.* Allais reported changes in a paraconical pendulum during the 1954 and 1959 eclipses. Saxl & Allen (1971) reported a large torsion-pendulum change, and Wang et al. (2000) reported a gravimeter dip of $(7.0 \pm 2.7)\times10^{-8}$ m/s² and suggested gravity shielding. Then replications:

- Kuusela (1991), torsion pendulum: any relative period change below $4.3\times10^{-6}$ at 90% confidence.
- Luo Jun et al. (1991), torsion pendulum: below $5.2\times10^{-5}$ at 90% confidence.
- Salva (2011), automated Foucault pendulum: no effect within error, and it would have detected 0.3°/hour.

Fair disclosure: the picture isn't perfectly clean. Other groups have reported positive signals, and Duif (2004) argued that the proposed mundane explanations don't quantitatively explain every old report. Van Flandern & Yang (2003) proposed an air-density explanation for the gravimeter result, which Duif disputes. So what caused each old anomaly is partly **Open**, while the effect itself isn't accepted, because it doesn't reproduce.

Two scale checks on it:

- The ordinary tide from Sun and Moon at such a site is about 160 µGal, varying by about 120 µGal in the 6 h around an eclipse. Wang's 7 µGal, at $2.6\sigma$, is about 6% of that. It has to be pulled out from a large changing signal, plus air and ground effects.
- If the Moon screened the Sun's pull, 7 µGal needs a screening fraction of $1.2\times10^{-5}$. Lunar laser ranging limits Majorana-type shielding to $h < 10^{-22}$ m²/kg (Eckhardt 1990), which for the Moon's column is about $10^{-12}$. Simple shielding is too weak by a factor of about $10^7$. This rules out only that specific model of shielding.

**Check:** [`code/verify_eclipse.py`](code/verify_eclipse.py) gets the bending three ways: Newtonian orbit, null geodesic, and tracing a ray through the refractive index. The first two use the orbit equation numerically; the third uses a different equation and has to agree.

```text
GM/c^2 for the Sun = 1476.6 m, Sun radius = 6.957e+08 m, ratio 2.123e-06

--- Three ways of bending a ray that grazes the Sun
  Newtonian particle at c: 0.8756 arcsec
  Schwarzschild null geodesic: 1.7512 arcsec
[PASS] Newton gives 2GM/(c^2 b) = 0.875 arcsec
[PASS] GR gives 4GM/(c^2 b) = 1.75 arcsec
[PASS] GR is twice Newton: ratio 2.00000
  Ray through an index with gamma = 0 (time curvature only): 0.8756 arcsec
  Ray through an index with gamma = 1 (time + space curvature): 1.7512 arcsec
[PASS] index route with gamma = 1 matches the geodesic
[PASS] gamma = 0 matches the Newtonian half: space curvature provides the other half

--- 1919 and the 1979 reanalysis (Kennefick 2009 table, Dyson et al. 1920)
  measurement                                value    sigma from Einstein   sigma from Newton
  Sobral 4-inch, 1919 reduction              1.98+-0.18  +1.27                +6.14
  Sobral 4-inch, 1979 re-reduction           1.90+-0.11  +1.35                +9.31
  Principe, 1919 reduction                   1.61+-0.30  -0.47                +2.45
  Sobral astrographic, 1979 re-reduction     1.55+-0.34  -0.59                +1.98
  Sobral astrographic, 1919 reduction: 0.93, no error stated; dropped by the team for a suspected
  scale change. Its reanalysis moves to 1.55, as the table above shows.
[PASS] Sobral 4-inch excludes Newton's half value, both reductions, over 5 sigma
[PASS] Principe alone sits within 3 sigma of both predictions
  Modern radio VLBI: gamma = 0.99983 (Shapiro et al. 2004), i.e. 1.75105 arcsec at the limb
  Cassini: gamma - 1 = (2.1 +- 2.3)e-5 (Bertotti et al. 2003)

--- Newtonian tide around an eclipse against the reported anomaly
  Sun + Moon vertical tide at a site under both at mid-eclipse: 163.5 uGal, range over +-3 h 119.3 uGal (1 uGal = 1e-8 m/s^2)
  Wang et al. 2000 anomaly: (7.0 +- 2.7) uGal, 2.6 sigma, 5.9% of that tide's 6 h range
[PASS] reported anomaly is a small fraction of the tide that must be removed

--- Moon as a gravity shield?
  Sun's pull at Earth 5.930e-03 m/s^2; a 7 uGal effect needs shielding fraction 1.18e-05
  Lunar laser ranging limit h < 1e-22 m^2/kg with Moon column 1.16e+10 kg/m^2: fraction < 1.16e-12
  So simple Majorana-type shielding would be too weak by a factor 1.0e+07
[PASS] shielding of that size is excluded by over 1e5 in that specific model
```

The significance values treat quoted errors as Gaussian and statistical, which is optimistic for 1919 plates. The tide uses a simplified site with the Sun and Moon both in the equatorial plane, and so it's only a scale estimate, not a prediction for a specific eclipse.

**Limits:** The first-order bending formula is enough for the Sun, where $GM/c^2b \approx 2\times10^{-6}$. Data come from secondary reports of the original papers. Some replication details and the list of positive claims come from Wikipedia's Allais effect article, so they are weaker than journal abstracts. That page also says NASA's 1999 coordinated tests were never published, which I could not confirm elsewhere.

**Sources:** [Kennefick, Physics Today 62(3), 37 (2009)](https://aip.brightspotcdn.com/PTO.v62.i3.37_1.online.pdf); [Kennefick, arXiv:0709.0685](https://arxiv.org/abs/0709.0685); [Principe value, British Astronomical Association](https://britastro.org/?p=213065); [Einstein Online on VLBI](https://www.einstein-online.info/en/spotlight/measuring-light-deflection-with-radio-telescopes/); [Cassini and VLBI summary, arXiv:0912.3421](https://arxiv.org/abs/0912.3421); [Helium, Solar eclipse of 18 August 1868](https://www.wikipedia.org/wiki/Solar_eclipse_of_18_August_1868); [Coronium (Wikipedia)](https://en.wikipedia.org/wiki/Coronium); [Coster et al., GRL 44, 12041 (2017)](https://munin.uit.no/handle/10037/13255); [Allais effect (Wikipedia)](https://en.wikipedia.org/wiki/Allais_effect); [Saxl & Allen, Erwin Saxl (Wikipedia)](https://en.wikipedia.org/wiki/Erwin_Saxl); [Wang et al., PRD 62, 041101 (2000)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.62.041101); [Kuusela, PRD 43, 2041 (1991)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.43.2041); [Luo Jun et al., PRD 44, 2611 (1991)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.44.2611); [Salva, PRD 83, 067302 (2011)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.83.067302); [Duif, arXiv:gr-qc/0408023](https://arxiv.org/abs/gr-qc/0408023); [Eckhardt, PRD 42, 2144 (1990)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.42.2144).

---

## General knowledge used without an online confirmation

These are standard values or facts that I did not confirm from a source I could open:

- L01: comparator resolution 0.1 µg and bench scale 0.1 g; skin absorption 60%, skin heat capacity; all flashlight figures (chosen assumptions).
- L02: Io's period; that Rømer didn't compute $c$ himself; the 35 km Michelson baseline; the 400 rev/s, 20 m Foucault-style numbers are illustrative only.
- L04: Compton's 1923 paper and the gauge-symmetry argument for $m_\gamma = 0$.
- L10: GRB 090510 redshift $z = 0.903$; 31 GeV is from the paper's abstract.
- L11: Einstein's teenage thought experiment; the 0.99975 figure from Alväger et al.
- L12: $H_0 = 67.4$ km/s/Mpc; tachyons being hypothetical with no detection; the idea that expansion is not local motion.
- L13: raindrop speed 6.5 m/s at 2 mm; rainfall rate is an assumption.
- L14: cesium's reactivity; the numbers for $g$ at equator and pole; compensated pendulums.
- L15: copper's Debye temperature, electronic coefficient and bulk modulus; emissivity and wire values; adsorption energy; BEC density and atom number (order of magnitude); cryogenic silicon cavities; nuclear demagnetisation.
- L16: attributing $1.61'' \pm 0.30''$ to Principe alone (the page I could open gives the number without saying so); the corona-heating problem as an active topic.
