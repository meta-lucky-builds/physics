# Antimatter and the hypothetical engine

These 21 questions are one design discussion. They start with "is there any antimatter?" and end with a closed-cycle engine that annihilates electron-positron pairs, takes work from the photons with moving mirrors, and turns the photons back into pairs.

This is a paper analysis of a cycle that can't deliver net work. It stays at the level of conservation laws, cross sections and thermodynamic limits.

**The result in one paragraph.** The student's own step A18 is correct: a closed cycle gives no net work. It follows from energy conservation alone, before any material problem appears. Real losses make the result strictly negative, and the only way to approach zero is to run infinitely slowly. Annihilation can't be a net energy source in a loop, because the pair must be bought back at full price. As an open-loop fuel it works, but it's a very expensive energy store, with about $10^{-9}$ to $10^{-8}$ of the input electricity becoming stored annihilation energy. The picokelvin sink doesn't change this: once you pay to keep it cold, its temperature cancels out of the efficiency.

**Things the student got right that are easy to miss.** Pair production is the real reverse process. Taking work directly from photons, without first turning them into heat, really can beat a 1500 K heat engine. The photon energy drop in A14 is real, and A18 reaches the correct conclusion without help.

**Things this section corrects that the design notes state as fact.**

- Reversing magnetic polarity does not make a trap hold the other charge (A13).
- Crossed fields of the Wien filter type do not sort charge at the balance speed (A17).
- Colliding photons head-on have no momentum problem for a nucleus to solve (A16).
- The efficiency formula in A20 isn't an open-loop efficiency.
- Maxwell's demon isn't the obstacle (A21).

Scripts are in [`code/`](code/). Their full output is at the [end of this page](#full-script-output). All energies are per pair unless stated. Sources are [listed below](#sources).

## Answers

### A01: Is there antimatter left, and can we make it by reversing annihilation?

**Status:** Established + Derived + Code

**Short answer:** No bulk antimatter is known to exist, only traces: cosmic-ray positrons and antiprotons, positrons from beta-plus decay, and a Galaxy-wide 511 keV annihilation glow. Your reverse process is real and is called pair production, and we detect the products routinely. It needs at least $2m_ec^2 = 1.022$ MeV of energy in the centre-of-momentum frame, and it can't happen with one free photon alone.

**Starting from:** Special relativity (energy and momentum conservation, invariant mass), plus the measurements cited below.

**Derivation:**

*What exists.*

- Cosmic rays: AMS-02 recorded $3.49\times10^5$ antiprotons against $2.42\times10^9$ protons, a raw event ratio of $1.4\times10^{-4}$ (not a corrected flux ratio) [[1]](#sources). Its positron fraction falls from 1 to about 8 GeV, then rises and stops increasing near $275\pm32$ GeV [[2]](#sources). The origin of that excess isn't settled.
- Near Earth: PAMELA reported antiprotons trapped in Earth's magnetic field [[3]](#sources).
- Beta-plus decay: proton-rich nuclei emit positrons. Decay of $^{26}$Al supplies about $0.3\times10^{43}$ positrons per second in the Galaxy [[6]](#sources), and hospital PET scans and positron labs use such isotopes ($^{18}$F and $^{22}$Na, general knowledge).
- Lightning: photonuclear reactions in one thunderstorm event produced positrons, seen through a 511 keV line lasting about a minute [[5]](#sources).
- The Galaxy: INTEGRAL/SPI sees about $10^{-3}$ photons cm$^{-2}$ s$^{-1}$ at 511 keV, which implies roughly $2\times10^{43}$ positrons annihilating per second. The figure depends on the assumed shape of the emission, and one model gives $3.7\times10^{43}$. About 97% of them annihilate through positronium [[6]](#sources).
- Bulk: Cohen, De Rújula and Glashow argued that annihilation at domain boundaries would have been seen, which excludes a matter-antimatter symmetric universe [[7]](#sources).

*The reverse process.* Define $s = (\sum E)^2 - |\sum \vec p\,c|^2$, which all frames agree on. A final state with an $e^+e^-$ pair has $\sqrt{s}\ge 2m_ec^2$.

For two photons at angle $\theta$:

$$s = 2E_1E_2(1-\cos\theta) \ge 4m_e^2c^4$$

So head-on photons need $E_1E_2\ge (511\text{ keV})^2$. Since $E_1+E_2\ge 2\sqrt{E_1E_2}$, their total energy is at least 1.022 MeV in any frame. For a photon on a nucleus of mass $M$, $s = M^2c^4+2EMc^2\ge(M+2m_e)^2c^4$, which gives $E\ge 2m_ec^2(1+m_e/M)$.

*What the experiments showed.* Be careful here, because "light into matter" headlines blur the details.

- **SLAC E144 (1997)** saw $106\pm14$ positrons when a 46.6 GeV electron beam met terawatt 527 nm laser pulses [[8]](#sources). A laser photon backscattered off an electron into a high-energy gamma, which then collided with several laser photons at once. So this is multiphoton (nonlinear) Breit-Wheeler, with the energy traceable to the electron beam. A follow-up reports a yield consistent with up to about 5 laser photons [[9]](#sources).
- **RHIC STAR (2021)** measured 6085 exclusive $e^+e^-$ pairs from ultraperipheral Au+Au collisions at $\sqrt{s_{NN}}=200$ GeV [[10]](#sources). The photons are the boosted fields of the passing nuclei. The paper says they behave as real photons through the equivalent-photon picture. Strictly they are nearly on-shell virtual photons from a field, and the authors report no direct measurement of their virtuality. This isn't two beams of free light colliding.
- **Two free real photons:** the 2014 collider proposal states that pure Breit-Wheeler had not been seen in the laboratory [[11]](#sources). I did not find a later direct observation, but my search wasn't exhaustive.

**Check:** [`code/verify_thresholds.py`](code/verify_thresholds.py) finds the photon-on-nucleus threshold by solving the conservation laws with a root finder, and samples random photon pairs:

```text
   proton         : k = 1.022554499 MeV, excess above 2 m_e c^2 = 556.597 eV (2m^2/M = 556.597 eV)
   Au-197 nucleus : k = 1.022000748 MeV, excess above 2 m_e c^2 = 2.847 eV (2m^2/M = 2.847 eV)
   865510 of 2000000 random photon pairs are above threshold; least E1+E2 among them = 2.03363 m c^2
[PASS] no above-threshold photon pair has E1 + E2 < 2 m c^2
```

The test can't show that QED governs pair creation at all. It only shows what the conservation laws allow.

**Limits:** Bound final states have a slightly lower threshold, because positronium is bound by 6.8 eV. See A16. Antimatter in distant domains is constrained by the argument in [7], not by direct search.

### A02: What if we produce antimatter at near light speed?

**Status:** Derived + Code + Established

**Short answer:** In a fixed-target accelerator it already comes out fast, and that's a cost rather than a benefit. At threshold, $p+p\to p+p+p+\bar p$ needs a beam kinetic energy of $6m_pc^2$ (5.63 GeV), and all four products leave at $0.866c$. The speed adds to the energy bill, and a fast antiproton must be slowed by a factor of tens of thousands in energy before a trap can hold it.

**Starting from:** Invariant mass; baryon-number and charge conservation, which require a $\bar p$ to be made with an extra $p$.

**Derivation:**

For beam energy $E$ on a target at rest, $s = (E+m c^2)^2-(pc)^2 = 2m^2c^4+2Emc^2$. The final state needs $\sqrt s\ge 4mc^2$, with equality when all products are at rest in the centre-of-momentum frame. Then

$$2m^2c^4+2Emc^2 = 16m^2c^4 \implies E = 7mc^2,\quad T = 6mc^2$$

Total lab energy is $8mc^2$ and total rest mass is $4mc^2$, so the products share $\gamma=2$, $\beta=\sqrt3/2$. Of the $6mc^2$ paid, only $2mc^2$ becomes new rest mass, and the remaining $4mc^2$ is kinetic energy. In a symmetric collider, each beam needs total energy $2mc^2$, so kinetic energy $2mc^2$ in total, which is why colliders are efficient at threshold.

More generally, a pair made with Lorentz factor $\gamma$ in the lab costs at least $2\gamma mc^2$. Annihilation in flight returns exactly that, but as a forward photon of $\gamma(1+\beta)mc^2$ and a backward one of $\gamma(1-\beta)mc^2$. So the kinetic energy is not lost in principle, but it arrives as harder photons along one direction, which any mirror handles worse.

*Why fast antimatter is hard to keep.*

- CERN's AD target takes $1.5\times10^{13}$ protons at 26 GeV/c and collects antiprotons near 3.5 GeV/c [[12]](#sources). The AD and ELENA then slow them to 5.3 MeV and 100 keV. Degrader foils alone trap well below 1% [[13]](#sources).
- A static trap can't capture a particle by conservative forces alone, because a trajectory that enters must also be able to leave. Capture requires removing energy, and that means cooling or degrading.
- Antiprotons at these energies also collide with their surroundings and annihilate, because anything solid is an antiproton's fuel.

**Check:** [`code/verify_thresholds.py`](code/verify_thresholds.py) doesn't use the invariant. It minimises the beam energy over all 12 final momentum components subject to energy and momentum conservation:

```text
   beam energy 7.000000 m c^2, beam kinetic 6.000000 m c^2, total lab energy 8.000000 m c^2
   final-state speeds at the minimum: [0.86603 0.86603 0.86603 0.86603]
   p [MeV/c]    KE [MeV]     beta    B rho [T m]   KE / 10 kV well
      3500.0    2685.311    0.9659       11.675      2.685e+05
       100.0       5.314    0.1060        0.334      5.314e+02
        13.7       0.100    0.0146        0.046      1.000e+01
```

The 10 kV well in the last column is an assumed round number for scale.

**Limits:** This is threshold kinematics with free protons. A proton bound in a nucleus lowers the threshold slightly, and real production runs far above threshold because yield is negligible near it.

### A03: Can annihilation drive an engine directly, and is it reversible?

**Status:** Derived + Code + Established

**Short answer:** It can drive an engine as a one-way fuel. Microscopically each annihilation has a time-reversed partner with an exactly related rate, so your instinct is right at that level. Macroscopically it's about as irreversible as a process can be. At 1500 K the equilibrium pair density is about $10^{-1.7\text{ million}}$ per cubic metre, so the reaction runs only one way. And even a perfect reversal would just return the energy: back and forth gives a battery with zero net work, not a source.

**Starting from:** Conservation laws; time-reversal invariance of QED; chemical equilibrium for a reaction whose photons have zero chemical potential.

**Derivation:**

*One photon can't do it.* A free photon has $E=|\vec p|c$, so $s=0$. A pair has $s\ge 4m_e^2c^4$. Since $s$ is the same in every frame and is conserved, the process $\gamma\to e^+e^-$ in vacuum is forbidden. So is $e^+e^-\to\gamma$. A nucleus or a second photon is required.

*Microscopic reversibility.* The squared amplitude for $e^+e^-\to\gamma\gamma$ equals that for $\gamma\gamma\to e^+e^-$. The two cross sections differ only in flux, phase space and a factor of 1/2 for identical final photons, which at centre-of-momentum pair speed $\beta$ gives

$$\frac{\sigma_{\gamma\gamma\to ee}}{\sigma_{ee\to\gamma\gamma}} = 2\beta^2$$

At low speed $\sigma_{ee\to\gamma\gamma}v\to\pi r_e^2c$ (without Coulomb focusing), so $\sigma_{\gamma\gamma\to ee}\to\pi r_e^2\beta$, which vanishes at threshold. This is the most useful single fact about the engine, and A16 and A18 use it.

*Macroscopic irreversibility.* Treat the reaction as chemistry. It is reversible only when the affinity vanishes, $\mu_++\mu_-=2\mu_\gamma=0$. For a dilute non-relativistic gas this requires

$$n_+n_-=n_{\rm eq}^2,\quad n_{\rm eq}=2\left(\frac{2\pi m_ek_BT}{h^2}\right)^{3/2}e^{-m_ec^2/k_BT}$$

At ordinary temperatures this density is absurdly small, so any real stock of positrons is far from equilibrium, the affinity is near 1 MeV per pair, and the reaction runs one way. Pairs only appear freely near $6\times10^9$ K, when $k_BT\sim m_ec^2$. A perfect reversal also requires the outgoing photons to return in the exact time-reversed state.

*A case that really does reverse.* If a perfect, fixed, spherical mirror surrounds a decaying para-positronium atom, the system is closed and unitary, and the two photons can return to the centre and re-form the state. This is general quantum mechanics (cavity revival), not something I found measured for positronium. But notice what it needs: a mirror that doesn't move, so no work is taken. Reversibility holds exactly when no work is extracted.

*What you got right about engines.* A combustion engine first makes heat, then loses part of that to the Carnot limit. Photons from annihilation have far lower entropy per joule than heat. A directed 511 keV photon in a known mode could in principle give up nearly all its energy to a moving mirror. For an open-cycle engine, the simple route is to absorb the photons (2.6 cm of lead attenuates 99% of them) and run a heat engine at 1500 K. That works, up to 80% with a 300 K sink. See A06.

Also, not every annihilation gives two photons. Spin-triplet positronium decays to three, and any two of those three are below threshold, so they can't simply collide back. The lifetimes are 125 ps and 142 ns [[14]](#sources).

**Check:** [`code/verify_breit_wheeler.py`](code/verify_breit_wheeler.py) integrates the tree-level QED squared amplitude numerically and compares it with both closed forms. [`code/verify_closed_cycle.py`](code/verify_closed_cycle.py) computes equilibrium and three-photon kinematics. [`code/verify_thresholds.py`](code/verify_thresholds.py) tests the one-photon case.

```text
    beta     sigma_BW [b]   amplitude [b]   sigma_ann [b]   amplitude [b]   ratio/2beta^2
    0.10       0.02519        0.02519         1.2595         1.2595   1.00000000
    0.70       0.17006        0.17006         0.1735         0.1735   1.00000000
   1000000 random pairs: min (E^2-P^2)/(mc^2)^2 = 4.000000 (a photon needs 0)
    1.50e+03    2.530e-07        -1.7169e+06
    6.00e+09    1.012e+00         3.5922e+01
   74966 sampled decays: largest s of any two photons / (2m)^2 = 0.9978
```

The amplitude comes from a textbook expression (Peskin and Schroeder eq. 5.105, general knowledge, not re-derived). Its agreement with two separately known closed forms is the test.

**Limits:** Radiative corrections and Coulomb focusing change rates at low energy. They don't change the conservation or equilibrium arguments. The non-relativistic $n_{\rm eq}$ is only rough near $6\times10^9$ K.

### A04: Combining antimatter with a 38 pK rubidium condensate

**Status:** Established + Derived + Code

**Short answer:** The temperature is real: a 2021 Bremen drop-tower experiment reached 38 pK of internal kinetic energy in about $10^5$ rubidium atoms. But it can't work as your cold reservoir. It's a tiny, isolated, free-falling cloud that isn't even a thermal equilibrium state, it holds $10^{-15}$ of the energy of one annihilation, and a 511 keV photon passes straight through it.

**Starting from:** Deppner et al. (cited), kinetic theory, photon attenuation.

**Derivation:**

Deppner et al. combined an interaction-driven collective mode with a magnetic lens in the Bremen drop tower (4.74 s free fall), reaching an internal kinetic energy of $\tfrac32k_B\cdot 38$ pK in a condensate of about $1.0\times10^5$ $^{87}$Rb atoms [[15]](#sources). An APS Viewpoint calls such lens figures effective temperatures, and reports expansion near 60 µm/s [[16]](#sources).

Three things matter for an engine.

1. **It isn't a reservoir in the Carnot sense.** Carnot's $T_c$ is the temperature of a large bath in equilibrium. A collimated condensate has a very small spread of velocities because it was shaped with a lens, not because it was thermalised at 38 pK. So using it as $T_c$ is already generous.
2. **It's tiny.** Even with the ideal-gas value, which overstates a condensate, $C=\tfrac32Nk_B=2.07\times10^{-18}$ J/K, so the thermal content is $CT_0=7.9\times10^{-29}$ J.
3. **It's decoupled.** The Compton optical depth for a 511 keV photon is about $10^{-15}$. The annihilation photons don't touch it, and any heat must come through a wall, which is much hotter than 38 pK in every real apparatus.

The drop-tower setting makes a further point. Your design says "after the antimatter is produced, then produce 38 pK". That ordering is correct in one way, since it's a preparation step, but the sink has to stay cold while heat flows into it. It doesn't. See A06.

**Check:** [`code/verify_carnot_sink.py`](code/verify_carnot_sink.py):

```text
   C = 3/2 N k_B = 2.071e-18 J/K, thermal content C T0 = 7.870e-29 J
   pair energy 2 m c^2 = 1.637e-13 J = 2.08e+15 times the thermal content
   one-axis rms speed at 38 pK: 60 um/s (the APS Viewpoint quotes about 60 um/s)
   Compton optical depth of the cloud for 511 keV (37 electrons per atom, radius 150 um): 1.5e-15
```

The 150 µm radius is my assumption for a cloud size, chosen from the reported expansion. The optical depth scales as $1/r^2$, and the conclusion survives any plausible value.

**Limits:** I treat the cloud as an ideal gas at 38 pK to give the idea its best chance. A careful treatment would be worse for the engine, because a condensate's heat capacity falls steeply below its transition temperature.

### A05: Take 1500 K as the hot side

**Status:** Derived + Code

**Short answer:** Then $T_c/T_h = 38\text{ pK}/1500\text{ K}=2.53\times10^{-14}$, so the Carnot limit is $1-2.53\times10^{-14}$. 1500 K is a sensible choice, since it's roughly what a very good solid absorber can survive, but A06 shows why this number isn't what the engine would achieve.

**Starting from:** The second law for a cyclic engine between two large baths, $\eta\le 1-T_c/T_h$.

**Derivation:**

The hot side is set by the material, not by the fuel. An annihilation photon has MeV energy, but once it has been absorbed and thermalised in a 1500 K block, only the block's temperature counts. Everything that makes the photon special has been thrown away in that step. So this design already gives up the photon's low entropy, which is the one thing that could beat a heat engine. With $T_h=1500$ K:

| Cold side | Carnot limit |
|---|---|
| 300 K ambient | 0.800 |
| 2.725 K, the CMB [[17]](#sources) | 0.99818 |
| 38 pK | $1-2.53\times10^{-14}$ |

**Check:** [`code/verify_carnot_sink.py`](code/verify_carnot_sink.py):

```text
   T_c/T_h = 2.5333e-14, so eta_C = 1 - 2.533e-14
   eta_C with T_c = CMB 2.725 K  : 1 - 1.817e-03 = 0.998183
   eta_C with T_c = 300 K ambient: 1 - 2.000e-01 = 0.800000
```

**Limits:** The 1500 K value is your choice of operating temperature, not a limit derived here. Refractory materials can go somewhat higher, which would change the numbers only modestly.

### A06: What is the maximum efficiency we can reach?

**Status:** Derived + Code + Estimate

**Short answer:** The formula gives $1-2.5\times10^{-14}$, but that isn't reachable. Paying a perfect refrigerator to keep the sink cold cancels $T_c$ exactly, and leaves $1-T_a/T_h = 80\%$ with a 300 K environment, or at most 99.8% in deep space against the CMB. For proton-antiproton fuel, multiply by about 0.56, because neutrinos carry the rest away.

**Starting from:** First and second laws. A reversible engine and a reversible refrigerator. Kinetic theory for a dilute gas. Stated assumptions: $T_a=300$ K, $N=10^5$, cloud radius 150 µm.

**Derivation:**

**1. The sink warms up.** Take a hot bath at fixed $T_h$ and a sink of heat capacity $C$ starting at $T_0$. Small reversible cycles have $dQ_c = dQ_h\,T/T_h$ and $dT=dQ_c/C$, so

$$T(Q_h)=T_0\,e^{Q_h/(CT_h)},\qquad W=Q_h-C(T-T_0)$$

The textbook figure holds only while $Q_h\ll CT_h = 3.1\times10^{-15}$ J, which is 19 keV, about 2% of one pair. The engine stops when $T=T_h$, at

$$Q_{\max}=CT_h\ln(T_h/T_0)=9.7\times10^{-14}\text{ J},\qquad \bar\eta = 1-\frac{1-T_0/T_h}{\ln(T_h/T_0)}=0.968$$

So your whole sink can't absorb the entropy from even one annihilation at 1500 K. It accepts 59% of a pair's energy, at an average efficiency of 96.8%, not $1-10^{-14}$. To stay near 38 pK for just one pair would take over $5\times10^6$ atoms. For 1 kg of annihilated mass it takes about $4\times10^{11}$ kg of rubidium at 38 pK, which isn't a meaningful object to build.

**2. Paying for the sink cancels $T_c$.** Heat $Q_c=Q_hT_c/T_h$ must leave the sink and reach the environment at $T_a$. A perfect refrigerator has COP $=T_c/(T_a-T_c)$ and needs work $W_f=Q_c(T_a-T_c)/T_c$. Then

$$\frac{W_{\rm engine}-W_f}{Q_h}=\frac{1}{Q_h}\left[Q_h-Q_c-Q_c\frac{T_a-T_c}{T_c}\right]=1-\frac{T_a}{T_h}$$

$T_c$ is gone. The result can't depend on the size of the cold sink, because the true cold reservoir was always the environment. The Carnot COP from 38 pK to 300 K is $1.3\times10^{-13}$, which is exactly the factor that takes back the gain.

If the fridge achieves only a fraction $f$ of Carnot, the engine nets nothing unless $f>(T_a-T_c)/(T_h-T_c)=0.2$. Real picokelvin apparatus is far below that. The least ideal work to cool this cloud from 300 K is $C[T_a\ln(T_a/T_0)-(T_a-T_0)]=1.8\times10^{-14}$ J. If the apparatus uses about 1 kW for 10 s (my assumption, order of magnitude only), that's $10^4$ J, about $10^{17}$ to $10^{18}$ times the ideal minimum.

**3. Heat can barely enter the gas.** Atoms at 38 pK move at about 96 µm/s. A wall only slightly warmer than the gas, at $2T_0$, gives heat flow

$$P_c=N\frac{3\bar v}{4r}\,2k_B(T_w-T_0)\approx5\times10^{-29}\text{ W}$$

An ideal engine dumping into this rate gives $P_c(T_h/T_c-1)\approx2\times10^{-15}$ W. Rejecting one pair's waste takes about 82 s, longer than the 4.74 s drop. Radiation is hopeless too: a 38 pK blackbody radiates $10^{-49}$ W/m² at a peak wavelength longer than Earth's radius.

**4. The honest table.**

| Claim | Value | Status |
|---|---|---|
| Carnot arithmetic | $1-2.5\times10^{-14}$ | Correct arithmetic, not achievable |
| Finite real sink, ignoring its cost | 0.968 at the limit of its capacity | Reversible idealisation |
| Sink cost included, 300 K world | 0.80 | Hard upper limit for heat-based engine |
| Sink cost included, CMB radiator | 0.998 | Upper limit, radiator must be huge |
| With $p\bar p$ neutrino loss | $\times\approx0.56$ | See A19 |
| Direct photon work, bypassing heat | can in principle exceed all of these | See A14, but needs impossible mirrors |

**Check:** [`code/verify_carnot_sink.py`](code/verify_carnot_sink.py) marches small cycles step by step, uses sympy for the cancellation, tracks a full prepare-run-recover ledger, and checks the wall heat rate with a Monte Carlo gas:

```text
   heat accepted 9.725e-14 J (59.4% of the pair energy) before T_sink = T_h
   mean efficiency over the stroke 0.9681
   closed form: Q_max = C T_h ln(T_h/T_0) = 9.725e-14 J, mean eta = 0.9681
   sympy: (W_engine - W_fridge)/Q_h = (-T_a + T_h)/T_h
     Q_h [J]    W_engine [J]   W_prep [J]   W_recovered [J]   net/Q_h   1-T_a/T_h
     1.0e-14   1.0000e-14   1.7829e-14      1.5829e-14     0.80000     0.800
   kinetic theory: P_c = N (3 v/4r) 2 k (T_wall - T_gas) = 5.05e-29 W
   Monte Carlo (200000 atoms, 12182 wall hits): P_c = 5.08e-29 W
   time to dump one pair's Carnot waste: 82 s; drop-tower free fall is 4.74 s
```

The lifetime ledger shows the same $1-T_a/T_h$ with a finite sink, from separate integrations of the engine and of the sink's exergy. The Monte Carlo can't test whether real atoms would stick to a surface, and in practice they'd be lost from the cloud on contact, which makes things worse.

**Limits:** Constant heat capacity up to 1500 K is a generous fiction for the warming sink. The 1 kW figure is an assumption. The heat-rate numbers scale with the cloud size, which is also assumed.

### A07: A cascade making antiparticles of the debris

**Status:** Derived + Code + Established

**Short answer:** It can't gain energy. Conservation laws force an antiparticle to be made together with its partner particle, so each step costs at least $2mc^2$ and returns at most $2mc^2$, and it leaves you with just as much debris as you began with. Real production costs about $10^8$ times more than the output. Debris with kinetic energy is worth harvesting directly, but as ordinary energy, with no antiparticles involved.

**Starting from:** Conservation of charge, lepton number and baryon number. Invariant mass. Sourced production efficiencies.

**Derivation:**

Let debris $X$ have mass $m$. To destroy it you need one $\bar X$.

1. Additive quantum numbers forbid making $\bar X$ alone. The reaction is at best $\ldots\to X+\bar X+\ldots$, with $s\ge(2mc^2+\ldots)^2$, so cost $\ge 2mc^2$. With a target the cost is above that: $2mc^2(1+m/M)$ on a nucleus, and $6m_pc^2$ for a proton on a proton.
2. $\bar X+X$ releases $2mc^2$ plus whatever kinetic energy they carried. That kinetic energy was already yours, and doesn't depend on the cascade.
3. Net gain of the step $= 2mc^2-\text{cost}\le0$. And the new partner $X$ replaces the destroyed one: net $X$ number never changes.

Any quantity conserved additively blocks the cascade the same way. If annihilation output is recycled, the total can't exceed what you began with. Input $E_0$ with production efficiency $\varepsilon$ decays as $\varepsilon^k$ across generations.

**Practical facts make this worse.** A pion decays in 26 ns [[18]](#sources), so it's $10^9$ lifetimes gone by the time a 100 s AD cycle finishes. Neutrinos barely interact: at about 1 MeV the cross section is near $G_F^2s/\pi\sim10^{-44}$ cm², so the mean free path in lead is tens of light years (order of magnitude). You can't catch them to annihilate them. Your detector idea also misses this: detection itself costs almost nothing, but it doesn't help, since what's limited is energy and time.

**Isotopes.** Beta-plus nuclei do give positrons without a pair-production bill, but the energy comes from the nuclear mass difference. That was paid for in a star or in the accelerator that made the isotope. So they're nuclear fuel, not free antimatter.

**One real idea nearby.** Antimatter-catalysed fusion is studied because small amounts of antiprotons might release much larger nuclear energy [[19]](#sources). There the energy comes from the fusion fuel, not from a cascade. Your instinct to look for leverage there is sound.

**Check:** [`code/verify_closed_cycle.py`](code/verify_closed_cycle.py) tracks $X$ and $\bar X$ counts as integers, with no result assigned by hand:

```text
   X    2 m c^2 [MeV]   least cost on a proton target [MeV]     gain     net [MeV]
   e          1.022                  1.023                1.022       -0.001
   mu       211.316                235.112              211.316      -23.796
   pi       279.140                320.663              279.140      -41.523
   p       1876.544               5629.632             1876.544    -3753.088
   eps = 1       pairs per generation [1000, 1000, 1000, 1000]..., final energy 1.000000e+09 MeV, debris X = 1000, leftover X-bar = 0
   eps = 1e-08   pairs per generation [9, 0, 0, 0]..., final energy 8.020001e+07 MeV, debris X = 1000, leftover X-bar = 0
```

This models bookkeeping, not collision dynamics. It can't test QED, but the argument doesn't depend on dynamics.

**Limits:** The $e$, $\mu$ and $\pi$ rows use the photon-on-proton threshold, as a simple lower bound for that route. Colliders can reach $2mc^2$ total. That makes the step break even at best, which doesn't change the conclusion.

### A08: "If physics says yes, engineers will find a way"

**Status:** Derived + Conceptual

**Short answer:** That logic is fine, and it has a good record. But for the part that matters here physics says no. It says yes to nearly every component and no to the claim that a closed loop produces net work. No engineer can fix a false antecedent.

**Starting from:** Conservation laws and the second law, treated as the best-tested claims we have. They aren't logically necessary truths, so this is the strongest available argument and not a proof about all possible future physics (see section 07).

**Derivation:** Going through each element:

| Element | Does physics allow it? | What says so | Who stops it |
|---|---|---|---|
| Making antimatter | Yes | Pair production, measured | Engineering: efficiency $10^{-9}$ to $10^{-8}$ |
| Storing it | Yes | Penning traps, measured | Engineering: density and quantity |
| Annihilating on demand | Yes | QED | Engineering: rates |
| Reflecting at grazing angles | Yes | Fresnel, multilayers near 0.65 MeV | Engineering, narrow range |
| Normal-incidence mirror at 511 keV | No ordinary material | $R\approx(\delta/2)^2\sim10^{-17}$ | Physics of matter, not conservation laws |
| Work from photons with moving mirrors | Yes | Doppler, measured radiation pressure | Mirrors |
| Pair creation from two photons | Yes | Breit-Wheeler | Very small rate |
| Sorting charges and re-trapping | Yes | Lorentz force; needs cooling | Engineering |
| Net work from a closed pair loop | **No** | First law, $s$ threshold | Nothing can fix this |
| Net work from one heat bath | **No** | Second law (Kelvin) | Nothing can fix this |
| Single photon to free pair in vacuum | **No** | Energy and momentum | Nothing can fix this |
| Massive parts at light speed | **No** | Relativity | Nothing can fix this |
| A free 'boost' from $p\bar p$ | **No** | Energy conservation | Burns fuel |

A closed-cycle engine could only beat zero if the first law were wrong for it, and then the discovery would matter far more than the engine. See A18 for the same argument with the numbers.

**Limits:** "No" here means no within all tested physics. It's a claim about evidence, not a claim that physics is finished.

### A09: First list all the problems in the model

**Status:** Derived

**Short answer:** They fall into three groups: things physics forbids, things physics allows but nobody has built, and statements in the design that are just wrong. The third group is the most useful one, because those errors are easy to fix. Each entry is supported in the answer named.

**Starting from:** The results of A01 to A21 of this section.

**Derivation:**

**(i) Forbidden by physics, regardless of technology**

1. Net work from the closed pair cycle (A18).
2. Getting a free boost from $p\bar p$ while calling it catalysis (A19).
3. $\gamma\to e^+e^-$ with one free photon in vacuum (A03).
4. Massive mirrors or particles reaching $c$ (A10).
5. A stock of positrons in true equilibrium with photons at any engineering temperature (A03).
6. Recapturing neutrinos, in any practical sense (A07, A19). This one is a matter of scale rather than a strict law, and I list it here because the scale is tens of light years of lead.
7. A cold sink that gains the Carnot number without paying for it (A06).

**(ii) Allowed, not built, and possibly not buildable**

1. A mirror with $1-R\lesssim10^{-14}$ at 511 keV, or any mirror near normal incidence (A14).
2. A closed gamma-ray storage ring (A15).
3. Colliding free gamma beams with useful optical depth, needing over 10 kJ per pulse in a micron spot (A16).
4. Antimatter production better than about $10^{-8}$ (A20).
5. Storage at useful quantity. Even the Brillouin limit stores less energy per volume than petrol (A21).
6. A sink coupled to the engine with useful power (A06).

**(iii) Wrong claims in the design as written**

1. "Opposite field polarity" as the way to hold the other charge (A13).
2. Magnetic confinement "ensures they meet at the centre" (A13).
3. A heavy nucleus or strong field being needed for head-on photons (A16).
4. Crossed fields sorting charge in "mass spectrometer style" (A17).
5. Usable photon energy of "1 to 2 GeV" from $p\bar p$ (A19). It's about 0.7 GeV.
6. The $\eta_{\max}$ formula as an open-loop efficiency (A20).
7. Capture being limited by a quantum measurement problem (A21 B).
8. Maxwell's demon being the essential obstacle (A21 D).
9. Perfect alignment being required for pair production (A21 C).
10. The 1 − pico efficiency as an achievable value (A11).

**Limits:** The grouping of item (i)6 is a judgement. The rest follow from the derivations cited.

### A10: What if we could produce light speeds in the engine?

**Status:** Derived + Code

**Short answer:** Massive parts can't reach $c$, and the photons already travel at $c$. What you can do is move a mirror away at high speed, so that the photon gives up more energy per bounce. The fraction per bounce is $2\beta/(1+\beta)$, which tends to 100% only as $\beta\to1$. That limit is unreachable, and even near it the photon takes very long to catch the mirror.

**Starting from:** Lorentz transformations. A perfect mirror reverses a photon in the mirror's rest frame.

**Derivation:**

*Massive parts.* Energy $\gamma Mc^2$ diverges as $\beta\to1$. A photon can't push a mirror to $c$ either: for a mirror at rest, energy-momentum conservation gives a final speed strictly below 1 for any $E/Mc^2$.

*The moving mirror.* A photon of energy $E$ meets a mirror receding at $\beta$. In the mirror frame, $E_m=\gamma(1-\beta)E$. Reflection reverses it, and transforming back gives

$$E'=\gamma(1-\beta)E_m=\gamma^2(1-\beta)^2E=\frac{1-\beta}{1+\beta}E$$

The work on the load per bounce is $E-E'=\frac{2\beta}{1+\beta}E$. This is where your intuition is right: faster recession means more extraction, and $\beta\to1$ means all of it.

*Why that limit doesn't deliver.* The photon takes $t=L_0/[c(1-\beta)]$ to catch the mirror, which diverges as $\beta\to1$. So you can take nearly everything in one bounce, in nearly infinite time. Also a mirror that recedes at fixed $\beta$ needs a load that holds its speed.

*A free mirror.* With no such load, the mirror recoils and $E'=E/(1+2E/Mc^2)$. The fraction extracted is $2\varepsilon/(1+2\varepsilon)$ with $\varepsilon=E/Mc^2$. For useful extraction in one bounce, the mirror's mass-energy has to be about the photon energy. A connection worth seeing: put $M=m_e$, $E=511$ keV, and this is exactly Compton backscatter to 170 keV. A single free electron *is* a mirror, one that gets carried off at $0.8c$ with two-thirds of the energy. That's why a real gamma-ray photon loses so much in matter.

**Check:** [`code/verify_optical_piston.py`](code/verify_optical_piston.py) does each reflection by Lorentz-transforming, and audits it against lab-frame conservation and a separate root-finding solve:

```text
    beta     E'/E (sim)   (1-b)/(1+b)    work/E
   0.500000   0.33333333   0.33333333 0.66666667
   0.990000   0.00502513   0.00502513 0.99497487
   mirror = one free electron: E' = 170.333 keV; Compton 180 deg formula 170.333 keV
   the electron leaves with beta = 0.8000, carrying 340.7 keV
   E/Mc^2 =    1e+02: fraction extracted in one bounce = 0.995025, mirror beta = 0.99995050
```

**Limits:** This treats the mirror as a single rigid object with one velocity. A real mirror has a finite sound speed, and at near-$c$ recession it couldn't stay rigid.

### A11: Reading it as 1 − pico efficiency

**Status:** Derived + Code

**Short answer:** As arithmetic your reading is about right: $1-2.5\times10^{-14}$, which is $1-0.025$ pico rather than $1-1$ pico, since 1500 K is not 1 K. As an efficiency you could reach, it's wrong. Carnot's formula assumes an infinite cold bath that costs nothing and accepts heat at any rate, and none of those hold.

**Starting from:** The Carnot formula and the results of A06.

**Derivation:** The form of the expression really is $\eta=1-T_c/T_h$, so a small $T_c$ does give a number just below 1, and you're right that it's not 100%. The ratio is $38\times10^{-12}/1500=2.53\times10^{-14}$. For the loss to be one pico, you'd need $T_h=38$ K.

The honest correction has three steps. A real sink warms, so the mean value is 0.968 at capacity. Paying for the sink takes it to $1-T_a/T_h$. And the power at the ideal value is near 2 fW for this cloud. In short, the number is the right answer to a question that this engine doesn't get to ask.

**Check:** [`code/verify_carnot_sink.py`](code/verify_carnot_sink.py). The table shows how quickly the textbook value goes away as the heat load grows on the $10^5$ atom sink:

```text
     Q_h [J]      T_final [K]     mean eta      1 - mean eta
     1.0e-18      3.801e-11     1.000000      2.534e-14
     1.0e-14      9.502e-10     1.000000      1.889e-13
     5.0e-14      3.715e-04     1.000000      1.539e-08
```

**Limits:** Same as A06.

### A12: The core concept, a closed reversible annihilation cycle

**Status:** Derived + Code

**Short answer:** A closed cycle returning to the same state has zero net energy change, so net work out equals net energy in. With no heat bath, the best possible net work is zero. A reversible cycle isn't a way to gain energy, it's a way to lose none. The idea is a well-defined thought experiment, and works as an energy store with no gain.

**Starting from:** The first law, and the second law in Kelvin's form: no cycle turns heat from one bath wholly into work.

**Derivation:** Over one cycle the stored pair returns to its initial state, so $\Delta U=0$ and $W_{\rm net}=Q_{\rm net}$.

- If the only energy exchanged is work, $W_{\rm net}=0$, and any loss gives $W_{\rm net}<0$.
- If there's a heat bath at one temperature, Kelvin says $W_{\rm net}\le0$.
- If there are baths at $T_h$ and $T_c$, $W_{\rm net}\le Q_h(1-T_c/T_h)$. In that case the pairs are just an exotic working fluid. The annihilation term cancels over the cycle, and the work is ordinary heat-engine output paid for by the bath. And thermally recreating pairs from a bath needs about $6\times10^9$ K.

The rest of this section walks through each step, because why it fails is instructive. See the [whole-cycle account](#whole-cycle-energy-account).

**Check:** The event-level run in [`code/verify_closed_cycle.py`](code/verify_closed_cycle.py) closes the first law per event and finds no positive net work. See A18 for output.

**Limits:** The argument assumes the first law holds across the whole cycle, including the rest energy of the particles. That's the best-tested assumption available.

### A13: Steps 1 and 2, storage and annihilation

**Status:** Derived + Code + Established

**Short answer:** Storing the species separately and getting 511 keV back-to-back photons from pairs at rest are both fine. Three statements are wrong. A Penning trap isn't a magnetic bottle. Reversing magnetic polarity doesn't change which charge is held; the electrostatic voltage does that. Magnetic fields don't steer particles to meet at the centre. Annihilation at laboratory densities is also very slow without Coulomb help.

**Starting from:** The Lorentz force $\vec F=q(\vec E+\vec v\times\vec B)$. Quoted measurements.

**Derivation:**

*What a Penning trap is.* A uniform $\vec B$ gives radial confinement. A static electric quadrupole $\phi\propto s(z^2-r^2/2)$ gives axial confinement [[20]](#sources). Earnshaw's theorem means the electric part alone can't hold anything in 3D: it's a well along $z$ and a hill along $r$. The axial motion has no $\vec B$ force at all, so axial stability depends only on the sign of $qs$. Flipping $\vec B$ only reverses the sense of the rotation. Radial stability requires $\omega_c^2>2\omega_z^2$, which depends on $B^2$ and not its sign.

*What a magnetic bottle is.* It uses field strength alone: a particle's magnetic moment is conserved, so it reflects off stronger field unless its pitch angle is inside the loss cone. Both charges reflect the same way, whichever way the field points. So a bottle sorts by pitch angle, not by charge, and could hold both species together, where they'd annihilate.

*So for storage.* Two Penning traps with **opposite electrode voltage** work. Same magnetic field, in fact, which is convenient. In one structure, a *nested* Penning trap holds opposite charges in adjoining axial wells. This has been done, with cold positrons and antiprotons stored together at 4.2 K [[21]](#sources).

*Meeting at the centre.* Nothing magnetic pulls them to a point. They meet because axial motion carries energised particles through each other in a nested well, or because an electric potential is changed to overlap them.

*The rate.* The low-energy rate coefficient is $\langle\sigma v\rangle\approx\pi r_e^2c=7.5\times10^{-21}$ m³/s. At $n=10^{15}$ m⁻³ the annihilation time is $1.3\times10^5$ s, about 37 hours. Coulomb attraction raises this by a Sommerfeld factor near 800 at 4 K, and positronium formation adds other channels, but "controlled mixing" isn't a fast power stroke.

*The photons.* Spin-singlet annihilation at rest gives two photons of 511 keV each, back to back. Your statement is right there. Triplet positronium gives three photons with a continuous spectrum. From the lifetimes, $(142\text{ ns}/125\text{ ps})/3\approx380$ for the spin-averaged $2\gamma:3\gamma$ ratio. The commonly quoted value is about 372 (general knowledge), and in the Galaxy the positronium fraction is about 97% [[6]](#sources).

**Check:** [`code/verify_traps.py`](code/verify_traps.py) integrates orbits with only the Lorentz force and classifies each case from the trajectory:

```text
     q    B     s     r_max     z_max     result
   +1  +1.0  +1       1.21      1.00   confined
   +1  +1.0  -1       1.19    200.21   ESCAPES
   +1  -1.0  +1       1.00      1.00   confined
   +1  -1.0  -1       1.00    200.21   ESCAPES
   -1  +1.0  +1       1.00    200.21   ESCAPES
   -1  +1.0  -1       1.00      1.00   confined
   -1  -1.0  +1       1.19    200.21   ESCAPES
   -1  -1.0  -1       1.21      1.00   confined
   pitch[deg]  +,B+    +,B-    -,B+    -,B-
         20   lost   lost   lost   lost 
         35   held   held   held   held 
   cold positive stays in z = [1.78, 2.20], cold negative in z = [-0.20, 0.20]
   positive given energy above the central hill passes z = 0: range [-2.93, 2.93]
```

The runs are single-particle, in scaled units, with no space charge. The nested-well test is axial only. The real radial holding comes from the uniform $B$ tested in the first table.

**Limits:** Space charge and rotating-wall effects matter at high density, and they're not simulated. The Sommerfeld numbers use a mean relative speed, and aren't a full thermal average.

### A14: Step 3, the optical piston and the perfect gamma mirror

**Status:** Derived + Code + Established

**Short answer:** The physics of the piston is right: photons do work on a receding mirror and redshift exactly as much, with energy conserved at each bounce. Your acknowledgement that the mirror is fictional is what saves the step. At 511 keV, real matter reflects about $10^{-17}$ of the light at normal incidence, and grazing reflection only works below about 0.16 mrad. And one photon needs about $10^{14}$ bounces to give its energy to 1 g mirrors.

**Starting from:** Lorentz transformations. Fresnel reflection with a free-electron refractive index. NIST attenuation coefficients. Sourced optics measurements.

**Derivation:**

*The piston.* Using $E'=E(1-\beta)/(1+\beta)$ from A10 with one mirror receding, consecutive bounces have lengths $L'=L(1+\beta)/(1-\beta)$, so

$$E\,L=\text{constant at each bounce}$$

This is exact for a constant-speed mirror, not just approximate. It's the same as the average force $E/L$ doing work $dW=E\,dL/L$. In the free-mirror case it's the adiabatic invariant, and total photon energy plus mirror kinetic energy stays at $E_0$.

*How many bounces.* For two free mirrors, each of rest energy $Mc^2$, nearly all of $E_0$ ends up as kinetic energy, so each moves at $v\to c\sqrt{E_0/Mc^2}$, and the bounce count scales as $\sqrt{Mc^2/E_0}$. The simulation gives $N\approx3\sqrt{Mc^2/E_0}$ to reach $L=100L_0$. For 1 g mirrors and one 511 keV photon: $N\approx10^{14}$. To survive those bounces with a $1/e$ chance, each reflection must lose less than about $10^{-14}$.

*What matter does at 511 keV.* The photon wavelength is 2.43 pm, smaller than atoms, and its energy equals $m_ec^2$.

- **Refractive index.** $n=1-\delta$ with $\delta=r_e\lambda^2n_e/2\pi$, which is $1.2\times10^{-8}$ for gold.
- **Normal incidence.** $R=(\delta/2)^2=3.8\times10^{-17}$.
- **Grazing.** Below the critical angle $\theta_c=\sqrt{2\delta}=0.157$ mrad (about 0.5 arcmin) gold reflects well in the calculation, with $R=0.988$ at 0.126 mrad. A real multilayer has reflected up to 0.65 MeV, which is the record I could confirm, on WC/SiC with 1 to 2 nm periods [[22]](#sources). I couldn't read the reflectivity values in that paper.
- **Bragg and Laue crystals.** Si(111) has a Bragg angle near 0.22° and a width near 0.1 arcsec. Measured Laue peak reflectivity is 26% to 31% near 600 keV for Cu and Au mosaics [[23]](#sources), and mosaic crystals are capped at 50% [[24]](#sources). These deflect by under a degree. A lens isn't a mirror.
- **True backscatter from a crystal.** This requires plane spacing $\lambda/2=1.2$ pm. No lattice has it. The thermal amplitude factor at that momentum transfer is about $10^{-328}$, and the 1s form factor is about $10^{-4}$.
- **Compton scattering.** Each scatter takes on average 34% of the energy, and backscatter takes two-thirds. This isn't reflection: it's the $M=m_e$ mirror of A10.

*Your two geometries conflict.* A sphere around the source receives radial photons at normal incidence, which is the worst case. Two facing mirrors catch only photons near their axis, and at grazing angles a cone of half-angle $\theta_c$ covers $6\times10^{-9}$ of the sphere. Photon direction is effectively random (A21 B), so the design needs a surface reflecting at every angle, and that's what's unavailable.

*How dense would a true metal mirror have to be?* A plasma reflects below its plasma frequency. For $\hbar\omega_p=511$ keV, $n_e\approx1.9\times10^{38}$ m⁻³, about $4\times10^7$ times solid gold. That's a mass density near $6\times10^{11}$ kg/m³, on the scale of a white dwarf core. The electrons there are relativistic, so this is only a scale. It shows why "fictional" is the right word.

**Check:** [`code/verify_optical_piston.py`](code/verify_optical_piston.py) and [`code/verify_gamma_mirror.py`](code/verify_gamma_mirror.py):

```text
    beta   bounces  E_final/E0   W_total/E0   sum-1      E L/(first E L)  time
   0.010      231     0.00985     0.990149  -1.1e-16   1.00000000    9950.46
   0.500        5     0.00412     0.995885  -1.1e-16   1.00000000     322.00
   E0/Mc^2    bounces(to L=100 L0)  E_f/E0     sum check   beta_mirror   sqrt(E0/M)   N*sqrt(eps)
      1e-08              29933     0.009998    1.5e-12    9.9504e-05  1.0000e-04     2.993
   fitted N ~ eps^-0.500, expected eps^-0.5
       1e-03 kg    9.93e+13              1.01e-14              9.05e-06 m/s
   Au        4.667e+30    1.232e-08       0.1570           0.540
   normal incidence: R = 3.795e-17, (delta/2)^2 = 3.795e-17
       0.1256     0.988079     0.987903
       0.2355     0.021287     0.019265
   Klein-Nishina at 511 keV: sigma = 0.2865 b per electron (closed form 0.2865 b), mean energy lost per scatter 34.4%
```

The piston simulation doesn't use the Doppler formula. The 1 g result extrapolates a law fitted between $10^{-4}$ and $10^{-8}$ by a further 20 orders of magnitude in $\varepsilon$, so its correct reading is "of order $10^{14}$".

**Limits:** Fresnel gives only the coherent specular part. The Debye-Waller number uses an assumed thermal displacement of 0.075 Å, and the exact exponent is not meaningful beyond "hopeless". The 0.3 nm roughness column uses a standard Nevot-Croce factor (general knowledge). The Bragg width uses an assumed form factor.

### A15: Step 4, the photon storage ring

**Status:** Derived + Code

**Short answer:** A storage ring for 511 keV photons can't be built from known reflection. Photons have no charge, so nothing bends them but reflection, and a full turn by grazing bounces on smooth gold loses the photon with survival about $10^{-78}$ per lap, independent of the ring's size. Crystals would be worse.

**Starting from:** Fresnel reflection with the complex index of gold, using NIST attenuation.

**Derivation:** A ring of bounces at grazing angle $\theta$ turns the photon by $2\theta$ each time, so a full turn takes $N=\pi/\theta$ bounces. For small $\theta$ the loss per bounce is proportional to $\theta$:

$$1-R\approx\frac{4\theta\beta_n}{\theta_c^3}$$

where $\beta_n$ is the imaginary part of the index. So the loss per lap has a constant limit:

$$N(1-R)\to\frac{4\pi\beta_n}{\theta_c^3}\approx180$$

Making more, gentler bounces doesn't help, because the extra bounces exactly offset the smaller loss per bounce. A larger ring doesn't help either. Survival per lap is $e^{-180}$.

Laue diffraction with Au(111) deflects by $2\theta_B\approx10$ mrad, so about 610 steps per turn. At the measured 29% that's $10^{-328}$ per lap, and this also ignores that a Laue crystal accepts a very narrow band and a lens has no way to steer the photon along a closed path.

Note a further catch in the design: this ring holds photons which must then be raised back above threshold to be useful, so A16's energy has to be paid at this point.

**Check:** [`code/verify_gamma_mirror.py`](code/verify_gamma_mirror.py) computes Fresnel reflectivity from the full complex index, and compares with the limit:

```text
    theta/theta_c      N       -ln(survival per lap)   predicted 4 pi beta/theta_c^3
         0.5000   4.00e+04           207.85               180.01
         0.0100   2.00e+06           180.01               180.01
   Survival per lap in the smooth limit: exp(-180.0) = 10^-78.2
   Au(111) at 511 keV: theta_B = 5.151 mrad, so 610 diffractions per 2 pi turn; 0.29^610 = 10^-328
```

**Limits:** The imaginary index uses total attenuation, which includes Compton scattering, and the true loss for an evanescent wave at a surface is approximately but not exactly this. Roughness would make it worse. The Laue estimate applies the 588 keV reflectivity at 511 keV.

### A16: Step 5, the reverse process

**Status:** Derived + Code + Established + Estimate

**Short answer:** The condition you wrote is right: centre-of-momentum energy at least 1.022 MeV. But a nucleus or field isn't needed for head-on photons, because they have no net momentum in that frame. A static field or nucleus gives momentum and no energy, so neither lets redshifted photons cross the threshold. And right at threshold the cross section is exactly zero, so the pair must carry away some kinetic energy.

**Starting from:** Invariant mass. The Breit-Wheeler cross section (closed form, checked against the QED amplitude).

**Derivation:**

*The momentum problem doesn't exist here.* Equal photons head-on have total momentum zero. A pair at rest also has zero. The constraint you're remembering is for **one** photon, A03. If you do use a nucleus with one photon, the threshold rises: $E\ge2m_ec^2(1+m_e/M)$, which adds 2.85 eV for gold and 557 eV for a proton. A nucleus costs energy, it doesn't provide it.

*Redshifted photons.* With $E_1E_2<m_e^2c^4$ head-on, no pair forms from those two photons alone. A nucleus's Coulomb field takes up momentum and adds a recoil cost, not energy. A field can only supply energy by giving up its own, and that energy then belongs in the ledger. Examples: the laser in SLAC E144 with several photons absorbed, or the Schwinger process at about $1.3\times10^{18}$ V/m, in which the field supplies $m_ec^2$ over a reduced Compton wavelength and is itself depleted.

*The cross section.* In the centre-of-momentum frame, with pair speed $\beta$:

$$\sigma_{\gamma\gamma}=\frac{3}{16}\sigma_T(1-\beta^2)\left[(3-\beta^4)\ln\frac{1+\beta}{1-\beta}-2\beta(2-\beta^2)\right]$$

This is the Breit-Wheeler result, in the form quoted from Gould and Schréder [[25]](#sources). It gives $\sigma\to\pi r_e^2\beta$ near threshold, a peak of $0.256\sigma_T=1.7\times10^{-29}$ m² at $\beta=0.70$, which is photons near 717 keV. Above that it falls. So a pair formed at exact threshold has zero rate, and useful rates need $\beta\sim0.7$ with photons of about 717 keV, where $2(\gamma-1)m_ec^2\approx 410$ keV goes into the pair's kinetic energy. That has to be recovered or it's lost.

*How likely a conversion is.* Two photons of wavelength $\lambda$ can't be focused below about $\lambda^2$, by diffraction. So $P=\sigma/A\lesssim\alpha^2\sim5\times10^{-5}$ per crossing. For two single photons this is a fundamental weakness of the coupling. In a cavity with reflectivity $R$ per bounce, the pair is formed before the photons are lost only with probability $\varepsilon=P/(P+1-R^2)$. At $R=0.99$ that's $9\times10^{-4}$. For bulk beams, optical depth 1 in a 1 µm spot needs over 20 kJ of 722 keV gammas per pulse.

*An edge case.* Because positronium is bound by 6.8 eV, two photons can resonantly form it just below the free threshold. Its width is about $\hbar/125$ ps, or 5 µeV. Perfectly matched, the resonance is strong (estimated peak cross section near $10^{-24}$ m², from a Breit-Wigner resonance formula, not derived here), but you then need 6.8 eV to free the pair, which puts back the energy you skipped. This is why the threshold in the ledger stays $2m_ec^2$.

**Check:** [`code/verify_breit_wheeler.py`](code/verify_breit_wheeler.py):

```text
   beta=1e-03: sigma_BW/(pi r_e^2 beta) = 1.000001
   peak: beta = 0.7013, sigma = 1.701e-29 m^2 = 0.2556 sigma_T, head-on photons of 717 keV each
   E_S = m^2 c^3/(e hbar) = 1.32e+18 V/m; work over one reduced Compton wavelength = 511 keV
   1 um radius spot, sigma at peak: tau = 1 needs N = 1.85e+17 photons per pulse = 2.14e+04 J per pulse of 722 keV gammas
         R         eps       photons lost per pair made
    0.990000000   9.09e-04       1.099e+03
```

The focus area $\lambda^2/\pi$ is an estimate of the best possible, so these numbers are good to a factor of a few.

**Limits:** The cross section is lowest order in $\alpha$. Polarised photons, stimulated effects and dense-field effects change the rates. The bound-state remark is an estimate.

### A17: Steps 6 and 7, sorting and repeating

**Status:** Derived + Code + Established

**Short answer:** Separating the charges is easy, and no demon is needed. But your method is the wrong one. Crossed E and B fields of the Wien filter type are charge-blind at the balance speed and make both charges drift the same way. Use E alone or B alone. Then you must cool the particles into the traps, which removes their kinetic energy as heat.

**Starting from:** The Lorentz force. Liouville's theorem for a static field.

**Derivation:**

In crossed fields, the force vanishes when $\vec E+\vec v\times\vec B=0$, and this doesn't involve $q$. So at $v=E/B$ both charges go straight. In long crossed fields, every charge's guiding centre drifts at $\vec E\times\vec B/B^2$ regardless of its sign. Two charges started together separate only by the gyration, at most $4(v-v_d)/\omega_c$. In a real filter, a short region can separate off-balance particles, but that's the $\vec B$ force, not the crossed geometry, and the balance condition picks out speed, not charge.

- **E field alone:** forces are opposite, so separation grows as $Et^2$ in scaled units. This is simple and works.
- **B field alone:** a magnetic spectrometer bends the two charges to opposite sides.

*Is it thermodynamically free?* One positron and one electron have a mixing entropy of at most $k_B\ln2$ each, about 0.04 eV of work per pair at 300 K. Compare 1.022 MeV. Separation is not the expensive part. The species can also be pulled apart without compression, each keeping its volume, exactly like Gibbs's semipermeable membranes, and a field works as such a membrane.

*Where you do pay.* A static trap can't capture without energy removal (A02), so the kinetic energy of a newly made pair must leave as heat unless recovered, in a slowdown that isn't free. The same particles aren't needed, you're right: electrons are identical, so "energy equivalents" is a meaningful phrase.

Attraction between the charges costs at most 6.8 eV to break if they form positronium, and this returns when they meet again.

**Check:** [`code/verify_traps.py`](code/verify_traps.py):

```text
   at the balance speed v = E/B: both go straight; separation after t=50: 5.42e-14
   off balance, v = 0.5: drift speed e+ 0.1002, e- 0.1002, E/B = 0.1
   separation stays bounded: max 1.60 over t = 2000; gyroradius-scale 2(v-E/B)/w_c x 2 = 1.60
   E field alone: separation after t=50 is 250.0, prediction E t^2 = 250.0
   B field alone, half turn: y = -1.000 for e+, +1.000 for e- (opposite bends)
```

**Limits:** Non-relativistic orbits, uniform fields. A relativistic pair doesn't change the charge blindness, because the condition for zero force contains no $q$.

### A18: The ideal energy flow

**Status:** Derived + Code

**Short answer:** Your conclusion is right: $W\le0$, and it follows from conservation alone. I'd make it stronger in three ways. It holds for any photon energies and angles, not just equal head-on ones. It's strictly negative at any finite rate, because the cross section vanishes exactly at the zero-loss point. Losses are large, not small: with realistic reflectivity the net per pair is about $-1430$ keV, because the fuel itself is lost.

**Starting from:** Energy and momentum conservation. The Breit-Wheeler rate from A16.

**Derivation:**

*Any geometry.* Start from a pair at rest with energy $2m_ec^2$. Work taken leaves photons of energies $E_1,E_2$ at angle $\theta$, with energy conservation $W=2m_ec^2-E_1-E_2$ in the lab frame where the mirrors started at rest. To re-form the pair, $E_1E_2(1-\cos\theta)/2\ge m_e^2c^4$. So $\sqrt{E_1E_2}\ge m_ec^2$, and $E_1+E_2\ge2\sqrt{E_1E_2}\ge2m_ec^2$. Therefore

$$W=2m_ec^2-(E_1+E_2)\le0$$

Equality needs $E_1=E_2=m_ec^2$ and head-on photons. Boosting the photons first with an inward mirror costs at least what it adds.

*Strictly negative.* Expanding to $E_f$ gives $W_{\rm out}=2(m_ec^2-E_f)$. Recompressing to $E_c=\gamma m_ec^2$ costs $W_{\rm in}=2(E_c-E_f)$. So $E_f$ cancels:

$$W_{\rm out}-W_{\rm in}=-2(\gamma-1)m_ec^2$$

That means the piston's work never contributes, however deep the expansion. The pair must carry $2(\gamma-1)m_ec^2$ of kinetic energy, and the pair's kinetic energy is recovered with fraction $\rho\le1$. With photons lost per crossing,

$$\langle W_{\rm net}\rangle=-(1-\varepsilon)\,2E_c-\varepsilon(1-\rho)\,2(\gamma-1)m_ec^2$$

Here $\varepsilon$ is the chance of converting before the photons are lost. A lost pair of photons costs $2E_c$: the recompression work already spent, plus a replacement pair. Now every route to zero closes:

- Taking $\gamma\to1$ makes $\sigma\to\pi r_e^2\beta\to0$, so the crossings per conversion diverge as $1/\beta$. The cycle time goes to infinity.
- $\varepsilon\to1$ needs $1-R\ll P$, and $P\propto\beta$ is falling at the same time.
- $\rho\to1$ needs reversible slowing of a just-made pair. A nucleus's recoil, if used, adds 2.85 eV (gold) or 557 eV (proton) of energy to be recovered.

$W=0$ is reached only at zero power with perfect everything. Any finite rate gives $W<0$. This matches the thermodynamic picture: reversible means infinitely slow, and any finite-rate process produces entropy.

*Carnot comparison.* With baths, $W\le Q_h(1-T_c/T_h)$ and the pair term cancels. For the mirror and photon part alone, the entropy of a 511 keV photon in 1 m³ is about 77 $k_B$, costing $T_aS\approx2$ eV at 300 K. So the second law's price for this closed cycle is about $10^{-6}$ of the energy. That's why the first law, not entropy, sets the headline: the entropy cost is small and real, while the threshold gap is exact.

**Check:** [`code/verify_closed_cycle.py`](code/verify_closed_cycle.py) runs a stochastic cycle, closing the first law on each event, and also compares with the analytic expectation above. [`code/verify_optical_piston.py`](code/verify_optical_piston.py) runs the expansion and recompression bounce by bounce. [`code/verify_thresholds.py`](code/verify_thresholds.py) samples the geometry claim.

```text
   R          beta_pair  recover   P per cross  made/total   net [keV]    ideal bound 0?
   0.29          0.700    0.50      1.78e-05   1.50e-05   -1431.0649     crossings/cycle 1.1e+00
   0.99          0.700    1.00      1.78e-05   9.00e-04   -1429.7953     crossings/cycle 5.0e+01
   0.999999      0.700    1.00      1.78e-05   8.99e-01    -145.1333     crossings/cycle 5.1e+04
   1             0.010    0.99      1.33e-07   1.00e+00      -0.0005     crossings/cycle 7.5e+06
   1             0.010    1.00      1.33e-07   1.00e+00       0.0000     crossings/cycle 7.5e+06
   E_f =  500.0 keV: net = -1430.0243 keV
   E_f =   10.0 keV: net = -1430.0243 keV
   beta_pair   KE per pair [keV]   P per crossing   crossings needed   mean cycle time [s]
     7.0e-01     4.0909e+02           1.78e-05        5.62e+04           3.75e-04
     1.0e-06     5.1104e-10           1.33e-11        7.51e+10           5.01e+02
   target 0.50: W_out = 0.503426, W_in = 0.503426, net = -0.000000
```

A caution about reading this test. The stochastic cycle is a bookkeeping model built on the same first-law assumption it checks, so it can't by itself prove the first law. What it adds is the size of the loss under each assumed $R$, $\beta$ and $\rho$. The independent test of the threshold is the random sampling in A01, and the independent test of the piston is the Lorentz simulation.

**Limits:** The cycle model treats the photons as a pair crossing once per round trip with a diffraction-limited chance $P$. A many-photon cavity differs. The rule $W\le0$ doesn't depend on this. Near threshold the times shown use a 1 m cavity.

### A19: The catalytic energy boost

**Status:** Derived + Code + Established

**Short answer:** This isn't catalysis. A catalyst comes out unchanged, and here the proton and antiproton are consumed. It's burning fuel with 1836 times the mass-energy, and if you remade the $p\bar p$, you'd lose at least the neutrino share, about 830 MeV. The usable budget is also smaller than stated: only about 0.7 GeV of 1.876 GeV appears as photons. One useful bonus: charged pions decay into about 1.5 $e^+e^-$ pairs per annihilation, so no separate pair regeneration is even needed.

**Starting from:** Energy conservation. Measured pion multiplicities in $\bar pp$ at rest: $\langle n_{\pi^\pm}\rangle=3.05\pm0.04$, $\langle n_{\pi^0}\rangle=1.93\pm0.12$, total $4.98\pm0.13$, with kaons in about 5% of events [[26]](#sources). PDG decay modes and masses [[18]](#sources).

**Derivation:**

*The real reaction is $p+\bar p\to$ pions*, then $\pi^0\to2\gamma$, $\pi^\pm\to\mu\nu$, $\mu\to e\nu\nu$. A charged pion at rest gives its $\nu_\mu$ 29.8 MeV. Its muon gives on average 38 MeV to the electron or positron, and the rest to two neutrinos. So neutrinos take about 72.5% of each charged pion's energy.

*The account.* The Monte Carlo below gives, for 1876.5 MeV:

| Channel | Energy per annihilation (MeV) | Fraction | Can it do work? |
|---|---|---|---|
| Neutrinos | about 830 | 44% (43 to 46%) | No, it escapes |
| Photons | about 730 | 39% (37 to 41%) | Only with a mirror or absorber |
| $e^\pm$ rest + kinetic | about 316 | 17% | Yes, then 511 keV pairs |
| **Upper bound on work** | **about 1046** | **56%** | Only without regeneration |

The $e^\pm$ line includes about 1.5 pairs of rest energy, which is a small part. At pion level, charged pions take about 61% of the energy and neutral pions about 39%. Photons average about 180 MeV, with about 4 per annihilation. So the claim of "1 to 2 GeV usable photon energy" is too large by a factor of 1.4 to 2.7.

*Why there's no boost.* Your statement that remaking $e^+e^-$ only costs 1.022 MeV is true, and it's 0.054% of the fuel. But it misses what got consumed: the net reaction is $p+\bar p\to$ work + heat + neutrinos, with no $p$ or $\bar p$ left. Baryon number allows it, because $p\bar p$ has $B=0$.

- **If you don't regenerate $p\bar p$:** this is simply fuel. Each cycle burns one antiproton. The "headroom" you describe is the antiproton's energy, not a gain.
- **If you do regenerate it:** you must pay at least 1876.5 MeV, and have at most 1046 MeV. So net $\le-830$ MeV, before any other loss, and at the fixed-target threshold $6m_pc^2=5630$ MeV the gap is much worse.

**Check:** [`code/verify_annihilation_budget.py`](code/verify_annihilation_budget.py) generates pions by n-body phase space, decays them with explicit four-vectors, and conserves four-momentum per event:

```text
   mode  <n+->  <n0>   charged pions  neutral pions   gamma    nu     e+-    sum    photons/ann  <E_gamma>
   A      3.10  1.97         0.631          0.369    0.369  0.458  0.173  1.0000      3.93      175.9 MeV
   B      2.94  2.06         0.590          0.410    0.410  0.427  0.164  1.0000      4.11      186.8 MeV
   fuel rest energy consumed                    1876.5
   lost to neutrinos                             830.1
   in photons                                    730.1
   in e+- (rest + kinetic)                       316.3
   upper bound on work                          1046.4  (55.8% of fuel)
   e+e- pairs the pion decays leave behind        1.51  per annihilation
   Close the p pbar loop too, at the bare threshold: net = 1046.4 - 1876.5 = -830.1 MeV
```

The two modes use different assumed multiplicity tables, chosen by me to reproduce the measured means. I did not use the measured channel table. The spread between them is a rough model uncertainty.

**Limits:** Phase space ignores resonances such as $\rho$ and $\omega$, which do matter, and kaon events (about 5%) are left out. Muon polarisation is ignored. I found no primary source giving these energy fractions directly, so the fractions here are this model's output, not a citation. The commonly quoted figure of "about half to neutrinos" fits this within the model spread.

### A20: The efficiency ceiling

**Status:** Derived + Code + Established + Estimate

**Short answer:** The formula $\eta_{\max}\approx1-m_ec^2/E_\gamma$ isn't derived from anything in the design, and it's not an open-loop efficiency. It gives 0 for $e^+e^-$ at rest, where nothing stops 100%. What it does correctly describe is a different quantity: the fraction of one photon's energy you can remove while keeping it above the 511 keV regeneration threshold. Your closed-loop statement of 0% is correct, and should read "at most 0%". True round-trip efficiency via real antimatter production is about $10^{-9}$ to $10^{-8}$.

**Starting from:** Energy accounting per unit of fuel. Sourced production figures.

**Derivation:**

*Test the formula.* For $e^+e^-$ at rest, $E_\gamma=m_ec^2$ gives $\eta=0$, yet open-loop fuel isn't remade, so no threshold applies and a perfect mirror could approach $\eta=1$. For $\bar pp$ with mean photon near 185 MeV, the formula gives 0.997, but 44% of the energy has already gone to neutrinos. So it fails in both directions.

*Where it could come from.* Suppose one photon of energy $E_\gamma$ gives up work while remaining above $m_ec^2$. The extractable fraction is $(E_\gamma-m_ec^2)/E_\gamma=1-m_ec^2/E_\gamma$. So the formula has a real meaning, but it's about the photon alone, in the A19 loop that burns $p\bar p$.

*The correct open-loop expression.* Let $E_i$ be the fuel energy, $f_u$ the fraction in channels that can do work, $f_c$ the fraction captured, and $E_{\rm left}/E_u$ the energy remaining after extraction:

$$\eta_{\rm open}=f_uf_c\left(1-\frac{E_{\rm left}}{E_u}\right)$$

For $e^+e^-$, $f_u=1$. For $\bar pp$, $f_u\approx0.56$. Using a heat engine in the last factor gives $1-T_a/T_h$.

*The price of the fuel.* Round-trip efficiency multiplies by production:

$$\eta_{\rm round}=\varepsilon_{\rm prod}\,\eta_{\rm open}$$

Estimating $\varepsilon_{\rm prod}$ from CERN's numbers:

- Input: $1.5\times10^{13}$ protons at 26 GeV/c per pulse [[12]](#sources), which is $6\times10^4$ J of beam.
- Output: about $3$ to $4\times10^7$ antiprotons per cycle [[13]](#sources), with annihilation energy $1.876$ GeV each.
- Beam-to-annihilation energy is $1.5$ to $2\times10^{-7}$. With a wall-plug efficiency of 5%, as assumed by LaPointe and Reddy, this is $0.8$ to $1\times10^{-8}$. Foils trapping below 1% reduce it further, under $10^{-10}$; ELENA aims for 10 to 100 times more captured.
- LaPointe and Reddy give 1 antiproton per $2.5\times10^6$ protons, or $3\times10^{-8}$ at beam level [[27]](#sources), so $1.5\times10^{-9}$ at the wall.

These agree on order $10^{-9}$ to $10^{-8}$. They differ because they use different eras, stages and definitions.

*As energy storage.* Antimatter stores $1.8\times10^{17}$ J per kg of antimatter (with an equal mass of matter). That's the highest known energy density, which is why it attracts interest in propulsion. As a battery, its round trip is at most about $6\times10^{-9}$ for $\bar pp$, against something near 0.9 for a lithium battery (general knowledge). At the AD's rate, one gram takes about $5\times10^{10}$ years. So it's the most compact known store and the worst known source.

**Check:** [`code/verify_annihilation_budget.py`](code/verify_annihilation_budget.py):

```text
        0.511            0.00000    e+e- at rest
      185.666            0.99725    mean pi0 photon
     p pbar  : f_usable <= 0.558, so eta <= 0.558 for any E_gamma
   3e+07 pbar/cycle: annihilation energy / beam energy = 1.5e-07; with 5% wall plug = 7.5e-09
   LaPointe & Reddy figure: 1 pbar per 2.5e6 protons, beam efficiency 2.9e-08 (their 3e-8)
   production 1e-08 x usable 0.56 (ideal conversion) = 5.6e-09
   AD rate 3.5e+05 pbar/s means 1 g takes 5.4e+10 years
```

**Limits:** Order of magnitude only. The 5% wall-plug figure is an assumption from one paper and I saw it through a search summary, because the report's PDF is a scan I couldn't read directly. The 3 to $4\times10^7$ count is delivered antiprotons after losses. This doesn't count facility power beyond the proton beam, which would lower it further. Positron routes have different efficiencies, which I didn't source.

### A21: The claimed critical challenges

**Status:** Derived + Code + Established

**Short answer:** (A) and (E) are right. (C) is right in substance, though perfect alignment isn't needed. (B) has the wrong reason: quantum measurement isn't the obstacle. (D) has a true fact and the wrong emphasis: the demon's cost is real but about a millionth of the energy, so it doesn't decide anything. The real showstopper is the first law, plus the missing mirror.

**Starting from:** The results above.

**Derivation:**

**(A) No material reflects 511 keV efficiently. Mostly correct.** At normal incidence, $R\sim10^{-17}$. Grazing reflection works below about 0.16 mrad, and Laue crystals manage 26 to 31% near 600 keV over small deflections. So it's true for any useful mirror, while limited optics do exist (A14).

**(B) Random direction and quantum measurement. Wrong reason.** Direction is random even for a pair at rest, not just "if particles have momentum". In its centre-of-momentum frame the emission is isotropic for an unpolarised pair. Moving particles only add Doppler shift and aberration (A02). And no knowledge is needed to catch a photon: a closed sphere catches every direction without measuring anything, and a 2.6 cm lead shell attenuates 99% of them. Quantum mechanics also isn't blocking this: a spherical mirror would simply act on the whole spherical wave. The real limit is that no surface reflects at those angles. The measurement problem doesn't enter.

**(C) Low yield, alignment and density. Correct in substance.** Yield is low: a best case of about $\alpha^2$ per crossing for two single photons, or over 20 kJ per pulse for bulk beams. But:

- Perfect alignment isn't needed. With 766 keV photons, the cross section at 120° crossing is 0.99 of its head-on value, and is even slightly larger at 150°. It vanishes only once $s$ is below threshold, below about 84° here.
- Higher energy doesn't always help. The cross section peaks near 717 keV and then falls. At 10 MeV it's 40 times smaller.

**(D) Entropy and Maxwell's demon. True in a small way, but not the obstacle.** Entropy does increase in real steps. Landauer's principle does set a minimum heat of $k_BT\ln2$ per erased bit, which has been experimentally confirmed [[28]](#sources). In numbers, at 300 K:

| Item | Size |
|---|---|
| Landauer, per bit | 0.018 eV |
| Mode entropy $T_aS$, per 511 keV photon in 1 m³ | about 2 eV (111 bits) |
| Mixing $e^+$ and $e^-$, per pair | about 0.04 eV |
| Rest energy of a pair | 1,022,000 eV |

So the demon's bill is about $10^{-6}$ of the energy at most. Sorting $e^+$ from $e^-$ needs no demon at all. A static field does it, because charge is a physical handle and no measurement is required (A17). A demon would be relevant for sorting by unknown microstate at a small energy scale. Even a perfectly free demon wouldn't help, because $W\le0$ comes from the first law.

**(E) Storage. Correct.** Positrons are stored in Penning-Malmberg traps. Numbers near $10^{10}$ at about 0.2 eV have been reported [[29]](#sources), and the Danielson et al. review covers long-term storage [[4]](#sources). That's only $1.6\times10^{-3}$ J of pair energy. Density is limited: the Brillouin limit at 5 T is $1.2\times10^{20}$ m⁻³, corresponding to $2\times10^7$ J/m³, well below petrol at about $3\times10^{10}$ J/m³ (general knowledge). I couldn't confirm a storage lifetime figure from the sources I read.

**Check:** [`code/verify_closed_cycle.py`](code/verify_closed_cycle.py), [`code/verify_breit_wheeler.py`](code/verify_breit_wheeler.py) and [`code/verify_traps.py`](code/verify_traps.py):

```text
   modes in 1 m^3 within dE/E=0.2%: 3.52e+33; S = k ln(modes) = 77 k = 111 bits
   T S at 300 K = 2.00 eV per photon = 3.9e-06 of its energy
   Landauer cost at 300 K = 0.0179 eV per bit
   2 k T ln 2 at 300 K = 0.036 eV of 1.022 MeV = 3.5e-08
            150           0.7236            1.010
            120           0.6383            0.990
             90           0.3333            0.536
             80           0.0000            0.000
          10000        0.9987       0.0041
   Brillouin limit at 5.0 T: n = 1.21e+20 1/m^3
   stored rest energy if all paired at that density: 1.99e+07 J/m^3 (petrol ~3e10 J/m^3, general knowledge)
```

The entropy figure depends only logarithmically on the assumed 1 m³ and 0.2% bandwidth, so it's insensitive to both.

**Limits:** The entropy count treats one photon in a set of modes, as a scale. The trap values are single-source numbers.

## Whole-cycle energy account

All values per stored pair, in keV. The "ideal limit" column assumes the fictional perfect mirror, perfect recovery, and $\beta\to0$. The "realistic best" column uses $R=0.99$ (far above anything real at these angles), $\beta=0.7$, with full recovery of pair kinetic energy, from the [run shown in A18](#a18-the-ideal-energy-flow).

| Step | Energy flow | Ideal limit | Realistic best | Law or source |
|---|---|---|---|---|
| Start: stored pair | rest energy | +1022 | +1022 | $E=mc^2$ |
| Annihilation | into photons | 1022 in photons | 1022 in photons (2γ) | Conservation; triplet Ps gives 3γ |
| Expansion on piston | $W_{\rm out}=2(511-E_f)$ | any, up to +1022 | any, up to +1022 | Doppler, A14 |
| Recompression | $W_{\rm in}=2(E_c-E_f)$ | equal to $W_{\rm out}$ | $W_{\rm out}$ + 409 | Threshold, A16, A18 |
| **Piston subtotal** | independent of $E_f$ | **0** | **−409** | First law |
| Pair created | photons → pair + KE | no KE, rate → 0 | +409 KE, 1 in 1100 succeed | Breit-Wheeler, A16 |
| KE recovered | from slowing pair | 0 | up to +409 | Assumed ideal |
| Photons lost in mirror | into heat | 0 | −1430 on most cycles | Cavity race, A18 |
| Replacement pair | buy new fuel | 0 | counted in the line above | First law, A20 |
| Recoil, if nucleus used | into nucleus | +2.85 eV of cost (Au) | | A16 |
| Sorting | field does it | about 0 | about 0.04 eV entropy price | A17 |
| Entropy of photon modes | at 300 K | about 2 eV per photon as minimum | | A18, A21 |
| **Net per pair** | | **0, at zero power** | **−1430** | **Always ≤ 0** |

And for comparison, using annihilation as fuel in an open cycle:

| Fuel use | Energy out per energy of fuel | Energy out per energy spent making it |
|---|---|---|
| $e^+e^-$, absorbed, 1500 K to 300 K | 0.80 | Not sourced for positrons |
| $\bar pp$, absorbed, 1500 K to 300 K | about 0.45 | about $4.5\times10^{-9}$ at $\varepsilon=10^{-8}$ |
| $\bar pp$, ideal direct conversion of all non-neutrino energy | at most 0.56 | about $6\times10^{-9}$ at $\varepsilon=10^{-8}$ |

The picokelvin sink appears in neither table because it cancels (A06). The $p\bar p$ boost isn't a line in the first table because it's a different, one-way fuel (A19).

## Sources

Every entry below was located online while writing this section, except where noted. Items marked "abstract only" were read through their abstract or a search summary, and I did not read the full text.

1. AMS Collaboration, [Antiproton Flux, Antiproton-to-Proton Flux Ratio, and Properties of Elementary Particle Fluxes](https://link.aps.org/doi/10.1103/PhysRevLett.117.091103), Phys. Rev. Lett. 117, 091103 (2016). Abstract only.
2. AMS Collaboration, [High Statistics Measurement of the Positron Fraction](https://arxiv.org/abs/1612.09579) and CERN's [summary of the measurements](https://home.cern/latest-measurements-from-the-ams-experiment-unveil-new-territories-in-the-flux-of-cosmic-rays/). Abstract only.
3. Adriani et al. (PAMELA), [The discovery of geomagnetically trapped cosmic ray antiprotons](https://arxiv.org/abs/1107.4882), Astrophys. J. Lett. 737, L29 (2011). Abstract only.
4. Danielson, Dubin, Greaves and Surko, [Plasma and trap-based techniques for science with positrons](https://link.aps.org/accepted/10.1103/RevModPhys.87.247), Rev. Mod. Phys. 87, 247 (2015). Table of contents and abstract only.
5. Enoto et al., [Photonuclear reactions triggered by lightning discharge](https://arxiv.org/abs/1711.08044), Nature 551, 481 (2017). Abstract only.
6. Prantzos et al., [The 511 keV emission from positron annihilation in the Galaxy](https://arxiv.org/abs/1009.4620), Rev. Mod. Phys. 83, 1001 (2011).
7. Cohen, De Rújula and Glashow, [A Matter-Antimatter Universe?](https://arxiv.org/abs/astro-ph/9707087), Astrophys. J. 495, 539 (1998). Abstract only.
8. Burke et al., [Positron Production in Multiphoton Light-by-Light Scattering](https://www.lle.rochester.edu/media/timeline/documents/1997/PhysRevLett.79.1626.pdf), Phys. Rev. Lett. 79, 1626 (1997). Abstract and summary.
9. Melissinos, [The Spontaneous Breakdown of the Vacuum](https://arxiv.org/abs/hep-ph/9805507), arXiv:hep-ph/9805507. Summary only.
10. STAR Collaboration, [Measurement of $e^+e^-$ Momentum and Angular Distributions from Linearly Polarized Photon Collisions](https://arxiv.org/abs/1910.12400), Phys. Rev. Lett. 127, 052302 (2021).
11. Pike, Mackenroth, Hill and Rose, A photon-photon collider in a vacuum hohlraum, Nature Photonics (2014), doi:10.1038/nphoton.2014.95. Through [Imperial College's summary](https://www.imperial.ac.uk/news/149087/scientists-discover-turn-light-into-matter/) and other news reports.
12. Torregrosa Martin et al., [Towards the last stages of the CERN's AD-target area consolidation project](https://proceedings.jacow.org/ipac2021/papers/WEPAB366.pdf), IPAC2021 WEPAB366.
13. Oelert, [The ELENA project at CERN](https://arxiv.org/abs/1501.05728), arXiv:1501.05728.
14. Gajos, [Investigations of fundamental symmetries with the electron-positron systems](https://arxiv.org/abs/1805.06009), arXiv:1805.06009, for the lifetimes of 125 ps and 142 ns. Summary only.
15. Deppner et al., [Collective-Mode Enhanced Matter-Wave Optics](https://physics.aps.org/featured-article-pdf/10.1103/PhysRevLett.127.100401), Phys. Rev. Lett. 127, 100401 (2021).
16. APS Physics Viewpoint on Deppner et al., [Physics 14, 119](https://physics.aps.org/articles/v14/119) (2021).
17. Fixsen, [The Temperature of the Cosmic Microwave Background](https://arxiv.org/abs/0911.1955), Astrophys. J. 707, 916 (2009).
18. Particle Data Group, [pion listing](https://pdglive.lbl.gov/DataBlock.action?node=S008T). Pion lifetime confirmed here. The muon lifetime of 2.197 µs is from general knowledge and I didn't confirm it online.
19. Gerrish and Schmidt, [Antimatter Production for Near-Term Propulsion Applications](https://ntrs.nasa.gov/citations/19990080056), NASA (1999). Abstract only.
20. Blaum, Novikov and Werth, [Penning traps as a versatile tool for precise experiments in fundamental physics](https://arxiv.org/abs/0909.1095), Contemp. Phys. 51, 149 (2010). Summary only. The trap theory traces to Brown and Gabrielse, Rev. Mod. Phys. 58, 233 (1986), which I did not read.
21. Gabrielse et al., [The ingredients of cold antihydrogen: Simultaneous confinement of antiprotons and positrons at 4 K](https://cfp.physics.northwestern.edu/documents/Ingredients.pdf), Phys. Lett. B 455, 311 (1999).
22. Brejnholt et al., Demonstration of multilayer reflective optics at photon energies above 0.6 MeV, Opt. Express 22, 15364 (2014), doi:10.1364/OE.22.015364, through [Phys.org's report](https://phys.org/news/2014-09-energy-multilayer-coated-mirrors.html) and the publisher's summary.
23. Barrière et al., [Experimental and theoretical study of the diffraction properties of various crystals for the realization of a soft gamma-ray Laue lens](https://arxiv.org/abs/0907.0458), J. Appl. Cryst. 42, 834 (2009).
24. Frontera and von Ballmoos, [Laue gamma-ray lenses for space astrophysics: status and prospects](https://arxiv.org/abs/1007.4308), arXiv:1007.4308.
25. Breit and Wheeler, Phys. Rev. 46, 1087 (1934); Gould and Schréder, Phys. Rev. 155, 1404 (1967). The formula as cited in [arXiv:0807.3135](https://arxiv.org/abs/0807.3135). I did not read the originals. The formula is verified numerically here against the QED amplitude.
26. Klempt, Batty and Richard, [The antinucleon-nucleon interaction at low energy: annihilation dynamics](https://arxiv.org/abs/hep-ex/0501020), Phys. Rep. 413, 197 (2005), eq. 5.3.
27. LaPointe and Reddy, [Antimatter Production at a Potential Boundary](https://ntrs.nasa.gov/citations/20010091016), NASA/CR-2001-211116 (2001). The efficiency figures were seen through a search summary. The report's body is a scan I couldn't read directly. The report proposes a speculative production method that I don't rely on.
28. Bérut et al., [Experimental verification of Landauer's principle](https://www.physics.rutgers.edu/grad/677/Physics_677_2023_files/Berut_Lutz_Nature2012.pdf), Nature 483, 187 (2012); Landauer, IBM J. Res. Dev. 5, 183 (1961), which is general knowledge.
29. Surko group, [Phys. Plasmas 13, 123502](https://marvin.byu.edu/beryllium7/references/SourceDesign/PhysPlasmas-13-123502.pdf) (2006). Summary only.

NIST mass attenuation coefficients are from the tables for [gold](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z79.html) and [lead](https://physics.nist.gov/PhysRefData/XrayMassCoef/ElemTab/z82.html).

## Full script output

Unedited output of every script from the last run. Each script exits non-zero if a check fails.

<details><summary><code>code/verify_thresholds.py</code></summary>

```text
1. Fixed-target threshold from direct minimisation (units of the particle mass)
   beam energy 7.000000 m c^2, beam kinetic 6.000000 m c^2, total lab energy 8.000000 m c^2
   final-state speeds at the minimum: [0.86603 0.86603 0.86603 0.86603]
[PASS] p p -> p p p pbar threshold kinetic energy is 6 m_p c^2  (5.630 GeV)
[PASS] total lab energy is 8 m c^2 with beam 7 m c^2 + target 1 m c^2
[PASS] all four final particles share beta = sqrt(3)/2 at threshold  (beta=0.86603)

2. Photon on a nucleus: solve the conservation laws by root finding (no s formula)
   proton         : k = 1.022554499 MeV, excess above 2 m_e c^2 = 556.597 eV (2m^2/M = 556.597 eV)
[PASS] recoil excess on proton matches 2 m^2 c^2 / M
   Au-197 nucleus : k = 1.022000748 MeV, excess above 2 m_e c^2 = 2.847 eV (2m^2/M = 2.847 eV)
[PASS] recoil excess on Au-197 nucleus matches 2 m^2 c^2 / M

3. A single free photon cannot become a free pair (A03)
   1000000 random pairs: min (E^2-P^2)/(mc^2)^2 = 4.000000 (a photon needs 0)
[PASS] every pair has invariant mass >= 2m, a photon has 0
   min (E - P c)/(m c^2) over samples: 8.800e-03  (always > 0, so E = Pc never holds)
[PASS] E - Pc stays strictly positive

4. Two photons: threshold on s means total energy >= 2 m c^2 in every lab frame (A16, A18)
   865510 of 2000000 random photon pairs are above threshold; least E1+E2 among them = 2.03363 m c^2
[PASS] no above-threshold photon pair has E1 + E2 < 2 m c^2
[PASS] sub-threshold redshifted head-on pair fails  (E1=E2=0.9 m c^2)
[PASS] asymmetric head-on pair passes when E1 E2 >= m^2 c^4  (E1=0.5, E2=2.0 m c^2, total 2.5 m c^2)

5. Kinetic energy adds to the bill (A02): pair made moving with Lorentz factor g
     beta      gamma    energy/2mc^2   forward photon/mc^2   backward photon/mc^2
     0.000     1.000        1.000              1.000                1.0000
[PASS] in-flight annihilation returns exactly 2 g m c^2 at beta=0
     0.500     1.155        1.155              1.732                0.5774
[PASS] in-flight annihilation returns exactly 2 g m c^2 at beta=0.5
     0.866     2.000        2.000              3.732                0.2680
[PASS] in-flight annihilation returns exactly 2 g m c^2 at beta=0.866
     0.990     7.089        7.089             14.107                0.0709
[PASS] in-flight annihilation returns exactly 2 g m c^2 at beta=0.99
     0.999    22.366       22.366             44.710                0.0224
[PASS] in-flight annihilation returns exactly 2 g m c^2 at beta=0.999

6. How hard a fast antiproton is to hold (A02)
   p [MeV/c]    KE [MeV]     beta    B rho [T m]   KE / 10 kV well
      3500.0    2685.311    0.9659       11.675      2.685e+05
       100.0       5.314    0.1060        0.334      5.314e+02
        13.7       0.100    0.0146        0.046      1.000e+01
[PASS] 100 MeV/c antiproton has the AD's quoted 5.3 MeV  (5.314 MeV)

ALL CHECKS PASSED
```

</details>

<details><summary><code>code/verify_carnot_sink.py</code></summary>

```text
1. The Carnot number (A06, A11)
   T_c/T_h = 2.5333e-14, so eta_C = 1 - 2.533e-14
[PASS] 1 - eta_C is 2.5e-14, about 0.025 pico rather than 1 pico
   eta_C with T_c = 38 pK        : 1 - 2.533e-14 = 1.000000
   eta_C with T_c = CMB 2.725 K  : 1 - 1.817e-03 = 0.998183
   eta_C with T_c = 300 K ambient: 1 - 2.000e-01 = 0.800000

2. Size of the reservoir (ideal-gas value, generous for a condensate)
   C = 3/2 N k_B = 2.071e-18 J/K, thermal content C T0 = 7.870e-29 J
   pair energy 2 m c^2 = 1.637e-13 J = 2.08e+15 times the thermal content
   Carnot waste per pair at fixed T_c: 4.148e-27 J = 52.7 times the thermal content
[PASS] the waste heat of one pair exceeds the sink's whole thermal content
   one-axis rms speed at 38 pK: 60 um/s (the APS Viewpoint quotes about 60 um/s)
[PASS] 38 pK corresponds to the reported 60 um/s expansion
   Compton optical depth of the cloud for 511 keV (37 electrons per atom, radius 150 um): 1.5e-15
[PASS] a 511 keV photon goes straight through the cloud

3. Engine stepped against the warming sink, N = 1e5 atoms, offered Q_h = 2 m c^2
   heat accepted 9.725e-14 J (59.4% of the pair energy) before T_sink = T_h
   work 9.415e-14 J, rejected 3.107e-15 J, final sink temperature 1500.1 K
   mean efficiency over the stroke 0.9681
   closed form: Q_max = C T_h ln(T_h/T_0) = 9.725e-14 J, mean eta = 0.9681
[PASS] stepped heat limit matches C T_h ln(T_h/T_0)
[PASS] stepped mean efficiency matches 1 - 1/ln(T_h/T_0)
[PASS] first law holds in the stepped run
[PASS] one pair cannot be fully absorbed by this sink  (0.59 of the pair)
   instantaneous eta falls from 1 - 2.5e-14 to 0 during the stroke.

   Efficiency against heat handled, same sink
     Q_h [J]      T_final [K]     mean eta      1 - mean eta
     1.0e-22      3.800e-11     1.000000      2.533e-14
     1.0e-18      3.801e-11     1.000000      2.534e-14
     1.0e-16      3.924e-11     1.000000      2.575e-14
     1.0e-15      5.243e-11     1.000000      2.989e-14
     1.0e-14      9.502e-10     1.000000      1.889e-13
     5.0e-14      3.715e-04     1.000000      1.539e-08
[PASS] small heat loads recover the textbook Carnot loss  (loss 2.534e-14 at Q_h=1e-18 J)

   atoms needed so one pair only doubles T_c: 5.27e+06
   Rb needed at 38 pK to sink the waste of 1 kg of annihilated mass the same way: 4.18e+11 kg

4. Paying for the cold sink: T_c cancels exactly (A06)
   sympy: (W_engine - W_fridge)/Q_h = (-T_a + T_h)/T_h
[PASS] net efficiency simplifies to 1 - T_a/T_h with T_c absent
   fridge must exceed this fraction of its Carnot COP to net anything: (-T_a + T_c)/(T_c - T_h)
   numerically 0.2000 for T_a = 300 K, T_h = 1500 K
[PASS] break-even fridge quality is 20% of Carnot
   Carnot COP from 38 pK to 300 K = 1.267e-13

5. Whole-life ledger with a finite sink: prepare, run, return to ambient
   ideal least work to cool the sink from 300 K to 38 pK: 1.783e-14 J = 111 keV
     Q_h [J]    W_engine [J]   W_prep [J]   W_recovered [J]   net/Q_h   1-T_a/T_h
     1.0e-16   1.0000e-16   1.7829e-14      1.7809e-14     0.80000     0.800
[PASS] lifetime net efficiency is 1 - T_a/T_h at Q_h=1e-16
     1.0e-15   1.0000e-15   1.7829e-14      1.7629e-14     0.80000     0.800
[PASS] lifetime net efficiency is 1 - T_a/T_h at Q_h=1e-15
     1.0e-14   1.0000e-14   1.7829e-14      1.5829e-14     0.80000     0.800
[PASS] lifetime net efficiency is 1 - T_a/T_h at Q_h=1e-14
     5.0e-14   5.0000e-14   1.7829e-14      7.8294e-15     0.80000     0.800
[PASS] lifetime net efficiency is 1 - T_a/T_h at Q_h=5e-14
     9.0e-14   8.9699e-14   1.7829e-14      1.3016e-16     0.80000     0.800
[PASS] lifetime net efficiency is 1 - T_a/T_h at Q_h=9e-14
   assumed real apparatus input 1 kW x 10 s = 1e+04 J per shot, which is 6e+17 times the ideal 1.8e-14 J

6. Heat-rejection rate into a dilute gas: Monte Carlo against kinetic theory
   mean speed 96 um/s at 38 pK; cloud radius assumed 150 um
   kinetic theory: P_c = N (3 v/4r) 2 k (T_wall - T_gas) = 5.05e-29 W
   Monte Carlo (200000 atoms, 12182 wall hits): P_c = 5.08e-29 W
[PASS] Monte Carlo heat rate agrees with kinetic theory within 10%  (ratio 1.007)
   work this allows with an ideal engine: P_c (T_h/T_c - 1) = 1.99e-15 W
   time to dump one pair's Carnot waste: 82 s; drop-tower free fall is 4.74 s
[PASS] the waste of one pair takes longer than the drop
   gas total thermal energy / rate = 1.56 s, the sink's own turnover time

7. Radiative cooling of the 38 pK side
   blackbody flux at 38 pK: 1.18e-49 W/m^2; over the cloud 3.34e-56 W
   Wien peak wavelength: 7.63e+07 m, compared with Earth's radius 6.4e6 m
[PASS] blackbody wavelength at 38 pK is larger than the Earth

ALL CHECKS PASSED
```

</details>

<details><summary><code>code/verify_optical_piston.py</code></summary>

```text
1. Single reflection off a mirror receding at fixed speed (the load holds its speed)
    beta     E'/E (sim)   (1-b)/(1+b)    work/E
   0.000001   0.99999800   0.99999800 0.00000200
[PASS] Lorentz reflection gives (1-b)/(1+b) at beta=1e-06
   0.010000   0.98019802   0.98019802 0.01980198
[PASS] Lorentz reflection gives (1-b)/(1+b) at beta=0.01
   0.100000   0.81818182   0.81818182 0.18181818
[PASS] Lorentz reflection gives (1-b)/(1+b) at beta=0.1
   0.500000   0.33333333   0.33333333 0.66666667
[PASS] Lorentz reflection gives (1-b)/(1+b) at beta=0.5
   0.900000   0.05263158   0.05263158 0.94736842
[PASS] Lorentz reflection gives (1-b)/(1+b) at beta=0.9
   0.990000   0.00502513   0.00502513 0.99497487
[PASS] Lorentz reflection gives (1-b)/(1+b) at beta=0.99
   0.999999   0.00000050   0.00000050 0.99999950
[PASS] Lorentz reflection gives (1-b)/(1+b) at beta=0.999999
[PASS] an approaching mirror blueshifts, so work is paid in  (E'/E=1.2222)
[PASS] no finite beta < 1 removes all the energy

2. Audit of one free-mirror reflection against lab-frame conservation laws
   energy before 5.061862178479, after 5.061862178479
   momentum before 2.612372435696, after 2.612372435696
[PASS] energy conserved
[PASS] momentum conserved
[PASS] direct lab-frame solve agrees with the transform  (E'=0.638367177)

3. One bounce off a free mirror at rest: E' = E/(1 + 2E/Mc^2), and Compton as a check
   mirror = one free electron: E' = 170.333 keV; Compton 180 deg formula 170.333 keV
   the electron leaves with beta = 0.8000, carrying 340.7 keV
[PASS] a free electron used as a mirror reproduces Compton backscatter
[PASS] backscatter leaves one third of 511 keV
   E/Mc^2 =    1e-06: fraction extracted in one bounce = 0.000002, mirror beta = 0.00000200
[PASS] mirror stays below c at E/Mc^2=1e-06
   E/Mc^2 =    1e-02: fraction extracted in one bounce = 0.019608, mirror beta = 0.01980004
[PASS] mirror stays below c at E/Mc^2=1e-02
   E/Mc^2 =    1e+00: fraction extracted in one bounce = 0.666667, mirror beta = 0.80000000
[PASS] mirror stays below c at E/Mc^2=1e+00
   E/Mc^2 =    1e+02: fraction extracted in one bounce = 0.995025, mirror beta = 0.99995050
[PASS] mirror stays below c at E/Mc^2=1e+02
   E/Mc^2 =    1e+06: fraction extracted in one bounce = 1.000000, mirror beta = 1.00000000
[PASS] mirror stays below c at E/Mc^2=1e+06

4. Optical piston with a clamped mirror receding at fixed beta: event-driven run
   (left mirror fixed at x=0, right mirror starts at L0=1, c=1)
    beta   bounces  E_final/E0   W_total/E0   sum-1      E L/(first E L)  time
   0.001     2303     0.00999     0.990008   6.7e-16   1.00000000   98983.19
[PASS] W + E_photon = E0 at beta=0.001
[PASS] E L is exactly invariant at bounce instants at beta=0.001
   0.010      231     0.00985     0.990149  -1.1e-16   1.00000000    9950.46
[PASS] W + E_photon = E0 at beta=0.01
[PASS] E L is exactly invariant at bounce instants at beta=0.01
   0.100       23     0.00990     0.990102   2.2e-16   1.00000000     908.46
[PASS] W + E_photon = E0 at beta=0.1
[PASS] E L is exactly invariant at bounce instants at beta=0.1
   0.500        5     0.00412     0.995885  -1.1e-16   1.00000000     322.00
[PASS] W + E_photon = E0 at beta=0.5
[PASS] E L is exactly invariant at bounce instants at beta=0.5

5. Two free mirrors with no load: photon energy flows into mirror kinetic energy
   E0/Mc^2    bounces(to L=100 L0)  E_f/E0     sum check   beta_mirror   sqrt(E0/M)   N*sqrt(eps)
      1e-04                300     0.009868   -1.3e-15    9.9995e-03  1.0000e-02     3.000
      1e-05                948     0.009909    6.7e-16    3.1515e-03  3.1623e-03     2.998
      1e-06               2994     0.009985    1.0e-14    9.9549e-04  1.0000e-03     2.994
      1e-07               9466     0.009996    1.5e-13    3.1469e-04  3.1623e-04     2.993
      1e-08              29933     0.009998    1.5e-12    9.9504e-05  1.0000e-04     2.993
[PASS] energy ledger closes at eps=1e-04  (err 1.3e-15)
[PASS] adiabatic invariant E L holds at eps=1e-04  (E L=0.9901)
[PASS] mirror speed tends to sqrt(E0/M) at eps=1e-04
[PASS] energy ledger closes at eps=1e-05  (err 6.7e-16)
[PASS] adiabatic invariant E L holds at eps=1e-05  (E L=0.9969)
[PASS] mirror speed tends to sqrt(E0/M) at eps=1e-05
[PASS] energy ledger closes at eps=1e-06  (err 1.0e-14)
[PASS] adiabatic invariant E L holds at eps=1e-06  (E L=0.9990)
[PASS] mirror speed tends to sqrt(E0/M) at eps=1e-06
[PASS] energy ledger closes at eps=1e-07  (err 1.5e-13)
[PASS] adiabatic invariant E L holds at eps=1e-07  (E L=0.9997)
[PASS] mirror speed tends to sqrt(E0/M) at eps=1e-07
[PASS] energy ledger closes at eps=1e-08  (err 1.5e-12)
[PASS] adiabatic invariant E L holds at eps=1e-08  (E L=0.9999)
[PASS] mirror speed tends to sqrt(E0/M) at eps=1e-08
   fitted N ~ eps^-0.500, expected eps^-0.5
[PASS] bounce count scales as sqrt(M/E0)

6. Extrapolating to a real pair of mirrors, one 511 keV photon (c restored)
   measured coefficient: N = 3.00 * sqrt(Mc^2/E0) to reach L = 100 L0
   mirror mass    bounces N      needed 1-R (for survival 1/e)   final mirror speed
       1e-03 kg    9.93e+13              1.01e-14              9.05e-06 m/s
       1e-09 kg    9.93e+10              1.01e-11              9.05e-03 m/s
       1e-15 kg    9.93e+07              1.01e-08              9.05e+00 m/s
[PASS] a 1 g mirror pair needs over 1e13 lossless reflections per photon  (N=9.9e+13)

7. Closed cycle in the piston alone: expand, then recompress to threshold (A18)
   E_final/mc^2 per photon: W_out and the W_in needed to climb back (per photon, units mc^2)
   target 0.99: W_out = 0.019802, W_in = 0.019802, net = -0.000000
[PASS] returning to threshold costs at least the work taken (target 0.99)
   target 0.50: W_out = 0.503426, W_in = 0.503426, net = -0.000000
[PASS] returning to threshold costs at least the work taken (target 0.5)
   target 0.10: W_out = 0.901734, W_in = 0.901734, net = -0.000000
[PASS] returning to threshold costs at least the work taken (target 0.1)

ALL CHECKS PASSED
```

</details>

<details><summary><code>code/verify_breit_wheeler.py</code></summary>

```text
1. Three independent routes to the cross sections
    beta     sigma_BW [b]   amplitude [b]   sigma_ann [b]   amplitude [b]   ratio/2beta^2
    0.01       0.00249        0.00249        12.4746        12.4746   1.00000000
[PASS] closed BW formula equals amplitude integral at beta=0.01
[PASS] Dirac formula equals amplitude integral at beta=0.01
[PASS] detailed balance sigma_pair/sigma_ann = 2 beta^2 at beta=0.01
    0.10       0.02519        0.02519         1.2595         1.2595   1.00000000
[PASS] closed BW formula equals amplitude integral at beta=0.1
[PASS] Dirac formula equals amplitude integral at beta=0.1
[PASS] detailed balance sigma_pair/sigma_ann = 2 beta^2 at beta=0.1
    0.30       0.08015        0.08015         0.4453         0.4453   1.00000000
[PASS] closed BW formula equals amplitude integral at beta=0.3
[PASS] Dirac formula equals amplitude integral at beta=0.3
[PASS] detailed balance sigma_pair/sigma_ann = 2 beta^2 at beta=0.3
    0.50       0.13819        0.13819         0.2764         0.2764   1.00000000
[PASS] closed BW formula equals amplitude integral at beta=0.5
[PASS] Dirac formula equals amplitude integral at beta=0.5
[PASS] detailed balance sigma_pair/sigma_ann = 2 beta^2 at beta=0.5
    0.70       0.17006        0.17006         0.1735         0.1735   1.00000000
[PASS] closed BW formula equals amplitude integral at beta=0.7
[PASS] Dirac formula equals amplitude integral at beta=0.7
[PASS] detailed balance sigma_pair/sigma_ann = 2 beta^2 at beta=0.7
    0.90       0.11280        0.11280         0.0696         0.0696   1.00000000
[PASS] closed BW formula equals amplitude integral at beta=0.9
[PASS] Dirac formula equals amplitude integral at beta=0.9
[PASS] detailed balance sigma_pair/sigma_ann = 2 beta^2 at beta=0.9
    0.99       0.02178        0.02178         0.0111         0.0111   1.00000000
[PASS] closed BW formula equals amplitude integral at beta=0.99
[PASS] Dirac formula equals amplitude integral at beta=0.99
[PASS] detailed balance sigma_pair/sigma_ann = 2 beta^2 at beta=0.99

2. Threshold behaviour
   beta=1e-01: sigma_BW/(pi r_e^2 beta) = 1.009760
   beta=1e-02: sigma_BW/(pi r_e^2 beta) = 1.000100
   beta=1e-03: sigma_BW/(pi r_e^2 beta) = 1.000001
   beta=1e-04: sigma_BW/(pi r_e^2 beta) = 1.000000
[PASS] sigma_BW -> pi r_e^2 beta, which is zero at threshold
[PASS] sigma_BW vanishes at beta=0
   slow annihilation rate coefficient sigma v = 7.479e-21 m^3/s, pi r_e^2 c = 7.479e-21 m^3/s (no Coulomb focusing)
[PASS] sigma v -> pi r_e^2 c at low speed
   peak: beta = 0.7013, sigma = 1.701e-29 m^2 = 0.2556 sigma_T, head-on photons of 717 keV each
[PASS] peak is about 0.26 sigma_T near beta 0.7
   E_gamma head-on [keV]    beta      sigma [b]
            511        0.0020       0.0005
            512        0.0625       0.0157
            520        0.1853       0.0477
            600        0.5241       0.1443
            722        0.7065       0.1700
           1000        0.8596       0.1359
           3000        0.9854       0.0292
          10000        0.9987       0.0041

3. Alignment: must the photons be perfectly head-on? (A21 C)
   each photon 1.5 m c^2 = 766 keV
   crossing angle [deg]    beta_cm    sigma/sigma(head-on)
            180           0.7454            1.000
            150           0.7236            1.010
            120           0.6383            0.990
             90           0.3333            0.536
             80           0.0000            0.000
[PASS] 90 degree crossing still works for 1.5 m c^2 photons

4. Chance that a single photon pair converts per crossing, best possible focus
   (assumes Gaussian waist w0 = lambda/pi, area A = lambda^2/pi; beam geometry is an estimate)
    E_gamma[keV]   beta      lambda[pm]   P = sigma/A
          512   0.0625        2.422     8.39e-07
          520   0.1853        2.384     2.63e-06
          600   0.5241        2.066     1.06e-05
          722   0.7065        1.717     1.81e-05
         1000   0.8596        1.240     2.78e-05
         3000   0.9854        0.413     5.38e-05
   largest P in this list 5.38e-05; alpha^2 = 5.33e-05
[PASS] single-pair conversion is below 1e-3 per crossing even at best focus

5. Conversion before the photons are lost: eps = P/(P + 1 - R^2) for a cavity
   (P at 722 keV, the cross-section peak; R is reflectivity per bounce, assumed)
         R         eps       photons lost per pair made
    0.290000000   1.98e-05       5.056e+04
    0.990000000   9.09e-04       1.099e+03
    0.999999000   9.01e-01       1.104e-01
    0.999999999   1.00e+00       1.104e-04
[PASS] even R=0.999999 converts only a part of the photons

6. Bulk photon beams: optical depth tau = N sigma / A
   1 um radius spot, sigma at peak: tau = 1 needs N = 1.85e+17 photons per pulse = 2.14e+04 J per pulse of 722 keV gammas
[PASS] optical depth 1 needs >1 kJ of MeV gammas in a micron spot

7. Static fields supply momentum, not energy; a pair from the field alone needs Schwinger's scale
   E_S = m^2 c^3/(e hbar) = 1.32e+18 V/m; work over one reduced Compton wavelength = 511 keV
[PASS] work over a reduced Compton wavelength at E_S equals m c^2

ALL CHECKS PASSED
```

</details>

<details><summary><code>code/verify_gamma_mirror.py</code></summary>

```text
511 keV photon: lambda = 2.4263 pm
NIST-interpolated mu/rho at 511 keV: Au 0.1485, Pb 0.1565 cm^2/g

1. Refractive decrement and critical angle (free-electron scattering)
   element   n_e [1/m^3]    delta       theta_c [mrad]   theta_c [arcmin]
   C         1.056e+30    2.788e-09       0.0747           0.257
   Si        6.994e+29    1.847e-09       0.0608           0.209
   W         4.678e+30    1.235e-08       0.1572           0.540
   Au        4.667e+30    1.232e-08       0.1570           0.540
   Pb        2.703e+30    7.136e-09       0.1195           0.411
[PASS] gold critical angle at 511 keV is below 0.2 mrad  (0.157 mrad)

2. Fresnel reflectivity for gold, exact complex index
   delta = 1.232e-08, beta = 5.541e-11, beta/delta = 4.497e-03
   theta [mrad]   R (smooth)   R (sigma = 0.3 nm, Nevot-Croce)
       0.0157     0.999096     0.999083
       0.0471     0.997175     0.997133
       0.0785     0.994821     0.994744
       0.1256     0.988079     0.987903
       0.1570     0.827461     0.825130
       0.2355     0.021287     0.019265
       0.4709     0.000867     0.000523
       1.5698     0.000006     0.000000
   normal incidence: R = 3.795e-17, (delta/2)^2 = 3.795e-17
[PASS] normal-incidence reflectivity follows (delta/2)^2
[PASS] no flat metal mirror reflects 511 keV at normal incidence above 1e-15
[PASS] above 10 theta_c gold reflects less than 1e-4
   Note: at this energy the real and imaginary parts come from incoherent processes too, so this is the coherent specular part only.

3. A closed sphere cannot use grazing incidence (A14)
   C: photons in a cone of half-angle theta_c around one axis: 1.39e-09 of 4 pi
   Au: photons in a cone of half-angle theta_c around one axis: 6.16e-09 of 4 pi
[PASS] grazing-capable solid angle is below 1e-7 of the sphere for Au
   Mirror placed around the source: radial photons meet a concentric sphere at normal incidence.

4. Photon storage ring from many grazing bounces (A15)
   Bending by 2 pi needs N = pi/theta bounces at grazing angle theta.
    theta/theta_c      N       -ln(survival per lap)   predicted 4 pi beta/theta_c^3
         0.5000   4.00e+04           207.85               180.01
         0.1000   2.00e+05           180.91               180.01
         0.0100   2.00e+06           180.01               180.01
         0.0010   2.00e+07           180.01               180.01
         0.0001   2.00e+08           180.01               180.01
   Survival per lap in the smooth limit: exp(-180.0) = 10^-78.2
[PASS] per-lap loss tends to 4 pi beta/theta_c^3, independent of ring size
[PASS] even a perfectly smooth gold gallery loses the photon in one lap

5. Crystal diffraction
   Si(111), d = 3.136 A: Bragg angle 3.869 mrad = 0.2217 deg
   Darwin width estimate 0.51 urad = 0.104 arcsec, deflection only 2 theta_B
[PASS] Si(111) Bragg angle at 511 keV is under 0.25 degree
[PASS] Si(111) accepts under 1 arcsec
   true backscatter needs d = lambda/2 = 1.21 pm, far below any lattice spacing
   backscatter momentum q = 5.179e+12 1/m; 1s form factor per core electron 1.07e-04
   thermal amplitude factor exp(-q^2 u^2 / 2) = 10^-328
[PASS] backscatter Debye-Waller suppression is below 1e-100
   Measured Laue peak reflectivity: Au(111) 29% at 588 keV (Barriere et al. 2009)
   Au(111) at 511 keV: theta_B = 5.151 mrad, so 610 diffractions per 2 pi turn; 0.29^610 = 10^-328
[PASS] a Laue-crystal ring would lose the photon in a few steps

6. Compton scattering is what the photon actually does in matter
   Klein-Nishina at 511 keV: sigma = 0.2865 b per electron (closed form 0.2865 b), mean energy lost per scatter 34.4%
[PASS] numerical Klein-Nishina integral matches the closed form
   Compton share for Au: 0.0692 cm^2/g of NIST total 0.1485 cm^2/g (47%); the rest is mainly photoabsorption

7. Absorbing is easy: thickness of lead needed to stop 511 keV
   stop 90.0%: 1.30 cm of Pb (attenuation, first interaction)
   stop 99.0%: 2.59 cm of Pb (attenuation, first interaction)
   stop 99.9%: 3.89 cm of Pb (attenuation, first interaction)
[PASS] 3 cm of Pb attenuates > 99% of 511 keV

8. What electron density would make a true metal mirror at 511 keV?
   plasma frequency = gamma frequency at n_e = 1.89e+38 1/m^3 = 4.1e+07 times solid gold, mass density ~ 6.3e+11 kg/m^3 for Z/A = 1/2
[PASS] critical density is at white-dwarf-core scale, above 1e11 kg/m^3
   (At that density electrons are relativistic and degenerate, so this is only a scale.)

ALL CHECKS PASSED
```

</details>

<details><summary><code>code/verify_traps.py</code></summary>

```text
1. Penning trap: uniform B plus quadrupole E. Try every sign combination.
   Potential phi = s (w_z^2/2)(z^2 - r^2/2), so s=+1 is a well for positive charge.
     q    B     s     r_max     z_max     result
   +1  +1.0  +1       1.21      1.00   confined
   +1  +1.0  -1       1.19    200.21   ESCAPES
   +1  -1.0  +1       1.00      1.00   confined
   +1  -1.0  -1       1.00    200.21   ESCAPES
   -1  +1.0  +1       1.00    200.21   ESCAPES
   -1  +1.0  -1       1.00      1.00   confined
   -1  -1.0  +1       1.19    200.21   ESCAPES
   -1  -1.0  -1       1.21      1.00   confined
[PASS] each sign of charge is confined only when s matches its sign
[PASS] flipping B never changes which charge is confined
   weak field test: B = 0.1 is below sqrt(2) w_z = 0.141, escapes = True
[PASS] below the Penning stability limit w_c^2 > 2 w_z^2 even the right charge escapes
   B = 0.2 is above it, escapes = False
[PASS] just above the limit the correct charge is held
   B=+1 net rotation -5.03 rad, B=-1 +5.03 rad: polarity only reverses the sense of rotation
[PASS] reversing B reverses the magnetron rotation

2. Magnetic bottle: B alone, mirror ratio 5
   mirror ratio 5, predicted loss-cone half-angle 26.6 deg from the axis
   pitch[deg]  +,B+    +,B-    -,B+    -,B-
         20   lost   lost   lost   lost 
[PASS] pitch 20 deg behaves as the loss cone predicts
         35   held   held   held   held 
[PASS] pitch 35 deg behaves as the loss cone predicts
         60   held   held   held   held 
[PASS] pitch 60 deg behaves as the loss cone predicts
         80   held   held   held   held 
[PASS] pitch 80 deg behaves as the loss cone predicts
[PASS] a magnetic bottle treats both charges and both polarities alike

3. Nested well: one axial potential holding opposite charges in different places
   phi(z) = A[(z/a)^4 - 2 (z/a)^2]. Positive charge sits at z = +-a, negative near z = 0.
   cold positive stays in z = [1.78, 2.20], cold negative in z = [-0.20, 0.20]
[PASS] both species are axially confined at once
[PASS] cold species do not overlap
   positive given energy above the central hill passes z = 0: range [-2.93, 2.93]
[PASS] an energised positive crosses through the negative species' region

4. Crossed E and B fields (A17). Does a 'mass spectrometer style' filter sort charge?
   at the balance speed v = E/B: both go straight; separation after t=50: 5.42e-14
[PASS] at the balance speed a Wien filter does not separate e+ from e-
   off balance, v = 0.5: drift speed e+ 0.1002, e- 0.1002, E/B = 0.1
   separation stays bounded: max 1.60 over t = 2000; gyroradius-scale 2(v-E/B)/w_c x 2 = 1.60
[PASS] long crossed fields drift both charges the same way
[PASS] in long crossed fields the separation never grows
   E field alone: separation after t=50 is 250.0, prediction E t^2 = 250.0
[PASS] an E field alone separates the charges and keeps doing so
   B field alone, half turn: y = -1.000 for e+, +1.000 for e- (opposite bends)
[PASS] a B-only spectrometer bends the two charges to opposite sides

5. Stored density, energy and annihilation time (SI)
   rate coefficient sigma v = pi r_e^2 c = 7.479e-21 m^3/s
    n [1/m^3]     annihilation time (no Coulomb focusing)
       1e+15           1.34e+05 s
       1e+16           1.34e+04 s
       1e+20           1.34e+00 s
   Brillouin limit at 5.0 T: n = 1.21e+20 1/m^3
   stored rest energy if all paired at that density: 1.99e+07 J/m^3 (petrol ~3e10 J/m^3, general knowledge)
[PASS] even the Brillouin-limit positron plasma stores far less than chemical fuel per volume
   1e10 positrons (Danielson et al. storage scale): rest-energy of pairs 1.64e-03 J
   Coulomb (Sommerfeld) enhancement at 4 K relative speed: 782x
   Coulomb (Sommerfeld) enhancement at 300 K relative speed: 90x

ALL CHECKS PASSED
```

</details>

<details><summary><code>code/verify_annihilation_budget.py</code></summary>

```text
1. Monte Carlo of p pbar at rest, 1876.5 MeV
   mode A: charged 0/2/4/6 prongs (4/44/46/6 %), Poisson neutrals; mode B: 50% (2,3), 50% (4,1).
   mode  <n+->  <n0>   charged pions  neutral pions   gamma    nu     e+-    sum    photons/ann  <E_gamma>
   A      3.10  1.97         0.631          0.369    0.369  0.458  0.173  1.0000      3.93      175.9 MeV
[PASS] mode A: four-momentum conserved per event  (worst 7.3e-12 MeV)
[PASS] mode A: final classes sum to the initial energy (per event)
[PASS] mode A: <n+-> within 10% of 3.05  (3.10)
[PASS] mode A: <n0> within 20% of 1.93  (1.97)
   B      2.94  2.06         0.590          0.410    0.410  0.427  0.164  1.0000      4.11      186.8 MeV
[PASS] mode B: four-momentum conserved per event  (worst 3.2e-12 MeV)
[PASS] mode B: final classes sum to the initial energy (per event)
[PASS] mode B: <n+-> within 10% of 3.05  (2.94)
[PASS] mode B: <n0> within 20% of 1.93  (2.06)

   Rounded result: photons 39%, neutrinos 44%, electrons and positrons 17%.
[PASS] neutrinos take roughly 40-50% of the energy
[PASS] photons take well under the 1-2 GeV stated in A19  (730 MeV)

   Single decay check: one charged pion at rest (exact two-body plus mean Michel spectrum)
   E(nu_mu) = 29.79 MeV, <E(e)> = 38.4 MeV: neutrino share of a charged pion = 72.5%
[PASS] neutrino share of a charged pion is about 72%

2. A19 ledger per annihilation, using the Monte Carlo means (MeV)
   fuel rest energy consumed                    1876.5
   lost to neutrinos                             830.1
   in photons                                    730.1
   in e+- (rest + kinetic)                       316.3
   upper bound on work                          1046.4  (55.8% of fuel)
   e+e- pairs the pion decays leave behind        1.51  per annihilation
   cost to remake one e+e- pair                  1.022  (0.054% of fuel)

   Close the p pbar loop too, at the bare threshold: net = 1046.4 - 1876.5 = -830.1 MeV
[PASS] closing the proton loop is strictly negative by at least the neutrino energy
   Do not close it: work comes from burning one antiproton per cycle. That is fuel.

3. A20: the claimed formula against a real energy accounting
   E_gamma [MeV]   claimed 1 - m_e c^2/E    meaning
        0.511            0.00000    e+e- at rest
        1.000            0.48900    1 MeV photon
      185.666            0.99725    mean pi0 photon
[PASS] the claimed formula gives 0 for e+e- at rest, where 100% is not excluded by conservation
   Derived open-loop accounting: eta = f_usable * (1 - E_final/E_initial) * capture
     e+e-    : f_usable = 1, so eta -> 1 as E_final -> 0 (no threshold matters if fuel is not remade)
     p pbar  : f_usable <= 0.558, so eta <= 0.558 for any E_gamma
[PASS] open-loop p pbar ceiling is below the 1 the formula implies

4. Real production cost (sourced inputs, order of magnitude)
   beam energy per pulse: 1.5e+13 x 25.08 GeV = 6.0e+04 J
   3e+07 pbar/cycle: annihilation energy / beam energy = 1.5e-07; with 5% wall plug = 7.5e-09
   4e+07 pbar/cycle: annihilation energy / beam energy = 2.0e-07; with 5% wall plug = 1.0e-08
   LaPointe & Reddy figure: 1 pbar per 2.5e6 protons, beam efficiency 2.9e-08 (their 3e-8)
   wall-plug to delivered beam ~ 9e-09; foil capture below 1% makes it below 9e-11
[PASS] delivered-beam wall-plug efficiency is in the 1e-9 to 1e-7 range  (8.7e-09)

5. Round-trip energy storage: efficiency of put-in energy back as work
   production 1e-09 x usable 0.56 (ideal conversion) = 5.6e-10
   production 1e-08 x usable 0.56 (ideal conversion) = 5.6e-09
   Fuel energy density: 2 x 9e16 J per kg of antimatter = 1.8e17 J/kg
   making 1 g of antiprotons: rest energy of the pair mass 1.8e+14 J; wall-plug at 1e-8: 1.8e+22 J
   AD rate 3.5e+05 pbar/s means 1 g takes 5.4e+10 years
[PASS] one gram at the AD rate takes longer than the age of the universe

ALL CHECKS PASSED
```

</details>

<details><summary><code>code/verify_closed_cycle.py</code></summary>

```text
1. Stochastic cycle, per stored pair. Net = work out - work in - replacement energy (keV)
   R          beta_pair  recover   P per cross  made/total   net [keV]    ideal bound 0?
   0.29          0.700    0.50      1.78e-05   1.50e-05   -1431.0649     crossings/cycle 1.1e+00
[PASS] first law holds per event at R=0.29, beta=0.7
[PASS] net work <= 0 at R=0.29, beta=0.7, recover=0.5
[PASS] net work strictly negative at R=0.29, beta=0.7  (-1431.065 keV)
[PASS] Monte Carlo net matches analytic expectation at R=0.29, beta=0.7  (expected -1431.0595 keV)
   0.99          0.700    0.50      1.78e-05   8.65e-04   -1430.0223     crossings/cycle 5.0e+01
[PASS] first law holds per event at R=0.99, beta=0.7
[PASS] net work <= 0 at R=0.99, beta=0.7, recover=0.5
[PASS] net work strictly negative at R=0.99, beta=0.7  (-1430.022 keV)
[PASS] Monte Carlo net matches analytic expectation at R=0.99, beta=0.7  (expected -1429.9875 keV)
   0.99          0.700    1.00      1.78e-05   9.00e-04   -1429.7953     crossings/cycle 5.0e+01
[PASS] first law holds per event at R=0.99, beta=0.7
[PASS] net work <= 0 at R=0.99, beta=0.7, recover=1.0
[PASS] net work strictly negative at R=0.99, beta=0.7  (-1429.795 keV)
[PASS] Monte Carlo net matches analytic expectation at R=0.99, beta=0.7  (expected -1429.8047 keV)
   0.999999      0.700    1.00      1.78e-05   8.99e-01    -145.1333     crossings/cycle 5.1e+04
[PASS] first law holds per event at R=0.999999, beta=0.7
[PASS] net work <= 0 at R=0.999999, beta=0.7, recover=1.0
[PASS] net work strictly negative at R=0.999999, beta=0.7  (-145.133 keV)
[PASS] Monte Carlo net matches analytic expectation at R=0.999999, beta=0.7  (expected -144.5896 keV)
   0.999999999    0.700    1.00      1.78e-05   1.00e+00      -0.1646     crossings/cycle 5.6e+04
[PASS] first law holds per event at R=0.999999999, beta=0.7
[PASS] net work <= 0 at R=0.999999999, beta=0.7, recover=1.0
[PASS] net work strictly negative at R=0.999999999, beta=0.7  (-0.165 keV)
[PASS] Monte Carlo net matches analytic expectation at R=0.999999999, beta=0.7  (expected -0.1608 keV)
   0.999999999    0.100    1.00      1.36e-06   9.99e-01      -1.5304     crossings/cycle 7.4e+05
[PASS] first law holds per event at R=0.999999999, beta=0.1
[PASS] net work <= 0 at R=0.999999999, beta=0.1, recover=1.0
[PASS] net work strictly negative at R=0.999999999, beta=0.1  (-1.530 keV)
[PASS] Monte Carlo net matches analytic expectation at R=0.999999999, beta=0.1  (expected -1.5107 keV)
   0.99999999999    0.010    1.00      1.33e-07   1.00e+00      -0.1278     crossings/cycle 7.5e+06
[PASS] first law holds per event at R=0.99999999999, beta=0.01
[PASS] net work <= 0 at R=0.99999999999, beta=0.01, recover=1.0
[PASS] net work strictly negative at R=0.99999999999, beta=0.01  (-0.128 keV)
[PASS] Monte Carlo net matches analytic expectation at R=0.99999999999, beta=0.01  (expected -0.1535 keV)
   1             0.010    0.99      1.33e-07   1.00e+00      -0.0005     crossings/cycle 7.5e+06
[PASS] first law holds per event at R=1, beta=0.01
[PASS] net work <= 0 at R=1, beta=0.01, recover=0.99
[PASS] net work strictly negative at R=1, beta=0.01  (-0.001 keV)
[PASS] Monte Carlo net matches analytic expectation at R=1, beta=0.01  (expected -0.0005 keV)
   1             0.010    1.00      1.33e-07   1.00e+00       0.0000     crossings/cycle 7.5e+06
[PASS] first law holds per event at R=1, beta=0.01
[PASS] net work <= 0 at R=1, beta=0.01, recover=1.0
[PASS] Monte Carlo net matches analytic expectation at R=1, beta=0.01  (expected -0.0000 keV)

   Sweep of expansion depth E_f with R=0.99, beta=0.7, recover=1: net must not depend on E_f
   E_f =  500.0 keV: net = -1430.0243 keV
   E_f =  255.5 keV: net = -1430.0243 keV
   E_f =  100.0 keV: net = -1430.0243 keV
   E_f =   10.0 keV: net = -1430.0243 keV
[PASS] work taken by the piston is exactly repaid: net is independent of expansion depth

2. Rate price of the zero-loss limit: conversion P shrinks with the pair speed
   Assumes cavity length 1 m. Time for one conversion even with R = 1.
   beta_pair   KE per pair [keV]   P per crossing   crossings needed   mean cycle time [s]
     7.0e-01     4.0909e+02           1.78e-05        5.62e+04           3.75e-04
     1.0e-01     5.1486e+00           1.36e-06        7.36e+05           4.91e-03
     1.0e-02     5.1104e-02           1.33e-07        7.51e+06           5.01e-02
     1.0e-03     5.1100e-04           1.33e-08        7.51e+07           5.01e-01
     1.0e-04     5.1100e-06           1.33e-09        7.51e+08           5.01e+00
     1.0e-06     5.1104e-10           1.33e-11        7.51e+10           5.01e+02
[PASS] P tends to 0 as the pair kinetic energy tends to 0

3. Nuclear recoil in A16: gamma + nucleus. Budget per pair at threshold (eV)
   proton : recoil excess = 556.596 eV of 1021998 eV
[PASS] recoil on proton makes the cycle strictly negative when not recovered  (net -556.596 eV)
   Au     : recoil excess = 2.846 eV of 1021998 eV
[PASS] recoil on Au makes the cycle strictly negative when not recovered  (net -2.846 eV)

4. A07: cascade with charges tracked as integers
   X    2 m c^2 [MeV]   least cost on a proton target [MeV]     gain     net [MeV]
   e          1.022                  1.023                1.022       -0.001
[PASS] cascade step net is not positive for e
   mu       211.316                235.112              211.316      -23.796
[PASS] cascade step net is not positive for mu
   pi       279.140                320.663              279.140      -41.523
[PASS] cascade step net is not positive for pi
   p       1876.544               5629.632             1876.544    -3753.088
[PASS] cascade step net is not positive for p

   Recycling all output into new antiparticles, electrons as debris, 1e9 MeV to start
   eps = 1       pairs per generation [1000, 1000, 1000, 1000]..., final energy 1.000000e+09 MeV, debris X = 1000, leftover X-bar = 0
[PASS] energy never grows at eps=1.0
[PASS] net lepton number X - Xbar is unchanged at eps=1.0
[PASS] debris count is still 1000 after 8 generations at eps=1.0
   eps = 0.5     pairs per generation [1000, 1000, 1000, 1000]..., final energy 9.999918e+08 MeV, debris X = 1000, leftover X-bar = 0
[PASS] energy never grows at eps=0.5
[PASS] net lepton number X - Xbar is unchanged at eps=0.5
[PASS] debris count is still 1000 after 8 generations at eps=0.5
   eps = 1e-08   pairs per generation [9, 0, 0, 0]..., final energy 8.020001e+07 MeV, debris X = 1000, leftover X-bar = 0
[PASS] energy never grows at eps=1e-08
[PASS] net lepton number X - Xbar is unchanged at eps=1e-08
[PASS] debris count is still 1000 after 8 generations at eps=1e-08

   Lifetimes compared with a typical AD cycle of ~100 s
   pi0   tau = 8.40e-17 s; a 100 s production cycle is 1.2e+18 lifetimes
   pi+-  tau = 2.60e-08 s; a 100 s production cycle is 3.8e+09 lifetimes
   mu    tau = 2.20e-06 s; a 100 s production cycle is 4.6e+07 lifetimes
[PASS] a charged pion is gone 1e9 lifetimes before a new antipion could be made
   neutrino at 1 MeV CM energy: sigma ~ G_F^2 s/pi = 1.8e-44 cm^2; mean free path in lead ~ 22 light years
[PASS] a ~MeV neutrino passes through light years of lead

5. A03: reaction equilibrium, Kelvin and Carnot
   Reversible means zero affinity: n+ n- = n_eq^2 with n_eq = 2 (2 pi m k T/h^2)^(3/2) exp(-mc^2/kT)
    T [K]        kT/mc^2     log10(n_eq [1/m^3])
    3.80e-11    6.408e-21        -6.7772e+19
    2.73e+00    4.595e-10        -9.4507e+08
    3.00e+02    5.059e-08        -8.5844e+06
    1.50e+03    2.530e-07        -1.7169e+06
    1.00e+08    1.686e-02         7.9307e+00
    1.00e+09    1.686e-01         3.2609e+01
    6.00e+09    1.012e+00         3.5922e+01
    1.00e+10    1.686e+00         3.6426e+01
[PASS] at 1500 K pair equilibrium density is below 10^-1e6 per m^3  (10^-1.717e+06)
   at 6e9 K (non-relativistic form, rough): n_eq/n_gamma ~ 0.19: pairs appear freely
[PASS] thermal pair creation only matters near 6e9 K

   Entropy of a photon in a cavity, compared with its energy (A21 D)
   modes in 1 m^3 within dE/E=0.2%: 3.52e+33; S = k ln(modes) = 77 k = 111 bits
   T S at 300 K = 2.00 eV per photon = 3.9e-06 of its energy
   Landauer cost at 300 K = 0.0179 eV per bit
[PASS] the entropy price of knowing a 511 keV photon's mode is below 1e-5 of its energy

   Mixing entropy of one e+ and one e-: k ln 2 each at most
   2 k T ln 2 at 300 K = 0.036 eV of 1.022 MeV = 3.5e-08

   Kelvin test with a real heat bath: W <= Q_h (1 - T_c/T_h), with no pair term at all
   T_h=1500 K, T_c=300 K: eta_C = 1 - 2.000e-01
   T_h=1500 K, T_c=2.73 K: eta_C = 1 - 1.817e-03
   T_h=1500 K, T_c=3.8e-11 K: eta_C = 1 - 2.533e-14
   Over a closed cycle the pair stock returns to start, so the annihilation term is zero.
   Any positive W must be ordinary heat-engine output paid for by Q_h from the bath.

6. Three-photon decay of ortho-positronium cannot be undone by two-photon collision
   74966 sampled decays: largest s of any two photons / (2m)^2 = 0.9978
[PASS] for three photons, any two are below the pair threshold
   Ratio of lifetimes: 142 ns / 125 ps = 1136, so with 3:1 spin weight a spin-averaged free 2-gamma:3-gamma ratio near 379:1 (standard value about 372:1).

ALL CHECKS PASSED
```

</details>

