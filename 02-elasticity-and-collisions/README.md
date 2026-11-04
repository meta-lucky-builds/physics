# Elasticity, bouncing and collisions

Questions E01 to E09 share one chain of reasoning. Atoms sit in potential wells, and a small push in a well gives a linear restoring force. A body bounces if the energy of impact is stored in those wells, returned before it leaks away, and returned as upward motion of the centre of mass rather than as spin, vibration or waves in whatever it hit. Each answer below is one link of that chain.

Six scripts test the claims. Every one simulates the system or solves it by a separate route, then compares with the derivation. Run them from the repo root with `python 02-elasticity-and-collisions/code/verify_spring.py` and so on. Printed output here is copied from a real run and shortened only by leaving out lines.

Many inputs are assumed values for ordinary objects: bed mass, foam permeability, ball stiffness. They are marked where they appear. Only what those inputs do to the answer has been tested, not their exact size. Sources are listed at the end, with a note on how far each was checked.

### E01: What makes a spring behave the way it behaves?

**Status:** Derived + Code

**Short answer:** Compressing a coil spring twists its wire very slightly, which shifts the atoms by a tiny amount in their potential wells. The wells push them back in proportion to the shift, so the spring pushes back in proportion to how far you moved it. It overshoots the rest position because the moving mass still has momentum there, and that is why it compresses and then elongates in turn.

**Starting from:** Newton's laws. A pair potential between atoms with one minimum, repulsive close in (Pauli exclusion and Coulomb repulsion) and attractive further out. A small-strain linear elastic wire with shear modulus $G$. A close-coiled helix.

**Derivation:**

*Step 1, atoms.* Near the minimum of any smooth potential $U(r)$, with $x=r-r_0$,

$$U(r)=U(r_0)+\tfrac12 U''(r_0)x^2+\tfrac16 U'''(r_0)x^3+\dots$$

The linear term vanishes because $r_0$ is a minimum. So $F=-dU/dx=-kx+\dots$ with $k=U''(r_0)$. Hooke's law is therefore not a special property of springs. It is what any smooth minimum looks like up close. For two standard model potentials:

- Lennard-Jones, $U=4\epsilon[(\sigma/r)^{12}-(\sigma/r)^6]$: $r_0=2^{1/6}\sigma$, $k=36\cdot 2^{2/3}\,\epsilon/\sigma^2=57.15\,\epsilon/\sigma^2$.
- Morse, $U=D(1-e^{-a(r-r_0)})^2$: $k=2Da^2$, $U'''(r_0)=-6Da^3$.

$U'''<0$ in both, so the force law is lopsided. Compressing is stiffer than stretching. The law is only linear for small strain.

*Step 2, from a bond to a wire.* A wire is a three-dimensional network of such bonds, and at small strain the network inherits their linearity. Its response is summed up by a Young's modulus $E$ and a shear modulus $G$. The size estimate is $E\sim k/r_0$, derived and tested in [E04](#e04-what-is-elasticity-and-what-does-it-depend-on).

*Step 3, from a wire to a coil.* Push along the coil's axis with force $F$. Cut the wire anywhere: equilibrium requires a force $F$ across the cut and a moment $F\cdot D/2$ about the wire's axis, because the cut is $D/2$ from the line of the force. For a close-coiled helix the wire is nearly horizontal, so that moment is a twist of the wire, not a bend. The twist energy per unit length of a round wire is $T^2/(2GJ)$ with $J=\pi d^4/32$. Over the wire length $\pi D n$ that gives

$$U=\frac{(FD/2)^2}{2GJ}\,\pi D n=\frac{4F^2D^3n}{Gd^4}.$$

Energy is also $U=\tfrac12 Fx$ for a linear system, so $x=dU/dF$ gives

$$x=\frac{8FD^3n}{Gd^4}\quad\Rightarrow\quad k=\frac{Gd^4}{8D^3n}.$$

So a coil spring works mainly by torsion of the wire, not by stretching it. That also explains its strength. For the spring below, a force of 54 N shortens it by 30 mm, half its 60 mm length. The peak shear stress is $8FD/(\pi d^3)=344$ MPa, a shear strain of 344 MPa / 72 GPa = 0.48%. The coil turns a strain of about half a percent into a 50% change of length. Half a percent lies inside the nearly linear region of the table below, which is why a coil obeys Hooke's law over a long travel.

*Step 4, why it oscillates.* Let a mass $m$ sit on the spring. At rest position the force is zero but the speed is greatest, so $m\ddot x=-kx$ carries it through to the other side. Energy moves between $\tfrac12 kx^2$ and $\tfrac12 mv^2$ with period $T=2\pi\sqrt{m/k}$. If the spring's own mass $m_s$ matters, the period is close to $2\pi\sqrt{(m+m_s/3)/k}$. [E07](#e07-will-a-ball-between-two-opposed-springs-keep-going) tests that rule.

**Check:** [`code/verify_spring.py`](code/verify_spring.py) tests each link without using the formula under test. A relaxed chain of 20 atoms is solved for equilibrium at imposed strain. A 3D finite-element beam model of the real helix, with no assumption about torsion, gives the stiffness and splits the strain energy.

```text
     chain secant stiffness per bond vs strain (eps/sig^2):
       strain -0.020:   70.69   departure from k: +23.71 %
       strain -0.010:   63.52   departure from k: +11.15 %
       strain -0.005:   60.24   departure from k:  +5.41 %
       strain -0.002:   58.36   departure from k:  +2.13 %
       strain -0.001:   57.75   departure from k:  +1.06 %
       strain +0.001:   56.55   departure from k:  -1.04 %
       strain +0.002:   55.96   departure from k:  -2.07 %
       strain +0.005:   54.23   departure from k:  -5.10 %
       strain +0.010:   51.48   departure from k:  -9.91 %
       strain +0.020:   46.44   departure from k: -18.73 %
     chain, even part of the secant at 0.1% strain: 57.150 (U'' = 57.146); odd part gives U''' = -1069 (exact -1069)
  [PASS] relaxed chain reproduces U''(r0) to 0.1%
  [PASS] relaxed chain reproduces U'''(r0) to 2%
     at 5% strain: stretch 34.4, compression 98.3 (eps/sig^2)
  [PASS] anharmonic: compression is stiffer than stretch at 5% strain
4. Coil spring: finite-element helix vs k = G d^4 / (8 D^3 n)
     wire d = 2.0 mm, coil D = 20 mm, n = 10, pitch = 6 mm
        12 elements/turn: k_FE =  1910.89 N/m, k_FE/k_formula = 1.0616
        24 elements/turn: k_FE =  1823.60 N/m, k_FE/k_formula = 1.0131
        48 elements/turn: k_FE =  1802.45 N/m, k_FE/k_formula = 1.0014
        96 elements/turn: k_FE =  1797.21 N/m, k_FE/k_formula = 0.9984
     k_formula = 1800.00 N/m
     strain-energy shares: torsion 99.27 %, bending 0.73 %, stretch 0.002 %
  [PASS] wire torsion carries over 98% of the energy at small pitch
  [PASS] FE stiffness within 1% of the formula at small pitch
     FE scaling: 2d -> x16.00 (formula 16), 2D -> x0.1252 (formula 0.125), 2n -> x0.500 (formula 0.5)
  [PASS] FE scaling exponents d^4, D^-3, n^-1
     steep coils: exact helix energy result vs FE (formula only holds at small pitch)
       pitch angle  5 deg: k_FE   1797.9, helix energy   1796.2, simple formula   1800.0 N/m; torsion share 99.4 %
       pitch angle 15 deg: k_FE   1766.7, helix energy   1765.0, simple formula   1800.0 N/m; torsion share 94.7 %
       pitch angle 30 deg: k_FE   1652.4, helix energy   1650.8, simple formula   1800.0 N/m; torsion share 79.4 %
       pitch angle 45 deg: k_FE   1433.7, helix energy   1432.4, simple formula   1800.0 N/m; torsion share 56.2 %
  [PASS] FE agrees with the exact small-strain helix energy at all pitch angles
```

The chain's even and odd parts recover $U''$ and $U'''$ separately. The FE stiffness converges to the formula as the mesh refines, and matches a helix energy that keeps every pitch-angle term.

**Limits:** $k=Gd^4/(8D^3n)$ holds for small pitch and a coil much wider than its wire. The FE model leaves out transverse shear and wire curvature. An energy estimate of the shear term gives a correction of $1/(2\kappa C^2)$ with spring index $C=D/d$ and $\kappa\approx 0.9$, about 0.6% here. Springs with closed or ground ends have fewer active turns than total turns. Once the coils touch the law ends, and past the yield point the wire keeps a permanent set. The Lennard-Jones form describes noble-gas atoms well and metals poorly, and is used here for its shape. Metallic bonding is not a sum of pair potentials, but the Taylor argument needs only a smooth minimum.

### E02: Why does any object bounce? Why does a ball bounce but not a sponge?

**Status:** Derived + Code + Estimate

**Short answer:** An object bounces when impact energy is stored as elastic strain and comes back as upward motion before it is lost. A rubber ball returns most of it. A sponge loses far more for four separate reasons: viscous air flow through its pores and hysteresis in its skeleton, its low mass against air drag, landing tilted so energy goes to rotation, and contact time close to its relaxation time. The sponge does rebound a little, and the model gives about 6 cm from a 1 m drop.

**Starting from:** Energy conservation. Rigid-body impulse laws. A contact force that can only push: $F=\max(0,\,k\delta+c\dot\delta)$. Darcy's law for air in a porous solid. The numbers for the foam, the ball and the pen are assumed typical values, listed below, not measurements.

**Derivation:**

*Why a floor can push you up at all.* A floor does no work on a body that isn't moving through it. It can only stop the bottom of the ball, and the top keeps coming, so the ball squashes. Impact kinetic energy turns into strain energy. The body then springs back and lifts itself off. So a bounce needs elastic storage. A coefficient of restitution $e$ = (speed out)/(speed in) measures how much of it returned, and the rebound height is $e^2h$ if air doesn't matter. $e$ belongs to the pair of body and surface, not to the body alone.

*Damping.* Take a mass on a spring-damper with $\zeta=c/(2\sqrt{km})$ and $s=\sqrt{1-\zeta^2}$. The usual formula $e=\exp(-\pi\zeta/s)$ is wrong for a real contact, as Hunt and Crossley pointed out, because it makes the surfaces pull on each other just before they part. Requiring $F\ge 0$, the body leaves when $F=0$, which is where $\ddot\delta=0$. With $\delta\propto e^{-\zeta\omega t}\sin(s\omega t)$ the lift-off phase $\phi$ satisfies $\tan\phi=2\zeta s/(2\zeta^2-1)$ on $(0,\pi)$, and

$$e=-e^{-\zeta\phi/s}\left(\cos\phi-\frac{\zeta}{s}\sin\phi\right).$$

The body leaves while still compressed, then recovers slowly in the air. That matches a sponge that visibly re-expands after leaving the floor.

*Why foam loses so much.* A sponge is mostly air. Squeeze it and the air has to get out. Take strain rate $v/L$ and escape path $L$. Air speed in the pores rises linearly from zero at the closed end, Darcy's law $dp/dz=\mu u/\kappa$ gives mean pressure $\mu vL/(3\kappa)$, and this acts as a damper:

$$\zeta_{\text{air}}=\frac{\mu L}{6\kappa\sqrt{E\rho}},\qquad \tau_{\text{poro}}\sim\frac{\mu L^2}{\kappa E}.$$

Gent and Rusch attributed the large damping of open-cell foams to this air movement, with permeability as the governing property. A smaller pore size or a larger block means more loss. A closed-cell "sponge ball" keeps its air in, which then acts as a gas spring with no such loss, so a ball made of foam can bounce where a kitchen sponge does not. This may be why your "ball made of sponge" behaves differently from a duster, though from here I can't tell which foam yours was.

*Contact time against relaxation time.* Elastic and viscous behaviour are not fixed properties, as the same material can do both. Compare the relaxation time $\tau$ of the material with the contact time $t_c$. A standard linear solid (spring in parallel with a spring and dashpot in series) shows it clearly: when $\tau\ll t_c$ the dashpot has relaxed and only the soft spring acts, when $\tau\gg t_c$ the dashpot is locked and both springs act, and in both cases the bounce is nearly lossless. Loss is greatest when $\tau\sim t_c$. The foam block below has $\tau_{\text{poro}}$ = 1.8 ms and contact 3.2 ms, near that lossy band.

A material with no long-term spring (a putty) has $F$ obeying $m\ddot F+(m/\tau)\dot F+kF=0$. $F$ reaches zero only if $\tau\omega>1/2$. Below that the body sinks in and never leaves. So a putty should rest on the floor under a slow push and bounce in a fast impact, the behaviour commonly described for silly putty (general knowledge, unverified here).

*Air.* Terminal speed scales as $\sqrt{m/(\rho_{\text{air}}C_dA)}$. A very light sponge doesn't arrive fast, and its small rebound energy is soon spent against drag.

*Rotation.* See [E03](#e03-what-does-bouncing-depend-on-is-it-the-area-of-contact). A sphere's contact normal passes through its centre, so none of the rebound goes into spin. A tilted slab can turn most of it into spin.

**Check:** [`code/verify_bounce.py`](code/verify_bounce.py) integrates the contact, audits the energy and compares with the closed form.

```text
1. Spring-damper that cannot pull: e versus damping ratio
     zeta    e(sim)   e(closed form)  naive exp(-pi z/s)  rebound height e^2   energy audit   body still compressed at lift-off (units v0/w)
     0.001   0.9969    0.9969         0.9969             0.9937         +1.1e-10     0.002
     0.010   0.9693    0.9693         0.9691             0.9395         -6.4e-11     0.019
     0.030   0.9117    0.9117         0.9100             0.8311         -7.1e-11     0.055
     0.100   0.7441    0.7441         0.7292             0.5537         +1.3e-11     0.149
     0.200   0.5717    0.5717         0.5266             0.3269         -2.4e-12     0.229
     0.300   0.4510    0.4510         0.3723             0.2034         +3.4e-12     0.271
     0.500   0.2984    0.2984         0.1630             0.0891         -2.5e-12     0.298
     0.700   0.2103    0.2103         0.0460             0.0442         -1.5e-12     0.294
     0.900   0.1553    0.1553         0.0015             0.0241         -2.1e-11     0.280
       1.0   0.1353   (overdamped)                         0.0183         -9.7e-12     0.271
       2.0   0.0478   (overdamped)                         0.0023         -1.4e-12     0.191
       5.0   0.0093   (overdamped)                         0.0001         -1.4e-12     0.093
  [PASS] simulated e matches the closed form for zeta < 1
  [PASS] energy audit closes to 1e-6
  [PASS] result converged against integrator tolerance  (|de| = 1.2e-08)
  [PASS] e falls monotonically with damping
3. Contact time against relaxation time (standard linear solid, k_inf = k1 = 1)
     tau/t_contact   e        | putty-like solid (k_inf = 0): tau*w   e
           0.0003   0.998   |                        0.0010   0.000
           0.0010   0.995   |                        0.0032   0.000
           0.0032   0.985   |                        0.0100   0.000
           0.0102   0.953   |                        0.0316   0.000
           0.0330   0.868   |                        0.1000   0.000
           0.1144   0.715   |                        0.3162   0.000
           0.4259   0.712   |                        1.0000   0.163
           1.4136   0.859   |                        3.1623   0.605
           4.4982   0.948   |                       10.0000   0.854
          14.2342   0.983   |                       31.6228   0.952
          45.0155   0.994   |                      100.0000   0.984
         142.3524   0.998   |                      316.2278   0.995
         450.1581   0.999   |                     1000.0000   0.998
     minimum e = 0.712 at tau/t_contact = 0.43
  [PASS] elastic in both limits, lossy when contact time ~ relaxation time
  [PASS] putty-like solid: no lift-off when relaxation is quick, bounces when contact is fast  (e 0.00 -> 1.00)
     putty threshold, derived tau*w = 1/2: e(0.45) = 0.000, e(0.55) = 0.001
  [PASS] putty lift-off threshold at tau*w = 1/2
4. Drop from 1 m in air (inputs are assumed typical values, see README)
   foam air-flow damping estimate, E = 30 kPa, rho = 30 kg/m^3, L = air escape path:
     permeability 1e-09 m^2, L 0.5 cm: zeta_air 0.016, poroelastic relaxation 0.015 ms
     permeability 1e-09 m^2, L 3.0 cm: zeta_air 0.095, poroelastic relaxation 0.540 ms
     permeability 3e-10 m^2, L 0.5 cm: zeta_air 0.053, poroelastic relaxation 0.050 ms
     permeability 3e-10 m^2, L 3.0 cm: zeta_air 0.316, poroelastic relaxation 1.800 ms
     permeability 1e-10 m^2, L 0.5 cm: zeta_air 0.158, poroelastic relaxation 0.150 ms
     permeability 1e-10 m^2, L 3.0 cm: zeta_air 0.949, poroelastic relaxation 5.400 ms
     object                            m(g)     k(N/m)   zeta   v_in   e_contact  contact(ms)  rebound (cm)  no-drag e^2 h (cm)
     rubber ball, r 2.5 cm             71.99     200000  0.030   4.41    0.911       1.89       81.8        83.1
     pen on its tip (rigid plastic)    10.00     100000  0.060   4.42    0.834       1.00       69.0        69.5
     sponge block flat, 10x6x3 cm       5.40       6000  0.416   3.19    0.348       3.19        6.0        12.4
     pen-shaped sponge on its tip       0.33         17  0.153   4.12    0.638      14.18       33.6        41.6
     pen-shaped sponge falling flat     0.33       4200  0.153   1.85    0.644       0.89        6.0        41.6
  [PASS] rubber ball rebounds far higher than the flat sponge
  [PASS] same pen geometry: rigid pen rebounds higher than pen-shaped sponge on its tip
```

Inputs for the table: rubber ball $k$ = 200 kN/m, $\zeta$ = 0.03. Pen tip $k$ = 100 kN/m, $\zeta$ = 0.06. Foam $E$ = 30 kPa, $\rho$ = 30 kg/m³, permeability $3\times10^{-10}$ m², skeleton hysteresis $\zeta$ = 0.1. All are assumed. For a sense of scale, the ITF requires a tennis ball dropped from 254 cm to rebound 135 to 147 cm, which is $e\approx$ 0.73 to 0.76 and $\zeta\approx$ 0.1 in this model.

**Limits:** A single spring-damper is a cartoon of a body. Real rubber and foam are nonlinear and rate-dependent, and Cross measured distinct hysteresis curves for different balls. The foam numbers set the size of the rebound, which moves by a factor of several across the permeability range printed above, so the model supports "much lower" and not a precise 6 cm. Plastic deformation (a clay ball, a dented can), adhesion and sticking (a wet sponge) are further ways to lose a bounce that aren't modelled. No drop test of a real sponge was done.

### E03: What does bouncing depend on? Is it the area of contact?

**Status:** Derived + Code + Estimate

**Short answer:** Contact area is not the deciding factor, and the pen against the pen-shaped sponge proves it, because that comparison keeps the area and changes only the material. A pen bounces because it is hard and low-loss at the tip. It bounces only sometimes because a long thin body landing at an angle puts energy into spin. You were right about one thing: in a viscous material a smaller contact loses somewhat less, so area has a real but secondary effect.

**Starting from:** The contact model of E02. Rigid-body impulse with Newton's restitution at the contact point. A frictionless floor for the rigid-body part.

**Derivation:**

*The area test.* Hold the mass and the material fixed and scale the contact area by $A$. Then $k\propto A$. Two loss laws bracket real materials:

- Viscous (Kelvin) solid: $c\propto A$ too, so $\zeta=c/(2\sqrt{km})\propto\sqrt A$. Small contacts are less damped, in the direction you guessed.
- Loss proportional to elastic stress (Hunt and Crossley's law, $F=K\delta^n(1+\tfrac32\alpha\dot\delta)$): rescaling time removes $K$, so $e$ is independent of area, and also of the contact shape $n$.

A large contact can still bounce well: a basketball flattens over a wide patch and rebounds well.

*The geometry that does matter.* Let a rigid body of mass $m$ and radius of gyration $\rho_g$ hit the floor at a point whose horizontal distance from the centre of mass is $a$, without initial spin. Apply a vertical impulse $J$. Then $\Delta v=J/m$, $\Delta\omega=Ja/(m\rho_g^2)$, and the contact point speed is $v+\omega a$. Newton's rule gives

$$J=\frac{(1+e)\,m\,v_0}{1+a^2/\rho_g^2},\qquad v_{\text{COM}}'=v_0\left[\frac{1+e}{1+a^2/\rho_g^2}-1\right].$$

The centre of mass rises only if $e>a^2/\rho_g^2$. A sphere has $a=0$, so it always rises. A thin rod at angle $\theta$ from flat has $a^2/\rho_g^2=3\cos^2\theta$, so it rises only when $e>3\cos^2\theta$. With $e$ = 0.8 that needs $\theta>59°$. Hit nearly flat, on one end, a pen keeps moving down and starts spinning.

Spin is not lost for good. The other end then swings into the floor, and that impact can turn spin back into upward motion, as the pen's clatter shows. Stoianovici and Hurmuzlu's experiments on dropped steel bars attributed the varying restitution mainly to bar vibration and multiple impacts.

*Why your pen bounces and another does not.* The simulated rise is very sensitive to $e$ at the tip. Median rise falls from about 35% of the drop height at $e$ = 0.7, to 10% at 0.5, to none at 0.3. So a hard-tipped pen and a pen with a rubber grip or a loose refill can sit on opposite sides of that drop.

*Why the pen-shaped sponge fails.* It has the same shape, so the same rotation problem, but four additions:

1. Lower $e$, from the air and hysteresis of E02.
2. About 0.33 g of mass (assumed density), so drag dominates when it falls flat: it arrives at 1.85 m/s after 1 m, where free fall gives 4.43 m/s.
3. It buckles when it lands on its tip. The elastic peak force is about 41 times its Euler load of 7.4 mN, so it folds instead of rebounding along its axis. The straight-column rebound of 33.6 cm in the table of E02 is therefore an overestimate. A solid plastic rod of the same size has an Euler load near 740 N.
4. Lower stiffness means a long contact, 14 ms against 1 ms, leaving time for the body to tip over.

**Check:** [`code/verify_bounce.py`](code/verify_bounce.py), parts 2 and 5. The first-impact formula is tested against an event-driven simulation of a rigid rectangle, with sequential corner impacts and exact flights. The Monte Carlo drops 300 bodies from 0.5 m at random orientation with some initial spin. "COM rise" means the highest point the centre of mass reaches above its height at first contact, in any later flight.

```text
2. Contact-area test
   (a) Fixed material and mass; contact area scales both k and c (Kelvin solid), zeta ~ sqrt(A):
       area x0.01   zeta 0.0030 e 0.991, contact time x10.00
       area x0.1    zeta 0.0095 e 0.971, contact time x3.16
       area x1      zeta 0.0300 e 0.912, contact time x1.00
       area x10     zeta 0.0949 e 0.755, contact time x0.32
       area x100    zeta 0.3000 e 0.451, contact time x0.10
   (b) Fixed material, loss proportional to elastic stress (Hunt-Crossley); area sets K and n:
       K        1  n 1.0 (1 = flat pad, 1.5 = Hertz ball): e 0.9090, contact time 3.144
       K      100  n 1.0 (1 = flat pad, 1.5 = Hertz ball): e 0.9090, contact time 0.314
       K    10000  n 1.0 (1 = flat pad, 1.5 = Hertz ball): e 0.9090, contact time 0.031
       K        1  n 1.5 (1 = flat pad, 1.5 = Hertz ball): e 0.9090, contact time 3.252
       K      100  n 1.5 (1 = flat pad, 1.5 = Hertz ball): e 0.9090, contact time 0.515
       K    10000  n 1.5 (1 = flat pad, 1.5 = Hertz ball): e 0.9090, contact time 0.082
  [PASS] (b) e does not depend on area when loss follows elastic stress  (spread 3.9e-08)
  [PASS] (a) smaller contact bounces better in a Kelvin solid (the student's direction)  (e 0.991 vs 0.451)
     foam pen on its tip: Euler load 7.4 mN, own weight 3.2 mN, elastic peak force 307 mN = 41 x Euler load
  [PASS] foam pen must buckle on landing, so the straight-column rebound above overstates it
     solid plastic pen for comparison: Euler load 742 N
     pen-shaped sponge falling flat hits at only 1.85 m/s after 1 m; free fall would give 4.43 m/s

5. Rigid-body drop from 0.5 m at random orientation, frictionless floor
     thin rod tilted  0 deg from flat, e = 0.8: COM velocity after first hit -0.550 v, derived -0.550 v, rotation takes 61 % of the energy
     thin rod tilted 10 deg from flat, e = 0.8: COM velocity after first hit -0.540 v, derived -0.540 v, rotation takes 62 % of the energy
     thin rod tilted 30 deg from flat, e = 0.8: COM velocity after first hit -0.446 v, derived -0.446 v, rotation takes 69 % of the energy
     thin rod tilted 60 deg from flat, e = 0.8: COM velocity after first hit +0.029 v, derived +0.029 v, rotation takes 79 % of the energy
     thin rod tilted 90 deg from flat, e = 0.8: COM velocity after first hit +0.800 v, derived +0.800 v, rotation takes 0 % of the energy
  [PASS] simulated first impact matches the derived rod formula
     impact-finder grid check, 60 drops: median rise 0.1646 m vs 0.1646 m on a 5x finer grid
  [PASS] Monte Carlo median not sensitive to the impact-finder grid
     case                                      median COM rise  P(rise > 10% of drop)  P(first hit sends COM up)  median rotation share
     pen, 14 x 1 cm, e_tip = 0.7                   34.6 %          98.0 %               33.0 %               59.4 %
     pen, 14 x 1 cm, e_tip = 0.5                    9.8 %          46.7 %               34.3 %               45.6 %
     pen-shaped sponge, e = 0.5 (no drag)          10.2 %          51.7 %               36.0 %               45.9 %
     pen-shaped sponge, e = 0.3 (no drag)           0.0 %           0.0 %               20.3 %               35.4 %
     sponge block 10 x 3 cm, e = 0.5               12.8 %          85.0 %               39.3 %               46.7 %
     sponge block 10 x 3 cm, e = 0.3                1.6 %           0.0 %               39.0 %               34.3 %
     ball with e = 0.7 for comparison: COM rise always 49 % of the drop height, rotation share 0 % (the normal passes through the centre)
  [PASS] pen with e = 0.7 rises less than a ball of the same e
  [PASS] the COM-rise chance grows with tip e
```

The rows "pen, $e$ = 0.5" and "pen-shaped sponge, $e$ = 0.5" are the same physics run with different random numbers, so their differences show the sampling noise, a few percent.

**Limits:** The rigid-body part ignores friction, which would change the spin. It also ignores flexible modes of the body and treats $e$ as constant whatever the angle, though Stoianovici and Hurmuzlu found it is not constant for all bars. Random uniform orientation is an assumption about how pens are dropped. Tip values of $e$ are assumed, not measured for real pens.

### E04: What is elasticity and what does it depend on?

**Status:** Derived + Code + Established

**Short answer:** Elasticity is the ability of a body to return to its original shape and give back the work done on it when the load is removed. It is a different property from stiffness. Stiffness (Young's modulus) says how much force a given strain takes, and depends mainly on bond strength and spacing. So "stiffer" does not mean "more elastic" in this sense: steel is about 100,000 times stiffer than soft rubber, but rubber recovers from far larger strains.

**Starting from:** The bond potential of E01. Statistical mechanics for rubber: free energy $F=U-TS$. The measured values cited below.

**Derivation:**

*Three separate properties.*

| Property | Question it answers | Set mainly by |
|---|---|---|
| Stiffness: $E$, $G$ | How hard is it to strain? | Bond curvature, atomic spacing, and for foams the density and geometry |
| Elastic limit | How far before damage stays? | Dislocations, cracks and flaws in metals and glass; chain extensibility in rubber |
| Resilience, or loss | How much work comes back? | Internal friction, rate and temperature (E02) |

A common school usage says steel is "more elastic" than rubber because its modulus is larger. That is a statement about stiffness only. It's a convention worth recognising in exams, and it should not be confused with the everyday and physical meaning of reversibility.

*Stiffness from bonds.* A simple cubic crystal has one bond of stiffness $k$ per area $r_0^2$, so stress = $k\,\Delta r/r_0^2$ and strain = $\Delta r/r_0$, giving $E\approx k/r_0$. That suggests why E is so similar in magnitude across dense solids: bond stiffness over a few ångströms. A foam has far lower stiffness because most of its volume is air and its struts bend.

Solid argon provides a real test. Summing the Lennard-Jones potential over an fcc lattice gives $B\approx 75\,\epsilon/\sigma^3$ when the second derivative of energy is taken numerically. Using $\epsilon,\sigma$ values from pair-potential fits gives 3.1 to 3.5 GPa, against 2.88 ± 0.04 GPa measured by Wilkins (his value is for $T=0$, $P=0$; he notes disagreement with sound-velocity values). That is within about 20% from a two-parameter potential, with zero-point motion ignored.

*Limits from bonds and from flaws.* A perfect chain under a dead load has no equilibrium past the inflection of $U(r)$: $r=(26/7)^{1/6}\sigma$, a force of $2.40\,\epsilon/\sigma$, 10.9% strain. Real solids yield far earlier because dislocations and cracks concentrate stress (general knowledge). The music wire listed by MakeItFrom has yield 2050 MPa and $E$ = 190 GPa, an elastic strain of about 1.1%, which is high for a metal. Past that limit atoms change neighbours, and returning the load no longer returns the shape.

*Rubber is a different mechanism.* A rubber chain is a random coil, and stretching it reduces the number of its possible shapes. For $N$ links of length $b$, the end-to-end probability is Gaussian, $P(x)\propto e^{-3x^2/(2Nb^2)}$, so the free energy $-k_BT\ln P$ gives a spring

$$k_{\text{chain}}=\frac{3k_BT}{Nb^2}.$$

It is a spring made of entropy, and is stiffer when hot, the opposite of a metal. A network with $n$ chains per volume has $G\approx nk_BT$ (standard affine result, general knowledge). Assuming $\rho$ = 900 kg/m³ and 5 kg/mol between crosslinks gives $G$ = 0.45 MPa, $E$ = 1.35 MPa, the right order for soft rubber. This explains how rubber combines low stiffness with huge recoverable strain. Its chains uncoil without breaking bonds.

*What elasticity depends on, in summary:* bond type and spacing, crystal structure, temperature, loading rate compared with relaxation time, the defects present, and for structures their geometry. A coil spring has a low stiffness in newtons per metre while its steel has a high modulus. Those are properties of a shape and of a material respectively.

**Check:** [`code/verify_spring.py`](code/verify_spring.py), parts 2, 3 and 5.

```text
2. Elastic limit of a perfect chain, by searching for equilibrium under a dead load
     largest load with an equilibrium: 2.3964 eps/sig, at inflection predicted 2.3964 eps/sig, strain 10.9 %
  [PASS] breaking load equals the force at the potential's inflection point
3. Argon bulk modulus from a Lennard-Jones fcc lattice sum (static, T = 0)
     nearest neighbour 1.0906 sig, B = 74.76 eps/sig^3 (independent numerical second derivative)
     eps/k=119.8 K, sig=3.405 A (widely used): B = 3.13 GPa, measured 2.88 GPa, ratio 1.09
     eps/k=124 K, sig=3.42 A: B = 3.20 GPa, measured 2.88 GPa, ratio 1.11
     eps/k=125.7 K, sig=3.345 A (White 1999): B = 3.47 GPa, measured 2.88 GPa, ratio 1.20
  [PASS] gas/solid fitted pair potentials give argon's B within a factor 1.5
5. Entropic stiffness of a random chain (rubber), from sampled statistics
     fitted k = 0.03015 kT/b^2 for N = 100 links, Gaussian theory 3/N = 0.03000
  [PASS] entropic spring constant 3 kT / (N b^2) within 3%
     affine network G = n k_B T with n = 1.08e+26 /m^3 at 300 K: G = 0.45 MPa, E = 3G = 1.35 MPa
```

**Limits:** The argon test is one substance and uses parameters found in secondary sources. The value 119.8 K, 3.405 Å is from general knowledge. The test of $E\approx k/r_0$ is only for order of magnitude. The freely jointed chain is the simplest model of a polymer. Dislocation theory is not derived here. "100,000 times stiffer" compares 190 GPa with the order of 1 MPa.

### E05: Why do we bounce when we jump on a trampoline?

**Status:** Derived + Code + Estimate

**Short answer:** You fall, the bed and its springs stretch and store your energy, and they throw you back up, in about 0.4 s of contact. A passive body would bounce a little lower each time. You go higher because you straighten your legs against the bed while it is pushing hard on you, adding work $\int F\,dl$, and you retract your legs in the air, where it costs no force.

**Starting from:** Newton's laws. A bed of effective mass $M_b$ and stiffness $K$, a person whose centre of mass sits a controllable leg length $l(t)$ above the foot, and a one-sided foot contact. $K$ = 5000 N/m is the mid-range value used by Srinivasan et al. for their passive trampoline model. $M_b$ = 5 kg, leg stiffness 20 kN/m, a 0.3 m leg stroke and the bed loss are assumed.

**Derivation:**

*Passive bounce.* With a massless bed, $m\ddot x=-Kx-mg$ during contact, about the sagged equilibrium $-mg/K$. Starting from $x=0$ at speed $v$, the time in contact is

$$t_c=\frac{\pi+2\arctan\!\big(g/(\omega v)\big)}{\omega},\qquad \omega=\sqrt{K/m}.$$

Contact lasts a little over half a period, and only weakly depends on height. With $K$ alone that is 0.375 s from 3 m. Adding the leg spring in series gives 0.423 s. Reports of about 0.4 s of contact in trampolining are consistent (Blajer and Czaplicki, as cited by Burke). The peak force from height $h$ follows from energy, $mg(h+d)=\tfrac12Kd^2$:

$$F_{\max}=mg+\sqrt{(mg)^2+2Kmgh},$$

which is 6.9 body weights from 2 m and 12.7 from 8 m for this linear bed.

*Coupled masses.* The person and bed are two masses joined by a one-sided contact. When the foot lands, the bed mass must be brought up to speed, which costs an inelastic share of about $M_b/(M_b+m)$ of the kinetic energy if the contact is lossy. This loss and the bed's own damping make each passive bounce lower.

*Gaining height.* The force on the body comes through the leg, so the rate of energy input is $F\,\dot l$, and

$$\Delta E=\oint F\,dl-\text{losses}.$$

In flight $F=0$, so moving the legs does no work and doesn't change the path of the centre of mass. Under load $F$ is large, so a 0.3 m extension can add over a kilojoule. Over one cycle the leg returns to the same length, but the force differs on the two strokes, so the net work isn't zero. That asymmetry is the pump. It works like pumping a swing, with the work done where the force is high. The stroke should straddle the bottom of the bed, where $F$ is largest. The same logic in reverse explains how gymnasts stop: shorten the legs while the bed pushes up, so $F\dot l<0$. Height then grows until the work per bounce equals the loss per bounce.

**Check:** [`code/verify_trampoline.py`](code/verify_trampoline.py), part A. Velocity Verlet, with damping handled semi-implicitly, 0.5 ms step, checked at half that step. Energy audit: $\Delta E$ equals leg work minus heat to about 2%, and with no damping total energy drifts by $8\times10^{-6}$ of $mgh$.

```text
Part A: person on a bed (M_p = 60 kg, M_b = 5 kg, K = 5000 N/m)
A1. Passive person, stiff legs, no bed damping: energy and contact time
     apexes (m): [2.308 1.571 1.137 0.761 0.497 0.349]
     first contact time 0.431 s; theory 0.375 s for the bed spring alone, 0.423 s for bed and leg springs in series
  [PASS] contact time matches bed and leg springs in series, to 5%  (ratio 1.019)
     audit: (Delta E + heat)/heat = -1.1e-02
  [PASS] energy audit closes for the passive run
     no damping of any kind, dt = 0.1 ms: total energy drift 8.4e-06 of mgh
  [PASS] total energy conserved by the integrator with no damping
     step-size check: first apex 2.3102 m (dt/2) vs 2.3082 m
  [PASS] first apex converged with step size
A2. What a passive person loses per bounce
     M_b (kg)  bed zeta   apexes (m) from 3 m
        0.5     0.00     [2.8  2.58 2.38 2.2  2.05 1.91]
        5.0     0.00     [2.31 1.57 1.14 0.76 0.5  0.35]
        5.0     0.01     [2.18 1.56 1.11 0.75 0.49 0.32]
        5.0     0.03     [1.94 1.22 0.75 0.43 0.22 0.1 ]
  [PASS] passive bounces die away with damping  (1.94 -> 0.02 m)
A3. Pumping with the legs, bed zeta = 0.03, 0.3 m stroke, no force limit
     apexes from a 0.5 m start: [1.24 1.98 2.67 3.28 3.82 4.28 4.67 4.99 5.26 5.5 ] ... [6.16 6.23 6.29]
     leg work 23.5 kJ in 17 bounces, mean 1385 J per bounce
  [PASS] pumping gains height  (1.24 -> 4.28 m)
     audit: (Delta E - leg work + heat)/leg work = -1.7e-02
  [PASS] energy gained equals leg work minus heat  (-1.7e-02)
     Height after 60 s, where leg work per bounce has come near the loss per bounce:
     bed zeta  stroke (m)  force limit (body weights)  apex (m)
      0.030     0.30           inf                     6.50
      0.030     0.30           4.5                     1.69
      0.010     0.30           inf                    11.17
      0.010     0.30           4.5                     1.98
      0.003     0.30           4.5                     2.01
      0.010     0.15           inf                     2.24
     peak bed force on a person, massless linear bed: 6.9 body weights from 2 m, 12.7 from 8 m
  [PASS] height falls as bed loss grows, all else equal
  [PASS] force limit lowers the steady height
A4. Timing: apex after 5 bounces from 1 m, by when the leg extension starts (bed zeta = 0.03)
     start (quarter periods after touchdown)  stroke midpoint  apex after 5 bounces (m)  leg work per bounce (J)
     bed bottom with stiff legs is 1.32 quarter periods after touchdown; the stroke lasts 1.16
       0.00                                   0.58             2.66                      829
       0.25                                   0.83             3.46                     1241
       0.50                                   1.08             4.13                     1411
       0.75                                   1.33             4.31                     1206
       1.00                                   1.58             3.95                     1252
       1.25                                   1.83             3.15                     1002
       1.50                                   2.08             2.49                      594
       1.75                                   2.33             1.49                      436
       2.00                                   2.58             0.79                      250
  [PASS] best stroke straddles the bottom of the bed (largest force)  (best midpoint 1.33, bottom 1.32)
  [PASS] timing changes the outcome by a large factor  (0.79 to 4.31 m)
     'killing' the bounce (shorten the leg as the bed pushes up), from 3 m: apexes [0.7], leg work -1.22 kJ
     for comparison, stiff legs: apexes [1.94 1.22 0.75 0.43]
  [PASS] shortening the leg under load removes energy
```

The steady height ranges from about 2 m to 11 m across plausible assumptions, so this model shows the mechanism and doesn't predict a height. The force limit matters most. A person able to push only at 4.5 body weights stalls near 2 m, because the bed force is above the limit for most of the stroke.

**Limits:** One vertical degree of freedom, a linear bed, a prescribed leg program, massless legs. Real beds stiffen or soften with depth, and a real body is a set of linked segments with muscle force-speed limits. Burke, Yeadon and Hiley use a more complete bed model (seen only in a search excerpt). Bed mass and damping are assumed and were not fitted to any measured trampoline. Srinivasan et al. show that two people on one bed can transfer nearly all energy between them in a bounce, which this single-person model leaves out.

### E06: If we throw a ball on a trampoline, does the ball bounce?

**Status:** Derived + Code + Estimate

**Short answer:** Yes, but usually lower than on a hard floor. The model gives about $e$ = 0.45 to 0.5 on the bed for a tennis ball, a superball or a basketball, compared with 0.75 to 0.85 on a floor, because about 70% of the ball's energy goes into waves on the bed and hasn't returned by the time the ball leaves. A heavy object gets that energy back, because it stays in contact long enough for the waves to return from the rim.

**Starting from:** A circular membrane under uniform tension $T$ and areal density $\sigma$ obeying the wave equation with a small damping. A ball with its own contact stiffness and damping, pressing on a patch of the membrane. Bed values are assumed: radius 1.5 m, $\sigma$ = 1 kg/m², and $T$ = 2031 N/m chosen to give $K$ = 5000 N/m under a 15 cm foot. Ball masses and stiffnesses are assumed typical values, with contact times of the right order.

**Derivation:**

*Springs in series.* If the bed were massless, it and the ball would act as two springs in series, each storing energy in inverse proportion to its stiffness. A soft bed would then hold nearly all of it, sparing the ball's own hysteresis. The effect is real and appears at large mass below. Engineers use the name "trampoline effect" for the same thing in bats and rackets (general knowledge, unverified here).

*But a light ball sees the bed's inertia.* A disturbance spreads over the bed at $c=\sqrt{T/\sigma}$ (45 m/s here), so in a time $t$ the ball can only be moving a mass of about $\sigma\pi(ct)^2$. The static stiffness of the membrane under a patch of radius $r$ is

$$K=\frac{2\pi T}{\ln(a/r)+\tfrac14},$$

but it takes a wave return time $t_{\text{ret}}=2a/c$ for the rim to be felt at all. Before then, energy put into the bed is carried away outward and is not available to the ball. Whether that energy comes back depends on a timing comparison:

$$\pi\sqrt{\frac{m}{K_s}}\ \text{ vs }\ \frac{2a}{c}\quad\Rightarrow\quad m^*=K_s\left(\frac{2a}{\pi c}\right)^2,$$

with $K_s$ the bed and ball contact in series. Below $m^*$ the ball is gone before the echo returns. Above $m^*$ it is still pressing when the echo arrives and takes the energy back. This is an impedance argument in its time-domain form: a light ball is poorly matched to a large bed, and the mismatch sends energy into waves and not into rebound.

*Why bounce is so similar for different balls.* On the bed the deformation is mostly in the bed, so the ball's own $e$ matters little. The rebound is set by what the membrane returns, and that depends only weakly (logarithmically) on the patch and the mass. A superball and a dead foam ball differ by a factor of 3 on a floor and by 1.3 on the bed.

*Rule of thumb.* A heavy body such as a person is far above $m^*$. A person with lossy legs doesn't bounce at all on a floor in this model, because gravity during the long contact cancels the weak rebound, yet the same person bounces well on the bed.

**Check:** [`code/verify_trampoline.py`](code/verify_trampoline.py), part B. An axisymmetric finite-volume wave equation with leapfrog time stepping. The solver's static stiffness and lowest mode are tested against the formula and the Bessel zero $j_{01}c/a$. Energy is conserved to 0.01% with no damping. The echo threshold is found by bisection on the mass at three bed radii, with no fitting.

```text
Part B: ball on a circular membrane, solved as a wave equation
     tension 2031 N/m, areal density 1.0 kg/m^2, radius 1.5 m (assumed), wave speed 45.1 m/s
     static stiffness under a 15 cm foot patch 4991 N/m (formula 2 pi T / (ln(a/r) + 1/4) = 5000)
  [PASS] numerical static stiffness matches the formula
     lowest mode 72.02 rad/s, Bessel-zero value 72.26 rad/s
  [PASS] membrane's lowest mode matches j01 c/a
     lossless, gravity-free run: final total energy / initial = 1.0001
  [PASS] wave solver conserves energy with no damping
     Named balls from 1 m. 'Floor' uses the same code with a rigid floor.
     ball                  m(kg)  k_c(kN/m)  e_floor  e_bed  contact: floor(ms)  bed(ms)  bed mass reached/m  KE left in bed at lift-off
     superball             0.050      490     0.852   0.452       0.96           9.9       0.12           76.4 %
     cricket ball          0.156      680     0.551   0.520       1.32          19.3       0.07           70.5 %
     tennis ball           0.058       23     0.753   0.463       4.65          11.0       2.37           72.6 %
     basketball            0.600       18     0.797   0.483      17.40          38.4       3.22           69.4 %
     dead foam ball        0.020        1     0.288   0.345      10.92          12.7      38.02           20.3 %
     medicine ball         5.000       50     0.270   0.902      25.12         114.1       0.81            9.6 %
     person, lossy legs   60.000       20     0.000   0.850    1911.17         437.3     388.48            0.5 %
  [PASS] gravity-free floor runs recover each ball's intended e to 0.01
     energy budget, rebound e^2 + bed share + heat in ball, each ball: [0.992 0.998 0.997 0.973 0.931 0.981 0.962]
  [PASS] rebound, bed and ball shares add to about 1 (bed damping and gravity in contact aside)
  [PASS] balls do bounce on the bed
  [PASS] a person with lossy legs bounces better on the bed than on a floor
  [PASS] a tennis ball and a superball bounce less on the bed than on a floor
     Mass scan: k_c = 20 kN/m, e_floor = 0.75, 5 cm patch
     m (kg)    e_floor   e_bed    contact (ms)   KE left in bed at lift-off
       0.01     0.759    0.252        3.6        72.0 %
       0.03     0.754    0.337        6.9        76.5 %
       0.10     0.751    0.432       14.2        73.3 %
       0.30     0.749    0.498       27.2        68.6 %
       1.00     0.745    0.547       55.1        64.2 %
       3.00     0.741    0.950      108.5         3.2 %
      10.00     0.733    0.962      193.8         1.9 %
      30.00     0.718    0.967      358.6         1.1 %
      60.00     0.704    0.969      319.5         0.9 %
     lowest bounce at m = 0.01 kg (e_bed = 0.252); highest at m = 60.0 kg (e_bed = 0.969)
  [PASS] in the heavy range, more mass does not lose the bed's rebound advantage
  [PASS] rebound rises from light masses to about 10 kg
     Is the jump set by the echo from the rim? Mass where e_bed crosses 0.75, three bed radii.
     Prediction: contact time pi sqrt(m/K_s) equals the wave return time 2a/c
     radius (m)  return time (ms)  predicted m* (kg)  simulated m* (kg)
       1.00          44.4              0.65              0.69
       1.50          66.6              1.33              1.45
       2.25          99.8              2.74              3.03
  [PASS] threshold mass within a factor 3 of the echo prediction at all radii
  [PASS] threshold mass grows with bed radius, as an echo requires  (0.69 -> 3.03 kg)
     Sensitivity of the tennis ball to the assumed bed density:
       areal density 0.5 kg/m^2: e_bed = 0.508
       areal density 1.0 kg/m^2: e_bed = 0.463
       areal density 2.0 kg/m^2: e_bed = 0.411
     Grid check on the tennis ball, dr = 1 cm vs 0.5 cm:
     e_bed 0.463 vs 0.459
  [PASS] tennis-ball result stable under grid refinement
```

The simulated $m^*$ tracks the prediction at each radius within about 10%. The floor value falls with mass in the scan because gravity acts during a long contact.

**Limits:** A real bed is rectangular, woven and hung from discrete springs, so the rim reflects only part of the wave and at spread-out times. The sharp jump in the scan will be smoother on real equipment, and no real trampoline was measured. The tennis ball result moves from 0.51 to 0.41 as the assumed bed density changes from 0.5 to 2 kg/m². Ball contact stiffness is taken as a linear spring. "Bed mass reached" in the table uses the floor contact time.

### E07: Will a ball between two opposed springs keep going?

**Status:** Derived + Code

**Short answer:** In the ideal model yes: with massless, lossless springs and no air the motion repeats forever. Any real damping makes it decay, and even with no damping at all, springs that have their own mass drain the ball's motion into their internal vibration. A ball exactly on the axis may also be pushed sideways, since compressed springs can buckle.

**Starting from:** Newton's laws. Two identical springs that can push but never pull. Motion along the axis, plus a sideways coordinate in part 4. Units $m=k=1$. Gravity is left out, and would just shift the centre by $mg/k$.

**Derivation:**

*The contact is one-sided, so the system is not automatically a linear oscillator.* Each spring acts only when compressed, $F=-k\max(0,\,\cdot)$. Three cases follow.

1. **Just touching at the centre.** For $x>0$ only the top spring acts and $F=-kx$. For $x<0$ only the bottom one acts and $F=-kx$. The two pieces join into one straight line, so this special case is exactly linear, with $T=2\pi\sqrt{m/k}$ whatever the amplitude. So your picture of one spring and then the other is right.
2. **A gap of half-width $g$.** The ball coasts between springs. If $A$ is the furthest point reached, the coasting speed is $\omega(A-g)$, and $T=2\pi/\omega+4g/[\omega(A-g)]$. The period now depends on amplitude.
3. **Both preloaded by $p$.** For $|x|<p$ both act, so the stiffness is $2k$. Beyond $p$ one has lost contact, so the stiffness is $k$. The period is a nonlinear function of amplitude, given by $T=4\int_0^A dx/\sqrt{2(V(A)-V(x))/m}$.

In all three, energy is conserved, so $x(t)$ is periodic forever.

*Decay.* Three damping laws give three different decays, from a work-per-cycle argument:

- A dashpot inside each spring: $A\propto e^{-\zeta\omega t}$. The spring can lag behind the ball when it recovers. It lets go when its force reaches zero, while still compressed, and that stored energy turns into heat. Separation is real, for 2% to 21% of the time in the cases below, and yet the decay rate stays within about 3% of $\zeta\omega$.
- Dry friction $F$: amplitude falls by $2F/k$ per half cycle, linearly, and the ball stops for good in finite time once $|x|\le F/k$.
- Air drag $bv|v|$: energy lost per cycle is $\tfrac83 bA^3\omega^2$, so $dA/dt=-\tfrac{4}{3\pi}(b\omega/m)A^2$ and $A\propto 1/t$. It never reaches exactly zero in this model.

*Without any damping, the ball's motion still isn't kept.* A real spring has mass and a large number of internal modes. When the ball leaves, the spring is still vibrating. Energy stays in the whole system, but it moves out of the single ball motion you were picturing. The first-cycle share depends on spring mass over ball mass: 0.99, 0.95 and about 0.7 for ratios of 0.01, 0.1 and 1. Internal friction would then turn the vibrations into heat.

*Going sideways.* The exact on-axis path is a solution, but is it stable? Let the springs act through hinged ends a distance $L$ from the ball, with preload force $P$ and a lateral stiffness $k_s$ each. For a sideways offset $y$ the spring length is $r\approx L+y^2/(2L)$, so

$$V(y)\approx \text{const}+\left(k_s-\frac{P}{L}\right)y^2.$$

If $k_s<P/L$, the on-axis position is an energy maximum and the ball is thrown out sideways, exponentially at rate $\sqrt{2(P/L-k_s)/m}$. This is buckling, in the form Haringx treated for helical springs. Perfect centring doesn't rescue it, since any offset grows. A small enough $P/L$ is safe, and that depends on the spring's slenderness. The FE table below gives $k_s/k$ from 1.2 for a squat coil down to 0.06 for a slender one.

Oscillation adds a second route. With touching springs and no preload, the compression varies with time, so

$$\ddot y+\left[k_s-\frac{k|z|}{L_0-|z|}\right]\frac{y}{m}=0,\quad z=A\cos\omega t.$$

The bracket repeats at $2\omega$, which is a Hill or Mathieu equation. So the sideways motion can grow by parametric resonance, when the lateral frequency is near the bounce frequency. Floquet analysis maps those bands, and the full 2D simulation confirms the growth rate. At amplitude $A/L_0$ = 0.3 there are two main bands. One is $k_s/k$ below about 0.26, which slender coils (length/diameter of about 4 or more in the table) fall into. The other is a narrow resonance at 1.16 to 1.37, close to where a squat coil sits. Between them the on-axis motion is stable. A lateral damping ratio of 0.02 removes little of this.

A practical remark: a small ball can drop inside a coil, so it would need a plate on each spring end. A sphere between two flat plates then has nothing holding it on the axis except friction and the springs' lateral stiffness.

**Check:** [`code/verify_two_springs.py`](code/verify_two_springs.py). Symplectic integration for conservative runs, with each period compared to an independently computed value. Chains of 12 and 24 masses for the massive spring, with one-sided stiff contact. Energy audits on the dissipative runs.

```text
1. Massless springs, one-sided contact (no damping)
   a. touching: mean period 6.283179 (2 pi = 6.283185), energy drift 6.2e-06 over 100 periods
  [PASS] touching springs are exactly a linear oscillator
  [PASS] no amplitude loss over 100 periods
   b. gap of half-width g: T = 2 pi + 4 g / (A - g) (A = turning point), simulated:
      g = 0.05, A = 0.30: period 7.0832, predicted 7.0832, energy drift 6.3e-06
      g = 0.10, A = 0.30: period 8.2832, predicted 8.2832, energy drift 6.3e-06
      g = 0.10, A = 0.15: period 14.2832, predicted 14.2832, energy drift 6.4e-06
  [PASS] gap makes the period depend on amplitude, as derived
   c. preload p: stiffness 2k inside |x| < p, k outside. Period from the energy integral:
      p = 0.10, A = 0.05: period 4.4429, energy integral 4.4429
      p = 0.10, A = 0.20: period 4.9593, energy integral 4.9593
      p = 0.10, A = 0.40: period 5.4838, energy integral 5.4838
  [PASS] preloaded system is nonlinear: period from simulation equals the energy integral

2. Springs with their own mass (chains of 12 masses; contact is one-sided)
   spring mass/ball mass  first period  rule 2 pi sqrt(1 + ms/3)  energy drift   ball's energy share: cycle 1  late mean  late min
       0.01                   6.348         6.294                     3.3e-04        0.990          0.951   0.941
       0.10                   6.457         6.387                     6.4e-05        0.954          0.853   0.833
       1.00                   7.130         7.255                     1.9e-04        0.685          0.075   0.053
  [PASS] total energy conserved in the chain model (drift below 0.5%)
  [PASS] first period follows the one-third-of-spring-mass rule within 5% for ms <= 0.1
      mesh check, share after the first cycle and late-time mean:
        ms = 0.1, N = 24: first cycle 0.956 (N = 12: 0.954), late mean 0.915, drift 2.8e-05
        ms = 1.0, N = 24: first cycle 0.706 (N = 12: 0.685), late mean 0.109, drift 7.9e-05
  [PASS] loss in the first cycle is mesh-independent to 0.03
  [PASS] light spring keeps nearly all the energy in the ball
  [PASS] ms = 0.1: the ball slowly leaks energy into the spring (late share 0.6-0.95, not converged in detail)
  [PASS] ms = 1: most energy has left the ball within a few cycles, with none lost overall

3. Three kinds of damping, three decay laws (A = 0.3, touching springs, one-sided)
   a. dashpot in the spring: amplitude should fall as exp(-zeta t) if contact never broke
      zeta   fitted decay rate  rate/zeta  fraction of time a compressed spring is out of contact  energy audit
      0.005       0.00498       0.996           0.024                                          1.0002
      0.020       0.01973       0.987           0.101                                          1.0003
      0.050       0.04838       0.968           0.213                                          1.0003
  [PASS] exponential decay close to the linear rate, energy audit closes
   b. dry friction F = 0.01 k: amplitude drop per half cycle 0.0200 (derived 2F/k = 0.0200); stops at t = 47.1, derived about 15 half cycles = 47.1
  [PASS] dry friction: linear decay and a finite stopping time
   c. air drag F = b v|v|, b = 0.5: slope of 1/A vs t = 0.21221 (derived 4 b / (3 pi) = 0.21221)
  [PASS] air drag: amplitude falls as 1/t, never exactly zero

4. Sideways stability
   a. preload strain 0.111, no lateral stiffness: growth rate 0.4713, derived sqrt(2 P / (L m)) = 0.4714
  [PASS] ball exactly on the line is unstable under preload if the springs have no lateral stiffness
      threshold k_s/k = P/(k L) = 0.1111: at 0.9 of it y reaches 1.8e-01, at 1.1 of it y stays at 1.0e-06
  [PASS] buckling threshold k_s = P/L
   b. What is k_s for a real coil? Finite-element helix (music wire, d 2 mm, D 20 mm, pitch 6 mm), ends kept parallel
      coils  length/D   k_axial (N/m)  k_lateral (N/m)  lateral/axial
         4     1.20          4559            5542        1.216
         7     2.10          2605            1521        0.584
        10     3.00          1824             590        0.324
        15     4.50          1216             188        0.155
        25     7.50           729              42        0.058
  [PASS] slender springs have much less relative lateral stiffness (they buckle earlier)
   c. Touching, no preload, A/L0 = 0.3: lateral motion y'' + [k_s - |z|/(L0-|z|)] y = 0, z = A cos t
      k_s/k:       0.02   0.05   0.10   0.15   0.20   0.30   0.40   0.50   0.60   0.70   0.80   0.90
      growth/s:   0.477  0.445  0.384  0.312  0.217 -0.000 -0.000 -0.000 -0.000 -0.000 -0.000 -0.000
      k_s/k:       1.00   1.10   1.20   1.30   1.50   2.00   3.00   3.50   3.90   4.10   4.50   6.00
      growth/s:  -0.000 -0.000  0.041  0.041 -0.000 -0.000 -0.000 -0.000 -0.000 -0.000 -0.000 -0.000
      unstable bands of k_s/k, undamped: 0.02-0.27, 1.18-1.38
      with lateral damping ratio 0.02 (relative to bounce frequency): 8 % of the scanned range stays unstable (undamped 8 %)
      strongest band above k_s/k = 0.3: k_s/k = 1.23, growth 0.047
  [PASS] resonant band near k_s = k (lateral frequency = bounce frequency)
      full 2D simulation, 1e-7 sideways offset: at k_s/k = 1.23 growth 0.046 (Floquet 0.047); at k_s/k = 1.53 y stays below 1.0e-07
  [PASS] full simulation confirms Floquet growth rate within 25%
  [PASS] away from the bands the on-line motion stays stable

All checks passed.
```

**Limits:** The late-time energy share in the massive spring is not mesh-converged in detail (0.85 and 0.92 at $m_s/m$ = 0.1 with 12 and 24 masses, and it wanders in time), though the first-cycle loss is mesh-independent. I report only the coarse trend. The hinged spring with an assigned $k_s$ is a simplification, and a real coil couples axial and sideways motion through its twist. The lateral stiffness table holds the spring ends parallel. No physical apparatus was tested, and no buckling load of a real coil was computed.

### E08: Can we calculate bounce height from release velocity and angle alone?

**Status:** Derived + Code + Estimate

**Short answer:** Partly. The angle at which the ball strikes the ground is not an unknown: gravity fixes it from release speed, angle and height, so that part of your question has a definite answer. What you can't get from those alone is what the pitch does, mainly its restitution and how much it deforms, and also its friction and the ball's seam and spin. With release data only, the height at the batter has a 5 to 95% range of about 0.8 m for a short ball, wider than the 0.4 m separating chest from head, so you can't call it. Calibrated to the day's pitch, the range is about 0.3 to 0.45 m, which decides only when the ball is well inside one band.

**Starting from:** Projectile motion with drag. A rigid sphere with restitution $e$ on the normal speed and Coulomb friction $\mu$. Cricket dimensions from the MCC Laws: ball 155.9 to 163 g, circumference 22.4 to 22.9 cm, popping crease 1.22 m in front of the stumps. Pitch 20.12 m (22 yards, general knowledge). I take release above the bowler's popping crease and the batter at the other, 17.68 m away. Release height $H$ = 2.2 m and $C_d$ = 0.45 are assumed.

**Derivation:**

*Release to pitch, without air.* With speed $v$ at angle $\alpha$ below horizontal, the fall gives

$$v_y=\sqrt{v^2\sin^2\alpha+2gH},\quad v_x=v\cos\alpha,\quad t=\frac{v_y-v\sin\alpha}{g},\quad x_p=v_xt,$$

$$\tan\theta_{\text{strike}}=\frac{\sqrt{v^2\sin^2\alpha+2gH}}{v\cos\alpha}.$$

So the strike angle needs nothing but $v,\alpha,H$. For 40 m/s it is 10° to 15° across all sensible lengths. It agrees with the "average seam bowler" of James et al., as quoted by Howard, of 36 m/s at 12.5°. Drag at 40 m/s is about 11 m/s² for $C_d$ = 0.45, comparable to gravity, but acts for only a few tenths of a second, so it trims the speed by several percent.

*The bounce.* Normal: $v_y'=e\,v_y$. Tangential: let the bottom of the ball slip at $u=v_x+R\omega$, with $\omega$ the backspin. A friction impulse $J$ per unit mass changes $u$ by $-(1+1/\beta)J$, where $I=\beta mR^2$ and $\beta=2/5$ for a uniform sphere. Friction can stop the slip only if it has enough impulse, so

$$J=\min\!\left(\frac{|u|}{1+1/\beta},\ \mu(1+e)v_y\right),\quad v_x'=v_x-J.$$

Fast bowling at usual friction is usually in the slip regime. A ball that does grip leaves with $v_x'=\tfrac57v_x$ if it had no spin. Cross measured such grip in bounces, and found a horizontal restitution of 0.1 to 0.2 for most balls and 0.5 to 0.6 for a superball, so a uniform grip rule is itself a simplification.

*Height at the batter.* Over the remaining distance $d$, $t_2\approx d/v_x'$ and $h=v_y't_2-\tfrac12gt_2^2$. So height depends on $e$ directly and on length through $v_y$ and $d$. The maximum height $v_y'^2/2g$ may not be reached before the ball arrives, so "how high it bounces" and "how high it is at the bat" are different numbers.

*What the pitch adds.* Howard measured $e\approx$ 0.3 to 0.4 on lab-made pitches at about 6.8 m/s, falling suddenly when brittle cracking occurred. He reports James et al.'s finding that a rigid-body model predicts rebound speed well but "always underestimates" the rebound angle, because the pitch deforms and kicks the ball up. James et al. also found that no single pitch measurement reliably predicts pace or bounce. I represent that missing effect by rotating the rebound upward by a kick angle at constant speed, with an assumed range of 0 to 5°. This is a labelled stand-in, not a fitted law. At 9 m length, 1° of kick adds 0.14 m, so this unmeasured number matters as much as anything else.

*What else you can't get from release data.* The seam: Howard, citing Carré et al., reports more deformation for landing on the seam than across it. Spin, including topspin gained at the bounce, gives a downward Magnus force. Swing: Grimshaw et al. give a sideways 0.52 m at 85 mph for $C_S$ = 0.3 bouncing 15 m from release. That changes where the ball lands sideways, and hardly its height.

*Sensitivities.* The run below shows, at 9 m length: +0.05 of $e$ gives +0.23 m, +1° of kick gives +0.14 m, ±0.5° of release angle gives only ±0.04 m at fixed nominal length, and +1 m/s of speed gives almost nothing. Speed and angle mostly decide where the ball lands, and the pitch decides how high it comes up from there.

**Check:** [`code/verify_cricket.py`](code/verify_cricket.py). A vectorised RK4, checked against adaptive DOP853 and by halving the step. The impact model is checked on 20000 random impacts, none of which gains energy, and against the 5/7 limit. As an outside comparison, the rigid-impact model is run on the published delivery, and it reproduces the range of heights (0.15 to 0.40 m) and pace (0.76 to 0.90) that Howard shows for that same kind of model. This confirms the code, not the physics, since Howard states the rigid model is not accurate for height.

```text
Ball: m = 0.159 kg, R = 3.61 cm (MCC Law 4). Batter at x = 17.68 m from release, stumps at 18.90 m.

1. Release to pitch, no air: strike angle from speed, angle and height
     v0 (m/s)  angle  H (m)   x_pitch(m)  strike angle(deg) sim   closed form
        40      3.0  2.20      19.555        9.792              9.792   (x closed 19.555)
        40      7.0  2.20      13.381       11.642             11.642   (x closed 13.381)
        30      3.0  2.10      15.388       12.437             12.437   (x closed 15.388)
        25     -2.0  2.00      18.334       14.203             14.203   (x closed 18.334)
        40      0.0  2.20      26.793        9.326              9.326   (x closed 26.793)
  [PASS] strike angle is fixed by v0, alpha and H under gravity alone

2. Integrator checks with drag (Cd = 0.45)
     RK4 x_pitch 15.72495 m, scipy DOP853 15.72496 m; vy -6.99107 vs -6.99107 m/s
  [PASS] vectorised RK4 agrees with an adaptive integrator
     step-size: x_pitch 15.722227 (1 ms) vs 15.722232 (0.25 ms)
  [PASS] converged with step size
3b. Against the published 'average seam bowler' of James et al. (2005), as used by Howard (2012):
     36 m/s at 12.5 deg at the pitch, bouncing about 5 m in front of the batter, no drag after the bounce
       e = 0.30, mu = 0.3: height at batter 0.245 m, pace (speed after/before) 0.894
       e = 0.35, mu = 0.3: height at batter 0.306 m, pace (speed after/before) 0.892
       e = 0.40, mu = 0.3: height at batter 0.368 m, pace (speed after/before) 0.890
       e = 0.30, mu = 0.6: height at batter 0.257 m, pace (speed after/before) 0.810
       e = 0.40, mu = 0.6: height at batter 0.395 m, pace (speed after/before) 0.799
     Howard's Figures 32-33 span heights 0.15-0.40 m and pace 0.76-0.90, and he states the rigid model underestimates the rebound angle
  [PASS] rigid model reproduces the height and pace range published for the same rigid model
     release 39 m/s, H = 2.2 m, alpha 7.49 deg gives strike 12.5 deg, arrival speed 36.2 m/s, length 6.34 m from stumps
  [PASS] release-to-pitch model gives the published arrival speed and angle at a sensible length

4. A 40 m/s delivery, H = 2.2 m, e = 0.40, mu = 0.35, Cd = 0.45. Height at the batter with no kick and with a 2.5 and 5 deg kick
     alpha(deg)  x_pitch(m)  length from stumps(m)  strike(deg)  v_y in(m/s)  rebound(deg)  v_x out  h, no kick(m)   h, 2.5 deg   h, 5 deg
        1.0        22.91        -4.01               10.40        6.23        4.61     30.90     -0.421*       -0.652*     -0.885*
        2.0        20.83        -1.93               10.41        6.33        4.62     31.34     -0.253*       -0.392*     -0.532*
        3.0        18.93        -0.03               10.54        6.49        4.68     31.70     -0.102*       -0.157*     -0.213*
        4.0        17.24         1.66               10.78        6.71        4.80     31.98      0.036         0.056       0.076
        5.0        15.72         3.18               11.11        6.99        4.97     32.17      0.152         0.238       0.325
        6.0        14.39         4.51               11.54        7.32        5.18     32.28      0.247         0.392       0.538
        7.0        13.21         5.69               12.04        7.70        5.44     32.32      0.330         0.527       0.726
        8.0        12.16         6.74               12.61        8.11        5.74     32.30      0.407         0.650       0.895
       10.0        10.43         8.47               13.90        9.04        6.43     32.09      0.558         0.877       1.199
       12.0         9.06         9.84               15.35       10.05        7.23     31.69      0.716         1.096       1.479
     * pitched beyond the batter: a full toss, the height shown is then not a bounce
5. Sensitivity for a short ball: alpha set so it lands 9 m from the stumps at 40 m/s
     alpha = 10.72 deg, strike 14.40 deg, v_y in 9.39 m/s, height at batter 0.613 m
     input               step              change in height at batter (m)
     release speed       +1 m/s            +0.001
     release angle       +0.5 deg          +0.039
     release angle       -0.5 deg          -0.039
     release height      +0.1 m            -0.007
     restitution e       +0.05             +0.234
     friction mu         +0.10             +0.013
     drag Cd             -0.10             +0.008
     lift Cl (backspin)  +0.10             -0.111
     spin                +10 Hz            +0.000
     deformation kick    +1 deg            +0.137
  [PASS] height sensitive to e by more than 0.1 m per 0.05
  [PASS] height sensitive to the unmodelled deformation by more than 0.1 m per degree
6. Uncertainty on height at the batter, Monte Carlo, N = 5000, four target lengths
     A: v0, angle, height known exactly; pitch and ball unknown over the ranges below
     B: as A, plus realistic release measurement error (v0 +/- 1 m/s, angle +/- 0.5 deg, H +/- 0.05 m, normal 1 sigma)
     C: pitch calibrated on the day (e 0.40 +/- 0.02, mu 0.35 +/- 0.05, kick 2.5 +/- 0.5 deg), delivery-to-delivery e scatter +/- 0.03, release exact
     D: as C, plus release measurement error as in B
     ranges for A, B: e 0.30-0.50, mu 0.20-0.60, kick 0-5 deg, Cd 0.30-0.55, Cl 0-0.1, spin 5-15 Hz (uniform, assumed)
     length(m)  case   median h(m)   5%..95% (m)   width(m)  P(below chest) P(chest 1.1-1.5)  P(head 1.5-1.9)  P(over 1.9)
        7.0      A      0.627      0.36 ..  0.90    0.54       1.00           0.00           0.00         0.00
        7.0      B      0.630      0.35 ..  0.94    0.59       1.00           0.00           0.00         0.00
        7.0      C      0.624      0.52 ..  0.73    0.21       1.00           0.00           0.00         0.00
        7.0      D      0.619      0.49 ..  0.77    0.29       1.00           0.00           0.00         0.00
        9.0      A      0.894      0.51 ..  1.30    0.79       0.78           0.22           0.00         0.00
        9.0      B      0.897      0.51 ..  1.32    0.81       0.77           0.22           0.01         0.00
        9.0      C      0.893      0.75 ..  1.04    0.30       0.99           0.01           0.00         0.00
        9.0      D      0.889      0.72 ..  1.07    0.35       0.98           0.02           0.00         0.00
       11.0      A      1.265      0.73 ..  1.83    1.10       0.32           0.43           0.23         0.03
       11.0      B      1.271      0.73 ..  1.84    1.11       0.31           0.43           0.22         0.03
       11.0      C      1.260      1.05 ..  1.47    0.42       0.11           0.86           0.03         0.00
       11.0      D      1.256      1.03 ..  1.49    0.46       0.13           0.83           0.04         0.00
       13.0      A      1.907      1.15 ..  2.70    1.55       0.04           0.18           0.28         0.51
       13.0      B      1.904      1.14 ..  2.70    1.56       0.04           0.17           0.29         0.50
       13.0      C      1.883      1.58 ..  2.19    0.62       0.00           0.02           0.51         0.46
       13.0      D      1.878      1.56 ..  2.22    0.66       0.00           0.02           0.52         0.46
  [PASS] with only release data the height at the batter spreads over more than one 0.4 m band  (5-95% width 0.79 m)
  [PASS] calibrating the pitch narrows the spread by at least a factor 2 if the release is exact  (0.79 -> 0.30 m)
  [PASS] release error alone keeps the calibrated forecast from a clean chest-or-head call
7. Release height at the same length, 6.82 m from the stumps, 37 m/s, 2.5 deg kick
     H = 1.9 m: alpha  6.32 deg, strike 11.67 deg, height at batter 0.588 m, apex 0.859 m
     H = 2.1 m: alpha  7.25 deg, strike 12.59 deg, height at batter 0.632 m, apex 0.946 m
     H = 2.3 m: alpha  8.17 deg, strike 13.50 deg, height at batter 0.676 m, apex 1.037 m
     H = 2.5 m: alpha  9.08 deg, strike 14.41 deg, height at batter 0.721 m, apex 1.132 m
     model: 0.22 m of bounce per m of release height at this length
  [PASS] taller release gives higher bounce at the same length

8. Speed at the same 9 m length, and a spinner
     v0 = 25 m/s: alpha  7.92 deg, strike 17.19 deg, v_y out  3.70 m/s, height at batter 0.665 m, apex 0.698 m
     v0 = 30 m/s: alpha  9.32 deg, strike 15.81 deg, v_y out  4.15 m/s, height at batter 0.816 m, apex 0.879 m
     v0 = 35 m/s: alpha 10.16 deg, strike 14.96 deg, v_y out  4.64 m/s, height at batter 0.902 m, apex 1.097 m
     v0 = 40 m/s: alpha 10.72 deg, strike 14.40 deg, v_y out  5.15 m/s, height at batter 0.955 m, apex 1.351 m
     v0 = 45 m/s: alpha 11.09 deg, strike 14.02 deg, v_y out  5.67 m/s, height at batter 0.991 m, apex 1.640 m
     spinner 22 m/s, 30 Hz topspin, mu 0.6, 4 m length: regime grip, v_x 19.8 -> 16.1 m/s, height at batter 0.27 m
```

*Honest error bar for head or chest.* Taking chest as 1.1 to 1.5 m and head as 1.5 to 1.9 m for a 1.8 m batter in stance (assumed):

| Length from stumps, 40 m/s | Release data only (A) | Pitch calibrated, exact release (C) | Calibrated with release error (D) |
|---|---|---|---|
| 9 m | 0.51 to 1.30 m | 0.75 to 1.04 m | 0.72 to 1.07 m |
| 11 m | 0.73 to 1.83 m | 1.05 to 1.47 m | 1.03 to 1.49 m |
| 13 m | 1.15 to 2.70 m | 1.58 to 2.19 m | 1.56 to 2.22 m |

So from release speed and angle alone you can say "this is short, so it will come up high" and you can't say "chest, not head". With a pitch calibrated on the day you often can, unless it lies near the boundary between the bands. Those absolute heights are likely to be low, since the 0 to 5° kick range is my assumption and I found no ball-tracking heights to compare with. The width is the more reliable output.

**Limits:** The ranges for $e$, $\mu$, kick and $C_d$ are assumptions informed by the sources, not measured distributions. Howard's $e$ values come from a final-year undergraduate project using lab-made pitches, at one impact speed, so real match pitches may differ. $e$ varies with impact speed, and $\mu$ with moisture and grass. I treat a real ball as a uniform sphere. The Monte Carlo doesn't include wear or uneven bounce from cracks: Howard's model suggests bounce height varying by up to 22.8% on old dry pitches, which would widen case C. I could not verify a measured drag coefficient for a cricket ball, nor any measured height at the batter.

### E09: In a room with bouncy walls, will a ball keep having infinite collisions?

**Status:** Derived + Code

**Short answer:** It depends on $e$ and on gravity. With perfectly elastic walls ($e=1$) yes: collisions go on forever at a steady rate. With $e<1$ under gravity, the ideal equations give infinitely many bounces in a *finite* time (a Zeno sequence), at which point the model stops being valid. A real ball makes a finite number of bounces and comes to rest near that same time. With $e<1$ and no gravity, the ideal ball still hits the walls forever, but ever more slowly.

**Starting from:** Newton's laws. Restitution $e$ applied to the normal velocity, with instantaneous impacts in parts 1, 4 and 5, and a spring-damper contact in part 2. A frictionless floor and no air unless stated.

**Derivation:**

*Under gravity, $e<1$.* A drop from $h_0$ takes $t_0=\sqrt{2h_0/g}$. Each bounce multiplies the speed by $e$, so the $n$-th rebound rises to $h_n=e^{2n}h_0$ and spends $2e^nt_0$ in the air. Summing the geometric series,

$$T=t_0+2t_0\sum_{n=1}^{\infty}e^n=t_0\,\frac{1+e}{1-e},\qquad D=h_0+2h_0\sum_{n=1}^\infty e^{2n}=h_0\frac{1+e^2}{1-e^2}.$$

Infinitely many bounces, in finite time and over finite distance. With $e=0.8$, $h_0$ = 1 m: $T$ = 4.06 s, $D$ = 4.56 m. There is no mathematical paradox, since an infinite sum can be finite. But the ball then has zero speed at $t=T$, and the equations don't define what comes after. The rule "reverse the velocity at contact" has to be replaced by "stays on the floor". Control theorists call this behaviour Zeno, after the same paradox.

*What happens in reality.* The model assumes contact takes no time. A real contact lasts about $\tau=\pi\sqrt{m/k}$, while the flight time $2e^nv_0/g$ keeps shrinking. They cross at

$$n^*=\frac{\ln\big(2v_0/(g\tau)\big)}{\ln(1/e)},$$

after which the ball no longer leaves the floor but oscillates on it, damping out at the static sag $mg/k$. Falcon et al. saw this in experiments: a finite number of bounces, with the bead ending in oscillation on the deformed surface when flight time becomes comparable to contact time. They also found the effective $e$ roughly constant over most bounces and then falling to zero near the end, and the simulation below shows the same trend. The rest time is close to the ideal Zeno time, because all the later bounces add little.

Other limits come later: height drops below an atom's size after about 52 bounces (for $e$ = 0.8, $h_0$ = 1 m), but contact time ends the sequence first, at about 26.

*Under gravity, $e=1$.* Periodic motion for ever. The bounce interval stays at $2t_0$, so the count grows without limit, but only as time grows without limit.

*No gravity, $e=1$.* A straight line of steady speed between walls, at $v_0/L$ collisions per second, for ever.

*No gravity, $e<1$.* This is a different case from Zeno. Speed falls as $e^n$, and each crossing takes longer, $L/(v_0e^n)$. The gaps grow, so the sum diverges. There is no finite-time accumulation, and the ball hits walls for ever, in the ideal model, with

$$t_n=\frac{L}{v_0}\,\frac{e^{-n}-1}{1/e-1},\qquad v\to\frac{eL}{(1-e)\,t}.$$

So gravity is what makes the collisions pile up, because it brings the ball back quickly however slow it is. Without it, a slow ball takes a long time to arrive.

*Add air.* With linear drag $b$, a ball that keeps speed only by bouncing with $e=1$ can travel only $mv_0/b$ in total, so the count is finite, about $mv_0/(bL)$. Quadratic drag alone gives $v\propto1/t$ and a distance growing as $\ln t$, so the ball would continue, more and more slowly. At low speed real air drag is linear, so a real ball stops, and at that point only thermal motion is left. A 50 g ball in still air would take about $3\times10^7$ years to wander 1 m by Brownian motion (Stokes estimate), which I don't count as bouncing.

*The room itself.* In a box with frictionless walls the three axes don't affect each other, so each can be solved alone. With gravity and $e<1$, the vertical motion stops first. It collapses at a finite time, after a finite number of ceiling hits and a Zeno sequence on the floor. The horizontal motion is unaffected by that in this idealisation, and the ball slides into the walls at ever longer gaps. In a real room, friction at each floor contact and rolling resistance stop the sideways motion too.

So the real answer to your question: real materials have $e<1$, so a real ball stops. No material returns all its energy, so the infinite case is an idealisation.

**Check:** [`code/verify_room.py`](code/verify_room.py). Exact flights between events, without using the series formula, compared to the series. A phase-by-phase DOP853 integration of the finite contact. A separate event simulation of each axis of the room.

```text
1. Point ball under gravity, h0 = 1.0 m, e = 0.8
   t0 = 0.4516 s; series limit T = t0 (1+e)/(1-e) = 4.064407 s; event-driven time at 131 bounces 4.064407 s
   bounce n    time (event-driven)   time (series)    height after bounce (m)   e^(2n) h0
       1           0.451601          0.451601      6.400000e-01          6.400000e-01
       2           1.174162          1.174162      4.096000e-01          4.096000e-01
       3           1.752211          1.752211      2.621440e-01          2.621440e-01
       5           2.584601          2.584601      1.073742e-01          1.073742e-01
      10           3.579504          3.579504      1.152922e-02          1.152922e-02
  [PASS] event-driven bounce times equal the series
  [PASS] heights equal e^(2n) h0
  [PASS] simulation accumulates at the series limit (infinitely many bounces, finite time)  (131 bounces before v < 1e-12 m/s)
   total path length also finite: h0 (1+e^2)/(1-e^2) = 4.5556 m
   e = 1: gap between bounces constant at 0.903202 s (2 t0 = 0.903202), 1001 bounces take 903.65 s
  [PASS] e = 1: bounce time grows without limit

2. Same ball with a real contact: spring-damper, tension-free, gravity on during contact
     m (kg)  k (N/m)    contact time (ms)   e_target   bounces   rest at t (s)   Zeno T (s)   sag mg/k (mm), final depth (mm)
     0.050     200000      1.571          0.80        26        4.0659       4.0644     0.002, 0.002
     0.050      20000      4.967          0.80        21        4.0666       4.0644     0.025, 0.025
     0.050     200000      1.571          0.90        56        8.5875       8.5804     0.002, 0.002
     0.050     200000      1.571          0.60        11        1.8055       1.8064     0.002, 0.002
     0.005    2000000      0.157          0.95       161       17.6156      17.6124     0.000, 0.000
   crossover when flight time = contact time: n* = ln(2 v0 / (g tau)) / ln(1/e) = 28.5, simulated bounces 26
  [PASS] finite number of bounces with a real contact
  [PASS] bounce count is close to the flight-time = contact-time estimate
  [PASS] it comes to rest near the Zeno time
   per-bounce restitution: first 0.7999, middle 0.7974, last three [0.751 0.718 0.63 ]
  [PASS] effective e constant for most bounces, falls in the last ones (gravity during contact)

3. Other limits on the series (e = 0.8, h0 = 1 m, m = 50 g)
   height below one atomic size (1e-10 m) after n = 52 bounces
   bounce energy below k_B T at 300 K after n = 104 bounces; that hop height would be 8.4e-21 m, far below an atom, so the point is never reached
   50 g, 2.5 cm ball in still air, no gravity: Stokes diffusion D = 4.9e-16 m^2/s, time to wander 1 m = 3.2e+07 years

4. No gravity: a ball between two walls L = 4 m apart, v0 = 5 m/s
   e = 1: 1250 collisions in 1000 s; rate 1.250/s, v0/L = 1.250/s
  [PASS] e = 1: collision count grows linearly forever
   e = 0.8: n-th collision time (L/v0) (e^-n - 1)/(1/e - 1). Simulated vs series at n = 10, 50, 100:
      n =  10: 2.6602e+01 s, 2.6602e+01 s
      n =  50: 2.2420e+05 s, 2.2420e+05 s
      n = 100: 1.5709e+10 s, 1.5709e+10 s
  [PASS] e < 1 without gravity: gaps grow geometrically, collisions go on with no finite-time limit
   speed against time: v t (1-e)/(e L) at n = 100 is 1.0000 (derived value 1 as n grows)
  [PASS] speed falls as 1/t
   adding air: linear drag b, and the distance it allows
      linear drag, m/b = 50 s, e = 1: total reach v0 m/b = 250 m, so 62 collisions, derived floor(250/4) = 62
  [PASS] linear drag: finite collisions even with e = 1

5. Room with gravity, 4 x 3 x 2.5 m, start (1, 1, 1) m with v = (3, 2, 6) m/s, frictionless, axes uncoupled
   e = 1.0: by t = 100 s, side walls 142, ceiling 100, floor 100
  [PASS] e = 1: collisions keep accumulating at a steady rate  (count ratio over double time 2.00)
      vertical energy per kg conserved by construction 27.81 J/kg; ceiling reached since 2.84 m > 2.5 m
   e = 0.8: by t = 100 s, side walls 25, ceiling 1, floor 133
      vertical motion stops at T = 6.8377 s after 133 floor hits (speed cut-off 1e-12 m/s)
      side-wall hits before that time 6; floor hits are infinitely many in the ideal limit, side-wall hits are not
      with a frictionless floor the ball then slides: side-wall hit no. 10, 30, 50 at t = 44, 4.3e+03, 3.74e+05 s
  [PASS] vertical collapse happens in finite time
  [PASS] side-wall hits stay finite before the vertical collapse

All checks passed.
```

**Limits:** Constant $e$ is an idealisation: real $e$ depends on speed. The room run assumes frictionless walls, no spin and no corner hits, and a real corner impact is not defined by a single-axis rule. The spring-damper contact is linear, whereas a ball on a floor is closer to Hertzian. The Falcon comparison is qualitative, since their bead and surface were not modelled.

## Sources

Confirmed to exist by a search record or by opening the document. "Abstract" means only the abstract or a catalogue record was read, not the full text.

- Blajer and Czaplicki on trampoline contact: figure of 0.4 s seen only as cited in D. J. Burke, *The mechanics of the contact phase in trampolining* ([record](https://www.semanticscholar.org/paper/The-mechanics-of-the-contact-phase-in-trampolining-Burke/66f26c6e4cbc7ba5a5a8033ebb1cbba44c329b83)), through a search excerpt.
- Carré, Baker, Newell and Haake, "The dynamic behaviour of cricket balls during impact and variations due to grass and soil type", *Sports Engineering* 2 (1999) 145 ([link](https://link.springer.com/article/10.1007/BF02844050)). Abstract, background only. Carré et al. (2004) is cited only through Howard.
- R. Cross, "The bounce of a ball", *Am. J. Phys.* 67 (1999) 222 ([ADS](https://ui.adsabs.harvard.edu/abs/1999AmJPh..67..222C/abstract)). Abstract. Horizontal restitution figures from his [bounce page](https://www.physics.usyd.edu.au/~cross/BOUNCE.htm).
- E. Falcon, C. Laroche, S. Fauve and C. Coste, "Behavior of one inelastic ball bouncing repeatedly off the ground", *Eur. Phys. J. B* 3 (1998) 45 ([link](https://link.springer.com/article/10.1007/s100510050283)). Abstract.
- A. N. Gent and K. C. Rusch, "Viscoelastic behavior of open cell foams", *Rubber Chem. Technol.* 39 (1966) 389 ([record](https://scilit.com/publications/71a8948802223bb644791265fe218390)). Abstract.
- S. D. Grimshaw, A. Briggs and N. R. Atkins, "A review and reassessment of the aerodynamics of cricket ball swing", *Flow* 4 (2024) E6, [doi:10.1017/flo.2024.4](https://doi.org/10.1017/flo.2024.4). Read for the swing figure.
- J. A. Haringx, "On highly compressible helical springs and rubber rods, and their application for vibration-free mountings", *Philips Res. Rep.* (1948). Seen only through citing papers, and its theory is mentioned here only in outline.
- R. Howard, *Cricket Ball Bounce*, fourth-year undergraduate project, Cambridge University Engineering Department, 2011/2012 ([PDF](https://www-geo.eng.cam.ac.uk/sites/default/files/howard.pdf)). Read. It is a student report, which is the most detailed source I could open. It is also my route to James et al. and Carré et al. (2004).
- K. H. Hunt and F. R. E. Crossley, "Coefficient of restitution interpreted as damping in vibroimpact", *J. Appl. Mech.* 42 (1975) 440, [doi:10.1115/1.3423596](https://doi.org/10.1115/1.3423596). Abstract.
- ITF, ball rebound test ([technical booklet](https://www.itftennis.com/media/4420/2023-technical-booklet.pdf)), seen through a search excerpt. The values of $e$ are my arithmetic.
- D. M. James, M. J. Carré and S. J. Haake, "Predicting the playing character of cricket pitches", *Sports Engineering* (2005) ([record](https://shura.shu.ac.uk/650/)). Seen through an abstract excerpt and through Howard.
- MCC Laws of Cricket: [the ball](https://www.lords.org/mcc/the-laws/the-ball), [the creases](https://www.lords.org/mcc/the-laws/the-creases). Seen through search excerpts.
- MakeItFrom, [ASTM A228 music wire](https://makeitfrom.com/material-properties/ASTM-A228-SWP-A-K08500-Music-Wire): $E$ = 190 GPa, $G$ = 72 GPa, yield 2050 MPa. A materials database, not a primary source.
- M. Srinivasan, Y. Wang and A. Sheets, "People bouncing on trampolines", *PLoS ONE* 8 (2013) e78645 ([link](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0078645)). Read.
- D. Stoianovici and Y. Hurmuzlu, "A critical study of the applicability of rigid-body collision theory", *J. Appl. Mech.* 63 (1996) 307, [doi:10.1115/1.2788865](https://doi.org/10.1115/1.2788865). Abstract.
- R. W. Wilkins, *Low temperature bulk modulus of solid argon*, PhD thesis, University of Illinois, 1973 ([record](https://ideals.illinois.edu/items/30928)). Abstract: 28.8 ± 0.4 kbar.
- Argon pair-potential values as quoted: 125.7 K, 3.345 Å attributed to J. A. White (1999) in a [LAMMPS forum post](https://matsci.org/t/incorrect-well-depth-for-lj-cut-potential/33857); 124 K, 3.42 Å in a [University of Ulm exercise](https://www.uni-ulm.de/fileadmin/website_uni_ulm/nawi.inst.250/Lehre/SS_2010/ModSim/Uebung02_ModSim_SS2010.pdf). Weak sources, used only to show the spread.

From general knowledge and not verified here: the pitch length of 22 yards, 119.8 K and 3.405 Å for argon, dislocations as the cause of yield, the affine rubber result $G=nk_BT$, silly putty behaviour.
