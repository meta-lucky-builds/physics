# 1. Forces and rigid bodies

Ten questions about pushing, friction, balance and spin. Almost all of them reduce to three
ideas: force components (normal and tangential), torque about a pivot, and conservation of
momentum and angular momentum. Where a claim can be tested, a script in [`code/`](code/)
simulates the system rather than re-evaluating the formula. Measured values carry a source. Values
marked "assumed" or "general knowledge" were not confirmed against a source and can be changed in
the script.

Run any script from the repo root, for example `python 01-forces-and-rigid-bodies/code/verify_tipping.py`.
[`code/prem.py`](code/prem.py) is a helper holding the Earth density model used by F09 and F10.

### F01: Why glass doors are push from one side and pull from the other

**Status:** Derived + Code + Established

**Short answer:** It isn't a separate decision for each side. A door with a frame stop swings one
way only, so the same motion is automatically a push from one side and a pull from the other. The
real design decision is which way it swings, and for exits serving many people, building codes
generally require it to swing in the direction of escape. Handles and plates then only tell you
what the door already does.

**Starting from:** Rigid-body rotation about a fixed hinge, $\tau = rF\sin\alpha = I\ddot\theta$;
momentum and angular momentum balance for the impact on a stop. For the design part, the code texts
as reported by the sources linked below.

**Derivation:**

*Part 1: mechanics.*

1. **One force, two names.** Put the hinge at the origin with the closed door along $+x$. Let a stop
   on the $-y$ side block rotation that way. Only a force towards $+y$ gives torque in the allowed
   direction, $\tau = r F_y > 0$. Someone standing on the $-y$ side experiences a force towards $+y$
   as a push, because it points away from them. Someone standing on the $+y$ side experiences that
   same force as a pull. So one-way swing alone gives push on one side and pull on the other. A
   double-acting door with no stop, such as a spring-pivot door, is push from both sides.
2. **Handle position.** The force needed against a resisting torque $\tau_0$ is $F = \tau_0/r$. This
   is why the handle sits at the free edge. At $r = 0.8$ m the force is one eighth of that at
   $r = 0.1$ m.
3. **Door stops.** A door with moment of inertia $I = mw^2/3$ about the hinge hits a stop at distance
   $d$ from the hinge with angular speed $\omega$. For an elastic stop, the angular impulse is
   $Jd = 2I\omega$, so $J = 2I\omega/d$. The centre of mass, at $w/2$, changes momentum by
   $m(w/2)(2\omega)$, which must equal $J + J_h$, where $J_h$ is the hinge's impulse. Therefore

   $$J_h = J\left(\frac{3d}{2w} - 1\right).$$

   The hinge feels no impulse if $d = 2w/3$, the centre of percussion. A stop close to the hinge takes a
   much larger impulse and pulls hard on the hinge screws. A wedge near the hinge is also a lever: a
   steady force $F$ at $r$ puts $Fr/d$ on it.

*Part 2: design and codes.* These are rules and history rather than physics, so they are cited.

- **United States, IBC.** The International Building Code requires doors to swing in the direction of
  egress travel where serving an occupant load of 50 or more, or a Group H (high-hazard) occupancy,
  at section 1010.1.2 in the 2015 and 2024 editions and 1008.1.2 in the 2009 and 2012 editions
  ([Allegion's iDigHardware summary](https://idighardware.com/2026/03/decoded-door-swing/)). I read
  this summary, not the code text. Adopted editions and local amendments vary.
- **United States, NFPA 101.** The same summary reports that the Life Safety Code requires this for
  high-hazard contents of any occupant load, for 50 or more occupants, and for doors into exit
  enclosures, with limited exceptions.
- **India.** A [guide to NBC door requirements](https://www.studiomatrx.org/guides/nbc-door-requirements-india)
  and an [IISc fire-exit policy](https://olseh.iisc.ac.in/wp-content/uploads/2020/11/Policy-on-Fire-Exit-Emergency.pdf)
  point to outward-opening exit doors under the National Building Code 2016 Part 4. I could not check the official
  BIS text, so treat this one as unverified at clause level.
- **Why the rule exists.** At the Iroquois Theatre fire in Chicago in 1903, about 600 people died,
  and accounts describe locked or inward-opening exit doors that the crowd could not get open
  ([Wikipedia](https://en.wikipedia.org/wiki/Iroquois_Theatre_fire), which also records disagreement
  on details). The panic bar was developed in response
  ([Indianapolis encyclopedia](https://indyencyclopedia.org/panic-bar/)), and one
  [UL document](https://code-authorities.ul.com/wp-content/uploads/sites/2/2019/03/RI_BLST_AR_047_UL303PanicHardware_0319_en.pdf)
  gives a release force limit of 220 N.
- **The physics behind the rule.** Fruin reports bent steel railings showing crowd forces above 4500 N
  ([quoted on G. K. Still's site](https://gkstill.com/ExpertWitness/CrowdDisasters.html)). Spread over a
  0.9 m door, that is about 2025 N m about the hinge. Pulling such a door inward at the handle takes
  about 2500 N, 3.7 body weights, and the door also has to clear 0.64 m² of floor occupied by people.
  With an outward door, the same crowd force helps open it.
- **Below the threshold.** Small rooms usually open inward so that the door doesn't sweep into a
  corridor. This is general design practice, from general knowledge and unverified here.
- **Handles.** A flat plate signals pushing, and a bar signals pulling. Don Norman's
  [*The Design of Everyday Things*](https://en.wikipedia.org/wiki/The_Design_of_Everyday_Things) made
  this argument well known, and doors that get it wrong are now called
  ["Norman doors"](https://99percentinvisible.org/article/norman-doors-dont-know-whether-push-pull-blame-design/).
  Glass makes the handle matter more, because it hides the hinges and stops that would otherwise
  show the swing direction. That last point is my own reasoning.

**Check:** [`code/verify_door_mechanics.py`](code/verify_door_mechanics.py) treats the door as a free
rigid body held by stiff springs at the hinge and stop, so the hinge impulse is measured and not
assumed. An early version gave impulses 1.75 times too large. The cause was a bug in the script,
which used the hinge's moment of inertia for rotation about the centre of mass. The simulation exposed
it, and the corrected run below agrees with the formula. It also shows something the impulse formula
doesn't: with a linear stop, peak force is about the same for every $d$, $\omega\sqrt{kI} \approx 9.3$ kN,
and the extra impulse near the hinge arrives as a longer contact.

```text
person on -y side pushes (force towards +y): door turns 41.7 deg in 1 s
    person on +y side pulls (force towards +y): door turns 41.7 deg in 1 s
    person on +y side pushes (force towards -y): door turns 0.0 deg in 1 s
    person on -y side pulls (force towards -y): door turns 0.0 deg in 1 s
[PASS] same door: push from one side and pull from the other are one force direction

Force to overcome a 2 N m resisting torque, and time to reach 90 deg with 10 N extra:
    handle at 0.10 m: threshold  20.0 N, time 5.82 s
    handle at 0.45 m: threshold   4.4 N, time 2.75 s
    handle at 0.80 m: threshold   2.5 N, time 2.06 s
[PASS] same extra force: handle at 0.8 m opens the door faster than at 0.45 m and 0.1 m

Door swung into a stop at 2 rad/s (elastic stop, simulated free body on stiff hinge):
 d (m)   stop impulse (N s)  2 I om / d  hinge impulse (N s)  rigid prediction  peak stop force (N)
  0.10                432.2       432.0               -360.2            -360.2                 9261
[PASS] d=0.10: stop impulse matches 2 I omega / d within 2%
[PASS] d=0.10: hinge impulse matches J (3 d / 2 w - 1) within 3 N s
  0.30                144.0       144.0                -72.0             -72.0                 9282
[PASS] d=0.30: stop impulse matches 2 I omega / d within 2%
[PASS] d=0.30: hinge impulse matches J (3 d / 2 w - 1) within 3 N s
  0.45                 96.0        96.0                -24.0             -24.0                 9292
[PASS] d=0.45: stop impulse matches 2 I omega / d within 2%
[PASS] d=0.45: hinge impulse matches J (3 d / 2 w - 1) within 3 N s
  0.60                 72.0        72.0                 -0.0              -0.0                 9295
[PASS] d=0.60: stop impulse matches 2 I omega / d within 2%
[PASS] d=0.60: hinge impulse matches J (3 d / 2 w - 1) within 3 N s
  0.80                 54.0        54.0                 18.0              18.0                 9289
[PASS] d=0.80: stop impulse matches 2 I omega / d within 2%
[PASS] d=0.80: hinge impulse matches J (3 d / 2 w - 1) within 3 N s
  0.90                 48.0        48.0                 24.1              24.0                 9281
[PASS] d=0.90: stop impulse matches 2 I omega / d within 2%
[PASS] d=0.90: hinge impulse matches J (3 d / 2 w - 1) within 3 N s
[PASS] hinge impulse nearly vanishes for a stop at 2w/3 = 0.60 m: -0.00 N s vs -360.2 N s at d = 0.10 m
[PASS] a stop near the hinge takes a far larger impulse

Crowd force 4500 N spread across a 0.9 m door: moment 2025 N m about the hinge
    pull needed at a handle 0.8 m out to open it against the crowd: 2531 N = 3.7 x the weight of a 70 kg person
    floor area the door must sweep to open fully: 0.64 m^2
[PASS] an inward door against this crowd needs more than three body weights of pull

12 mm glass door 0.9 x 2.1 m: 57 kg; horizontal hinge reactions 139 N with hinges 1.8 m apart

All F01 checks passed.
```

**Limits:** The door is rigid and uniform, and stop stiffness (2 MN/m) is assumed. Door mass and glass
thickness are illustrative. Code statements come from secondary summaries and apply only where those
editions are adopted, so a real project needs the adopted code and the local authority.

### F02: Pressing down: does friction act upward?

**Status:** Derived + Code

**Short answer:** No. A force pressing straight down onto a floor produces no friction at all, and the
upward force you get is the normal force. Friction only ever acts along the surface. Pressing does
matter, though, because it raises the maximum friction available against sideways forces.

**Starting from:** Newton's laws and the Coulomb model of dry friction. The contact force is split into
a part along the surface normal ($N$) and a part along the surface (friction $f$), with
$|f| \le \mu N$ and $f$ opposing actual or impending sliding.

**Derivation:**

1. The rule "force left, friction right" is a special case. It holds because that force lies along the
   surface and tries to make the surfaces slide.
2. On a floor, a push $F$ straight down has no component along the surface. Balance along the surface
   reads $f = 0$. Balance along the normal reads $N = mg + F$. The floor pushes back upward, but this is
   the normal force.
3. Sideways pushes $P$ and a press $F$ together give $f = -P$ if $P \le \mu(mg + F)$. So pressing
   doesn't create friction but raises its limit.
4. Friction can act upward when the *surface is vertical*. Press a block against a wall and gravity
   tries to slide it along the wall, so friction acts upward with $f = mg$. A further downward force
   $F$ gives $f = mg + F$, up to $\mu N$. So your idea is correct for a wall, because there down is
   along the surface.
5. Combining these gives the friction cone. A force within angle $\arctan\mu$ of the normal can't cause sliding,
   however large it is, since both the push and the limit grow together. Outside the cone, even a small
   force slides it.
6. Microscopically, pressing a soft object does make it spread a little, so small local friction forces
   act inward in the contact patch. By symmetry they add to zero, which matches the net result above.

**Check:** [`code/verify_contact_friction.py`](code/verify_contact_friction.py) simulates a block with a
stiff normal contact and a tangential spring that slips at $\mu N$. The contact law isn't told which way
friction should point.

```text
Case 1: press straight down on a floor with 20 N
    N = 29.81 N up, friction = [0. 0.] N
[PASS] no friction at all; only the normal force goes up, N = m g + F
Case 2: press down 20 N and push right 5 N (limit mu N = 14.7 N)
    N = 29.81 N, friction = [-5.  0.] N, drift = 0.025 mm
[PASS] friction balances only the sideways part and points left
Case 3: pressing harder raises the limit: push right 8 N with 0 N press, then with 20 N press
    no press: speed 4.64 m/s; with press: speed 0.0000 m/s
[PASS] the press stops the slide by raising mu N
Case 4: slide-free block pressed against a wall with 30 N, gravity down
    N = 30.00 N, friction = [0.   9.81] N
[PASS] friction points up on a wall, equal to m g
Case 5: on that wall add a 5 N downward force
    N = 30.00 N, friction = [ 0.   14.81] N
[PASS] downward force on a vertical surface: friction does act upward
Case 6: friction cone. Weightless block, floor, mu = 0.5, cone half-angle 26.57 deg
       10 N at 25 deg from the normal: slip speed 0.0000 m/s
[PASS] inside the cone it never slides (10 N)
       10 N at 28 deg from the normal: slip speed 0.142 m/s
[PASS] just outside the cone it slides (10 N)
     1000 N at 25 deg from the normal: slip speed 0.0000 m/s
[PASS] inside the cone it never slides (1000 N)
     1000 N at 28 deg from the normal: slip speed 14.223 m/s
[PASS] just outside the cone it slides (1000 N)

All F02 checks passed.
```

**Limits:** The Coulomb model is an approximation (see F08). The test shows the model gives these
directions, not that the model is exact for real surfaces. Rolling contact and adhesion are left out.

### F03: Roads move and cars stay still: will that work?

**Status:** Derived + Code

**Short answer:** For the wheels and tyres, yes: a car held still on a road moving backward at 30 m/s
has the same wheel speed, tyre forces and rolling loss as a car driving at 30 m/s on a still road. But
a car that stays still relative to the ground goes nowhere, because your destination is fixed to the
ground. To travel, the car has to stay still relative to the *road* and be carried, and then it
meets the air just as it does now.

**Starting from:** Newton's laws, which take the same form in all inertial frames (Galilean invariance:
$x' = x - ut$, so accelerations and forces are unchanged). Tyre forces and rolling resistance depend on
relative motion between tyre and road.

**Derivation:**

1. The wheel-road contact only knows the slip speed $v - \omega R - v_\text{road}$ and the wheel's rotation
   rate. Add a constant $u$ to both $v$ and $v_\text{road}$ and nothing in the contact changes. So "car moves on
   still road" and "road moves under still car" are one physical situation seen from two frames.
   You were right about this.
2. The symmetry covers only what moves along with the road. Four things break it:
   - **Air.** Air is at rest relative to the ground, not the road. A car standing still on a treadmill
     has no drag, saving 9.7 kW at 30 m/s with the inputs in F04. A car carried by the road
     sees the full drag again.
   - **Destination.** Buildings don't move with the road, so standing still relative to the ground
     means no trip.
   - **Who supplies power.** Tyre loss is force times relative speed in either frame, 4.4 kW here. It is
     paid by the engine if the car holds itself in place, or by the belt motor if the car is tethered.
   - **Starting up.** While the belt accelerates, the frame isn't inertial. A free-wheeling car
     barely follows at first, because only the wheels' inertia couples to it:
     $a_\text{car}/a_\text{belt} = 1/(1 + MR^2/J)$, about 3%. With brakes locked it rides along up to the traction limit $\mu g$.
3. A surprise from the simulation: a free-wheeling car doesn't remain still. Rolling resistance acts as a
   steady force dragging it along with the belt until it balances air drag, at
   $v = \sqrt{C_{rr}Mg/(\tfrac12\rho C_dA)} \approx 20$ m/s, reached within minutes. To stay still it has to brake or drive.
4. Maximum traction is the same on a moving road, since it depends only on $\mu N$.
5. In a strict sense roads do move a little: a car reaching 30 m/s pushes the Earth back at about
   $7.5\times10^{-21}$ m/s, so momentum stays balanced.

**Check:** [`code/verify_moving_road.py`](code/verify_moving_road.py) integrates a car with wheel
dynamics and a tyre slip law, with road speed as an input. Identical torque programmes on a still road and on a
moving road agree to $10^{-7}$ m/s, which tests the invariance numerically. Output for this question:

```text
Same torque programme, road still vs road at -30 m/s, no air: final car speed relative to road 15.26 m/s; differences: speed 5.8e-08 m/s, wheel rate 2.0e-07 rad/s
[PASS] car motion relative to road is the same in both cases
Cruise at 30 m/s still road: v = 30.000; stationary on belt: v = 0.000 m/s, wheel 100.0 rad/s; power at wheels 4.41 kW
[PASS] the car stays still on the belt at wheel torque that holds cruise on a still road
With still air and the same torque: on the road the car slows to 24.8 m/s in 30 s; on the belt it stays at 0.000 m/s. Drag at 30 m/s = 324 N, 9.7 kW
[PASS] air breaks the symmetry
Free wheels, belt reaching 30 m/s at 2.0 m/s^2: car at 2.29 m/s, 7.6% of belt speed; wheel-inertia estimate without rolling resistance 2.9%, rolling resistance adds 1.47 m/s
[PASS] free-wheeling car barely follows the belt at first
Left free for 30 min with still air: car reaches 20.22 m/s (half way after 108 s); balance Crr M g = drag gives 20.22 m/s
[PASS] free car ends up carried at the rolling-resistance/drag balance
Locked wheels, belt 3.0 m/s^2 (gentle): car acceleration 3.00 m/s^2 (mu g = 8.83)
[PASS] car rides the belt below the traction limit
Locked wheels, belt 15.0 m/s^2 (hard): car acceleration 8.83 m/s^2 (mu g = 8.83)
[PASS] car slides at the traction limit above it
Wheelspin launch: acceleration 8.83 m/s^2 on still road, 8.83 m/s^2 on moving road
[PASS] maximum traction is the same on a moving road
```

**Limits:** Flat road, one-dimensional motion, tanh tyre law, constant $\mu$ and $C_{rr}$. Input values are listed
under F04. Tyres on a thin belt also behave differently from tyres on asphalt.

### F04: Can we design a road like a conveyor belt? (Lubricants, Teflon, central motors)

**Status:** Derived + Code + Estimate + Established

**Short answer:** Physically yes, and short versions exist: wind-tunnel rolling roads, airport-style walkways, the
1900 Paris moving walkway. As ordinary long roads they make poor sense, because you must move the whole
road even for one car. At normal traffic a steel-decked lane costs about 2.6 to 14 times the wheel power of
driving the cars, and much more at night, and it can't easily curve, merge or stop. Teflon doesn't help
because rollers already do better, and oil films lose more energy at road speed. A central motor is a good idea that
already works in electric railways, which move the vehicle and not the road.

**Starting from:** Force times speed for power, the F03 results, and the sized inputs below. All costs are
order-of-magnitude.

**Derivation:**

1. **Inputs.** Lane 10 km long, 3.7 m wide (assumed), at 30 m/s (assumed). Car $M = 1500$ kg (assumed),
   $C_dA = 0.6$ m² ([Wikipedia's definition](https://en.wikipedia.org/wiki/Automobile_drag_coefficient) with a
   typical $C_d \approx 0.3$ and about 2 m² frontal area), $C_{rr} = 0.01$ (midpoint of a
   [secondary source's](https://sportfactoryproshop.com/blog/tire-rolling-resistance) 0.007 to 0.014), air
   1.2 kg/m³. Roller drag factor $f = 0.02$, the DIN 22101 default for belt conveyors
   ([conveyor glossary](https://bison-convey-website.pages.dev/glossary/friction-factor)).
   Traffic of 2000 vehicles per hour per lane ([LA street design manual citing the HCM](https://completestreetdesignmanual.engineering.lacity.gov/e-300-roadway-design-controls-and-criteria/e-340-highway-street-capacity/e-341-uninterrupted-flow)).
   Steel 7850 kg/m³ and 250 MPa yield and concrete 2400 kg/m³ are general knowledge.
2. **Power of a belt.** The top and return strands together have mass $2m_LL$. Drag is
   $f(2m_LL + m_\text{cars})g$, and power is that times $v$, plus air drag for the carried cars.
   A car driven in the normal way needs $(C_{rr}Mg + \tfrac12\rho C_dAv^2)v = 14.1$ kW at the wheels.

   | Traffic | Cars on 10 km | Cars driven | 1 mm steel | 10 mm steel | 200 mm concrete |
   |---|---|---|---|---|---|
   | 2000 veh/h | 185 | 2.6 MW | 6.9 MW (2.6×) | 37.6 MW (14×) | 213 MW (81×) |
   | 18 000 veh/h (6 m spacing) | 1667 | 23.6 MW | 34.3 MW (1.5×) | 65.1 MW (2.8×) | 240 MW (10×) |
   | 200 veh/h (night) | 19 | 0.26 MW | 3.8 MW (14×) | 34.5 MW (132×) | 209 MW (800×) |

   The design principle here: move the smallest mass that can carry the load. A car weighs far less than
   the road section it would drag along.
3. **Credit where due.** Locked to the belt, cars can't brake unpredictably, so they could sit much closer.
   The 18 000 veh/h row assumes 6 m spacing, nine times normal capacity. Close spacing would also reduce
   drag, which the table leaves out. At saturation the gap narrows to under a factor of 2 for the
   lightest deck. Electrified belts also avoid engine losses, and the table compares only power at the wheels.
4. **Materials.** A 10 mm steel deck is 5800 t of moving steel per 10 km lane, before rollers
   and supports. The 1 mm deck is 581 t. Such a thin deck can only carry wheel loads with
   continuous support underneath, as wind-tunnel belts have.
5. **Teflon.** PTFE slides at $\mu \approx 0.05$ to 0.2 on steel
   ([Engineering ToolBox](https://engineeringtoolbox.com/friction-coefficients-d_778.html)). That is more than
   rollers' 0.02, giving 85 to 342 MW instead of 34 MW for the 10 mm deck, with 1.2 to 4.6 kW/m² of heat at the sliding
   interface. PTFE melts at 327 °C with a 260 °C working limit
   ([Piedmont Plastics sheet](https://www.piedmontplastics.com/resources/literatures/vital-polymers-ptfe-virgin)), and
   vendors warn of a low pressure-velocity limit for virgin PTFE, though I couldn't confirm a
   numeric limit. And it can't go on top of the road, where the car needs grip.
6. **Lubricants.** A fluid film has viscous drag $\eta v/h$, which grows with speed, and it doesn't care
   about the load. At 30 m/s: water at 1 mm gives 67 MW, water at 0.1 mm gives 666 MW, oil at 1 mm gives 6.7 GW.
   An air film at 0.1 mm gives 12 MW, the lowest here, though the pumping power and the wheels' concentrated
   loads aren't counted. Viscosities are general knowledge. Lubricants work well on small
   bearings, and sharing them across 74 000 m² of interface is where the numbers go wrong.
7. **Central motors.** One motor at the end has to pull the whole belt: 1140 kN of tension for the 10 mm deck, or 31 MPa,
   which steel handles. Steel pulling only itself reaches yield at $\sigma_y/(f\rho g) \approx 162$ km, so the tension
   limit is real for long roads and calls for distributed drives. The idea has precedent:
   San Francisco's cable cars have no motors on board; [one powerhouse pulls cables under the street](https://www.streetcar.org/wheels-motion/cable-cars-work/) at 9.5 mph, and cars grip the moving cable.
8. **Other problems that are not about energy.**
   - **Boarding.** Paris's 1900 [Rue de l'Avenir](https://en.wikipedia.org/wiki/Rue_de_l%27Avenir) used
     steps of 4.2 km/h. Reaching 108 km/h that way takes 26 parallel strips, so cars would need their
     own motors just to get on.
   - **Curves and merges.** A rigid deck doesn't turn, and a single speed applies to everyone.
   - **Thermal growth.** 10 km of steel changes length by 3.6 m over 30 K.
   - **Stopping.** At 0.3 g the 10 mm deck stops in 10 s, and 0.73 MWh must go somewhere, averaging 256 MW.
   - **Failure.** One jam stops every vehicle.
9. **What exists.** [Wind-tunnel rolling roads](https://mts.com/in/articles/automotive/case-study-windshear) run a 1 mm
   stainless belt, 3.2 m × 8.9 m, at 290 km/h. Long belts carry bulk freight, such as the
   [roughly 98 km Bou Craa system](https://science.nasa.gov/earth/earth-observatory/worlds-longest-conveyor-belt-system-150869/),
   though NASA gives 61 miles for it. So the speed and the length have each been achieved, but not together with a
   road's width and wheel loads.

**Check:** [`code/verify_moving_road.py`](code/verify_moving_road.py), sizing part:

```text
One car driven at 30 m/s: 14.1 kW at the wheels (drag part 9.7 kW)

Moving road L = 10 km, width 3.7 m, roller factor f = 0.02, top + return strand
            deck   kg/m  moving mass t   KE MWh  drag force kN  belt MW empty  stress MPa
      1 mm steel     29            581     0.07            114            3.4        30.8
     10 mm steel    290           5809     0.73           1140           34.2        30.8
 200 mm concrete   1776          35520     4.44           6969          209.1         nan

Comparison with the same traffic driven on an ordinary road:
                       traffic   cars  cars driven MW       1 mm steel      10 mm steel  200 mm concrete
    2000 veh/h (lane capacity)    185            2.62     6.9MW   2.6x    37.6MW  14.4x   212.5MW  81.2x
     18000 veh/h (6 m spacing)   1667           23.56    34.3MW   1.5x    65.1MW   2.8x   240.0MW  10.2x
             200 veh/h (night)     19            0.26     3.8MW  14.4x    34.5MW 131.9x   209.4MW 800.1x
[PASS] at normal lane capacity a 10 mm steel belt costs over 10 x the wheel power of driving: 14.4 x
    even the lightest deck at saturation costs 1.46 x driving (without any drafting benefit)

Support options under a 10 mm steel deck (contact pressure 0.77 kPa), belt force and power for both strands of the 10 km road, heat flux at the interface:
        rollers f = 0.02: shear     15.40 Pa, power     34.19 MW, heat flux    0.462 kW/m^2
          PTFE mu = 0.05: shear     38.50 Pa, power     85.48 MW, heat flux    1.155 kW/m^2
           PTFE mu = 0.2: shear    154.02 Pa, power    341.92 MW, heat flux    4.621 kW/m^2
         water film 1 mm: shear     30.00 Pa, power     66.60 MW, heat flux    0.900 kW/m^2
       water film 0.1 mm: shear    300.00 Pa, power    666.00 MW, heat flux    9.000 kW/m^2
           oil film 1 mm: shear   3000.00 Pa, power   6660.00 MW, heat flux   90.000 kW/m^2
         air film 0.1 mm: shear      5.40 Pa, power     11.99 MW, heat flux    0.162 kW/m^2
[PASS] a 1 mm oil film at 30 m/s is the worst option (viscous shear dominates)
    film pressure needed to lift the deck alone: 770 Pa; each car adds wheel loads of about 3.7 kN on a contact patch

Central drive on a 10 mm steel deck: tension 1140 kN, stress 30.8 MPa; yield-limited length for steel pulling itself on rollers = sigma_y / (f rho g) = 162 km
Thermal change in length: 12e-6 /K x 30 K x 10 km = 3.6 m
Emergency stop at 0.3 g from 30 m/s takes 10.2 s over 153 m; KE 0.73 MWh must leave as heat or regeneration, average 256 MW
Speed steps to reach 108 km/h at the 4.2 km/h Paris 1900 step: 26 parallel strips

Car reaching 30 m/s pushes the Earth back at 7.5e-21 m/s

All F03/F04 checks passed.
```

**Limits:** These are estimates. $f = 0.02$ was set for bulk conveyors, not road decks at 30 m/s, and a deck that carries
wheel loads could be worse. Costs in money aren't estimated, as I have no sourced figure. Deck designs are placeholders,
not structural designs, and fatigue, wear, weather and maintenance are left out.

### F05: Centre of mass of a donut

**Status:** Derived + Code

**Short answer:** It's exactly at the centre of the hole, with no donut there. The centre of mass is an average
position, not a piece of the object, so it doesn't have to be in the material. The hole is right where the COM point should be.

**Starting from:** The definition $\vec r_\text{cm} = \frac{1}{M}\int \vec r\,\rho\,dV$, and Newton's second and third laws.

**Derivation:**

1. **By symmetry.** A uniform torus is unchanged by a half-turn about its axis, and by reflection in its mid-plane.
   Those operations send $(x, y, z)$ to $(-x, -y, z)$ and to $(x, y, -z)$. The COM is defined by the shape alone, so it must
   stay put under each. The only point fixed by both is the origin, which is in the hole.
2. **By integration.** Use $x = (R + s\cos\phi)\cos\theta$, $y = (R + s\cos\phi)\sin\theta$, $z = s\sin\phi$ with
   volume element $s(R + s\cos\phi)\,ds\,d\phi\,d\theta$. The volume is $V = 2\pi^2Rr^2$. Then
   $\int x\,dV$ contains $\int_0^{2\pi}\cos\theta\,d\theta = 0$, $\int y\,dV$ contains $\int_0^{2\pi}\sin\theta\,d\theta = 0$, and
   $\int z\,dV$ contains $\int_0^{2\pi}\sin\phi(R + s\cos\phi)\,d\phi = 0$. So the COM is at the origin.
3. **It needn't lie in the body.** Since the COM is a weighted average of points, it always lies in the convex hull of the body,
   meaning the shape you'd get by shrink-wrapping it. For a convex body that means it's inside. A non-convex body can have it
   outside: a donut, a C shape (a half ring has it at $2R/\pi$ from the centre, off the wire), a boomerang.
4. **Why that empty point still means something.** Sum Newton's second law over all particles. Internal forces cancel in pairs,
   leaving $M\ddot{\vec r}_\text{cm} = \vec F_\text{ext}$. So a tossed donut's COM follows a parabola even though nothing is
   there, and the donut spins around that empty point. In uniform gravity the total torque about the COM is zero, so
   a donut hung on a nail settles with the hole's centre directly below the nail.
5. **What it doesn't mean.** The mass doesn't act as if it were all located there in every respect. The donut's gravity at its own
   centre is zero by symmetry, not infinite. Far away it looks like a point at the COM.

**Check:** [`code/verify_torus_com.py`](code/verify_torus_com.py) does the exact integral with sympy and a Monte Carlo integral over a
bounding box with no symmetry input, then tosses a lopsided spinning ring of particles and tracks the COM.

```text
exact volume V = 2*pi**2*R*r**2
    exact integrals of x, y, z dV = [0, 0, 0]
[PASS] volume is 2 pi^2 R r^2
[PASS] first moments vanish: COM at the geometric centre
    Monte Carlo: 979710 points, V = 177.817 cm^3 (exact 177.653), COM = [0.005 0.016 0.002] mm +- [0.03  0.03  0.008] mm
[PASS] Monte Carlo COM within 4 standard errors of the centre
[PASS] Monte Carlo volume within 0.2%
[PASS] the COM point itself is not in the material
[PASS] COM lies inside the convex hull of the donut
[PASS] half ring COM at 2R/pi, off the wire: 2.546 cm from the centre; the wire is 4.0 cm away at that angle
    iced donut COM = [0.01 0.01 1.28] mm
[PASS] icing shifts COM along the axis only, still not in material
[PASS] tossed lopsided ring: COM follows the parabola to under 1 micron: miss 7.1e-16 m
[PASS] no particle sits at the COM at the end of the toss: nearest particle 3.05 cm away
    ring on a nail, potential minimum at tilt 0.00 deg from straight down
[PASS] hanging ring rests with its empty centre directly below the nail
    field at centre / G M R^-2 = 8.4e-04; field at 3 m / (G M / 9) = 1.00014
[PASS] field at the donut's centre is near zero, far field looks like a point at the COM

All F05 checks passed.
```

**Limits:** Assumes a rigid or at least closed system for the COM motion, and uniform gravity for the torque statement. A non-uniform
donut, as with icing, moves the COM along the axis, and the script finds it still in the hole.

### F06: Pen in soil, deep foundations, heat, and a giant ball at the base

**Status:** Derived + Code + Estimate + Established

**Short answer:** Your pen idea is real physics, but tall towers mostly stand for a different reason: their huge weight
over a wide base beats the wind. Foundations mainly stop the ground from crushing or settling. Heat doesn't limit foundations at all: the
deepest are about 120 m, which is only about 3 °C warmer. A heavy ball only helps if its weight is low in the ball, because a uniform ball
with a tower on top falls over, and a ball able to hold a tower would need tens of tower masses and would sink into the ground.

**Starting from:** Moment balance about a pivot. Rankine passive earth pressure with Broms' resistance $3K_p\gamma zd$ for a pile, three times the Rankine value
([Pedone et al. 2026, describing Broms 1964](https://www.issmge.org/uploads/publications/1/141/600.pdf)). Rolling without slip, and
energy for stability. Burj Khalifa and other values are cited below. Soil values $\gamma = 16$ kN/m³ and $\phi = 30°$ are assumed.

**Derivation:**

1. **The pen.** A rigid pen of diameter $d$ pushed to depth $L$ turns about its tip when pushed, and soil resists. With
   resistance per unit depth $3K_p\gamma zd$ and lever arm $L - z$ about the tip:

   $$M = \int_0^L 3K_p\gamma zd\,(L - z)\,dz = \tfrac12K_p\gamma dL^3.$$

   So holding power grows as $L^3$. Doubling the depth gives eight times the capacity, which explains why a pen stands when it's deep and
   falls when it's shallow. Resistance grows with depth *and* the lever gets longer *and* more soil is engaged.
2. **Do towers work like that?** Partly. Burj Khalifa's piles run about 47 to 50 m below a 3.7 m raft ([Poulos & Bunce 2008](https://www.slideshare.net/slideshow/foundation-deisgn-of-burj-khalifa/16447928)),
   only about 6% of the 828 m height, compared to a third for the pen. The main stabiliser is weight. For a rigid
   block on a flat base of half-width $b$, tipping requires $M_\text{wind} > Wb$.
   Inputs: 500 000 t empty ([official fact sheet](https://www.burjkhalifa.ae/en/assets/FACT-SHEET.pdf), seen only as a search excerpt), $b = 35$ m (my estimate from the pile and raft
   layout), a 50 m/s wind, $C_d = 1.3$, and width tapering from 60 m to 10 m (all assumed). This gives
   $M_\text{wind} \approx 1.8\times10^{10}$ N m against $Wb \approx 1.7\times10^{11}$ N m, a factor of about 10. Poulos & Bunce quote about 5 for
   overturning, so my crude number has the right order.
3. **Scale.** Enlarge a tower by $s$. Wind force goes as area, so the wind moment goes as $s^3$. The weight moment goes as $s^4$. So toppling
   gets *easier* to resist as things grow. But pressure on the ground goes as $s$: weight $s^3$ over area $s^2$. The Burj's mean
   pressure on an estimated 4000 to 8000 m² is about 0.6 to 1.2 MPa. This is why it uses a piled raft socketed into weak rock
   (Poulos & Bunce): crushing and settlement govern here.
4. **Heat premise.** This needs correcting. Near-surface temperature rises about 25 to 30 °C per km ([Wikipedia](https://en.wikipedia.org/wiki/Geothermal_gradient)). The deepest
   building foundations are at about 120 m ([Guinness, Petronas Towers](https://www.guinnessworldrecords.com/world-records/102889-deepest-foundations-for-a-building)), so about 3 to 3.6 °C warmer. The [Kola borehole went 12 262 m](https://scienceblog.com/t-kola-superdeep-borehole-deepest-vertical-hole/), a hundred times deeper, and
   measured [about 180 °C at 10 km](https://www.science.ca/askascientist/viewquestion.php?qID=1404), still far below the melting point
   of rock, which is about 1000 °C from general knowledge (unverified). The practical limits are ground conditions, groundwater, construction and cost. You were right that the deep
   interior is hot, but foundations don't come close.
5. **Ball at the base, rigid analysis.** Consider a body with a spherical bottom of radius $R$, and COM at height $h_c$ above the lowest point, so
   $d = R - h_c$ below the sphere's centre. Rolling moves the centre sideways but keeps it at height $R$, so $U(\theta) = mg(R - d\cos\theta)$, and $U''(0) = mgd$. Therefore upright is stable if and only if $h_c < R$. This is the roly-poly toy.
6. **Uniform ball fails.** The ball's own COM is at its centre, so anything on top lifts the total COM above the centre, $h_c > R$, however heavy the ball. Such a tower falls. A bare
   uniform ball is exactly neutral. You need the mass low: with mass $m_b$ whose COM lies $\delta$ below the centre, stability needs $m_b\delta > M_t(R + H/2)$.
7. **Sizing.** Fill only the lower hemisphere ($m\delta = \rho\pi R^4/4$) under the same tower. Neutral stability needs $R = 80$ m of steel (17 tower masses) or 109 m of concrete (13). But a sphere tilts under *any*
   moment, $\theta = M/S$, whereas a flat base doesn't tilt at all until $M = Wb$. Limiting tilt to $H/500$ (a drift figure from general knowledge, unverified), or 1.7 m at the top, in the 50 m/s wind needs
   $R = 117$ m of steel (53 tower masses) or 159 m of concrete (40).
8. **The ground.** A sphere on a real surface has no point contact. With 2 to $2.7\times10^{11}$ N on ground bearing 10 MPa the contact radius is 80 to 92 m and it sinks 22 to 45 m. On 1 MPa ground it buries past its equator. Then it can't
   roll. Friction in the socket, $\sim\mu WR \approx 10^{13}$ N m with $\mu = 0.5$ (assumed), far exceeds the wind moment. So the ball ends up being an enormous embedded foundation, which is what you were trying to avoid, and far bigger than the piles. The answer to "will it stand" is: yes, if bottom-heavy enough, but it's a worse foundation.
9. **Where your idea does appear.** Low mass does stabilise real structures. Hywind Scotland's floating wind turbines have a [78 m draft](https://questfwe.com/wp-content/uploads/2018/07/Hywind-Scotland.pdf) with [dense ballast at the bottom](https://lkabminerals.com/wp-content/uploads/2025/03/Case-Study-Hywind-Floating-Offshore-Wind-_-Foundations-13-05EN-24-09.pdf). They float and may tilt, which a floating turbine can tolerate and an occupied tower can't. Taipei 101 hangs a [660 t mass](https://www.guinnessworldrecords.com/world-records/774616-largest-tuned-mass-damper) to damp sway, which is for motion, not for stopping toppling.

**Check:** [`code/verify_foundation.py`](code/verify_foundation.py) generates the rolling equation from the Lagrangian with sympy, integrates it, and bisects for the stability boundary. It finds $h_c$ within about 1% of $R$; the small excess is the finite test time, since growth is very slow near the boundary.

```text
Kp = 3.00. Pen d = 1 cm in soil, capacity about the toe:
    depth  0.5 cm:   0.030 mN m
    depth  1.0 cm:   0.240 mN m
    depth  2.0 cm:   1.920 mN m
    depth  4.0 cm:  15.360 mN m
[PASS] capacity scales as L^3 (doubling depth gives 8x)
[PASS] capacity equals (1/2) Kp gamma d L^3
    depth 0.5 cm: load   8.24 mN m vs capacity    0.03 mN m -> falls
    depth 2.0 cm: load   8.36 mN m vs capacity    1.92 mN m -> falls
    depth 5.0 cm: load   8.62 mN m vs capacity   30.00 mN m -> stands
    depth needed for that load: about 3.2 cm
[PASS] shallow pen falls, deep pen stands, in the model

Geometric scaling by s: wind moment ~ s^3, weight moment ~ s^4, bearing pressure ~ s.

Tower: W = 4.90e+09 N. 50 m/s wind: F = 5.65e+07 N, M = 1.78e+10 N m. Restoring W b = 1.72e+11 N m. Factor 9.6 (Poulos & Bunce quote about 5).
[PASS] gravity moment beats wind overturning by the same order as the published factor
    mean bearing pressure on 4000 m^2: 1.23 MPa
    mean bearing pressure on 8000 m^2: 0.61 MPa
    Burj Khalifa piles: 50 m -> 1.2 to 1.5 K warmer at 25-30 K/km
    Petronas foundations: 120 m -> 3.0 to 3.6 K warmer at 25-30 K/km
    Kola borehole: 12262 m -> 306.6 to 367.9 K warmer at 25-30 K/km
[PASS] deepest building foundation warms by under 4 K at 30 K/km
    Petronas depth is 1.0% of the Kola depth

Body on a sphere: stability boundary by bisection (simulated). h_c = R - d is the COM height.
    R =   1.0 m: boundary at d = -0.0071 R, that is h_c = 1.0071 R
[PASS] upright is stable only for h_c < R (R = 1.0)
    R =  50.0 m: boundary at d = -0.0071 R, that is h_c = 1.0071 R
[PASS] upright is stable only for h_c < R (R = 50.0)

Uniform ball of radius 100 m, ball mass 100 x tower, tower H = 828 m, on the ball:
    COM is 5.09 m above the ball's centre
[PASS] a uniform ball with a tower on top falls over in the simulation
[PASS] a bare uniform ball (d = 0) has zero angular acceleration at any tilt

Bottom-heavy hemisphere under a 500 000 t, 828 m tower, top half empty:
    steel: neutral at R = 80 m (17 x tower mass); tilt H/500 at R = 117 m, ball 2.67e+10 kg = 53 x tower
        ground bearing 1 MPa: contact radius 291 m > R, ball buries past its equator
[PASS] steel: ball cannot rest on a point at 1 MPa
        ground bearing 10 MPa: contact radius 92 m, sinks 45 m
[PASS] steel: ball cannot rest on a point at 10 MPa
    concrete: neutral at R = 109 m (13 x tower mass); tilt H/500 at R = 159 m, ball 2.01e+10 kg = 40 x tower
        ground bearing 1 MPa: contact radius 254 m > R, ball buries past its equator
[PASS] concrete: ball cannot rest on a point at 1 MPa
        ground bearing 10 MPa: contact radius 80 m, sinks 22 m
[PASS] concrete: ball cannot rest on a point at 10 MPa

Flat rigid base: no tilt until M = W b = 1.72e+11 N m. Sphere with restoring S: tilt = M / S for any M > 0.

All F06 checks passed.
```

**Limits:** Pen model ignores soil cohesion, so damp soil holds a shallow pen better than this. Tower geometry, wind and half-width are assumed, wind is static with no vortex shedding or resonance, and earthquakes are left out. I confirmed Broms' $3K_p\gamma zd$ from a paper describing it, not his original, and the $L^3$ moment is derived here. Taking the pivot at the tip is a simplification: Pedone et al. report a rotation point near 2/3 of the embedded length in their chalk models, which changes the capacity by a constant factor but keeps the $L^3$ law. Rocking is modelled as rigid rolling, and contact sinking is a simple area estimate rather than a soil analysis.

### F07: Pushing a 2 m and a 1 m person at the forehead

**Status:** Derived + Code

**Short answer:** With the same stance width and the same mass, the shorter person takes twice the steady force to tip: about 63 N against 32 N. But the reason is the height of the push, since that sets its lever arm about the feet, and not the height of the COM. COM height doesn't enter the steady force at all. It does matter for how far you can lean and how a quick shove plays out. If stance scales with height, the steady force is the same for both.

**Starting from:** Rigid-body torque about the edge of the base. The person is modelled as a rigid block, height $2h$ and depth $2b$, pivoting on the edge of the stance, with $b$ the horizontal distance from COM to that edge. The stance $b = 0.10$ m and forehead at 0.93 of height are assumed.

**Derivation:**

1. **Correcting the setup.** The COM is a point, so it has no height of 1 cm, and your 99.5 cm above and below are arithmetically right but not used. What counts is where the force acts and where the pivot is. A real standing COM is a little above half height, about 55% (general knowledge, unverified). Also, a 1 m person with the same proportions as a 60 kg, 2 m person would have mass $60/8 = 7.5$ kg, so a 60 kg, 1 m person must be much wider. That would raise $b$ and make them harder still to tip.
2. **Steady push.** A horizontal force $F$ at height $H$ turns the body about the rear edge with moment $FH$. Weight opposes it with moment $mgb$. Tipping begins when

   $$F_\text{tip} = \frac{mgb}{H}.$$

   The COM's height $h$ isn't in this formula, because gravity's line of action is vertical, so raising the COM doesn't change its horizontal lever arm.
3. **Your two people.** Same $b$: $F_\text{tip} \propto 1/H$, so 31.65 N and 63.29 N. Scaled stance, $b \propto$ height: 31.65 N each, since $b/H$ is the same. So your conclusion is right under the same-stance assumption, and depends on it.
4. **Sliding.** Feet slide at $\mu mg$, about 294 N with $\mu = 0.5$ (assumed), so tipping comes first when $b/H < \mu$, as here.
5. **Where COM height does matter.**
   - Tip-over angle: the body passes the point of no return at $\tan\theta = b/h$, 5.7° for tall and 11.3° for short.
   - Energy barrier: $mg(\sqrt{h^2 + b^2} - h) \approx mgb^2/(2h)$, or 2.9 J for tall and 5.8 J for short.
   - Quick shove: the angular impulse $JH$ must supply that energy, $(JH)^2/(2I) = E$. In the simulated 0.10 s shove, the same stance needs 13.4 N s against 20.1 N s. With scaled stance it reverses: 13.4 N s tall against 10.0 N s short. So for a shove, the answer depends on stance and on COM height.
   - Time to react: the fall time scale is $\sqrt{I/(mgL)}$, 0.37 s for tall against 0.26 s for short, so the taller person has longer to step.

**Check:** [`code/verify_tipping.py`](code/verify_tipping.py) simulates rocking about the edge and finds each threshold by bisection on whether the block falls.

```text
Sustained push held for 5 s; forehead at 0.93 of height; COM at 0.50 of height

case                                        simulated F (N)  m g b / H (N)
2 m person, same stance, b = 0.10 m                   31.65          31.65
[PASS] static threshold equals m g b / H (2 m person, same stance, b = 0.10 m)
1 m person, same stance, b = 0.10 m                   63.29          63.29
[PASS] static threshold equals m g b / H (1 m person, same stance, b = 0.10 m)
2 m person, stance scales with height, b = 0.05 s            31.65          31.65
[PASS] static threshold equals m g b / H (2 m person, stance scales with height, b = 0.05 s)
1 m person, stance scales with height, b = 0.05 s            31.65          31.65
[PASS] static threshold equals m g b / H (1 m person, stance scales with height, b = 0.05 s)
[PASS] same stance: the short person needs twice the force: ratio 2.000
[PASS] stance scaling with height: same force for both: ratio 1.000

COM height varied at fixed forehead height 0.93 m, b = 0.10 m:
    h = 0.35 m -> 63.29 N
    h = 0.50 m -> 63.29 N
    h = 0.55 m -> 63.29 N
    h = 0.70 m -> 63.29 N
[PASS] static threshold does not depend on COM height

Short shove lasting 0.10 s, threshold force and impulse:
    2 m person, b = 0.10 m: F =   133.6 N, impulse =   13.4 N s
    1 m person, b = 0.10 m: F =   200.9 N, impulse =   20.1 N s
    2 m person, b = 0.05 s: F =   133.6 N, impulse =   13.4 N s
    1 m person, b = 0.05 s: F =    99.6 N, impulse =   10.0 N s
    impulse ratio short/tall: same stance 1.50; scaled stance 0.75
[PASS] short shove, same stance: short person is harder
[PASS] short shove, scaled stance: tall person is harder

Geometry with b = 0.10 m:
    2 m person: tips past 5.7 deg, energy barrier 2.9 J, fall time scale 0.37 s
    1 m person: tips past 11.3 deg, energy barrier 5.8 J, fall time scale 0.26 s

Sliding needs mu m g = 294 N with mu = 0.5, far above the tipping force. Tipping comes first when b/H < mu.
[PASS] tipping precedes sliding for both people

All F07 checks passed.
```

**Limits:** People aren't rigid. A forehead push bends the neck and hips, ankles and hips produce active corrections ([Winter 1995, *Gait & Posture* 3, 193](https://edisciplinas.usp.br/pluginfile.php/5245117/mod_resource/content/1/winter_1996JNP.pdf) describes these strategies; the link is to Winter's 1996 follow-up), and people step. So these forces describe a passive, rigid person, the lower bound on what a real person resists. The block's inertia is a rough stand-in for a body, which changes the shove numbers but not the steady threshold.

### F08: Can the coefficient of friction exceed 1?

**Status:** Derived + Code + Established

**Short answer:** Yes. $\mu \le 1$ is not a law; it's just common for ordinary surfaces. Clean metals in vacuum exceed 1 and can weld solid, and rubber can exceed 1 too. Your intuition about molecular attraction is right, but polishing alone isn't what does it: removing the oxide and adsorbed films is. Perfect atomic flatness can even give nearly zero friction when the atoms of the two surfaces don't line up. So this one isn't simple, and it is well worth knowing.

**Starting from:** The macroscopic definition $\mu = F/N$. Bowden-Tabor adhesion, in which real contact occurs at plastically yielded junctions. Tabor's combined-stress yield condition. Periodic atomic potentials for the last part.

**Derivation:**

1. **Definition.** $\mu$ is the measured ratio of friction force to normal force. Nothing in the definition limits it. If $\mu > 1$, friction exceeds the normal force, which breaks no law. It also isn't a fundamental constant of a material pair. It is a summary of what happens at a contact.
2. **Why ordinary values are small.** Real contact area is only $A = N/H$, where $H$ is the hardness, since the tops of the bumps yield under the load. Sliding shears these junctions at stress $\tau$, so $F = \tau A$ and $\mu = \tau/H$. For a metal, $H \approx 3\sigma_y$ and $\tau \approx \sigma_y/\sqrt3$ give $\mu \approx 0.2$, a standard estimate. Surface contaminants weaken the interface further. This also explains why friction is proportional to load and doesn't depend on the apparent area, since load sets $A$.
3. **Why clean surfaces go far higher: junction growth.** With tangential force present, the junction yields according to $p^2 + \alpha s^2 = p_0^2$, where $p = N/A$ and $s = F/A$. So the area grows as tangential force grows. If the interface can carry a fraction $c$ of the bulk shear strength $s_0 = p_0/\sqrt\alpha$, sliding starts when $s = cs_0$, giving

   $$\mu = \frac{c}{\sqrt{\alpha(1 - c^2)}}.$$

   This diverges as $c \to 1$: a fully clean interface is as strong as the metal, so the junction keeps growing and there's no sliding, only seizure. [Tabor (1959)](https://www.sciencegate.app/document/10.1098/rspa.1959.0114) reports that an interface just 5% weaker than the bulk gives $\mu$ of order 1, and this formula reproduces that: 0.88 to 1.76 for $\alpha = 12$ to 3. Tiny traces of contamination therefore matter enormously.
4. **Adhesion breaks $F \propto N$.** If surfaces stick with no load, $F = \tau(A_0 + N/H)$, so $F/N$ diverges at small $N$. Then $\mu$ isn't a constant any more. So in this regime the coefficient itself stops being the useful idea, which is the proper form of your question.
5. **Your atomically flat surfaces.** van der Waals attraction between flat surfaces $D$ apart is $P = A_H/(6\pi D^3)$, about 200 MPa at 0.3 nm with $A_H = 10^{-19}$ J (general knowledge, unverified), which is 20 kN on 1 cm². So molecular attraction can indeed be huge. But three things stand in the way:
   - **Films.** In air, an oxide and adsorbed layer form quickly, so polishing in air doesn't give bare contact.
   - **Steps.** Cleaved surfaces can be atomically flat between steps, but steps stop full contact over large areas ([Buckley, NASA SP-277](https://ntrs.nasa.gov/api/citations/19720012801/downloads/19720012801.pdf)).
   - **Registry.** Atoms have to line up. Summing the sideways force over $N$ atoms on a periodic potential gives $N$ times the single-atom force if aligned. For incommensurate spacing the forces cancel, and the total stays bounded, so the force per atom falls about as $1/N$. Graphite shows this: rotated flakes have ultralow friction, called superlubricity ([Dienwiebel et al. 2004](https://link.aps.org/doi/10.1103/PhysRevLett.92.126101)). Soft enough lattices can distort and lock in anyway, called the Aubry transition (general knowledge, and shown in the simulation below).
6. **Measured values.**
   - 52100 steel on itself in vacuum goes from about 0.45 to nearly 5 as oxide wears away, then welds (Buckley, SP-277, fig. 23).
   - Diamond on itself goes from just over 0.1 in air to about 0.9 in vacuum (Bowden & Hanwell, in SP-277).
   - Pure metals sliding on themselves in high vacuum reached up to 6.5 ([NASA TN D-4824](https://ntrs.nasa.gov/api/citations/19690001426/downloads/19690001426.pdf)).
   - Buckley also reports complete seizure of clean elemental metals. Gold nanowires cold-weld within seconds on contact ([Lu et al. 2010](https://doi.org/10.1038/NNANO.2010.4)).
   - Ordinary tables list clean dry copper on copper at 1.6, silver on silver at 1.4, aluminium on aluminium at 1.05 to 1.35, and rubber on dry asphalt at 0.9 ([Engineering ToolBox](https://engineeringtoolbox.com/friction-coefficients-d_778.html)). Race tyres are reported at about 1.0 to 1.5 ([Building Speed](https://buildingspeed.org/2020/09/25/coefficient-of-friction-grip-nascar/)). Rubber's grip comes largely from a different mechanism, deformation and adhesion of a soft solid.

   The Galileo spacecraft antenna is often cited as a cold-welding failure, but the sources I could check disagree on the cause, so I don't rely on it.

**Check:** [`code/verify_friction_models.py`](code/verify_friction_models.py) solves junction growth numerically and compares with the formula and with Tabor's statement. It then tests aligned and misaligned contact, first as a rigid chain and then as an elastic chain relaxed under rising force. The force steps are coarse, a factor of 1.26: 0.005 means it slid at the lowest force tried, and 1.25 means it slid only at the top step.

```text
Junction growth, mu at gross sliding (solved numerically):
c = s_i/s0 alpha=3    alpha=9    alpha=12  
     0.200      0.118      0.068      0.059
     0.500      0.333      0.192      0.167
     0.800      0.770      0.444      0.385
     0.950      1.757      1.014      0.878
     0.990      4.052      2.339      2.026
     0.999     12.900      7.448      6.450
[PASS] numerical mu matches c / sqrt(alpha (1 - c^2)): worst fractional error 6.2e-13
[PASS] 5% weaker interface gives mu of order unity, as Tabor (1959) reports: [1.76 1.01 0.88]
[PASS] a nearly clean interface gives mu well above 1: 7.4

Adhesive contact: F = tau (A0 + N/H), tau = 200 MPa, A0 = 0.01 mm^2, H = 1 GPa
    N =    0.001 N: F =     2.000 N, F/N =  2000.200
    N =    0.100 N: F =     2.020 N, F/N =    20.200
    N =    1.000 N: F =     2.200 N, F/N =     2.200
    N =   10.000 N: F =     4.000 N, F/N =     0.400
    N =  100.000 N: F =    22.000 N, F/N =     0.220
    N = 1000.000 N: F =   202.000 N, F/N =     0.202
[PASS] F/N exceeds 1 at low load and settles to tau/H = 0.2 at high load

van der Waals pressure at D = 0.3 nm, A_H = 1e-19 J: 196 MPa; over 1 cm^2: 19.6 kN

Rigid chain, peak lateral force per atom (units of the single-atom force):
      N    aligned   golden-ratio spacing
     10     1.0000               0.057587
    100     1.0000               0.006213
   1000     1.0000               0.000114
  10000     1.0000               0.000094
[PASS] aligned contact: force per atom does not fall with size
[PASS] misaligned contact: force per atom falls roughly as 1/N: 0.0576 to 0.000094

Elastic chain, force per atom to slide / single-atom force, N = 55 on 34 periods (a near the golden ratio):
    spring K =  0.02: 0.993
    spring K =  0.05: 0.993
    spring K =  0.30: 0.250
    spring K =  1.00: 0.005
    spring K =  3.00: 0.005
    aligned (a = 1), K = 0.30: 1.250
[PASS] aligned chain needs the full single-atom force
[PASS] stiff misaligned chain (K = 3) slides almost freely
[PASS] soft misaligned chain pins: the Aubry transition

All F08 checks passed.
```

**Limits:** The value of $\alpha$ is empirical, and I couldn't read Tabor's full paper, so the formula is derived here from the yield condition rather than quoted. The junction model is a plasticity picture, not atomic. The chain models are one-dimensional, overdamped and at zero temperature, and they show mechanisms, not the $\mu$ of any material.

### F09: Can an asteroid speed up the Earth's spin?

**Status:** Derived + Code + Estimate + Established

**Short answer:** Yes, your cricket-ball picture is right. An off-centre hit changes the spin: a hit moving with the rotation speeds it up, and one against it slows it. But the effect is tiny: a Chicxulub-size impact changes the day by at most a few milliseconds, comparable to a century or so of tidal slowing. The actual Chicxulub trajectory had a westward part, so it probably slowed the Earth by about 1 to 4 ms.

**Starting from:** Conservation of angular momentum, $\vec L = I\vec\omega + m\vec r\times\vec v$. Earth values from the [NASA fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html): mass $5.9722\times10^{24}$ kg, $I/MR^2 = 0.3308$, sidereal period 23.9345 h. Impact sizes from the literature below.

**Derivation:**

1. **Cricket ball.** A tangential impulse $J$ at the edge of a ball of radius $r$ changes spin by $\Delta\omega = Jr/I$. With $I = \tfrac25mr^2$, $m = 0.16$ kg, $r = 3.6$ cm (general knowledge) and $J = 0.5$ N s, $\Delta\omega \approx 217$ rad/s. A hit through the centre, with no offset, changes no spin, unless the projectile sticks and adds inertia. The simulation shows all of these: from 2 rad/s, a sticking central hit gives 1.68 rad/s, a hit against the spin gives 1.26 rad/s, and a hit with the spin gives 2.11 rad/s.
2. **Earth.** Its spin angular momentum is $L_E = I\omega = 5.85\times10^{33}$ kg m²/s. An impactor at speed $v$ and angle $\gamma$ above the horizon, heading along the equator, carries $mvR\cos\gamma$ about the axis. Since $L_E$ is fixed, the day changes by

   $$\Delta T \approx -T\frac{mvR\cos\gamma}{L_E}.$$

   More generally only the part of $\vec r\times\vec v$ along the axis counts, and the rest tips the axis slightly. Added mass adds a small lengthening $\Delta I/I$ even for a vertical hit.
3. **Real numbers.** Density 2630 kg/m³ and sizes from [Schulte et al. 2010](https://doi.org/10.1126/science.1177265) (about 10 km) and the simulations in [Collins et al. 2020](https://doi.org/10.1038/s41467-020-15269-x) (17 km at 12 km/s; 20 km/s is the typical asteroid speed), plus the top of the range in [Durand-Manterola & Cordero-Tercero 2014](https://arxiv.org/abs/1403.6391). Best-case geometry, 45° at the equator:

   | Scenario | Mass | Energy | $\Delta T$ heading east | $\Delta T$ heading west |
   |---|---|---|---|---|
   | 10 km, 20 km/s | $1.4\times10^{15}$ kg | $2.8\times10^{23}$ J | −1.8 ms | +1.9 ms |
   | 17 km, 12 km/s | $6.8\times10^{15}$ kg | $4.9\times10^{23}$ J | −5.1 ms | +5.7 ms |
   | $4.6\times10^{17}$ kg, 15.9 km/s | $4.6\times10^{17}$ kg | $5.8\times10^{25}$ J | −465 ms | +506 ms |

   The third is the upper end of one published range, well above the usual size, so read it as a bound.
4. **Actual trajectory.** Collins et al. infer a steep impact, 45 to 60°, from the northeast. I take the heading as exact southwest and the modern site at 21.38° N. This trajectory moves westward, against the spin, so $\Delta T$ is about +0.9 to +3.8 ms, with the axis shifting about 4 to 10 milliarcseconds and Earth's own velocity changing by 5 to 14 µm/s.
5. **Scale.** Tides lengthen the day by about 1.8 ms per century ([Stephenson et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5247521/)), so Chicxulub is worth about 50 to 200 years of tidal slowing. NASA's estimate for the 2004 Sumatra earthquake is −2.68 µs ([release](https://www.eurekalert.org/pub_releases/2005-01/nsfc-nde011105.php)), and later model runs reported in the press give different values of the same order.
6. **Off-centre effect and energy.** The vertical hit gives no torque and only +0.06 ms from added inertia. The share of energy going to spin is $\omega\Delta L/E = 2\omega R\cos\gamma/v$. Earth's surface moves at about 465 m/s at the equator and the impactor at 20 km/s, so this is only 3.3% for the first scenario. The rest becomes heat, shock, ejecta and earthquakes. Rotation is the least of the problem. Impact speeds lie between about 11.2 km/s, the escape speed, and about 73 km/s, so these inputs can't change much.
7. **How big would it need to be?** Shortening the day by 1 hour at 20 km/s, 45°, takes $2.7\times10^{21}$ kg, a body about 1250 km across, with $5.4\times10^{29}$ J, more than the whole spin energy of the Earth ($2.1\times10^{29}$ J). It would melt the surface, so the rigid-Earth formula wouldn't be reliable there.
8. **Where it did matter.** A giant impact can set the spin. The present Earth-Moon system holds $3.46\times10^{34}$ kg m²/s, 5.9 times Earth's spin, and a Mars-mass body (6.4 × 10²³ kg, general knowledge) at 10 km/s in a grazing hit carries about that much. This is the standard picture of the [Moon's origin (Canup & Asphaug 2001)](https://doi.org/10.1038/35089010). So your idea really did act once in Earth's history, on a vastly larger scale.

**Check:** [`code/verify_earth_impact.py`](code/verify_earth_impact.py) builds a disc from particles and springs, spins it, and fires a projectile that sticks. New spin is measured by tracking the whole body's orientation and compared with angular momentum bookkeeping. Early versions disagreed by 3%, and the reason was physical: a projectile bonded to one particle swung freely, so the body wasn't rigid. With three bonds it agrees. The Earth part is a vector calculation with the inertia tensor.

```text
Spinning particle disc (10 kg, radius 1 m, 2 rad/s anticlockwise), hit by a 1 kg sticking projectile at 4 m/s
Positive offset b means the projectile passes above the centre (y > 0), moving +x, a clockwise push.
  b (m)  tracked spin rad/s  L_before/I_after   |L_f-L_0|/|L_0|
   0.00              1.6841            1.6841           5.3e-15
[PASS] spin rate from particle tracking equals L_before/I_after for b=0.0
[PASS] total angular momentum conserved for b=0.0
   0.70              1.2614            1.2614           2.4e-16
[PASS] spin rate from particle tracking equals L_before/I_after for b=0.7
[PASS] total angular momentum conserved for b=0.7
  -0.70              2.1125            2.1125           8.2e-15
[PASS] spin rate from particle tracking equals L_before/I_after for b=-0.7
[PASS] total angular momentum conserved for b=-0.7
[PASS] hit against the spin slows it
[PASS] hit along the spin speeds it up

Earth: I = 8.019e+37 kg m^2, omega = 7.2921e-05 rad/s, L = 5.847e+33 kg m^2/s, rotational KE = 2.13e+29 J
[PASS] PREM moment of inertia agrees with the NASA factor within 0.2%
    core holds 11.4% of the moment of inertia (PREM), mantle+crust 88.6%

Best case for spin-up: 45 deg impact on the equator, heading east. Then the worst case, heading west:
                                  scenario    m (kg)     E (J)  E (Mt TNT)  L/L_earth  dT east (ms)  dT west (ms)
        10 km, 20 km/s (Schulte 2010 size)  1.38e+15  2.75e+23    6.58e+07   2.12e-08        -1.768         1.888
        17 km, 12 km/s (Collins 2020 sim.)  6.77e+15  4.87e+23    1.16e+08   6.25e-08        -5.094         5.684
    D-M 2014 largest, 4.6e17 kg, 15.9 km/s  4.60e+17  5.81e+25    1.39e+10   5.63e-06      -465.454       505.584
[PASS] east-heading impact shortens the day, west-heading lengthens it
    estimate dT = -T L_imp / L_earth = -1.828 ms for the first scenario
[PASS] vector result matches the simple estimate within 10%
    head-on (vertical) impact: L_imp = 1.1e+10; dT = +0.0601 ms, from added mass alone
[PASS] vertical impact gives no spin torque; only a small lengthening from added inertia

Chicxulub site with the Collins et al. trajectory (heading SW, modern latitude):
    10 km, 20 km/s (Schulte 2010 size), 45 deg: dT = +1.256 ms, axis shift 4 mas, Earth dv = 4.6 micron/s
    10 km, 20 km/s (Schulte 2010 size), 60 deg: dT = +0.903 ms, axis shift 4 mas, Earth dv = 4.6 micron/s
    17 km, 12 km/s (Collins 2020 sim.), 45 deg: dT = +3.804 ms, axis shift 10 mas, Earth dv = 13.6 micron/s
    17 km, 12 km/s (Collins 2020 sim.), 60 deg: dT = +2.765 ms, axis shift 7 mas, Earth dv = 13.6 micron/s
[PASS] this trajectory has a westward part, so it slowed the spin

Energy for the first scenario: impact 2.75e+23 J; spin energy change omega L = 9.05e+21 J = 3.3% of it; ratio equals 2 (omega R) cos(angle) / v
[PASS] spin energy share is small and equals 2 omega R cos(angle)/v

For scale: tidal slowing +1.8 ms per century (Stephenson et al. 2016);
2004 Sumatra earthquake about -2.68 microseconds (NASA).
To shorten the day by 1 hour at 20 km/s, 45 deg: m = 2.70e+21 kg, diameter 1252 km, energy 5.4e+29 J = 1963721 x the first scenario
[PASS] that impactor is over 1000 km across and carries more energy than the Earth's spin
Moon's orbital L = 2.88e+34, Earth's spin L = 5.85e+33, total 3.46e+34 kg m^2/s; total / spin = 5.9
Mars-mass impactor (6.4e23 kg, general knowledge) at 10 km/s, grazing: L = 4.1e+34 kg m^2/s, 1.2 x the present system
[PASS] a Mars-sized grazing impact carries Earth-Moon scale angular momentum

Impact speed bounds: escape speed 11.2 km/s; solar-escape 42.1 + orbital 29.8 km/s head-on combined with Earth's well: 73 km/s

All F09 checks passed.
```

**Limits:** Rigid Earth and a sticking impactor, so escaping ejecta and its angular momentum are left out. Redistributed mass also matters: moving $10^{16}$ kg (assumed) from 21° latitude to a uniform blanket changes $I$ by a fraction of about $10^{-9}$, or about −87 µs. The core holds 11.4% of $I$, so if it lags at first, the mantle's change is up to about 13% larger until they couple. The heading and modern latitude are simplifications, and the size and speed of the impactor remain uncertain. The Schulte value is the hypothesis size given in that paper.

### F10: Drilling straight down: do we come out the other side?

**Status:** Derived + Code + Established

**Short answer:** Yes. A straight line that enters a ball must leave it, and through the centre it leaves at the antipode. From India that is in the eastern Pacific, so your hole would open on the seabed. The catch is that "keep going down" stops working at the centre, because down reverses there: past it you're drilling up, and you need a straight line rather than a plumb line.

**Starting from:** Euclidean geometry for a ball. Newton's gravity and the shell theorem, derived and integrated here. The PREM density model ([Dziewonski & Anderson 1981](https://resolver.caltech.edu/CaltechAUTHORS:20150218-110100613)), with polynomials [cross-checked in an independent code](https://geodynamics.org/pipermail/cig-commits/2009-April/009331.html) and against Earth's total mass. As you asked, heat, pressure and rock flow are ignored.

**Derivation:**

1. **Geometry.** For a line $\vec p + t\vec u$ from surface point $\vec p$ with unit direction $\vec u$, $|\vec p + t\vec u|^2 = R^2$ has roots $t = 0$ and $t = -2\vec p\cdot\vec u$. If the line goes inward, the second root is positive, so it exits. Through the centre, $\vec u = -\vec p/R$, $t = 2R$ and the exit is $-\vec p$, 12 742 km away. A real hole is short by comparison: the Kola borehole reached [12 262 m](https://scienceblog.com/t-kola-superdeep-borehole-deepest-vertical-hole/), or 0.096% of the diameter.
2. **Gravity inside.** Integrating the attraction over a thin shell gives zero at any point inside and $Gm/r^2$ outside (the shell theorem). So at radius $r$ only the mass inside $r$ pulls: $g(r) = Gm(r)/r^2$. For uniform density $g = g_0r/R$, which falls to zero at the centre. Real Earth has a dense core, so $g$ actually *rises* to 10.69 m/s² at the core boundary before falling to zero at the centre.
3. **Which way is down.** The pull points to the centre, so after the centre, down is back the way you came. A rule of "always drill down" ends at the centre, and getting through needs an inertial or surveyed straight line. A plumb line also doesn't aim exactly at the centre on an oblate Earth: it's off by up to 0.19° at 45° latitude, a 21 km miss at the centre.
4. **The rod.** A rod reaching the centre carries stress equal to its density times the gravity potential drop, $\rho\int g\,dr$. For steel that's 386 GPa, about 1500 times its yield strength, and steel holds its own weight only to about 3.2 km at surface $g$. So the "super rod" assumption is doing real work. Past the centre, the far part is pulled back towards the centre, so that section is in tension. The real Earth's central pressure is 364 GPa by this integration, which you asked to ignore.
5. **Falling through.** With no air, you fall to the far side and stop exactly at the opposite rim, then fall back. For uniform density it's simple harmonic motion with one-way time $\pi\sqrt{R/g_0}$, 42 min. With PREM it's shorter, 38 min, because $g$ stays nearly constant for much of the way down, as [Klotz (2015)](https://ar5iv.arxiv.org/html/1308.1342) showed. Peak speed is 9.9 km/s at the centre with PREM.
6. **Rotation.** A dropped object starts with the surface's eastward speed, but deeper rock moves more slowly, so it drifts east and hits the wall. Sideways shift is $\tfrac13\omega gt^3$. In a 5 m radius hole at the equator, that's after 27.6 s, at only 3.7 km depth. The hole's radius would have to be about 300 km to reach the centre without touching the wall. A pole-to-pole hole avoids it, and only that one. And an open hole at the seabed would let the ocean in, so yours needs sealing.

**Check:** [`code/verify_earth_tunnel.py`](code/verify_earth_tunnel.py) integrates the shell directly, integrates PREM and compares with measured values that weren't used to build it (mass, inertia factor), then simulates the drop. My fall times are 42 min 10 s, 37 min 58 s and 38 min 11 s against Klotz's 42:12, 38:00 and 38:11. The 2 s gaps in the first two come from the choice of $R$ and $g_0$.

```text
[PASS] line through interior exits on the far surface (20000 random lines): worst radius error fraction 1.0e-15
    antipode of central India (20.0N, 78.0E) is (20.0S, 102.0W), in the eastern Pacific
[PASS] shell theorem: zero force inside a shell: max |F| = 1.1e-15 (G m = 1, R = 1)
[PASS] shell theorem: outside it acts like a point mass: F r^2 / (G m) = [1.0, 1.0, 1.0]
[PASS] PREM mass matches NASA 5.9722e24 kg to 0.1%: M = 5.9732e+24 kg
[PASS] PREM I/(M R^2) matches NASA 0.3308 to 0.1%: factor = 0.3308
[PASS] PREM surface gravity near 9.8 m/s^2: g(R) = 9.822 m/s^2
    gravity peaks at r = 3480 km with g = 10.69 m/s^2
    hydrostatic central pressure = 364 GPa
    uniform density : 42 min 10.2 s, peak speed 7.91 km/s
    constant gravity: 37 min 58.0 s, peak speed 11.19 km/s
    PREM            : 38 min 10.7 s, peak speed 9.92 km/s
[PASS] uniform-density time is pi sqrt(R/g) and near Klotz's 42 min 12 s
[PASS] constant-gravity time near Klotz's 38 min 0 s
[PASS] PREM time near Klotz's 38 min 11 s
    5 m radius hole at the equator: wall reached after 27.6 s at depth 3.73 km
[PASS] wall-hit time matches the Coriolis estimate (1/3) omega g t^3 within 3%: estimate 27.6 s
    sideways offset by the time the depth reaches the Earth's centre: 299 km
[PASS] offset at the centre is hundreds of km
    gravitational potential drop surface to centre: 4.919e+07 J/kg
    steel rod stress at the end of a surface-to-centre rod: 386 GPa = 1545 x yield; self-supporting length at surface g: 3.2 km
[PASS] peak speed from energy equals the integrated peak speed: sqrt(2 dphi) = 9.92 km/s
    plumb line vs geocentric direction: up to 0.192 deg at latitude 45.1; straight-line miss at the centre 21.4 km
[PASS] plumb-line angle near 0.19 deg

All F10 checks passed.
```

**Limits:** Spherically symmetric Earth in the fall times, no air, equatorial hole for the Coriolis case. Klotz estimates rotation changes fall times by about 4 s. Steel yield and density are general knowledge. Real drilling is stopped long before this by temperature, pressure and rock flow, which the question sets aside. No material could do this; the answer is about geometry and gravity.
