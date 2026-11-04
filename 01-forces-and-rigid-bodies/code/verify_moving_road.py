"""F03, F04: a car on a moving road.

A one-dimensional car (body, driven wheels, tyre slip law, rolling resistance, optional air
drag) is integrated with the road speed as an input, so frame equivalence is tested rather than
assumed. The second half sizes a moving road from stated inputs.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
fails = []
g = 9.81


def check(name, ok, detail=""):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))
    if not ok:
        fails.append(name)


M = 1500.0      # kg, car
J = 4 * 1.0     # kg m^2, four wheels, assumed 1 kg m^2 each
Rw = 0.30       # m, wheel radius
mu = 0.9        # dry tyre on asphalt, Engineering ToolBox
Crr = 0.01      # rolling resistance coefficient
CdA = 0.60      # m^2, drag area
rho_air = 1.2   # kg/m^3
s0 = 0.05       # m/s, slip speed scale of the tanh tyre law
w0 = 0.5        # rad/s, smoothing of the rolling resistance sign
N = M * g


def car_rhs(t, y, road_speed, torque, air_on, wheels_locked=False, air_speed=0.0):
    """y = (v, omega). omega > 0 is forward rolling. All speeds are in the ground frame."""
    v, om = y
    vr = road_speed(t)
    if wheels_locked:
        om = 0.0
    slip = v - om * Rw - vr
    F_tyre = -mu * N * np.tanh(slip / s0)
    T_rr = Crr * N * Rw * np.tanh(om / w0)
    rel_air = v - air_speed
    F_air = -0.5 * rho_air * CdA * abs(rel_air) * rel_air if air_on else 0.0
    dv = (F_tyre + F_air) / M
    dom = 0.0 if wheels_locked else (torque(t) - F_tyre * Rw - T_rr) / J
    return [dv, dom]


def run(road_speed, torque, y0, T, air_on=False, locked=False, air_speed=0.0, n=601):
    ts = np.linspace(0, T, n)
    sol = solve_ivp(car_rhs, (0, T), y0, args=(road_speed, torque, air_on, locked, air_speed),
                    method="LSODA", t_eval=ts, rtol=1e-9, atol=1e-8)
    return ts, sol.y


# 1. Galilean invariance, no air. A torque programme on a still road and on a road at -U.
U = 30.0
prog = lambda t: 900.0 * np.exp(-t / 10.0)   # N m at the wheels
ts, A = run(lambda t: 0.0, prog, [0.0, 0.0], 40.0)
ts, B = run(lambda t: -U, prog, [-U, 0.0], 40.0)
err_v = np.max(np.abs((B[0] + U) - A[0]))
err_w = np.max(np.abs(B[1] - A[1]))
print(f"Same torque programme, road still vs road at -{U:.0f} m/s, no air: final car speed relative "
      f"to road {A[0, -1]:.2f} m/s; differences: speed {err_v:.1e} m/s, wheel rate {err_w:.1e} rad/s")
check("car motion relative to road is the same in both cases", err_v < 1e-3 and err_w < 1e-3)

# 2. Cruise on still road vs stationary car on a belt, no air.
T_hold = Crr * N * Rw
ts, C = run(lambda t: 0.0, lambda t: T_hold, [U, U / Rw], 20.0)
ts, D = run(lambda t: -U, lambda t: T_hold, [0.0, U / Rw], 20.0)
print(f"Cruise at {U:.0f} m/s still road: v = {C[0, -1]:.3f}; stationary on belt: v = {D[0, -1]:.3f} m/s, "
      f"wheel {D[1, -1]:.1f} rad/s; power at wheels {T_hold * D[1, -1] / 1e3:.2f} kW")
check("the car stays still on the belt at wheel torque that holds cruise on a still road",
      abs(D[0, -1]) < 1e-3 and abs(C[0, -1] - U) < 1e-3)

# 3. Now air at rest relative to the ground: invariance breaks.
ts, E = run(lambda t: 0.0, lambda t: T_hold, [U, U / Rw], 30.0, air_on=True)
ts, Fb = run(lambda t: -U, lambda t: T_hold, [0.0, U / Rw], 30.0, air_on=True)
F_drag = 0.5 * rho_air * CdA * U ** 2
print(f"With still air and the same torque: on the road the car slows to {E[0, -1]:.1f} m/s in 30 s; on the "
      f"belt it stays at {Fb[0, -1]:.3f} m/s. Drag at {U:.0f} m/s = {F_drag:.0f} N, "
      f"{F_drag * U / 1e3:.1f} kW")
check("air breaks the symmetry", E[0, -1] < U - 3 and abs(Fb[0, -1]) < 1e-3)

# 4. Belt starts under a free-wheeling car (neutral, no brakes).
ab = 2.0
belt = lambda t: min(ab * t, U)
ts, G1 = run(belt, lambda t: 0.0, [0.0, 0.0], U / ab, n=301)
frac = G1[0, -1] / U
pred = 1.0 / (1.0 + M * Rw ** 2 / J)
print(f"Free wheels, belt reaching {U:.0f} m/s at {ab} m/s^2: car at {G1[0, -1]:.2f} m/s, "
      f"{100 * frac:.1f}% of belt speed; wheel-inertia estimate without rolling resistance {100 * pred:.1f}%, "
      f"rolling resistance adds {Crr * g * U / ab:.2f} m/s")
check("free-wheeling car barely follows the belt at first",
      abs(G1[0, -1] - (pred * U + (Crr * g * (U / ab)) * (1 - pred))) < 0.1)

# 5. Long-run drift: rolling resistance pulls the car along until air drag balances it.
ts, G2 = run(belt, lambda t: 0.0, [0.0, 0.0], 1800.0, air_on=True, n=1801)
v_term = np.sqrt(Crr * M * g / (0.5 * rho_air * CdA))
t_half = ts[np.argmax(G2[0] > 0.5 * v_term)]
print(f"Left free for 30 min with still air: car reaches {G2[0, -1]:.2f} m/s (half way after {t_half:.0f} s); "
      f"balance Crr M g = drag gives {v_term:.2f} m/s")
check("free car ends up carried at the rolling-resistance/drag balance", abs(G2[0, -1] / v_term - 1) < 0.01)

# 6. Locked wheels on a belt: conveyor limit mu g.
for a_b, label in ((3.0, "gentle"), (15.0, "hard")):
    bb = lambda t, a_b=a_b: min(a_b * t, U)
    t_end = U / a_b
    ts, H = run(bb, lambda t: 0.0, [0.0, 0.0], t_end, locked=True, n=301)
    a_car = np.polyfit(ts[len(ts) // 2:], H[0][len(ts) // 2:], 1)[0]
    print(f"Locked wheels, belt {a_b} m/s^2 ({label}): car acceleration {a_car:.2f} m/s^2 (mu g = {mu * g:.2f})")
    if a_b < mu * g:
        check("car rides the belt below the traction limit", abs(a_car - a_b) < 0.05)
    else:
        check("car slides at the traction limit above it", abs(a_car - mu * g) < 0.05)

# 7. Traction is frame independent: wheel torque far above the grip gives mu g in both cases.
big = lambda t: 5000.0
ts, P1 = run(lambda t: 0.0, big, [0.0, 0.0], 1.0)
ts, P2 = run(lambda t: -U, big, [-U, 0.0], 1.0)
a1 = np.polyfit(ts[300:], P1[0][300:], 1)[0]
a2 = np.polyfit(ts[300:], P2[0][300:], 1)[0]
print(f"Wheelspin launch: acceleration {a1:.2f} m/s^2 on still road, {a2:.2f} m/s^2 on moving road")
check("maximum traction is the same on a moving road", abs(a1 - a2) < 1e-3 and abs(a1 - mu * g) < 0.05)

# 8. Sizing a 10 km one-lane moving road at 30 m/s. Input values are listed in the README.
L = 10e3
width = 3.7
v = U
f_roller = 0.02         # DIN 22101 default artificial friction factor
decks = {"1 mm steel": 7850 * 0.001 * width, "10 mm steel": 7850 * 0.010 * width,
         "200 mm concrete": 2400 * 0.20 * width}
P_car_drive = (Crr * M * g + 0.5 * rho_air * CdA * v ** 2) * v
P_air = 0.5 * rho_air * CdA * v ** 3
print(f"\nOne car driven at {v:.0f} m/s: {P_car_drive / 1e3:.1f} kW at the wheels (drag part {P_air / 1e3:.1f} kW)")
print(f"\nMoving road L = {L / 1e3:.0f} km, width {width} m, roller factor f = {f_roller}, top + return strand")
hdr = f"{'deck':>16} {'kg/m':>6} {'moving mass t':>14} {'KE MWh':>8} {'drag force kN':>14} {'belt MW empty':>14} {'stress MPa':>11}"
print(hdr)
rows = {}
for name, ml in decks.items():
    mass = 2 * ml * L
    ke = 0.5 * mass * v ** 2
    force = f_roller * mass * g
    area = ml / 7850 if "steel" in name else None
    stress = force / area / 1e6 if area else float("nan")
    rows[name] = (ml, mass, ke, force)
    print(f"{name:>16} {ml:6.0f} {mass / 1e3:14.0f} {ke / 3.6e9:8.2f} {force / 1e3:14.0f} {force * v / 1e6:14.1f} {stress:11.1f}")
print("\nComparison with the same traffic driven on an ordinary road:")
print(f"{'traffic':>30} {'cars':>6} {'cars driven MW':>15} " + " ".join(f"{n:>16}" for n in decks))
for flow, label in ((2000, "2000 veh/h (lane capacity)"), (18000, "18000 veh/h (6 m spacing)"), (200, "200 veh/h (night)")):
    ncars = flow / 3600 * L / v
    drive = ncars * P_car_drive
    cells = []
    for name, (ml, mass, ke, force) in rows.items():
        p = f_roller * (mass + ncars * M) * g * v + ncars * P_air
        cells.append(f"{p / 1e6:7.1f}MW {p / drive:5.1f}x")
    print(f"{label:>30} {ncars:6.0f} {drive / 1e6:15.2f} " + " ".join(f"{c:>16}" for c in cells))

# Check the belt beats driving only near saturation and for a very light deck.
ncars_sat = 18000 / 3600 * L / v
p_sat = f_roller * (rows["1 mm steel"][1] + ncars_sat * M) * g * v + ncars_sat * P_air
ncars_cap = 2000 / 3600 * L / v
p_cap = f_roller * (rows["10 mm steel"][1] + ncars_cap * M) * g * v + ncars_cap * P_air
check("at normal lane capacity a 10 mm steel belt costs over 10 x the wheel power of driving",
      p_cap / (ncars_cap * P_car_drive) > 10, f"{p_cap / (ncars_cap * P_car_drive):.1f} x")
print(f"    even the lightest deck at saturation costs {p_sat / (ncars_sat * P_car_drive):.2f} x driving "
      "(without any drafting benefit)")

# Lubrication options per unit area under a 10 mm steel deck: heat flux = shear stress x speed.
ml = rows["10 mm steel"][0]
pressure = ml * g / width
print(f"\nSupport options under a 10 mm steel deck (contact pressure {pressure / 1e3:.2f} kPa), belt force and "
      f"power for both strands of the {L / 1e3:.0f} km road, heat flux at the interface:")
A_int = 2 * width * L
opts = {
    "rollers f = 0.02": (f_roller * pressure, None),
    "PTFE mu = 0.05": (0.05 * pressure, None),
    "PTFE mu = 0.2": (0.20 * pressure, None),
    "water film 1 mm": (1e-3 * v / 1e-3, None),
    "water film 0.1 mm": (1e-3 * v / 1e-4, None),
    "oil film 1 mm": (0.1 * v / 1e-3, None),
    "air film 0.1 mm": (1.8e-5 * v / 1e-4, None),
}
for name, (tau, _) in opts.items():
    q = tau * v
    print(f"    {name:>20}: shear {tau:9.2f} Pa, power {tau * A_int * v / 1e6:9.2f} MW, heat flux {q / 1e3:8.3f} kW/m^2")
check("a 1 mm oil film at 30 m/s is the worst option (viscous shear dominates)",
      0.1 * v / 1e-3 > 0.2 * pressure)

# Air film feasibility: the deck's weight must be held by film pressure over the whole area.
print(f"    film pressure needed to lift the deck alone: {pressure:.0f} Pa; each car adds wheel loads of about "
      f"{M * g / 4 / 1e3:.1f} kN on a contact patch")

# Central motor: tension, steel area, thermal length change.
force = rows["10 mm steel"][3]
area = rows["10 mm steel"][0] / 7850
sigma_y = 250e6
print(f"\nCentral drive on a 10 mm steel deck: tension {force / 1e3:.0f} kN, stress {force / area / 1e6:.1f} MPa; "
      f"yield-limited length for steel pulling itself on rollers = sigma_y / (f rho g) "
      f"= {sigma_y / (f_roller * 7850 * g) / 1e3:.0f} km")
print(f"Thermal change in length: 12e-6 /K x 30 K x {L / 1e3:.0f} km = {12e-6 * 30 * L:.1f} m")
ke = rows["10 mm steel"][2]
print(f"Emergency stop at 0.3 g from {v:.0f} m/s takes {v / (0.3 * g):.1f} s over {v ** 2 / (2 * 0.3 * g):.0f} m; "
      f"KE {ke / 3.6e9:.2f} MWh must leave as heat or regeneration, average {ke / (v / (0.3 * g)) / 1e6:.0f} MW")
print(f"Speed steps to reach {v * 3.6:.0f} km/h at the 4.2 km/h Paris 1900 step: {v * 3.6 / 4.2:.0f} parallel strips")

# Earth recoil, for the 'road moves' intuition.
print(f"\nCar reaching {v:.0f} m/s pushes the Earth back at {M * v / 5.9722e24:.1e} m/s")

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All F03/F04 checks passed.")
