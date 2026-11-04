"""F01: door mechanics.

Checks: (a) one-way swing makes the same door push from one side and pull from the other,
(b) opening force against handle position, (c) a stop's load and the hinge's load, found by
simulating a free rigid door held by stiff springs at the hinge and stop, (d) a crowd against
an inward-opening door.
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


w, m = 0.9, 40.0
I = m * w ** 2 / 3

# (a) The door lies along +x from a hinge at the origin. Its stop sits on the -y side, so it can
# only rotate towards +y (positive z rotation). Each person pulls or pushes a handle at r.
r = 0.8


def swing_angle(force_y, T=1.0):
    """Hinged door with a stop: rotation allowed only for theta >= 0."""
    def rhs(t, y):
        th, om = y
        tau = r * force_y * np.cos(th)
        if th <= 0 and tau <= 0:
            return [0.0, 0.0]
        return [om, tau / I]
    sol = solve_ivp(rhs, (0, T), [0.0, 0.0], max_step=0.01)
    return sol.y[0, -1]


F = 20.0
cases = {
    "person on -y side pushes (force towards +y)": +F,
    "person on +y side pulls (force towards +y)": +F,
    "person on +y side pushes (force towards -y)": -F,
    "person on -y side pulls (force towards -y)": -F,
}
for name, fy in cases.items():
    print(f"    {name}: door turns {np.degrees(swing_angle(fy)):.1f} deg in 1 s")
check("same door: push from one side and pull from the other are one force direction",
      swing_angle(F) > 0.1 and swing_angle(-F) == 0.0)

# (b) Handle position. A fixed 2 N m closer torque needs F = tau / r; time to 90 deg simulated.
tau_c = 2.0
print("\nForce to overcome a 2 N m resisting torque, and time to reach 90 deg with 10 N extra:")
times = {}
for rr in (0.1, 0.45, 0.8):
    f_need = tau_c / rr

    def rhs(t, y, rr=rr, f=f_need + 10.0):
        return [y[1], (rr * f - tau_c) / I]

    def done(t, y):
        return y[0] - np.pi / 2
    done.terminal = True
    s = solve_ivp(rhs, (0, 60), [0, 0], events=done, max_step=0.01)
    times[rr] = s.t_events[0][0]
    print(f"    handle at {rr:.2f} m: threshold {f_need:5.1f} N, time {times[rr]:.2f} s")
check("same extra force: handle at 0.8 m opens the door faster than at 0.45 m and 0.1 m",
      times[0.8] < times[0.45] < times[0.1])

# (c) Door swings shut into a stop at distance d from the hinge. The door is a free rigid body:
# state x, y, theta plus velocities. Stiff springs stand in for the hinge pin and the stop.
k_h, c_h = 2e8, 2e5
k_s = 2e6  # stop rubber, N/m
I_cm = m * w ** 2 / 12  # about the door's own centre of mass


def hit_stop(d, om0=2.0, th0=0.05):
    def state_rhs(t, s):
        x, y, th, vx, vy, om = s
        cs, sn = np.cos(th), np.sin(th)
        hx, hy = x - 0.5 * w * cs, y - 0.5 * w * sn           # hinge point on the door
        hvx, hvy = vx + 0.5 * w * sn * om, vy - 0.5 * w * cs * om
        Fh = np.array([-k_h * hx - c_h * hvx, -k_h * hy - c_h * hvy])
        py = y + (d - 0.5 * w) * sn                            # door point over the stop
        Fs = k_s * (-py) if py < 0 else 0.0
        arm_h = np.array([-0.5 * w * cs, -0.5 * w * sn])
        arm_s_x = (d - 0.5 * w) * cs
        torque = arm_h[0] * Fh[1] - arm_h[1] * Fh[0] + arm_s_x * Fs
        return [vx, vy, om, Fh[0] / m, (Fh[1] + Fs) / m, torque / I_cm], Fh, Fs

    def rhs(t, s):
        return state_rhs(t, s)[0]

    x0, y0 = 0.5 * w * np.cos(th0), 0.5 * w * np.sin(th0)
    v0x = 0.5 * w * np.sin(th0) * om0
    v0y = -0.5 * w * np.cos(th0) * om0
    ts = np.linspace(0, 0.25, 5001)
    sol = solve_ivp(rhs, (0, 0.25), [x0, y0, th0, v0x, v0y, -om0], method="LSODA", t_eval=ts,
                    rtol=1e-8, atol=1e-9)
    Fh_y = np.zeros_like(ts)
    Fs = np.zeros_like(ts)
    for i, s in enumerate(sol.y.T):
        _, fh, fs = state_rhs(ts[i], s)
        Fh_y[i], Fs[i] = fh[1], fs
    contact = Fs > 0
    J_stop = np.trapezoid(Fs, ts)
    # Hinge impulse during the contact, relative to the small steady value before it.
    J_hinge = np.trapezoid(Fh_y[contact], ts[contact]) if contact.any() else 0.0
    return J_stop, J_hinge, Fs.max(), sol.y[5, -1]


print("\nDoor swung into a stop at 2 rad/s (elastic stop, simulated free body on stiff hinge):")
print(f"{'d (m)':>6} {'stop impulse (N s)':>20} {'2 I om / d':>11} {'hinge impulse (N s)':>20} {'rigid prediction':>17} {'peak stop force (N)':>20}")
res = {}
for d in (0.10, 0.30, 0.45, 0.60, 0.80, 0.90):
    J, Jh, Fmax, om_end = hit_stop(d)
    res[d] = (J, Jh, Fmax)
    pred = 2 * I * 2.0 / d
    pred_h = J * (3 * d / (2 * w) - 1)
    print(f"{d:6.2f} {J:20.1f} {pred:11.1f} {Jh:20.1f} {pred_h:17.1f} {Fmax:20.0f}")
    check(f"d={d:.2f}: stop impulse matches 2 I omega / d within 2%", abs(J / pred - 1) < 0.02)
    check(f"d={d:.2f}: hinge impulse matches J (3 d / 2 w - 1) within 3 N s", abs(Jh - pred_h) < 3.0)
check("hinge impulse nearly vanishes for a stop at 2w/3 = 0.60 m", abs(res[0.60][1]) < 3.0,
      f"{res[0.60][1]:.2f} N s vs {res[0.10][1]:.1f} N s at d = 0.10 m")
check("a stop near the hinge takes a far larger impulse", res[0.10][0] / res[0.90][0] > 8)

# (d) A crowd against the door. Input: Fruin's crowd force of more than 4500 N.
F_crowd = 4500.0
tau_crowd = F_crowd * w / 2  # spread uniformly across the door width
F_handle = tau_crowd / 0.8
print(f"\nCrowd force {F_crowd:.0f} N spread across a {w} m door: moment {tau_crowd:.0f} N m about the hinge")
print(f"    pull needed at a handle 0.8 m out to open it against the crowd: {F_handle:.0f} N "
      f"= {F_handle / (70 * g):.1f} x the weight of a 70 kg person")
print(f"    floor area the door must sweep to open fully: {np.pi * w ** 2 / 4:.2f} m^2")
check("an inward door against this crowd needs more than three body weights of pull", F_handle > 3 * 70 * g)

# Glass door weight and hinge couple.
t_glass, h_door, rho_glass, hinge_gap = 0.012, 2.1, 2500.0, 1.8
m_glass = t_glass * h_door * w * rho_glass
print(f"\n12 mm glass door {w} x {h_door} m: {m_glass:.0f} kg; horizontal hinge reactions "
      f"{m_glass * g * w / 2 / hinge_gap:.0f} N with hinges {hinge_gap} m apart")

print()
if fails:
    print("FAILED:", fails)
    sys.exit(1)
print("All F01 checks passed.")
