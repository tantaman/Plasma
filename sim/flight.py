"""
Free-flight simulation of the vehicle in MuJoCo (rigid-body dynamics + ground contact).

What MuJoCo integrates: the 6-DOF airframe built from the mass budget in
results/design.json, four thrust actuators (one per HV zone), gravity, and contact of
the landing feet with the floor.

What this file adds on top, every physics step:
  * thruster model  - stage voltage -> thrust, from the FEM corona solution + duct model
                      (results/design.json: thrust_map), with converter and duct-flow lags
  * electrical model - HV converter efficiency, LiPo pack (OCV curve, internal resistance,
                      state of charge), low-voltage behaviour
  * aerodynamics    - momentum ("ram") drag of the through-flow, bluff-body drag, damping
  * flight controller - cascaded position / attitude loops and a quadrotor-style mixer,
                      running on noisy state estimates at 100 Hz

Run:  python sim/flight.py            -> results/flight.{json,png,mp4}, results/flight_*.png
"""
from __future__ import annotations

import json
import os
import pathlib
import sys

os.environ.setdefault("MUJOCO_GL", "egl")

import mujoco
import numpy as np

import design as D

ROOT = pathlib.Path(__file__).resolve().parent.parent
RES = ROOT / "results"
G = 9.80665
RHO = 1.20


# --------------------------------------------------------------------------- model
def build_mjcf(rep: dict, payload: float) -> str:
    b = rep["design"]
    a = rep["assumptions"]
    cellg = rep["cell"]
    W, L = b["W"], b["L"]
    h = W / 2
    n = b["n"]
    S = cellg["d"] + cellg["c_c"] + cellg["c_d"]
    leg = 0.12
    m_pod = b["m_batt"] + b["m_conv"] + a["m_avionics"]
    # frame mass split over 12 box edges + rails (lumped into top/bottom edges)
    m_edge_h = (b["m_frame"] + a["m_misc"]) * 0.82 / 8
    m_edge_v = (b["m_frame"] + a["m_misc"]) * 0.18 / 4
    m_wall = b["m_skin"] / 4
    m_grille = b["m_strips"] / (n + 1)
    zt, zb = L / 2, -L / 2
    g = []
    # frame edges
    for z, tag in ((zt, "t"), (zb, "b")):
        for k, (x0, y0, x1, y1) in enumerate(((-h, -h, h, -h), (h, -h, h, h), (h, h, -h, h), (-h, h, -h, -h))):
            g.append(f'<geom name="e{tag}{k}" type="capsule" fromto="{x0} {y0} {z} {x1} {y1} {z}" size="0.004" '
                     f'mass="{m_edge_h:.6f}" rgba="0.12 0.12 0.12 1" contype="0" conaffinity="0"/>')
    for k, (x, y) in enumerate(((-h, -h), (h, -h), (h, h), (-h, h))):
        g.append(f'<geom name="ev{k}" type="capsule" fromto="{x} {y} {zb - leg} {x} {y} {zt}" size="0.004" '
                 f'mass="{m_edge_v:.6f}" rgba="0.12 0.12 0.12 1" contype="0" conaffinity="0"/>')
        g.append(f'<geom name="foot{k}" type="sphere" pos="{x} {y} {zb - leg}" size="0.008" mass="0.0003" '
                 f'rgba="0.9 0.5 0.1 1" friction="0.8 0.01 0.001"/>')
    # duct walls (film) and inlet lip
    for k, (x, y, sx, sy) in enumerate(((0, -h, h, 0.0005), (0, h, h, 0.0005), (-h, 0, 0.0005, h), (h, 0, 0.0005, h))):
        g.append(f'<geom name="wall{k}" type="box" pos="{x} {y} 0" size="{sx} {sy} {L/2}" mass="{m_wall:.6f}" '
                 f'rgba="0.75 0.85 0.95 0.22" contype="0" conaffinity="0"/>')
    for k, (x0, y0, x1, y1) in enumerate(((-h, -h, h, -h), (h, -h, h, h), (h, h, -h, h), (-h, h, -h, -h))):
        ox, oy = (0.0, -0.03) if k == 0 else (0.03, 0.0) if k == 1 else (0.0, 0.03) if k == 2 else (-0.03, 0.0)
        g.append(f'<geom name="lip{k}" type="capsule" fromto="{x0+ox} {y0+oy} {zt} {x1+ox} {y1+oy} {zt}" size="0.03" '
                 f'mass="1e-5" rgba="0.95 0.95 0.9 0.55" contype="0" conaffinity="0"/>')
    # electrode grilles: one mass-carrying slab each, plus visual strips
    z_first = zt - 0.06
    pol = ("0.85 0.25 0.2 0.9", "0.2 0.4 0.85 0.9")
    for i in range(n + 1):
        z = z_first - i * S
        g.append(f'<geom name="grille{i}" type="box" pos="0 0 {z - 0.5*(cellg["c_c"]+cellg["c_d"])}" '
                 f'size="{0.48*W} {0.48*W} {0.5*(cellg["c_c"]+cellg["c_d"])}" mass="{m_grille:.6f}" '
                 f'rgba="0 0 0 0" contype="0" conaffinity="0"/>')
        for j, y in enumerate(np.linspace(-0.465 * W, 0.465 * W, 16)):
            g.append(f'<geom type="box" pos="0 {y:.4f} {z - 0.5*(cellg["c_c"]+cellg["c_d"])}" '
                     f'size="{0.48*W} 0.0008 {0.5*(cellg["c_c"]+cellg["c_d"])}" mass="1e-7" rgba="{pol[i % 2]}" '
                     f'contype="0" conaffinity="0"/>')
    # central pod (battery + 4 HV modules + avionics) and payload
    g.append(f'<geom name="pod" type="box" pos="0 0 {zt + 0.015}" size="0.055 0.045 0.03" mass="{m_pod:.6f}" '
             f'rgba="0.15 0.15 0.18 1" contype="0" conaffinity="0"/>')
    g.append(f'<geom name="mast" type="capsule" fromto="0 0 {zb} 0 0 {zt + 0.04}" size="0.003" mass="1e-5" '
             f'rgba="0.12 0.12 0.12 1" contype="0" conaffinity="0"/>')
    if payload > 0:
        g.append(f'<geom name="payload" type="box" pos="0 0 {zb - 0.02}" size="0.02 0.02 0.015" mass="{payload:.6f}" '
                 f'rgba="0.1 0.7 0.3 1" contype="0" conaffinity="0"/>')
    sites, acts = [], []
    q = W / 4
    kz = 0.012  # N*m of yaw torque per N of zone thrust (canted exit tabs), sign alternates like rotor torques
    for k, (sx, sy, sgn) in enumerate(((1, 1, 1), (-1, 1, -1), (-1, -1, 1), (1, -1, -1))):
        sites.append(f'<site name="thr{k}" pos="{sx*q} {sy*q} 0" size="0.01" rgba="1 0.8 0 0.6"/>')
        acts.append(f'<motor name="zone{k}" site="thr{k}" gear="0 0 1 0 0 {sgn*kz}" ctrllimited="true" ctrlrange="0 5"/>')
    sites.append(f'<site name="inlet" pos="0 0 {zt}" size="0.005"/>')
    z0 = L / 2 + leg + 0.008
    return f"""
<mujoco model="ion_lantern">
  <option timestep="0.002" gravity="0 0 -{G}" integrator="implicitfast" density="0" viscosity="0"/>
  <visual><global offwidth="1280" offheight="720"/><quality shadowsize="4096"/>
    <headlight ambient="0.45 0.45 0.45" diffuse="0.5 0.5 0.5"/></visual>
  <asset>
    <texture name="grid" type="2d" builtin="checker" rgb1="0.86 0.87 0.88" rgb2="0.78 0.79 0.80" width="256" height="256"/>
    <material name="grid" texture="grid" texrepeat="12 12" reflectance="0.05"/>
    <texture type="skybox" builtin="gradient" rgb1="0.93 0.95 0.98" rgb2="0.75 0.80 0.88" width="64" height="256"/>
  </asset>
  <worldbody>
    <light pos="1.5 -2 4" dir="-0.3 0.4 -1" diffuse="0.6 0.6 0.6" castshadow="true"/>
    <geom name="floor" type="plane" size="6 6 0.1" material="grid" friction="0.8 0.01 0.001"/>
    <geom type="cylinder" pos="1.9 1.2 0.9" size="0.01 0.9" rgba="0.3 0.3 0.3 1" contype="0" conaffinity="0"/>
    <geom type="box" pos="1.9 1.2 1.0" size="0.06 0.002 0.002" rgba="0.9 0.1 0.1 1" contype="0" conaffinity="0"/>
    <body name="craft" pos="0 0 {z0}">
      <freejoint name="root"/>
      {' '.join(g)}
      {' '.join(sites)}
    </body>
    <camera name="chase" pos="2.6 -3.0 1.6" xyaxes="0.75 0.66 0 -0.2 0.23 0.95" fovy="42"/>
  </worldbody>
  <actuator>{' '.join(acts)}</actuator>
</mujoco>
"""


# --------------------------------------------------------------------------- subsystems
class Thruster:
    """Four HV zones.  Voltage -> thrust from the FEM/duct thrust map."""

    def __init__(self, rep, thrust_scale=1.0):
        tm = rep["thrust_map"]
        self.V = np.array(tm["V"])
        self.Tb = np.array(tm["T_body"]) / 4 * thrust_scale  # body force per zone (instant)
        self.T = np.array(tm["T"]) / 4 * thrust_scale        # steady ducted thrust per zone
        self.P = np.array(tm["P"]) / 4                       # electrical power per zone
        self.vd = np.array(tm["v_duct"])
        self.v_on = self.V[0]
        self.v_max = D.V_MAX
        self.f_max = float(np.interp(self.v_max, self.V, self.T))
        self.volts = np.full(4, 0.0)
        self.aug = np.zeros(4)  # duct-suction part of the thrust (lags the body force)
        self.tau_v = 0.015   # s, converter + multiplier response
        self.tau_f = 0.15    # s, duct through-flow build-up (~ duct length / flow speed)

    def volts_for(self, f):
        """Inverse map: wanted steady thrust per zone -> stage voltage command."""
        f = np.clip(f, 0.0, self.f_max)
        return np.where(f <= 1e-4, 0.0, np.interp(f, self.T, self.V))

    def step(self, v_cmd, dt, v_limit):
        v_cmd = np.clip(v_cmd, 0.0, min(self.v_max, v_limit))
        self.volts += (v_cmd - self.volts) * (1 - np.exp(-dt / self.tau_v))
        lit = self.volts > self.v_on
        body = np.where(lit, np.interp(self.volts, self.V, self.Tb), 0.0)
        steady = np.where(lit, np.interp(self.volts, self.V, self.T), 0.0)
        self.aug += ((steady - body) - self.aug) * (1 - np.exp(-dt / self.tau_f))
        p_el = np.where(lit, np.interp(self.volts, self.V, self.P), 0.0)
        v_duct = float(np.interp(self.volts.mean(), self.V, self.vd)) if lit.any() else 0.0
        return body + self.aug, p_el, v_duct


class Battery:
    """LiPo pack: open-circuit voltage vs state of charge, ohmic drop, coulomb counting."""

    SOC = np.array([0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00])
    OCV = np.array([3.27, 3.61, 3.69, 3.73, 3.77, 3.80, 3.84, 3.87, 3.94, 4.02, 4.11, 4.20])

    def __init__(self, spec):
        self.n = spec["cells"]
        self.cap = spec["capacity_ah"]
        self.R = spec["cells"] * spec["r_cell"] + spec["r_leads"]
        self.soc = 1.0
        self.v = self.ocv()
        self.i = 0.0

    def ocv(self):
        return self.n * float(np.interp(self.soc, self.SOC, self.OCV))

    def max_power(self):
        return self.ocv() ** 2 / (4 * self.R)

    def step(self, p_load, dt):
        e = self.ocv()
        p = min(p_load, 0.98 * self.max_power())
        self.i = (e - np.sqrt(e * e - 4 * self.R * p)) / (2 * self.R)
        self.v = e - self.i * self.R
        self.soc = max(0.0, self.soc - self.i * dt / 3600 / self.cap)
        return p


class Controller:
    """Cascaded hover controller + mixer (the same structure a quadrotor flight stack uses)."""

    def __init__(self, mass, inertia, arm, kz, f_max, mdot_hover):
        self.m, self.I, self.arm, self.kz, self.f_max = mass, inertia, arm, kz, f_max
        self.mdot = mdot_hover  # estimate of the through-flow, for ram-drag feed-forward
        self.iz = 0.0
        self.ixy = np.zeros(2)

    def update(self, t, dt, pos, vel, rpy, omega, ref_pos, ref_vel, ref_yaw, on_ground):
        # --- altitude -> collective thrust
        ez, evz = ref_pos[2] - pos[2], ref_vel[2] - vel[2]
        if not on_ground:
            self.iz = np.clip(self.iz + ez * dt, -1.0, 1.0)
        acc_z = 3.0 * ez + 3.0 * evz + 1.2 * self.iz
        tilt = np.cos(rpy[0]) * np.cos(rpy[1])
        T = self.m * (G + acc_z) / max(tilt, 0.8)
        # --- horizontal position -> tilt command (small angles), with integral against steady draughts
        exy, evxy = ref_pos[:2] - pos[:2], ref_vel[:2] - vel[:2]
        if not on_ground:
            self.ixy = np.clip(self.ixy + exy * dt, -1.0, 1.0)
        acc = 2.5 * exy + 1.0 * evxy + 0.6 * self.ixy + self.mdot / self.m * ref_vel[:2]
        cy, sy = np.cos(rpy[2]), np.sin(rpy[2])
        ax_b, ay_b = cy * acc[0] + sy * acc[1], -sy * acc[0] + cy * acc[1]
        lim = np.radians(7.0)
        pitch_ref = np.clip(ax_b / G, -lim, lim)
        roll_ref = np.clip(-ay_b / G, -lim, lim)
        # --- attitude -> torques
        tx = self.I[0] * (36.0 * (roll_ref - rpy[0]) - 9.0 * omega[0])
        ty = self.I[1] * (36.0 * (pitch_ref - rpy[1]) - 9.0 * omega[1])
        eyaw = (ref_yaw - rpy[2] + np.pi) % (2 * np.pi) - np.pi
        tz = self.I[2] * (2.5 * eyaw - 2.5 * omega[2])
        return self.mix(T, tx, ty, tz), dict(T=T, roll_ref=roll_ref, pitch_ref=pitch_ref)

    def mix(self, T, tx, ty, tz):
        """Zone order: (+x,+y), (-x,+y), (-x,-y), (+x,-y); yaw signs +,-,+,-."""
        a = self.arm
        sx = np.array([1, -1, -1, 1.0])
        sy = np.array([1, 1, -1, -1.0])
        sz = np.array([1, -1, 1, -1.0])
        T = np.clip(T, 0.0, 4 * self.f_max)
        base = T / 4
        d_rp = tx / (4 * a) * sy - ty / (4 * a) * sx
        d_y = tz / (4 * self.kz) * sz
        # authority: attitude first, then yaw, never beyond what the zones can deliver
        room = min(base, self.f_max - base)
        s = min(1.0, room / (np.abs(d_rp).max() + 1e-9)) if np.abs(d_rp).max() > room else 1.0
        d_rp = d_rp * s
        room_y = np.minimum(base + d_rp, self.f_max - base - d_rp).min()
        s = min(1.0, max(room_y, 0.0) / (np.abs(d_y).max() + 1e-9)) if np.abs(d_y).max() > room_y else 1.0
        return np.clip(base + d_rp + d_y * s, 0.0, self.f_max)


def smooth_ref(t, t0, t1, x0, x1):
    """Minimum-jerk transition from x0 (at t0) to x1 (at t1): value and rate."""
    if t <= t0:
        return x0, 0.0
    if t >= t1:
        return x1, 0.0
    s = (t - t0) / (t1 - t0)
    return x0 + (x1 - x0) * (10 * s**3 - 15 * s**4 + 6 * s**5), (x1 - x0) * (30 * s**2 - 60 * s**3 + 30 * s**4) / (t1 - t0)


# --------------------------------------------------------------------------- simulation
def simulate(rep, payload=None, thrust_scale=1.0, t_end=95.0, render=False, seed=1, gust=True, label="flight"):
    payload = rep["assumptions"]["m_payload"] if payload is None else payload
    model = mujoco.MjModel.from_xml_string(build_mjcf(rep, payload))
    data = mujoco.MjData(model)
    bid = model.body("craft").id
    mass = float(model.body_subtreemass[bid])
    Ri = np.zeros(9)
    mujoco.mju_quat2Mat(Ri, model.body_iquat[bid])
    Ri = Ri.reshape(3, 3)
    inertia = np.diag(Ri @ np.diag(model.body_inertia[bid]) @ Ri.T).copy()  # about body x, y, z
    b = rep["design"]
    W, L = b["W"], b["L"]
    A = W * W
    thr = Thruster(rep, thrust_scale)
    bat = Battery(rep["battery"])
    eff = rep["assumptions"]["conv_eff"]
    mdot_hover = RHO * A * float(np.interp(mass * G / 4, thr.T, thr.vd))
    ctl = Controller(mass, inertia, W / 4, 0.012, thr.f_max, mdot_hover)
    rng = np.random.default_rng(seed)
    dt = model.opt.timestep
    n_ctrl = int(round(0.01 / dt))
    z_ground = L / 2 + 0.12 + 0.008
    hover_z = z_ground + 1.0
    inlet = model.site("inlet").id

    # mission: spool-up, climb 1 m, hover, translate 0.6 m, yaw 45 deg, hover until the battery says land
    t_land = None
    wind = np.zeros(3)
    log = {k: [] for k in ("t", "x", "y", "z", "vz", "roll", "pitch", "yaw", "zref", "xref", "T", "W", "volts",
                           "p_el", "p_batt", "v_batt", "i_batt", "soc", "wind", "f")}
    frames = []
    renderer = mujoco.Renderer(model, 540, 960) if render else None
    f_cmd = np.zeros(4)
    landed_at = None
    max_tilt = 0.0
    step = 0
    airborne = False
    t_liftoff = None
    while data.time < t_end:
        t = data.time
        pos = data.qpos[:3].copy()
        quat = data.qpos[3:7].copy()
        R = np.zeros(9)
        mujoco.mju_quat2Mat(R, quat)
        R = R.reshape(3, 3)
        vel = data.qvel[:3].copy()           # world frame
        omega = data.qvel[3:6].copy()        # body frame
        roll = np.arctan2(R[2, 1], R[2, 2])
        pitch = -np.arcsin(np.clip(R[2, 0], -1, 1))
        yaw = np.arctan2(R[1, 0], R[0, 0])
        on_ground = pos[2] < z_ground + 0.004
        if not on_ground and not airborne and pos[2] > z_ground + 0.02:
            airborne, t_liftoff = True, t
        # ---- battery management: land when 22 % charge remains or the pack sags below 3.3 V/cell
        if t_land is None and airborne and (bat.soc < 0.22 or bat.v < 3.3 * bat.n):
            t_land = t
        if t_land is None and t > t_end - 14:
            t_land = t  # (only reached in the short robustness runs)
        # ---- reference trajectory
        if t_land is None:
            zr, vzr = smooth_ref(t, 2.0, 10.0, z_ground - 0.02, hover_z)
            xr, vxr = smooth_ref(t, 18.0, 28.0, 0.0, 0.6)
            yawr, _ = smooth_ref(t, 32.0, 40.0, 0.0, np.radians(45))
        else:
            if "z_at_land" not in log:
                log["z_at_land"] = zr
            zr, vzr = smooth_ref(t, t_land, t_land + 9.0, log["z_at_land"], z_ground - 0.05)
            vxr = 0.0
        # ---- controller at 100 Hz on noisy estimates
        if step % n_ctrl == 0:
            est_pos = pos + rng.normal(0, 0.004, 3)
            est_vel = vel + rng.normal(0, 0.01, 3)
            est_rpy = np.array([roll, pitch, yaw]) + rng.normal(0, np.radians(0.15), 3)
            est_om = omega + rng.normal(0, np.radians(0.5), 3)
            if t < 1.0 or (landed_at is not None):
                f_cmd = np.zeros(4)
            else:
                f_cmd, _ = ctl.update(t, 0.01, est_pos, est_vel, est_rpy, est_om, np.array([xr, 0.0, zr]),
                                      np.array([vxr, 0.0, vzr]), yawr, on_ground)
                if t < 2.0:  # spool-up: light the corona at low power before the climb
                    f_cmd = np.minimum(f_cmd, 0.3 * mass * G / 4)
        if t_land is not None and on_ground and t > t_land + 6.0 and landed_at is None:
            landed_at = t
        # ---- thruster + electrical
        # the converters cannot hold full voltage once the pack sags (input power limited)
        p_avail = 0.98 * bat.max_power() * eff
        v_cmd = thr.volts_for(f_cmd)
        f, p_el, v_duct = thr.step(v_cmd, dt, v_limit=D.V_MAX)
        if p_el.sum() > p_avail:
            f, p_el = f * p_avail / p_el.sum(), p_el * p_avail / p_el.sum()
        f = f * (1 + rng.normal(0, 0.01, 4))
        p_batt = bat.step(p_el.sum() / eff + 3.0, dt)
        data.ctrl[:] = np.clip(f, 0.0, 5.0)
        # ---- aerodynamics (applied at the inlet centre, in world frame)
        if gust:
            wind[:2] += (-wind[:2] / 4.0) * dt + 0.05 * np.sqrt(2 * dt / 4.0) * rng.normal(size=2)  # ~0.05 m/s rms draughts
            wind_now = wind + (np.array([0.0, 0.18, 0.0]) if 46.0 < t < 50.0 else 0.0)                # one deliberate gust
        else:
            wind_now = wind
        p_in = data.site_xpos[inlet].copy()
        v_pt = vel + R @ np.cross(omega, R.T @ (p_in - data.xipos[bid]))
        v_rel = v_pt - wind_now
        axis = R[:, 2]
        v_ax = float(v_rel @ axis)
        v_lat = v_rel - v_ax * axis
        mdot = RHO * A * v_duct
        F = -mdot * v_lat                                   # ram drag of the captured stream tube
        F += -mdot * v_ax * axis                            # thrust lapse with axial speed
        F += -0.5 * RHO * 1.1 * (W * L) * np.linalg.norm(v_lat) * v_lat      # bluff-body side drag
        F += -0.5 * RHO * 1.2 * A * 0.15 * abs(v_ax) * v_ax * axis           # drag of grilles when unpowered
        tau = -R @ (np.array([0.04, 0.04, 0.01]) * (0.3 + mdot) * omega)     # rotational damping
        data.qfrc_applied[:] = 0.0
        mujoco.mj_applyFT(model, data, F, tau, p_in, bid, data.qfrc_applied)
        mujoco.mj_step(model, data)
        step += 1
        max_tilt = max(max_tilt, np.degrees(np.arccos(np.clip(R[2, 2], -1, 1)))) if airborne else max_tilt
        if step % 10 == 0:
            for k, v in (("t", t), ("x", pos[0]), ("y", pos[1]), ("z", pos[2] - z_ground), ("vz", vel[2]),
                         ("roll", np.degrees(roll)), ("pitch", np.degrees(pitch)), ("yaw", np.degrees(yaw)),
                         ("zref", max(zr - z_ground, 0.0)), ("xref", xr), ("T", f.sum()), ("W", mass * G),
                         ("volts", thr.volts.copy()), ("p_el", p_el.sum()), ("p_batt", p_batt), ("v_batt", bat.v),
                         ("i_batt", bat.i), ("soc", bat.soc), ("wind", wind_now[:2].copy()), ("f", f.copy())):
                log[k].append(v)
        if render and step % 20 == 0:  # 25 fps
            cam = mujoco.MjvCamera()
            cam.type = mujoco.mjtCamera.mjCAMERA_FREE
            cam.lookat[:] = [0.3, 0.0, 0.7]
            cam.distance, cam.azimuth, cam.elevation = 3.4, 128 + 0.25 * t, -12
            renderer.update_scene(data, cam)
            frames.append(renderer.render().copy())
        if landed_at is not None and t > landed_at + 2.0:
            break
    log.pop("z_at_land", None)
    out = {k: np.array(v) for k, v in log.items()}
    air = out["z"] > 0.02
    hov = (out["t"] > 11) & (out["t"] < (t_land or t_end)) & air
    summary = dict(
        label=label, mass_g=mass * 1e3, payload_g=payload * 1e3, thrust_scale=thrust_scale,
        max_thrust_N=4 * thr.f_max, thrust_to_weight=4 * thr.f_max / (mass * G),
        inertia_kgm2=inertia.tolist(), lifted_off=bool(airborne), t_liftoff_s=t_liftoff,
        max_altitude_m=float(out["z"].max()), airborne_time_s=float(air.sum() * 10 * dt),
        hover_power_el_W=float(out["p_el"][hov].mean()) if hov.any() else None,
        hover_power_batt_W=float(out["p_batt"][hov].mean()) if hov.any() else None,
        hover_voltage_kV=float(out["volts"][hov].mean() / 1e3) if hov.any() else None,
        hover_current_A=float(out["i_batt"][hov].mean()) if hov.any() else None,
        hover_c_rate=float(out["i_batt"][hov].mean() / bat.cap) if hov.any() else None,
        pos_error_rms_m=float(np.sqrt(np.mean((out["x"][hov] - out["xref"][hov]) ** 2 + out["y"][hov] ** 2))) if hov.any() else None,
        alt_error_rms_m=float(np.sqrt(np.mean((out["z"][hov] - out["zref"][hov]) ** 2))) if hov.any() else None,
        max_tilt_deg=float(max_tilt), min_pack_voltage_V=float(out["v_batt"].min()),
        soc_at_landing=float(out["soc"][-1]), t_land_command_s=t_land, t_touchdown_s=landed_at,
        final_yaw_deg=float(out["yaw"][-1]),
    )
    return out, summary, frames


def plot(out, s, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = out["t"]
    fig, ax = plt.subplots(3, 2, figsize=(12.5, 9.5), sharex=True)
    a = ax[0, 0]
    a.plot(t, out["zref"], "k--", lw=1, label="altitude reference")
    a.plot(t, out["z"], color="C0", label="altitude")
    a.plot(t, out["xref"], color="0.5", ls="--", lw=1)
    a.plot(t, out["x"], color="C1", label="x")
    a.plot(t, out["y"], color="C2", label="y")
    a.set(ylabel="position [m]", title="Position (feet above floor)")
    a.legend(frameon=False, ncol=4, fontsize=8)
    a = ax[0, 1]
    a.plot(t, out["roll"], label="roll")
    a.plot(t, out["pitch"], label="pitch")
    a.plot(t, out["yaw"] / 10, label="yaw / 10")
    a.set(ylabel="angle [deg]", title="Attitude")
    a.legend(frameon=False, ncol=3, fontsize=8)
    a = ax[1, 0]
    a.plot(t, out["T"], color="C3", label="total thrust")
    a.plot(t, out["W"], "k:", label="weight")
    a.axhline(s["max_thrust_N"], color="C3", lw=0.8, ls="--")
    a.text(t[-1], s["max_thrust_N"], " max", va="center", fontsize=8, color="C3")
    a.set(ylabel="force [N]", title="Thrust vs weight")
    a.legend(frameon=False, fontsize=8)
    a = ax[1, 1]
    a.plot(t, out["volts"] / 1e3, lw=0.9)
    a.set(ylabel="stage voltage [kV]", title="HV zone voltages (4 zones)")
    a = ax[2, 0]
    a.plot(t, out["p_batt"], color="C0", label="battery power")
    a.plot(t, out["p_el"], color="C4", label="corona (electrical) power")
    a.set(ylabel="power [W]", xlabel="time [s]", title="Power")
    a.legend(frameon=False, fontsize=8)
    a = ax[2, 1]
    a.plot(t, out["v_batt"], color="C1", label="pack voltage [V]")
    a.plot(t, out["i_batt"], color="C2", label="pack current [A]")
    a.set(xlabel="time [s]", title="Battery")
    b = a.twinx()
    b.plot(t, out["soc"] * 100, color="k", ls="--", lw=1, label="state of charge [%]")
    b.set_ylim(0, 105)
    h1, l1 = a.get_legend_handles_labels()
    h2, l2 = b.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="center left")
    fig.suptitle(f"MuJoCo free flight: {s['mass_g']:.0f} g all-up incl. {s['payload_g']:.0f} g payload, "
                 f"T/W = {s['thrust_to_weight']:.2f}", y=0.995)
    fig.tight_layout()
    fig.savefig(path, dpi=130)


def main():
    rep = json.loads((RES / "design.json").read_text())
    render = "--no-render" not in sys.argv
    out, s, frames = simulate(rep, render=render)
    plot(out, s, RES / "flight.png")
    results = {"nominal": s}
    print(json.dumps(s, indent=2))
    if frames:
        import imageio.v2 as imageio

        imageio.mimsave(RES / "flight.mp4", frames, fps=25, quality=6, macro_block_size=None)
        for name, tt in (("liftoff", 6.0), ("hover", 15.0), ("translate", 24.0), ("landing", s["t_land_command_s"] + 5.0)):
            imageio.imwrite(RES / f"flight_{name}.png", frames[min(int(tt * 25), len(frames) - 1)])
    # robustness: what the same vehicle does with more payload, or if the thrusters under-deliver
    for label, kw in (("payload_50g", dict(payload=0.050)), ("payload_0g", dict(payload=0.0)),
                      ("thrust_minus_10pct", dict(thrust_scale=0.90)),
                      ("thrust_minus_20pct", dict(thrust_scale=0.80))):
        _, r, _ = simulate(rep, t_end=40.0, label=label, **kw)
        results[label] = r
        print(f"{label:20s} T/W={r['thrust_to_weight']:.3f} lifted_off={r['lifted_off']} max_alt={r['max_altitude_m']:.2f} m "
              f"hover_P={r['hover_power_batt_W']}")
    (RES / "flight.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
