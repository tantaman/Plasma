"""
System model of the "wet emitter" variant: the corona wires are replaced by wetted
fibres that electrospray charged droplets, everything else (duct, strips, staging,
alternating rails) stays as in sizing.py.

The stage physics comes from carrier.py (droplet drift + convection + outflow); the
duct, structure, battery and converter models are the ones in sizing.py, with these
additions, all carried as explicit assumptions in WetAssumptions:

  * a share x_ion of the supply current is carried by gas ions (dry spots, discharge
    at cone tips), more of it in the negative-polarity stages.  Ions are 10x faster
    than droplets, so they add current but almost no space charge and no thrust;
  * charge blown through a collector row is lost and costs k_leak times its share of
    the thrust (it neutralises fresh carriers of the next, opposite-polarity stage);
  * every HV module has a fixed mass on top of the part that scales with power;
  * liquid: what is sprayed during the flight (with margin) plus what the fibres hold.

Run:  python sim/wet.py            -> results/design_wet.json, results/wet_*.png, results/flight_wet.{json,png}
      (--no-flight skips the MuJoCo flight; --replot redraws the figures from the saved report)
"""
from __future__ import annotations

import functools
import json
import pathlib
from dataclasses import asdict, dataclass, replace

import numpy as np

import sizing as S
from carrier import CarrierProblem, mobility_band, spray_onset_field, stage_curve_carrier
from corona import MU_ION, StageCurve, mesh_cell, peek_field

G = S.G
RHO_AIR = S.RHO_AIR
CACHE = S.CACHE
RES = pathlib.Path(__file__).resolve().parent.parent / "results"
U_GRID = (1.0, 2.5, 4.0)        # duct velocities at which stage curves are tabulated [m/s]
SPLIT = 0.9                      # outflow cut at 90 % of the collector chord


@dataclass(frozen=True)
class WetAssumptions(S.Assumptions):
    # --- carriers (the range tried in the sensitivity study is given in brackets)
    Z: float = 2.0e-5            # m^2/(V s) droplet mobility [1.2e-5 .. 4e-5]; 50-300 nm drops at 50-100 % of the
    #                              Rayleigh limit come out at 1.3-3.5e-5 (carrier.mobility_band); ions are 2e-4
    q_over_m: float = 100.0      # C/kg charge per mass of the spray [30 .. 300]; cone-jet laws give 75-750 for K = 0.03-0.3 S/m
    gamma: float = 0.040         # N/m surface tension [0.035 .. 0.048]: EMI-Im 0.035, tri/tetraethylene glycol 0.045
    r_fibre: float = 40e-6       # m, radius of the wetted emitter (fibre or tooth tip)
    x_ion_pos: float = 0.25      # share of supply current carried by gas ions, positive stages [0 .. 0.6]
    x_ion_neg: float = 0.50      # ... negative stages, where discharge sets in earlier [0.25 .. 0.9 = corona takes over]
    k_leak: float = 2.0          # thrust lost per unit of charge blown through a collector row [1 .. 3]
    # --- override of sizing.Assumptions
    k_ehd: float = 0.80          # derate on computed body force: discrete sites, uneven wetting [0.65 .. 0.9]
    # --- hardware
    conv_fixed: float = 0.008    # kg per HV module that does not scale with power (controller, leads, potting)
    n_modules: int = 4
    wet_per_m: float = 0.20e-3   # kg per m of strip: wettable emitter, absorbent nose wrap, feeder threads, liquid held [0.10 .. 0.35]
    m_liquid_fixed: float = 0.004    # kg: priming losses, catch pads on the last grille
    liquid_margin: float = 2.0       # carried / consumed
    hover_time: float = 120.0        # s the battery and the liquid are sized for


# ------------------------------------------------------------------ stage curves
def _curve_file(tag: str, d: float, r_e: float, Z: float, u0: float, E_on: float) -> pathlib.Path:
    return CACHE / f"{tag}_d{d*1e3:.1f}_r{r_e*1e6:.1f}_Z{Z*1e6:.1f}_u{u0:.1f}_E{E_on/1e6:.2f}.json"


@functools.lru_cache(maxsize=None)
def _problem(d: float, r_e: float) -> CarrierProblem:
    g = S.cell_geometry(d, r_e)
    return CarrierProblem(mesh_cell(**g, n_wire=18, h_max=g["pitch"] / 60, h_coll=g["t_c"] / 6, split=SPLIT))


@functools.lru_cache(maxsize=None)
def stage_carrier(d: float, r_e: float, Z: float, u0: float, E_on: float, tag: str = "wet") -> StageCurve:
    """Cached stage curve for carriers of mobility Z in a duct flow u0 (outflow model)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    f = _curve_file(tag, d, r_e, Z, u0, E_on)
    if f.exists():
        return StageCurve.from_dict(json.loads(f.read_text()))
    c = stage_curve_carrier(_problem(d, r_e), E_on, Z, u0, n_pts=12, v_max_over_onset=3.6)
    f.write_text(json.dumps(c.to_dict()))
    return c


def stage_at(d: float, r_e: float, Z: float, E_on: float, V: float, u: float, tag: str = "wet"):
    """(current [A/m], body force [N/m], peak collector field, leak fraction) at voltage V and duct velocity u."""
    rows = []
    for u0 in U_GRID:
        c = stage_carrier(d, r_e, Z, u0, E_on, tag)
        I, F, Ec = (float(x) for x in c.at(V))
        lk = float(np.interp(V, c.V, c.meta["leak"]))
        rows.append((I, F, Ec, lk))
    rows = np.array(rows)
    return tuple(float(np.interp(u, U_GRID, rows[:, k])) for k in range(4))


def onset_voltage(d: float, r_e: float, E_on: float) -> float:
    return stage_carrier(d, r_e, 1.5e-5, U_GRID[1], E_on).V_onset if E_on else 0.0


# ---------------------------------------------------------------------- system
def evaluate(d: float, V: float, n: int, W: float, a: S.Assumptions = WetAssumptions(), battery: dict | None = None,
             carriers: str = "droplets", rails: str = "alternating") -> dict:
    """Thrust, power and mass of an n-stage craft.

    carriers="droplets": wet emitters, `a` must be a WetAssumptions.
    carriers="ions":     the dry corona stack with the same outflow stage model and the plain
                         sizing.Assumptions (a like-for-like reference for the frozen baseline).
    rails="alternating": grilles on +V/2, -V/2, +V/2 ... as in the baseline; every second stage emits
                         negative carriers.
    rails="monotonic":   grille potentials fall by V per stage, so every stage emits positive carriers;
                         the supply has to make n*V.
    """
    wet = carriers == "droplets"
    if wet:
        r_e, Z, tag = a.r_fibre, a.Z, "wet"
        E_on = spray_onset_field(r_e, a.gamma)
    else:
        r_e, Z, tag = 25e-6, MU_ION, "ion"
        E_on = peek_field(r_e)
    g = S.cell_geometry(d, r_e)
    A = W * W
    n_neg = n // 2 if rails == "alternating" else 0
    f_pol = 1.0 if wet else ((n - n_neg) + n_neg * a.neg_F) / n      # no thrust bonus assumed for negative sprays
    v = 2.0
    for _ in range(40):
        I1, F1, Ecoll, leak = stage_at(d, r_e, Z, E_on, V, v, tag)
        # Ions that pass the cut stay between the collector faces and land there.  Droplets are blown
        # on into the next stage: with alternating rails they meet a reversed field and fresh carriers
        # of the other sign (penalty k_leak); with monotonic rails they are simply accelerated again
        # (no penalty, and no credit taken).
        lost = a.k_leak * leak if (wet and rails == "alternating") else 0.0
        dp = a.k_ehd * f_pol * F1 * max(1.0 - lost, 0.0) / g["pitch"] * a.fill
        KL = S.stage_loss(g, v, a)
        v_new = np.sqrt(max(2 * n * dp / (RHO_AIR * (1 + a.K0 + n * KL)), 1e-9))
        if abs(v_new - v) < 1e-5:
            break
        v = 0.5 * (v + v_new)
    aug = 2.0 / (1 + a.K0 + n * KL)
    T = n * dp * A * aug
    strip_len = a.fill * A / g["pitch"]
    I_carrier = I1 * strip_len                                   # carrier current of one stage
    if wet:
        I_pos, I_neg = I_carrier / (1.0 - a.x_ion_pos), I_carrier / (1.0 - a.x_ion_neg)
    else:
        I_pos, I_neg = I_carrier, I_carrier * a.neg_I
    P_el = ((n - n_neg) * I_pos + n_neg * I_neg) * V
    V_rr = V if rails == "alternating" else n * V
    # --- masses
    Sp = g["c_c"] + g["c_d"] + d
    L = n * Sp + 0.5 * (g["c_c"] + g["c_d"]) + 0.10
    m_strips = (n + 1) * strip_len * (S.strip_mass_per_m(g, a) + (a.wet_per_m if wet else 0.0))
    m_frame = (2 * 4 * W + 2 * 4 * W) * a.tube_g_per_m + (4 * L + L) * a.longeron_g_per_m
    m_skin = a.film_wall * (4 * W * L + np.pi * 0.035 * 4 * W) * 1.15 + 0.004
    p_conv = a.conv_kw_per_kg * 1e3 * min(1.0, (20e3 / V_rr) ** a.conv_v_exp)
    m_conv = P_el / p_conv + (a.n_modules * a.conv_fixed if wet else 0.0)
    P_batt = P_el / a.conv_eff
    e_batt = max(P_batt * a.hover_time / 3600 / a.batt_dod, P_batt / a.batt_max_c)
    m_batt = e_batt / a.batt_wh_per_kg
    if battery is not None:
        e_batt = battery["cells"] * 3.7 * battery["capacity_ah"]
        m_batt = battery["mass"]
    endurance = e_batt * a.batt_dod * 3600 / P_batt if P_batt else 0.0
    mdot_liq = n * I_carrier / a.q_over_m if wet else 0.0        # kg/s; no credit for re-use down the stack
    m_liquid = (a.liquid_margin * mdot_liq * max(a.hover_time, endurance) + a.m_liquid_fixed) if wet else 0.0
    m_struct = m_strips + m_frame + m_skin + a.m_misc
    m_empty = m_struct + m_conv + m_batt + a.m_avionics + m_liquid
    m_total = m_empty + a.m_payload
    return dict(d=d, V=V, n=n, W=W, L=L, A=A, T=T, P_el=P_el, P_batt=P_batt, aug=aug, KL=KL, v_duct=v,
                I_stage=I_pos, I_carrier=I_carrier, dp_stage=dp, theta=T / P_el if P_el > 0 else 0.0,
                thrust_density=T / A, E_avg=V / d, E_coll=Ecoll, V_rr=V_rr, leak=leak, rails=rails, carriers=carriers,
                V_onset=stage_carrier(d, r_e, Z, U_GRID[1], E_on, tag).V_onset, E_on=E_on,
                uA_per_m=I1 * 1e6, body_N_per_kW=F1 / (I1 * V) * 1e3 if I1 > 0 else 0.0,
                mdot_liq=mdot_liq, m_liquid=m_liquid,
                m_strips=m_strips, m_frame=m_frame, m_skin=m_skin, m_struct=m_struct, m_conv=m_conv,
                m_batt=m_batt, e_batt_wh=e_batt, c_rate=P_batt / e_batt if e_batt else 0.0, p_conv=p_conv,
                m_empty=m_empty, m_total=m_total, TW=T / (m_total * G) if m_total else 0.0, endurance=endurance)


def payload_at(r: dict, tw: float) -> float:
    """Payload [kg] this craft carries with thrust/weight = tw at full throttle."""
    return r["T"] / (G * tw) - r["m_empty"]


# ------------------------------------------------------------- frozen design point
# "Ion Lantern W": the baseline airframe (same duct, same 4 grilles of integrated strips, same gap
# and stage voltage) with wet emitters, monotonic rails and a power system sized for the new load.
D_GAP = 0.030
N_STAGES = 3
DUCT_W = 1.00
V_MAX = 20.0e3            # V across every stage; grilles at +30, +10, -10, -30 kV
RAILS = "monotonic"
M_PAYLOAD = 0.100         # kg design payload
BATTERY = dict(
    name="3S 450 mAh 75C LiPo (mass taken as 44 g; check the listing of the pack you buy)",
    cells=3, capacity_ah=0.45, mass=0.044,
    r_cell=0.0133,        # ohm per cell (6 mOhm*Ah / 0.45 Ah)
    r_leads=0.006,
)
TAU_SPRAY = 0.05          # s, assumed response of the spray to a voltage step (corona: 0.015 s)


def slab_check() -> list[dict]:
    """1-D check of drift + convection against the exact solution (carrier.slab_analytic)."""
    from carrier import mesh_slab, slab_analytic

    d, h, Z, E0 = 0.03, 0.006, 2.0e-5, 2.0e5
    prob = CarrierProblem(mesh_slab(d, h, n=80))
    rows = []
    for u in (0.0, 2.5, 10.0):
        for V in (8e3, 12e3, 16e3, 20e3):
            prob.w = u / (Z * V)
            lo, hi, i, dg = 0.0, None, 0.3 * 9 / 8 * 8.854e-12 * h / d**3, None
            for _ in range(60):
                _, _, dg, it = prob.solve_unit(i)
                if it < 0 or dg["e_em"] < E0 / V:
                    hi = i
                else:
                    lo = i
                i = 0.5 * (lo + hi) if hi else 2 * i
                if hi and (hi - lo) < 1e-6 * hi:
                    break
            J = Z * dg["i_em"] * V * V / h
            Vx, dpx = slab_analytic(d, E0, Z, u, J)
            rows.append(dict(u=u, V=V, J=J, dp=dg["fx"] * V * V / h, V_exact=float(Vx[0]), dp_exact=float(dpx[0])))
    return rows


def thrust_map(a: WetAssumptions, n_pts: int = 50) -> dict:
    v_on = evaluate(D_GAP, V_MAX, N_STAGES, DUCT_W, a, BATTERY, rails=RAILS)["V_onset"]
    V = np.linspace(v_on, 1.1 * V_MAX, n_pts)
    rows = [evaluate(D_GAP, v, N_STAGES, DUCT_W, a, BATTERY, rails=RAILS) for v in V]
    return dict(V=V.tolist(), T=[r["T"] for r in rows], P=[r["P_el"] for r in rows],
                T_body=[r["T"] / r["aug"] for r in rows], v_duct=[r["v_duct"] for r in rows])


def design_report() -> dict:
    import design as D

    frozen = json.loads((RES / "design.json").read_text())
    tw_ref = frozen["design"]["TW"]                   # the baseline's thrust margin: payloads are compared at this T/W
    a = WetAssumptions(m_avionics=D.M_AVIONICS, m_payload=M_PAYLOAD)
    dry_a = replace(S.Assumptions(), m_avionics=D.M_AVIONICS, m_payload=D.M_PAYLOAD)
    ev = lambda aa, **kw: evaluate(D_GAP, V_MAX, kw.pop("n", N_STAGES), DUCT_W, aa, **{"battery": BATTERY, "rails": RAILS, **kw})
    cap = lambda r: payload_at(r, tw_ref)
    mk2 = ev(a)
    dry = evaluate(D.D_GAP, D.V_MAX, D.N_STAGES, D.DUCT_W, dry_a, battery=D.BATTERY, carriers="ions")
    out = dict(assumptions=asdict(a), design=mk2, tw_ref=tw_ref, payload_capacity=cap(mk2),
               cell=S.cell_geometry(D_GAP, a.r_fibre), battery=BATTERY, rails=RAILS, tau_v=TAU_SPRAY,
               frozen_baseline=dict(T=frozen["design"]["T"], P_el=frozen["design"]["P_el"], TW=tw_ref,
                                    m_total=frozen["design"]["m_total"], payload=frozen["assumptions"]["m_payload"]),
               dry_outflow=dict(dry, payload_capacity=cap(dry)))
    # --- one-at-a-time sensitivity of the payload carried at the baseline's thrust margin
    sens = []
    for label, lo, hi in [
        ("droplet mobility Z (4e-5 / 1.2e-5 m²/Vs)", dict(Z=4e-5), dict(Z=1.2e-5)),
        ("ion share of current (0.6 / 0)", dict(x_ion_pos=0.6), dict(x_ion_pos=0.0)),
        ("force derate k_ehd (0.65 / 0.9)", dict(k_ehd=0.65), dict(k_ehd=0.9)),
        ("surface tension (0.048 / 0.035 N/m)", dict(gamma=0.048), dict(gamma=0.035)),
        ("charge per mass (30 / 300 C/kg)", dict(q_over_m=30.0), dict(q_over_m=300.0)),
        ("converter kW/kg (2.0 / 3.5)", dict(conv_kw_per_kg=2.0), dict(conv_kw_per_kg=3.5)),
        ("converter fixed mass per module (15 / 4 g)", dict(conv_fixed=0.015), dict(conv_fixed=0.004)),
        ("fixed duct loss K0 (0.45 / 0.15)", dict(K0=0.45), dict(K0=0.15)),
        ("stage loss (MIT estimate 0.14 / as modelled 0.09)", dict(K_support=0.07), dict()),
        ("wet hardware on strips (0.35 / 0.10 g/m)", dict(wet_per_m=0.35e-3), dict(wet_per_m=0.10e-3)),
    ]:
        sens.append(dict(label=label, lo=cap(ev(replace(a, **lo))), hi=cap(ev(replace(a, **hi)))))
    out["sensitivity"] = sens
    cons = replace(a, Z=3e-5, x_ion_pos=0.5, k_ehd=0.65, conv_kw_per_kg=2.0, K0=0.4, q_over_m=30.0, wet_per_m=0.3e-3)
    opt = replace(a, Z=1.4e-5, x_ion_pos=0.1, k_ehd=0.9, conv_kw_per_kg=3.5, K0=0.2)
    out["scenarios"] = {k: dict(payload_capacity=cap(r), TW_design_payload=r["TW"], P_el=r["P_el"], T=r["T"])
                        for k, r in (("conservative", ev(cons)), ("baseline", mk2), ("optimistic", ev(opt)))}
    # --- layouts and stage count (battery sized for hover_time in every case)
    lay = []
    for n in (2, 3, 4, 6):
        for name, aa, rails in (("monotonic rails, all stages positive", a, "monotonic"),
                                ("alternating rails, negative spray works", a, "alternating"),
                                ("alternating rails, negative stages fall back to corona", replace(a, x_ion_neg=0.9), "alternating")):
            r = evaluate(D_GAP, V_MAX, n, DUCT_W, aa, rails=rails)
            lay.append(dict(layout=name, n=n, T=r["T"], P_el=r["P_el"], V_supply=r["V_rr"], m_empty=r["m_empty"],
                            m_conv=r["m_conv"], m_batt=r["m_batt"], payload_capacity=cap(r)))
    out["layouts"] = lay
    # --- carrier mobility sweep at the design point
    mob = []
    for Z in (0.7e-5, 1.2e-5, 2.0e-5, 3.0e-5, 4.0e-5, 6.0e-5, 1.0e-4):
        r = ev(replace(a, Z=Z))
        mob.append(dict(Z=Z, body_N_per_kW=r["body_N_per_kW"], theta=r["theta"], P_el=r["P_el"], T=r["T"],
                        leak=r["leak"], uA_per_m=r["uA_per_m"], payload_capacity=cap(r)))
    out["mobility_sweep"] = mob
    out["mobility_band"] = mobility_band(a.gamma)     # 50-300 nm drops at 50-100 % of the Rayleigh limit
    # --- stage curves at the design gap, 2.5 m/s (for the figure and the bench-test table)
    E_on = spray_onset_field(a.r_fibre, a.gamma)
    cw = stage_carrier(D_GAP, a.r_fibre, a.Z, U_GRID[1], E_on, "wet")
    ci = stage_carrier(D_GAP, 25e-6, MU_ION, U_GRID[1], peek_field(25e-6), "ion")
    out["stage_curves"] = dict(droplets=cw.to_dict(), ions=ci.to_dict(), frozen=S.stage(D.D_GAP, D.R_WIRE).to_dict())
    out["validation_slab"] = slab_check()
    # --- thrust map and hover point for the flight simulation
    tm = thrust_map(a)
    out["thrust_map"] = tm
    W_N = mk2["m_total"] * G
    out["hover"] = dict(V=float(np.interp(W_N, tm["T"], tm["V"])), P_el=float(np.interp(W_N, tm["T"], tm["P"])), weight_N=W_N)
    out["hover"]["P_batt"] = out["hover"]["P_el"] / a.conv_eff
    out["hover"]["endurance_s"] = mk2["e_batt_wh"] * a.batt_dod * 3600 / out["hover"]["P_batt"]
    return out


# ------------------------------------------------------------------------ figures
INK, INK2, MUTED, AXIS, GRID, SURF = "#0b0b0b", "#52514e", "#898781", "#c3c2b7", "#e9e8e3", "#fcfcfb"
C_ION, C_DROP, C_THIRD = "#2a78d6", "#eb6834", "#1baf7a"      # categorical slots 1-3 (validated all-pairs)
PARTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]


def _style(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURF)
    ax.set_title(title, fontsize=10.5, color=INK, loc="left", pad=8)
    ax.set_xlabel(xlabel, fontsize=9, color=INK2)
    ax.set_ylabel(ylabel, fontsize=9, color=INK2)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=8.5, length=0)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(AXIS)


def make_figures(rep: dict, outdir: pathlib.Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    a = rep["assumptions"]
    # ============================================================ stage physics
    fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.4), facecolor=SURF)
    # (a) exact 1-D check
    a0 = ax[0, 0]
    rows = rep["validation_slab"]
    from carrier import slab_analytic

    for u, c in zip((0.0, 2.5, 10.0), (C_ION, C_DROP, C_THIRD)):
        pts = [r for r in rows if r["u"] == u]
        J = np.linspace(0.0, 1.15 * max(r["J"] for r in pts), 80)[1:]
        Vx, dpx = slab_analytic(0.03, 2.0e5, 2.0e-5, u, J)
        a0.plot(Vx / 1e3, dpx, color=c, lw=2, label=f"exact, air at {u:g} m/s")
        a0.plot([r["V"] / 1e3 for r in pts], [r["dp"] for r in pts], "o", ms=7, mfc=SURF, mec=c, mew=1.8)
    a0.plot([], [], "o", ms=7, mfc=SURF, mec=INK2, mew=1.8, label="this solver")
    err = max(abs(r["dp"] / r["dp_exact"] - 1) for r in rows)
    _style(a0, f"A. Check against the exact 1-D solution (worst error {100*err:.1f} %)", "voltage across a 30 mm slab [kV]",
           "pressure rise [Pa]")
    a0.set_xlim(6, 21)
    a0.legend(frameon=False, fontsize=8.5, labelcolor=INK2)
    # (b), (c) stage curves
    sc = rep["stage_curves"]
    V = np.linspace(8e3, 22e3, 60)
    for k, (key, lab, ylab) in enumerate((("I", "B. Supply current of one stage", "current per metre of emitter [µA/m]"),
                                          ("F", "C. Force on the air from one stage", "body force per metre of emitter [mN/m]"))):
        a0 = ax[0, 1] if k == 0 else ax[1, 0]
        scale = 1e6 if key == "I" else 1e3
        for name, c, lbl in (("ions", C_ION, "corona ions (dry wire)"), ("droplets", C_DROP, "sprayed droplets (wet fibre)")):
            y = np.interp(V, sc[name]["V"], sc[name][key], left=0.0) * scale
            a0.plot(V / 1e3, y, color=c, lw=2, label=lbl)
            a0.annotate(f"{y[np.searchsorted(V, 20e3)]:.0f}", (20, y[np.searchsorted(V, 20e3)]), textcoords="offset points",
                        xytext=(-4, 7), ha="right", fontsize=8.5, color=INK2)
        a0.axvline(20, color=AXIS, lw=1)
        _style(a0, lab, "stage voltage [kV]", ylab)
        a0.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper left")
    # (d) thrust per watt against mobility
    a0 = ax[1, 1]
    mob = rep["mobility_sweep"]
    Zs = np.array([m["Z"] for m in mob])
    a0.plot(Zs * 1e5, [m["body_N_per_kW"] for m in mob], color=C_DROP, lw=2, marker="o", ms=7, mfc=C_DROP, mec=SURF, mew=1.5,
            label="this stage, 20 kV across 30 mm")
    dry = rep["dry_outflow"]
    a0.plot([20.0], [dry["body_N_per_kW"]], "s", ms=8, mfc=C_ION, mec=SURF, mew=1.5, label="corona ions")
    zb = [z * 1e5 for z in rep["mobility_band"]]
    a0.axvspan(*zb, color=C_DROP, alpha=0.10, lw=0)
    a0.text((zb[0] * zb[1]) ** 0.5, 176, "50–300 nm drops at\n50–100 % of the\nRayleigh limit", fontsize=8.5, color=INK2, ha="center", va="top")
    a0.set_xscale("log")
    a0.set_xticks([0.7, 1, 2, 4, 10, 20], ["0.7", "1", "2", "4", "10", "20"])
    a0.set_ylim(0, 185)
    _style(a0, "D. Body force per electrical watt against carrier mobility", "carrier mobility [10⁻⁵ m²/(V·s)]", "N per kW")
    a0.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="center right")
    fig.tight_layout()
    fig.savefig(outdir / "wet_stage.png", dpi=140, facecolor=SURF)
    plt.close(fig)

    # ============================================================ system
    fig, ax = plt.subplots(2, 2, figsize=(12.5, 8.8), facecolor=SURF)
    b = rep["design"]
    fz = json.loads((RES / "design.json").read_text())
    fb, fa = fz["design"], fz["assumptions"]
    # (a) mass budgets
    a0 = ax[0, 0]
    names = ["structure", "HV converters", "battery", "avionics", "liquid", "payload"]
    crafts = [("Baseline as frozen\n(corona, ±10 kV)", [fb["m_struct"], fb["m_conv"], fb["m_batt"], fa["m_avionics"], 0.0, fa["m_payload"]], fb["T"]),
              ("Wet emitters\n(same airframe)", [b["m_struct"], b["m_conv"], b["m_batt"], a["m_avionics"], b["m_liquid"], a["m_payload"]], b["T"])]
    for row, (lab, parts, T) in enumerate(crafts):
        y = 1 - row
        left = 0.0
        for k, m in enumerate(parts):
            if m <= 0:
                continue
            a0.barh(y, m * 1e3 - 2.0, left=left + 1.0, height=0.34, color=PARTS[k], label=names[k] if row == 1 or k != 4 else None)
            if m * 1e3 > 18:
                a0.text(left + m * 500, y - 0.27, f"{m*1e3:.0f}", ha="center", va="top", fontsize=8.5, color=INK2)
            left += m * 1e3
        a0.plot([T / G * 1e3] * 2, [y - 0.24, y + 0.24], color=INK, lw=2)
        a0.text(T / G * 1e3 + 6, y, f"thrust {T/G*1e3:.0f} gf", va="center", fontsize=8.5, color=INK)
        a0.text(-8, y, lab, ha="right", va="center", fontsize=9, color=INK)
    h, l = a0.get_legend_handles_labels()
    seen = dict(zip(l, h))
    a0.legend([seen[n] for n in names if n in seen], [n for n in names if n in seen], frameon=False, fontsize=8.5, ncol=3,
              loc="upper center", bbox_to_anchor=(0.5, 1.0), labelcolor=INK2)
    _style(a0, "A. Mass budget against maximum thrust", "grams", "")
    a0.set_yticks([])
    a0.set_ylim(-0.6, 2.05)
    a0.set_xlim(0, 720)
    a0.grid(axis="y", visible=False)
    # (b) sensitivity
    a0 = ax[0, 1]
    sens = sorted(rep["sensitivity"], key=lambda s: abs(s["hi"] - s["lo"]))
    for k, s in enumerate(sens):
        a0.barh(k, (s["hi"] - s["lo"]) * 1e3, left=s["lo"] * 1e3, height=0.5, color=C_DROP)
        a0.text(min(s["lo"], s["hi"]) * 1e3 - 3, k, f"{s['lo']*1e3:.0f}", ha="right", va="center", fontsize=8, color=INK2)
        a0.text(max(s["lo"], s["hi"]) * 1e3 + 3, k, f"{s['hi']*1e3:.0f}", ha="left", va="center", fontsize=8, color=INK2)
    a0.set_yticks(range(len(sens)), [s["label"] for s in sens], fontsize=8.5)
    a0.axvline(rep["payload_capacity"] * 1e3, color=INK, lw=1.2)
    a0.axvline(rep["frozen_baseline"]["payload"] * 1e3, color=MUTED, lw=1.2)
    a0.text(rep["frozen_baseline"]["payload"] * 1e3 - 3, len(sens) - 0.3, "baseline 20 g", fontsize=8, color=INK2, ha="right")
    a0.text(rep["payload_capacity"] * 1e3 + 3, len(sens) - 0.3, f"nominal {rep['payload_capacity']*1e3:.0f} g", fontsize=8, color=INK, ha="left")
    _style(a0, f"B. Payload at thrust/weight {rep['tw_ref']:.2f}: one assumption at a time", "payload [g]", "")
    a0.tick_params(axis="y", colors=INK2)
    a0.grid(axis="y", visible=False)
    a0.set_ylim(-0.6, len(sens) + 0.2)
    lo = min(min(s["lo"], s["hi"]) for s in sens) * 1e3
    hi = max(max(s["lo"], s["hi"]) for s in sens) * 1e3
    a0.set_xlim(min(lo, 0) - 25, hi + 25)
    # (c) layouts
    a0 = ax[1, 0]
    lay = rep["layouts"]
    for name, c, mk in (("monotonic rails, all stages positive", C_DROP, "o"), ("alternating rails, negative spray works", C_ION, "s"),
                        ("alternating rails, negative stages fall back to corona", C_THIRD, "^")):
        pts = [r for r in lay if r["layout"] == name]
        a0.plot([r["n"] for r in pts], [r["payload_capacity"] * 1e3 for r in pts], color=c, lw=2, marker=mk, ms=7.5, mec=SURF, mew=1.5,
                label=name)
        if name.startswith("monotonic"):
            for r in pts:
                a0.annotate(f"{r['V_supply']/1e3:.0f} kV", (r["n"], r["payload_capacity"] * 1e3), textcoords="offset points",
                            xytext=(0, 9), ha="center", fontsize=8, color=INK2)
    a0.axhline(rep["frozen_baseline"]["payload"] * 1e3, color=MUTED, lw=1.2)
    a0.text(6.0, rep["frozen_baseline"]["payload"] * 1e3 + 8, "baseline 20 g", fontsize=8, color=INK2, ha="right")
    a0.axhline(0, color=AXIS, lw=1)
    _style(a0, "C. Rail layout and stage count (supply voltage marked; battery sized for 120 s)", "acceleration stages",
           "payload at the baseline's thrust margin [g]")
    a0.set_xticks([2, 3, 4, 6])
    a0.set_ylim(-320, 360)
    a0.legend(frameon=False, fontsize=8.5, labelcolor=INK2, loc="upper left")
    # (d) payload against mobility
    a0 = ax[1, 1]
    a0.plot(Zs * 1e5, [m["payload_capacity"] * 1e3 for m in mob], color=C_DROP, lw=2, marker="o", ms=7, mec=SURF, mew=1.5)
    a0.axvspan(*zb, color=C_DROP, alpha=0.10, lw=0)
    a0.axhline(rep["frozen_baseline"]["payload"] * 1e3, color=MUTED, lw=1.2)
    a0.text(9.8, rep["frozen_baseline"]["payload"] * 1e3 + 5, "baseline 20 g", fontsize=8, color=INK2, ha="right")
    a0.axhline(0, color=AXIS, lw=1)
    a0.set_xscale("log")
    a0.set_xticks([0.7, 1, 2, 4, 10], ["0.7", "1", "2", "4", "10"])
    _style(a0, "D. Payload against droplet mobility (shaded: expected range)", "droplet mobility [10⁻⁵ m²/(V·s)]", "payload [g]")
    fig.tight_layout()
    fig.savefig(outdir / "wet_sizing.png", dpi=140, facecolor=SURF)
    plt.close(fig)


def fly(rep: dict) -> dict:
    """Free flight of the wet-emitter craft in the MuJoCo model of flight.py."""
    import flight

    b = rep["design"]
    frep = dict(rep, design=dict(b, m_strips=b["m_strips"] + b["m_liquid"]))      # the liquid rides on the strips
    out, s, _ = flight.simulate(frep, t_end=260.0, label="wet_nominal")
    flight.plot(out, s, RES / "flight_wet.png")
    results = {"nominal": s}
    for label, kw in (("payload_150g", dict(payload=0.150)), ("payload_50g", dict(payload=0.050)),
                      ("thrust_minus_10pct", dict(thrust_scale=0.90)), ("thrust_minus_20pct", dict(thrust_scale=0.80))):
        _, r, _ = flight.simulate(frep, t_end=40.0, label=label, **kw)
        results[label] = r
    return results


if __name__ == "__main__":
    import sys

    if "--replot" in sys.argv:
        rep = json.loads((RES / "design_wet.json").read_text())
    else:
        rep = design_report()
        (RES / "design_wet.json").write_text(json.dumps(rep, indent=2))
    make_figures(rep, RES)
    b, d0 = rep["design"], rep["dry_outflow"]
    print(f"baseline as frozen : T={rep['frozen_baseline']['T']:.2f} N  P_el={rep['frozen_baseline']['P_el']:.0f} W  T/W={rep['tw_ref']:.3f}")
    print(f"baseline, outflow  : T={d0['T']:.2f} N  P_el={d0['P_el']:.0f} W  T/W={d0['TW']:.3f}  body {d0['body_N_per_kW']:.1f} N/kW  "
          f"payload at T/W {rep['tw_ref']:.2f}: {d0['payload_capacity']*1e3:.0f} g")
    print(f"wet emitters       : T={b['T']:.2f} N  P_el={b['P_el']:.0f} W  theta={b['theta']*1e3:.1f} N/kW  body {b['body_N_per_kW']:.1f} N/kW  "
          f"mass={b['m_total']*1e3:.0f} g  T/W={b['TW']:.3f}  payload at T/W {rep['tw_ref']:.2f}: {rep['payload_capacity']*1e3:.0f} g")
    print(f"  masses [g]: struct {b['m_struct']*1e3:.0f} (strips {b['m_strips']*1e3:.0f}) conv {b['m_conv']*1e3:.0f} batt {b['m_batt']*1e3:.0f} "
          f"liquid {b['m_liquid']*1e3:.1f}; {b['uA_per_m']:.0f} uA/m, {b['I_stage']*1e3:.2f} mA/stage, liquid {b['mdot_liq']*1e6:.0f} mg/s, "
          f"leak {b['leak']:.3f}, V_on {b['V_onset']/1e3:.1f} kV, E_coll {b['E_coll']/1e6:.2f} MV/m, v_duct {b['v_duct']:.2f} m/s")
    print("  hover:", {k: round(v, 2) for k, v in rep["hover"].items()})
    print("  scenarios:", {k: round(v["payload_capacity"] * 1e3) for k, v in rep["scenarios"].items()})
    for s_ in rep["sensitivity"]:
        print(f"   {s_['label']:52s} {s_['lo']*1e3:6.0f} .. {s_['hi']*1e3:6.0f} g")
    for r in rep["layouts"]:
        print(f"   n={r['n']} {r['layout']:56s} T={r['T']:.2f} P={r['P_el']:4.0f} W supply {r['V_supply']/1e3:3.0f} kV payload {r['payload_capacity']*1e3:5.0f} g")
    if "--no-flight" not in sys.argv:
        res = fly(rep)
        (RES / "flight_wet.json").write_text(json.dumps(res, indent=2))
        for k, r in res.items():
            print(f"   flight {k:20s} T/W={r['thrust_to_weight']:.3f} lifted_off={r['lifted_off']} airborne={r['airborne_time_s']:.0f} s "
                  f"hover_P_batt={r['hover_power_batt_W']}")
