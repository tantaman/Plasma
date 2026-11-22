"""
System model: electrode cell (FEM) -> duct aerodynamics -> masses -> thrust/weight.

Everything that is an assumption rather than a computed or measured number lives in
`Assumptions`, with its source, so the sensitivity of the conclusion to each one can
be checked (see sweep()/sensitivity()).
"""
from __future__ import annotations

import functools
import json
import pathlib
from dataclasses import asdict, dataclass, replace

import numpy as np

from corona import CoronaProblem, StageCurve, mesh_cell, peek_field

G = 9.80665
RHO_AIR = 1.20
NU_AIR = 1.5e-5
CACHE = pathlib.Path(__file__).resolve().parent.parent / "results" / "stage_curves"


@dataclass(frozen=True)
class Assumptions:
    # --- plasma / aerodynamics
    k_ehd: float = 0.90        # derate on FEM body force (unvalidated strip geometry, mobility spread)
    K0: float = 0.30           # fixed duct loss (inlet, exit, non-uniformity); calibrated to MIT single stage
    K_support: float = 0.02    # per-stage blockage of rails/threads
    fill: float = 0.90         # fraction of duct cross-section covered by active electrodes
    # --- battery (high-rate LiPo pouch pack, FPV racing class)
    batt_wh_per_kg: float = 125.0   # pack level incl. leads (Tattu R-Line 6S 1050 mAh: 23.3 Wh / 186 g)
    batt_max_c: float = 45.0        # sustained discharge rate allowed in sizing
    batt_dod: float = 0.80
    batt_mohm_ah: float = 6.0       # cell resistance x capacity [mOhm*Ah] (~6 mOhm for a 1 Ah race cell)
    cells: int = 6
    # --- HV converter (resonant inverter + transformer + Cockcroft-Walton)
    conv_kw_per_kg: float = 2.7       # short-duration rating up to 20 kV rail-to-rail (MIT: 2.7 at 20 kV, 1.5 kW, 40 s)
    conv_v_exp: float = 0.6           # above 20 kV: ~ (20 kV / V)^exp  (1.15-1.2 kW/kg demonstrated at 40 kV, continuous)
    conv_eff: float = 0.88
    # --- negative-corona stages (every second stage), relative to a positive stage at the same voltage
    neg_F: float = 1.10               # MIT patent FIG. 19: 0.80 N at 16 kV negative vs ~0.72 N positive
    neg_I: float = 1.30               # ... 6.2 mA vs ~4.7 mA
    # --- structure
    film_strip: float = 5.4e-3   # kg/m^2, 6 um metallised polypropylene capacitor film
    film_wall: float = 8.4e-3    # kg/m^2, 6 um PET (Mylar)
    foam: float = 35.0           # kg/m^3, XPS/Depron leading-edge stick
    wire_g_per_m: float = 0.038e-3   # kg/m, 50 um tungsten
    glue_per_m: float = 0.05e-3      # kg/m of strip
    tube_g_per_m: float = 2.2e-3     # kg/m, 1.5 x 0.7 mm pultruded carbon tube
    longeron_g_per_m: float = 3.8e-3  # kg/m, 2 x 1 mm carbon tube
    m_misc: float = 0.012        # threads, bus wires, joints, feet
    # --- fixed equipment
    m_avionics: float = 0.012    # flight controller, receiver, ToF/flow sensor, wiring
    m_payload: float = 0.020
    hover_time: float = 60.0     # s of hover the battery is sized for


def cell_geometry(d: float, r_e: float = 25e-6) -> dict:
    """Proportions of the integrated strip, scaled with the gap d."""
    return dict(d=d, pitch=d, r_e=r_e, t_c=max(0.0010, 0.05 * d), c_c=0.5 * d, c_d=0.6 * d, mode="cascade")


@functools.lru_cache(maxsize=None)
def stage(d: float, r_e: float = 25e-6) -> StageCurve:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"cascade_d{d*1e3:.1f}_r{r_e*1e6:.1f}.json"
    if f.exists():
        return StageCurve.from_dict(json.loads(f.read_text()))
    g = cell_geometry(d, r_e)
    m = mesh_cell(**g, n_wire=18, h_max=g["pitch"] / 60, h_coll=g["t_c"] / 6)
    c = stage_curve_safe(CoronaProblem(m), peek_field(r_e))
    f.write_text(json.dumps(c.to_dict()))
    return c


def stage_curve_safe(prob, E_c):
    from corona import stage_curve

    return stage_curve(prob, E_c, n_pts=16, v_max_over_onset=4.5)


def strip_mass_per_m(g: dict, a: Assumptions) -> float:
    film = a.film_strip * (2 * g["c_c"] + g["c_d"])
    stick = a.foam * g["t_c"] * 0.003
    return film + stick + a.wire_g_per_m + a.glue_per_m


def stage_loss(g: dict, v: float, a: Assumptions) -> float:
    """Per-stage pressure-loss coefficient (referred to the duct dynamic pressure)."""
    chord = g["c_c"] + g["c_d"]
    re_c = max(v * chord / NU_AIR, 50.0)
    cf = 1.328 / np.sqrt(re_c)                       # Blasius, each face
    k_fric = 2 * cf * chord / g["pitch"]
    k_nose = 0.35 * g["t_c"] / g["pitch"] / (1 - g["t_c"] / g["pitch"]) ** 2
    re_w = max(v * 2 * g["r_e"] / NU_AIR, 0.1)
    cd_w = 1.0 + 10.0 / re_w ** (2 / 3)              # cylinder at Re ~ 1-100
    k_wire = cd_w * 2 * g["r_e"] / g["pitch"]
    return k_fric + k_nose + k_wire + a.K_support


def evaluate(d: float, V: float, n: int, W: float, a: Assumptions = Assumptions(), r_e: float = 25e-6,
             battery: dict | None = None) -> dict:
    """Thrust, power and mass of an n-stage square-duct craft of side W at stage voltage V."""
    g = cell_geometry(d, r_e)
    cur = stage(d, r_e)
    I1, F1, Ecoll = (float(x) for x in cur.at(V))
    A = W * W
    f_pol = ((n - n // 2) + (n // 2) * a.neg_F) / n   # mean force factor over +/- stages
    dp = a.k_ehd * f_pol * F1 / g["pitch"] * a.fill  # mean EHD pressure rise per stage [Pa]
    # duct: T = 2 n dp A / (1 + K0 + n K_L), iterate on velocity because K_L depends on Re
    v = 2.0
    for _ in range(30):
        KL = stage_loss(g, v, a)
        v_new = np.sqrt(max(2 * n * dp / (RHO_AIR * (1 + a.K0 + n * KL)), 1e-9))
        if abs(v_new - v) < 1e-6:
            break
        v = 0.5 * (v + v_new)
    aug = 2.0 / (1 + a.K0 + n * KL)
    T = n * dp * A * aug
    strip_len = a.fill * A / g["pitch"]              # metres of strip per grille
    I_stage = I1 * strip_len                         # positive-stage current
    n_neg = n // 2
    P_el = (n - n_neg + n_neg * a.neg_I) * I_stage * V
    V_rr = V                                         # grilles alternate between the +V/2 and -V/2 rails
    # --- masses
    S = g["c_c"] + g["c_d"] + d
    L = n * S + 0.5 * (g["c_c"] + g["c_d"]) + 0.10   # stack + inlet lip + exit margin
    m_strips = (n + 1) * strip_len * strip_mass_per_m(g, a)
    m_frame = (2 * 4 * W + 2 * 4 * W) * a.tube_g_per_m + (4 * L + L) * a.longeron_g_per_m
    r_lip = 0.035
    m_skin = a.film_wall * (4 * W * L + np.pi * r_lip * 4 * W) * 1.15 + 0.004  # walls + lip + seams, formers
    p_conv = a.conv_kw_per_kg * 1e3 * min(1.0, (20e3 / V_rr) ** a.conv_v_exp)
    m_conv = P_el / p_conv
    P_batt = P_el / a.conv_eff
    e_need = P_batt * a.hover_time / 3600 / a.batt_dod            # Wh
    e_rate = P_batt / a.batt_max_c                                 # Wh needed to respect the C limit
    e_batt = max(e_need, e_rate)
    m_batt = e_batt / a.batt_wh_per_kg
    if battery is not None:  # a specific pack instead of a sized one
        e_batt = battery["cells"] * 3.7 * battery["capacity_ah"]
        m_batt = battery["mass"]
    m_struct = m_strips + m_frame + m_skin + a.m_misc
    m_total = m_struct + m_conv + m_batt + a.m_avionics + a.m_payload
    return dict(d=d, V=V, n=n, W=W, L=L, A=A, T=T, P_el=P_el, P_batt=P_batt, aug=aug, KL=KL, v_duct=v,
                I_stage=I_stage, dp_stage=dp, theta=T / P_el if P_el > 0 else 0.0, thrust_density=T / A,
                E_avg=V / d, E_coll=Ecoll, V_rr=V_rr, V_onset=cur.V_onset,
                m_strips=m_strips, m_frame=m_frame, m_skin=m_skin, m_struct=m_struct, m_conv=m_conv,
                m_batt=m_batt, e_batt_wh=e_batt, c_rate=P_batt / e_batt if e_batt else 0.0,
                p_conv=p_conv, m_total=m_total, TW=T / (m_total * G) if m_total else 0.0,
                endurance=e_batt * a.batt_dod * 3600 / P_batt if P_batt else 0.0)


def best(a: Assumptions = Assumptions(), ds=(0.020, 0.025, 0.030, 0.040), ns=(2, 3, 4, 5, 6),
         Ws=(0.6, 0.8, 1.0), e_max=7.5e5, r_e=25e-6):
    """Grid search for the highest thrust-to-weight."""
    out = None
    for d in ds:
        cur = stage(d, r_e)
        for V in np.linspace(1.15 * cur.V_onset, min(e_max * d, cur.V.max()), 40):
            for n in ns:
                for W in Ws:
                    r = evaluate(d, V, n, W, a, r_e)
                    if out is None or r["TW"] > out["TW"]:
                        out = r
    return out


def thrust_map(a: Assumptions = Assumptions(), n_pts: int = 60) -> dict:
    """Total static thrust and electrical power of the frozen design versus stage voltage."""
    import design as D

    cur = stage(D.D_GAP, D.R_WIRE)
    V = np.linspace(cur.V_onset, 1.1 * D.V_MAX, n_pts)
    rows = [evaluate(D.D_GAP, v, D.N_STAGES, D.DUCT_W, a, D.R_WIRE, D.BATTERY) for v in V]
    return dict(V=V.tolist(), T=[r["T"] for r in rows], P=[r["P_el"] for r in rows],
                T_body=[r["T"] / r["aug"] for r in rows], v_duct=[r["v_duct"] for r in rows])


def design_report(a: Assumptions = Assumptions()) -> dict:
    """Evaluate the frozen design, its mass budget and its sensitivity to every assumption."""
    import design as D

    a = replace(a, m_avionics=D.M_AVIONICS, m_payload=D.M_PAYLOAD)
    ev = lambda aa, **kw: evaluate(D.D_GAP, kw.get("V", D.V_MAX), D.N_STAGES, D.DUCT_W, aa, D.R_WIRE, D.BATTERY)
    base = ev(a)
    sens = []
    for label, lo, hi in [
        ("EHD force derate k_ehd (0.8 / 1.0)", dict(k_ehd=0.8), dict(k_ehd=1.0)),
        ("fixed duct loss K0 (0.45 / 0.15)", dict(K0=0.45), dict(K0=0.15)),
        ("converter kW/kg (2.0 / 3.5)", dict(conv_kw_per_kg=2.0), dict(conv_kw_per_kg=3.5)),
        ("negative-stage current x (1.45 / 1.15)", dict(neg_I=1.45), dict(neg_I=1.15)),
        ("strip film g/m2 (8.4 / 3.6)", dict(film_strip=8.4e-3), dict(film_strip=3.6e-3)),
        ("frame tubes (x1.5 / x0.75)", dict(tube_g_per_m=3.3e-3, longeron_g_per_m=5.7e-3),
         dict(tube_g_per_m=1.65e-3, longeron_g_per_m=2.85e-3)),
        ("active area fraction (0.85 / 0.93)", dict(fill=0.85), dict(fill=0.93)),
        ("payload (50 g / 0 g)", dict(m_payload=0.050), dict(m_payload=0.0)),
    ]:
        sens.append(dict(label=label, lo=ev(replace(a, **lo))["TW"], hi=ev(replace(a, **hi))["TW"]))
    cons = replace(a, k_ehd=0.8, conv_kw_per_kg=2.0, K0=0.4)
    opt = replace(a, k_ehd=1.0, conv_kw_per_kg=3.5, K0=0.2)
    out = dict(assumptions=asdict(a), design=base, sensitivity=sens, thrust_map=thrust_map(a),
               scenarios=dict(conservative=ev(cons)["TW"], baseline=base["TW"], optimistic=ev(opt)["TW"]),
               cell=cell_geometry(D.D_GAP, D.R_WIRE), battery=D.BATTERY)
    # hover point: voltage at which thrust equals weight
    tm = out["thrust_map"]
    W_N = base["m_total"] * G
    out["hover"] = dict(V=float(np.interp(W_N, tm["T"], tm["V"])), P_el=float(np.interp(W_N, tm["T"], tm["P"])),
                        weight_N=W_N)
    out["hover"]["P_batt"] = out["hover"]["P_el"] / a.conv_eff
    out["hover"]["endurance_s"] = base["e_batt_wh"] * a.batt_dod * 3600 / out["hover"]["P_batt"]
    return out


def make_figures(rep: dict, outdir: pathlib.Path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import design as D

    b = rep["design"]
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.6))
    # (a) thrust / power curve
    a0 = ax[0, 0]
    tm = rep["thrust_map"]
    V = np.array(tm["V"]) / 1e3
    a0.plot(V, tm["T"], color="C3", lw=2, label="thrust (ducted)")
    a0.plot(V, tm["T_body"], color="C3", lw=1, ls="--", label="EHD body force only")
    a0.axhline(rep["hover"]["weight_N"], color="0.3", lw=1, ls=":")
    a0.text(V[0], rep["hover"]["weight_N"] * 1.03, f"weight {b['m_total']*1e3:.0f} g", fontsize=9)
    a0.axvline(D.V_MAX / 1e3, color="0.6", lw=0.8)
    a0.set(xlabel="stage voltage [kV]", ylabel="thrust [N]", title="Thrust and electrical power vs stage voltage")
    a1 = a0.twinx()
    a1.plot(V, tm["P"], color="C0", lw=1.6, label="electrical power")
    a1.set_ylabel("electrical power [W]")
    h1, l1 = a0.get_legend_handles_labels()
    h2, l2 = a1.get_legend_handles_labels()
    a0.legend(h1 + h2, l1 + l2, frameon=False, loc="upper left", fontsize=9)
    # (b) mass budget vs thrust
    a0 = ax[0, 1]
    parts = [("electrode strips", b["m_strips"]), ("frame", b["m_frame"]), ("duct skin + lip", b["m_skin"]),
             ("threads, wiring, feet", rep["assumptions"]["m_misc"]), ("HV converters (4)", b["m_conv"]),
             ("battery", b["m_batt"]), ("avionics", rep["assumptions"]["m_avionics"]),
             ("payload", rep["assumptions"]["m_payload"])]
    left = 0.0
    cols = ["#8c6d31", "#bd9e39", "#e7cb94", "#cedb9c", "#3182bd", "#e6550d", "#756bb1", "#31a354"]
    for (name, m), c in zip(parts, cols):
        a0.barh(0, m * 1e3, left=left, color=c, edgecolor="w", label=f"{name}: {m*1e3:.0f} g")
        left += m * 1e3
    a0.text(left + 8, 0, f"{left:.0f} g", va="center", fontsize=9)
    a0.barh(-1, b["T"] / G * 1e3, color="C3", edgecolor="w")
    a0.text(b["T"] / G * 500, -1, f"max thrust {b['T']/G*1e3:.0f} gf  (T/W = {b['TW']:.2f})", ha="center", va="center",
            color="w", fontsize=9)
    a0.set(yticks=[], xlabel="grams", title="Mass budget against maximum thrust", ylim=(-1.6, 3.6))
    a0.legend(frameon=False, fontsize=8, ncol=2, loc="upper left")
    # (c) sensitivity
    a0 = ax[1, 0]
    sens = sorted(rep["sensitivity"], key=lambda s: abs(s["hi"] - s["lo"]))
    for k, s in enumerate(sens):
        a0.barh(k, s["hi"] - s["lo"], left=s["lo"], color="C0", alpha=0.75)
        a0.text(s["lo"] - 0.01, k, f"{s['lo']:.2f}", ha="right", va="center", fontsize=8)
        a0.text(s["hi"] + 0.01, k, f"{s['hi']:.2f}", ha="left", va="center", fontsize=8)
    a0.set_yticks(range(len(sens)), [s["label"] for s in sens], fontsize=8)
    a0.axvline(b["TW"], color="k", lw=1)
    a0.axvline(1.0, color="C3", lw=1, ls="--")
    lo = min(s["lo"] for s in sens)
    hi = max(s["hi"] for s in sens)
    a0.set(xlim=(lo - 0.1, hi + 0.1), xlabel="thrust / weight at full throttle",
           title=f"Sensitivity (baseline {b['TW']:.2f}; lift-off needs > 1)")
    # (d) design space
    a0 = ax[1, 1]
    A = Assumptions(**rep["assumptions"])
    for n, c in zip((2, 3, 4, 6), ("C0", "C1", "C2", "C4")):
        Ws = np.linspace(0.5, 1.4, 19)
        tw = []
        for W in Ws:
            cur = stage(D.D_GAP, D.R_WIRE)
            tw.append(max(evaluate(D.D_GAP, v, n, W, A, D.R_WIRE)["TW"]
                          for v in np.linspace(1.2 * cur.V_onset, 22.5e3, 25)))
        a0.plot(Ws, tw, color=c, label=f"{n} stages")
    a0.plot([D.DUCT_W], [b["TW"]], "k*", ms=12, label="design point")
    a0.axhline(1.0, color="C3", lw=1, ls="--")
    a0.set(xlabel="duct side W [m]", ylabel="best thrust / weight", title="Scale and stage count (battery sized per case)")
    a0.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "sizing.png", dpi=140)


if __name__ == "__main__":
    import sys

    out = pathlib.Path(__file__).resolve().parent.parent / "results"
    rep = design_report()
    (out / "design.json").write_text(json.dumps(rep, indent=2))
    make_figures(rep, out)
    b = rep["design"]
    print(f"design: T={b['T']:.2f} N  P_el={b['P_el']:.0f} W  theta={b['theta']*1e3:.1f} N/kW  aug={b['aug']:.2f} "
          f"v_duct={b['v_duct']:.2f} m/s  mass={b['m_total']*1e3:.0f} g  T/W={b['TW']:.3f}")
    print(f"  masses [g]: strips {b['m_strips']*1e3:.0f} frame {b['m_frame']*1e3:.0f} skin {b['m_skin']*1e3:.0f} "
          f"conv {b['m_conv']*1e3:.0f} batt {b['m_batt']*1e3:.0f}; L={b['L']:.3f} m; I_stage={b['I_stage']*1e3:.2f} mA; "
          f"E_coll={b['E_coll']/1e6:.2f} MV/m; C-rate {b['c_rate']:.0f}")
    print("  hover:", {k: round(v, 3) for k, v in rep["hover"].items()})
    print("  scenarios:", {k: round(v, 3) for k, v in rep["scenarios"].items()})
    for s in rep["sensitivity"]:
        print(f"   {s['label']:45s} {s['lo']:.3f} .. {s['hi']:.3f}")
    if "--search" in sys.argv:
        cases = {
            "baseline": Assumptions(),
            "conservative": Assumptions(k_ehd=0.8, conv_kw_per_kg=2.0, K0=0.4, batt_max_c=35),
            "optimistic": Assumptions(k_ehd=1.0, conv_kw_per_kg=3.5, conv_eff=0.9, batt_max_c=60, hover_time=40),
        }
        for name, a in cases.items():
            for r_e in (25e-6, 12.5e-6):
                b = best(a, r_e=r_e)
                print(f"{name:12s} r_e={r_e*1e6:4.1f}um: T/W={b['TW']:.3f} d={b['d']*1e3:.0f}mm V={b['V']/1e3:.1f}kV n={b['n']} "
                      f"W={b['W']:.1f} T={b['T']:.2f}N P={b['P_el']:.0f}W theta={b['theta']*1e3:.1f}N/kW aug={b['aug']:.2f} "
                      f"m={b['m_total']*1e3:.0f}g [struct {b['m_struct']*1e3:.0f} conv {b['m_conv']*1e3:.0f} batt {b['m_batt']*1e3:.0f}] "
                      f"C={b['c_rate']:.0f} t={b['endurance']:.0f}s E={b['E_avg']/1e5:.1f}e5")
