"""
Validation of the corona solver (sim/corona.py) against

  A. Townsend's exact solution for a coaxial wire-cylinder corona.
  B. Moreau et al., J. Phys. D 46 (2013) 475204: 25 um tungsten wire to a 12 mm
     cylinder, single pair in open air (current-voltage fits from their Table I,
     thrust from Table III).
  C. The single-stage ducted thruster in US patent 12,404,844 (Gomez-Vega, Brown,
     Xu, Barrett; MIT): 20 wire/plate pairs x 0.6 m, 50 um (44 AWG) tungsten wire,
     flat-plate collectors 25.4 x 0.8 mm, gap = pitch = 20 mm, 600 x 400 mm duct.
     Values were read from the patent figures (so +-5 %).
  D. Mesh convergence on the design cell.

Run:  python sim/validate_corona.py     -> results/corona_validation.{png,json}
"""
from __future__ import annotations

import json
import pathlib
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from corona import MU_ION, CoronaProblem, coaxial_analytic, mesh_cell, mesh_coaxial, peek_field, stage_curve

OUT = pathlib.Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(exist_ok=True)


def main():
    res = {}
    fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.5))

    # ---------------------------------------------------------------- A. coaxial
    r0, R = 50e-6, 0.03
    Ec = peek_field(r0)
    c = stage_curve(CoronaProblem(mesh_coaxial(r0, R, n_wire=16, h_max_frac=0.015)), Ec, v_max_over_onset=6)
    Vex = coaxial_analytic(r0, R, Ec, c.I)
    err = np.abs(c.V - Vex) / Vex
    res["A_coaxial"] = dict(max_rel_error=float(err.max()), n_points=len(c.V))
    a = ax[0, 0]
    a.plot(Vex / 1e3, c.I * 1e3, "-", color="0.4", label="exact (Townsend)")
    a.plot(c.V / 1e3, c.I * 1e3, "o", ms=4, color="C0", label="this solver")
    a.set(xlabel="voltage [kV]", ylabel="current [mA/m]",
          title=f"A. Coaxial corona, 100 µm wire in 60 mm tube\nmax error {100*err.max():.1f} %")
    a.legend(frameon=False)

    # ---------------------------------------------------------------- B. Moreau 2013
    fits = {20: (4.2e-12, 6.1e3), 30: (2.27e-12, 7.5e3), 40: (1.2e-12, 8.8e3)}  # C [A/V^2/m], V0 [V]
    vmax = {20: 17e3, 30: 22e3, 40: 25e3}  # upper end of the measured range
    r_e = 12.5e-6
    Ec = peek_field(r_e)
    a = ax[0, 1]
    res["B_moreau"] = {}
    for k, (dmm, (C, V0)) in enumerate(fits.items()):
        d = dmm * 1e-3
        m = mesh_cell(d=d, pitch=0.6, r_e=r_e, t_c=0.012, c_c=0.012, mode="isolated",
                      upstream=0.4, downstream=0.3, h_max=0.012, h_coll=0.0008)
        cur = stage_curve(CoronaProblem(m), Ec, n_pts=18, v_max_over_onset=9)
        V = np.linspace(V0, vmax[dmm], 40)
        a.plot(V / 1e3, C * V * (V - V0) * 1e3, "-", color=f"C{k}", lw=1.2, label=f"measured fit, d = {dmm} mm")
        Vc = min(vmax[dmm], cur.V.max())  # never extrapolate the model beyond its last solved point
        Vm = np.linspace(cur.V_onset, Vc, 40)
        Im, Fm, _ = cur.at(Vm)
        a.plot(Vm / 1e3, Im * 1e3, "--", color=f"C{k}", lw=1.8)
        res["B_moreau"][f"d{dmm}"] = dict(
            V_onset_model_kV=cur.V_onset / 1e3, V_onset_meas_kV=V0 / 1e3,
            I_model_mA_m=float(cur.at(Vc)[0] * 1e3), I_meas_mA_m=float(C * Vc * (Vc - V0) * 1e3), at_kV=Vc / 1e3,
            F_per_I_model_N_per_A=float(cur.at(Vc)[1] / cur.at(Vc)[0]), F_per_I_meas_N_per_A=33.0 * dmm / 10)
    a.plot([], [], "--", color="0.3", lw=1.8, label="this solver (dashed)")
    a.set(xlabel="voltage [kV]", ylabel="current [mA/m]",
          title="B. Moreau 2013: 25 µm wire → 12 mm cylinder")
    a.legend(frameon=False, fontsize=8)

    # ---------------------------------------------------------------- C. MIT ducted single stage
    meas_V = np.array([10e3, 15e3, 19e3])
    meas_I = np.array([1.0e-3, 3.9e-3, 7.5e-3])  # A, 12 m of electrode
    meas_T = np.array([0.15, 0.59, 1.14])  # N, ducted
    span = 20 * 0.6
    r_e = 25e-6
    m = mesh_cell(d=0.02, pitch=0.02, r_e=r_e, t_c=0.0008, c_c=0.0254, mode="isolated", upstream=0.16)
    cur = stage_curve(CoronaProblem(m), peek_field(r_e), v_max_over_onset=4)
    V = np.linspace(cur.V_onset, 20e3, 60)
    I, F, _ = cur.at(V)
    Im, Fm, _ = cur.at(meas_V)
    aug = meas_T / (Fm * span)
    res["C_mit_single_stage"] = dict(
        V_kV=(meas_V / 1e3).tolist(), I_meas_mA=(meas_I * 1e3).tolist(), I_model_mA=(Im * span * 1e3).tolist(),
        T_meas_N=meas_T.tolist(), F_ehd_model_N=(Fm * span).tolist(), thrust_over_model_body_force=aug.tolist(),
        note="thrust/body-force is the duct augmentation that duct.py is calibrated to")
    a = ax[1, 0]
    a.plot(V / 1e3, I * span * 1e3, "-", color="C0", label="current, this solver")
    a.plot(meas_V / 1e3, meas_I * 1e3, "s", color="C0", mfc="w", label="current, measured")
    a.set(xlabel="voltage [kV]", ylabel="current [mA]",
          title="C. MIT ducted stage (patent US 12,404,844)\n20 pairs × 0.6 m, gap = pitch = 20 mm")
    b = a.twinx()
    b.plot(V / 1e3, F * span, "-", color="C3", label="EHD body force, this solver")
    b.plot(V / 1e3, F * span * aug.mean(), ":", color="C3", label=f"× {aug.mean():.2f} (duct)")
    b.plot(meas_V / 1e3, meas_T, "o", color="C3", mfc="w", label="ducted thrust, measured")
    b.set_ylabel("force [N]")
    h1, l1 = a.get_legend_handles_labels()
    h2, l2 = b.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="upper left")

    # ---------------------------------------------------------------- D. mesh convergence
    from design import STAGE_V, cell  # design point

    CELL = cell()

    a = ax[1, 1]
    res["D_mesh"] = []
    for n_wire, div in ((12, 30), (16, 45), (22, 70), (32, 110)):
        t0 = time.time()
        m = mesh_cell(**CELL, n_wire=n_wire, h_max=CELL["pitch"] / div, h_coll=CELL["t_c"] / (0.12 * div))
        pr = CoronaProblem(m)
        cur = stage_curve(pr, peek_field(CELL["r_e"]), v_max_over_onset=4)
        I, F, _ = cur.at(STAGE_V)
        res["D_mesh"].append(dict(nodes=pr.n, I_mA_m=float(I * 1e3), F_mN_m=float(F * 1e3),
                                  V_onset_kV=cur.V_onset / 1e3, seconds=time.time() - t0))
        a.plot(pr.n, F * 1e3, "o", color="C2")
        a.annotate(f"{I*1e3:.3f} mA/m", (pr.n, F * 1e3), textcoords="offset points", xytext=(6, -12), fontsize=8)
    a.set(xscale="log", xlabel="mesh nodes", ylabel="EHD force [mN/m]",
          title=f"D. Mesh convergence, design cell at {STAGE_V/1e3:.1f} kV")
    Fs = [r["F_mN_m"] for r in res["D_mesh"]]
    a.set_ylim(0.9 * min(Fs), 1.1 * max(Fs))
    fig.tight_layout()
    fig.savefig(OUT / "corona_validation.png", dpi=140)
    (OUT / "corona_validation.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
