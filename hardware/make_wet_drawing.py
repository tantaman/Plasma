"""
Concept sheet for the wet-emitter variant (docs/WET_EMITTER.md), drawn from
results/design_wet.json:

  hardware/05_wet_emitter_concept.{svg,png}

Run:  python hardware/make_wet_drawing.py     (after  python sim/wet.py)
This is a concept sheet, not a build drawing: the emitter and the liquid handling
have not been tested.
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle

import make_drawings as M

REP = json.loads((M.HERE.parent / "results" / "design_wet.json").read_text())
B, A, CELL = REP["design"], REP["assumptions"], REP["cell"]
D = CELL["d"] * 1e3
PITCH = CELL["pitch"] * 1e3
CC, CD, TC = CELL["c_c"] * 1e3, CELL["c_d"] * 1e3, CELL["t_c"] * 1e3
CHORD = CC + CD
S = D + CHORD
N = B["n"]
WET, PAD, DROP = "#1b7f8c", "#9fd3c7", "#e07b39"
RAIL_COL = ["#a61e1e", "#d98880", "#7fa6d9", "#1f4e96"]


def rail_kv(i: int) -> float:
    """Potential of grille i (0 = top) for monotonic rails, centred on the pod."""
    return (B["V"] * (N / 2 - i)) / 1e3


def wet_strip(ax, y, zn, col, emit=True, catch=False, lw=1.0):
    M.strip(ax, y, zn, col, wire=False, lw=lw)
    # absorbent wrap over the nose
    th = np.linspace(0, np.pi, 24)
    r = TC / 2 + 0.7
    ax.plot(y + r * np.cos(th), zn - TC / 2 + r * np.sin(th), color=PAD, lw=3.2 * lw, solid_capstyle="butt", zorder=3)
    for sgn in (-1, 1):
        ax.plot([y + sgn * r, y + sgn * r], [zn - TC / 2, zn - 8.0], color=PAD, lw=3.2 * lw, solid_capstyle="butt", zorder=3)
    if emit:
        ax.add_patch(Circle((y, zn - CHORD), 0.55, fc=WET, ec=WET, zorder=5))
    if catch:
        ax.add_patch(Rectangle((y - 1.6, zn - CHORD - 1.2), 3.2, 2.4, fc=PAD, ec=WET, lw=0.6, zorder=5))


def spray(ax, y, z0, z1, rng, n=70):
    """Droplets fanning from an emitter at (y, z0) to the collector row whose noses are at z1."""
    for _ in range(n):
        s = rng.random() ** 0.8
        spread = (0.5 * PITCH) * (0.15 + 0.85 * s) * rng.uniform(-1, 1) * s**0.6
        ax.add_patch(Circle((y + spread, z0 - s * (z0 - z1) - rng.uniform(0, 2.0)), 0.22, fc=DROP, ec="none", alpha=0.75, zorder=2))


def concept():
    rng = np.random.default_rng(3)
    fig = plt.figure(figsize=(16, 10))
    # ------------------------------------------------------------ left: three columns, all grilles
    ax = fig.add_axes([0.03, 0.05, 0.52, 0.87])
    ys = (-PITCH, 0.0, PITCH)
    z_top = 0.0
    for i in range(N + 1):
        zn = z_top - i * S
        for y in ys:
            wet_strip(ax, y, zn, RAIL_COL[i], emit=i < N, catch=i == N)
            if i < N:
                spray(ax, y, zn - CHORD - 0.8, zn - S + 0.5, rng)
        ax.text(ys[0] - 44, zn - CHORD / 2, f"G{i+1}\n{rail_kv(i):+.0f} kV", ha="center", va="center", fontsize=11,
                color=RAIL_COL[i], weight="bold")
    # air flow
    for y in (-PITCH / 2, PITCH / 2):
        ax.add_patch(FancyArrowPatch((y, z_top + 22), (y, z_top + 6), arrowstyle="-|>", mutation_scale=14, color=M.GREY, lw=1.2))
    ax.text(0, z_top + 26, f"air, {B['v_duct']:.1f} m/s at full throttle", ha="center", fontsize=9.5, color=M.DARK)
    zb = z_top - N * S - CHORD
    for y in (-PITCH / 2, PITCH / 2):
        ax.add_patch(FancyArrowPatch((y, zb - 6), (y, zb - 22), arrowstyle="-|>", mutation_scale=14, color=M.GREY, lw=1.2))
    # dimensions
    M.dim(ax, (ys[0] - 4, z_top - CHORD), (ys[0] - 4, z_top - S), f"{D:.0f} gap", off=-12, horiz=False)
    M.dim(ax, (ys[0] - 4, z_top), (ys[0] - 4, z_top - CHORD), f"{CHORD:.0f} strip", off=-12, horiz=False)
    M.dim(ax, (ys[0], z_top + 3), (ys[1], z_top + 3), f"{PITCH:.0f} pitch", off=8)
    # notes
    M.note(ax, (ys[1] + 5.5, z_top - CHORD - 12), (ys[2] + 30, z_top - CHORD - 10),
           f"positive droplets, 0.1–0.3 µm,\ndrift ≈ {A['Z']*B['V']/CELL['d']:.0f} m/s\n(air ions: ≈ 130 m/s)", fs=9)
    M.note(ax, (ys[2] + 0.6, z_top - S - CHORD), (ys[2] + 30, z_top - S - CHORD + 4),
           "wetted emitter on the trailing edge,\nat the potential of its own strip", fs=9)
    M.note(ax, (ys[2] + 1.6, z_top - 2 * S - 0.3), (ys[2] + 30, z_top - 2 * S + 6),
           "absorbent wrap on every nose:\nwhat lands is soaked up, no beads\nleft standing in the field", fs=9)
    M.note(ax, (ys[2] + 1.8, zb + 0.2), (ys[2] + 30, zb - 2),
           "last grille: catch pads, no emitters;\na few bare points below it bleed off\nthe charge of droplets that escape", fs=9)
    ax.text(ys[0] - 44, zb - 20, f"{N} stages × {B['V']/1e3:.0f} kV,\n{B['V_rr']/1e3:.0f} kV end to end.\nEvery stage sprays\npositive.", ha="center",
            va="top", fontsize=9.5, color=M.DARK)
    ax.set_xlim(ys[0] - 62, ys[2] + 95)
    ax.set_ylim(zb - 40, z_top + 34)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Section through three columns of strips (flow downwards)", fontsize=11, loc="left")

    # ------------------------------------------------------------ right: one strip, enlarged
    ax = fig.add_axes([0.58, 0.44, 0.40, 0.48])
    k = 1.0
    wet_strip(ax, 0.0, 0.0, RAIL_COL[1], lw=1.6)
    # feeder thread and emitter teeth, drawn schematically beside the section
    ax.plot([0.35, 0.35], [-CC + 1.0, -CHORD], color=WET, lw=1.4, zorder=4)
    x0 = 9.0
    ax.add_patch(Rectangle((x0, -CHORD + 1.2), 22.0, CD - 1.2, fc="#eaf3fa", ec="#7fb2d9", lw=0.8))
    ax.add_patch(Rectangle((x0, -CC), 22.0, CC, fc="#f2dede", ec=RAIL_COL[1], lw=0.8))
    ax.add_patch(Rectangle((x0, -8.0), 22.0, 8.0, fc=PAD, ec=WET, lw=0.6, alpha=0.9))
    teeth = np.arange(x0, x0 + 22.0 + 1e-6, 2.0)
    pts = [(x0, -CHORD + 1.2)]
    for xa, xb in zip(teeth[:-1], teeth[1:]):
        pts += [((xa + xb) / 2, -CHORD - 0.4), (xb, -CHORD + 1.2)]
    ax.add_patch(Polygon(pts + [(x0 + 22.0, -CHORD + 2.4), (x0, -CHORD + 2.4)], closed=True, fc=WET, ec=WET, lw=0.5))
    for xf in (x0 + 4.0, x0 + 18.0):
        ax.plot([xf, xf], [-CC + 1.0, -CHORD + 2.2], color=WET, lw=1.2)
    for xa, xb in zip(teeth[:-1], teeth[1:]):
        xm = (xa + xb) / 2
        for j in range(5):
            ax.add_patch(Circle((xm + rng.uniform(-0.5, 0.5) * (j + 1) * 0.35, -CHORD - 1.0 - 0.9 * j), 0.13, fc=DROP, ec="none"))
    ax.text(x0 + 11.0, 2.6, "same strip seen face-on (22 mm of span)", ha="center", fontsize=8.5, color=M.GREY)
    M.note(ax, (x0 + 11.0, -4.0), (x0 + 28.0, -3.0), "nose wrap: 8 mm of lens tissue\nover the metallised film", fs=9)
    M.note(ax, (x0 + 15.0, -CC / 2 - 5.0), (x0 + 28.0, -11.5), "collector: metallised film,\nas in the dry strip", fs=9)
    M.note(ax, (x0 + 18.0, -CC - CD / 2), (x0 + 28.0, -22.0), "tail: clear film, kept dry.\nFeeder threads every ≈ 15 mm\ncarry liquid and the strip's potential\nto the emitter", fs=9)
    M.note(ax, (x0 + 20.0, -CHORD + 1.0), (x0 + 28.0, -33.5),
           f"emitter: porous ribbon cut to teeth at 2 mm,\ntip radius ≈ {A['r_fibre']*1e6:.0f} µm, pre-wetted.\n"
           f"Needs {B['uA_per_m']:.0f} µA per metre: ≈ {B['uA_per_m']*2:.0f} nA per tooth", fs=9)
    M.dim(ax, (-3.0, 0.0), (-3.0, -CC), f"{CC:.0f}", off=-5, horiz=False)
    M.dim(ax, (-3.0, -CC), (-3.0, -CHORD), f"{CD:.0f}", off=-5, horiz=False)
    ax.set_xlim(-14, x0 + 78)
    ax.set_ylim(-CHORD - 9, 6)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("One integrated strip with wet emitter (section and face view)", fontsize=11, loc="left")

    # ------------------------------------------------------------ numbers
    fz = REP["frozen_baseline"]
    txt = (f"Model, full throttle ({B['V']/1e3:.0f} kV per stage)        dry baseline      wet emitters\n"
           f"thrust                                     {fz['T']:5.2f} N         {B['T']:5.2f} N\n"
           f"electrical power                           {fz['P_el']:5.0f} W         {B['P_el']:5.0f} W\n"
           f"thrust per electrical power                {fz['T']/fz['P_el']*1e3:5.1f} N/kW      {B['theta']*1e3:5.1f} N/kW\n"
           f"all-up mass / payload                    {fz['m_total']*1e3:4.0f} / {fz['payload']*1e3:3.0f} g     {B['m_total']*1e3:4.0f} / {A['m_payload']*1e3:3.0f} g\n"
           f"liquid sprayed                                  —           {B['mdot_liq']*1e6:5.0f} mg/s")
    fig.text(0.585, 0.37, txt, fontsize=9.3, family="monospace", va="top", color=M.DARK)
    fig.text(0.585, 0.20, "Liquid: non-volatile, surface tension ≤ 0.045 N/m, conductivity 0.03–0.3 S/m\n"
             "(tetraethylene glycol with a dissolved salt, or an ionic liquid diluted in it).\n"
             "Unproven: emission at this density from a passive emitter; behaviour over minutes as the collectors load.",
             fontsize=9, va="top", color=M.DARK)
    fig.text(0.5, 0.985, "Ion Lantern W — wet-emitter concept", ha="center", va="top", fontsize=14, weight="bold")
    fig.text(0.995, 0.005, "Ion Lantern W  |  concept sheet  |  dimensions in mm  |  generated from results/design_wet.json",
             ha="right", va="bottom", fontsize=7, color=M.GREY)
    M.save(fig, "05_wet_emitter_concept")


if __name__ == "__main__":
    concept()
    print("concept sheet written to", M.HERE)
