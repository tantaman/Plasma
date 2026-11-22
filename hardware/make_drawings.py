"""
Generates the build drawings from the frozen design (results/design.json):

  hardware/01_general_arrangement.{svg,png}
  hardware/02_electrode_detail.{svg,png}
  hardware/03_hv_module_schematic.{svg,png}
  hardware/04_system_wiring.{svg,png}

Run:  python hardware/make_drawings.py
All dimensions on the drawings are millimetres unless noted.
"""
from __future__ import annotations

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle

HERE = pathlib.Path(__file__).resolve().parent
REP = json.loads((HERE.parent / "results" / "design.json").read_text())
B, CELL = REP["design"], REP["cell"]

RED, BLUE, DARK, FILM, FOAM, GREY = "#c0392b", "#2c5fa8", "#222222", "#bcd9ee", "#e8d9a8", "#8a8a8a"
W = B["W"] * 1e3
L = round(B["L"] * 1e3)
N = B["n"]
D = CELL["d"] * 1e3
PITCH = CELL["pitch"] * 1e3
CC, CD, TC = CELL["c_c"] * 1e3, CELL["c_d"] * 1e3, CELL["t_c"] * 1e3
CHORD = CC + CD
S = D + CHORD
LIP = 35.0
LEG = 120.0
TOP_MARGIN = 60.0
ROWS = int(round(0.96 * W / PITCH))
Y0 = -(ROWS - 1) / 2 * PITCH
HALF = 0.47 * W
VHALF = B["V"] / 2e3  # kV


def save(fig, name):
    for ext in ("svg", "png"):
        fig.savefig(HERE / f"{name}.{ext}", dpi=170, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def dim(ax, p0, p1, text, off=0.0, horiz=True, fs=8, tick=6.0):
    """Dimension line between two points, offset perpendicular by `off` (mm)."""
    (x0, y0), (x1, y1) = p0, p1
    if horiz:
        y = y0 + off
        ax.plot([x0, x0], [y0, y], color=GREY, lw=0.5)
        ax.plot([x1, x1], [y1, y], color=GREY, lw=0.5)
        ax.annotate("", (x1, y), (x0, y), arrowprops=dict(arrowstyle="<->", color="k", lw=0.7, shrinkA=0, shrinkB=0))
        ax.text((x0 + x1) / 2, y + tick * 0.6, text, ha="center", va="bottom", fontsize=fs)
    else:
        x = x0 + off
        ax.plot([x0, x], [y0, y0], color=GREY, lw=0.5)
        ax.plot([x1, x], [y1, y1], color=GREY, lw=0.5)
        ax.annotate("", (x, y1), (x, y0), arrowprops=dict(arrowstyle="<->", color="k", lw=0.7, shrinkA=0, shrinkB=0))
        ax.text(x - tick * 0.6, (y0 + y1) / 2, text, ha="right", va="center", fontsize=fs, rotation=90)


def note(ax, xy, xytext, text, fs=8, ha="left"):
    ax.annotate(text, xy, xytext, fontsize=fs, ha=ha, va="center",
                arrowprops=dict(arrowstyle="-", color="k", lw=0.6, shrinkA=2, shrinkB=1))


def title_block(fig, title, sheet):
    fig.text(0.5, 0.985, title, ha="center", va="top", fontsize=14, weight="bold")
    fig.text(0.995, 0.005, f"Ion Lantern  |  sheet {sheet}/4  |  dimensions in mm  |  generated from results/design.json",
             ha="right", va="bottom", fontsize=7, color=GREY)


# ------------------------------------------------------------------ 1. general arrangement
def general_arrangement():
    fig = plt.figure(figsize=(17, 9.2))
    ax = fig.add_axes([0.03, 0.05, 0.44, 0.86])
    h = W / 2
    # lip and duct
    ax.add_patch(FancyBboxPatch((-h - LIP, -h - LIP), W + 2 * LIP, W + 2 * LIP, boxstyle=f"round,pad=0,rounding_size={LIP}",
                                fc="#f3f1e7", ec=DARK, lw=1.0))
    ax.add_patch(Rectangle((-h, -h), W, W, fc="white", ec=DARK, lw=1.2))
    zc = {"A": (1, 1, "#fdecea"), "B": (-1, 1, "#eaf1fb"), "C": (-1, -1, "#fdecea"), "D": (1, -1, "#eaf1fb")}
    for nm, (sx, sy, col) in zc.items():
        ax.add_patch(Rectangle((min(0, sx * h), min(0, sy * h)), h, h, fc=col, ec="none", zorder=0.5))
        ax.text(sx * h * 0.5, sy * h * 0.5, f"ZONE {nm}", fontsize=15, color="#555555", ha="center", va="center",
                weight="bold", alpha=0.55, zorder=3)
    for i in range(ROWS):
        y = Y0 + i * PITCH
        for sx in (-1, 1):
            ax.plot([sx * 10, sx * (10 + HALF)], [y, y], color="#b06a2c", lw=0.9, zorder=1)
    for x in (-490, -245, 0, 245, 490):
        ax.plot([x, x], [-h, h], color=DARK, lw=1.4 if x == 0 else 1.0, ls=(0, (6, 3)), zorder=2)
    for sx in (-1, 1):
        for sy in (-1, 1):
            ax.add_patch(Circle((sx * h, sy * h), 7, fc=DARK, ec=DARK, zorder=4))
    ax.add_patch(Rectangle((-55, -45), 110, 90, fc="#2b2b30", ec="k", zorder=5))
    ax.text(0, 0, "POD", color="w", ha="center", va="center", fontsize=8, zorder=6)
    dim(ax, (-h, -h - LIP), (h, -h - LIP), f"{W:.0f} inside", off=-45)
    dim(ax, (-h - LIP, -h), (-h - LIP, h), f"{W + 2*LIP:.0f} over lip" if False else f"{W:.0f} inside", off=-45, horiz=False)
    dim(ax, (10, h + LIP), (10 + HALF, h + LIP), f"strip {HALF:.0f}", off=22)
    note(ax, (245, 330), (300, 610), "support rails: 1.5 mm carbon tube, 5 per level\n(top and bottom), strips hang from them", fs=8)
    note(ax, (-300, Y0 + 3 * PITCH), (-520, -640), f"{ROWS} strip rows per grille at {PITCH:.0f} pitch,\neach row split in two half-strips (one per zone)", fs=8)
    note(ax, (h, h), (380, -640), "corner longeron, 2 mm carbon tube,\nextends 120 below as landing leg", fs=8)
    note(ax, (-h - LIP * 0.5, 200), (-640, 640), f"inlet lip: R{LIP:.0f} half-round,\n6 µm film over foam formers", fs=8)
    ax.text(0, h + LIP + 75, "TOP VIEW  (looking into the inlet)", ha="center", fontsize=11, weight="bold")
    ax.set_xlim(-700, 700)
    ax.set_ylim(-700, 720)
    ax.set_aspect("equal")
    ax.axis("off")

    # ---- section
    ax = fig.add_axes([0.50, 0.05, 0.49, 0.86])
    ax.plot([-h, -h], [0, L], color="#3a8fd0", lw=1.6)
    ax.plot([h, h], [0, L], color="#3a8fd0", lw=1.6)
    for sx in (-1, 1):
        ax.add_patch(Circle((sx * (h + LIP), L), LIP, fc="#f3f1e7", ec=DARK, lw=1.0))
        ax.plot([sx * h, sx * h], [-LEG, L], color=DARK, lw=1.0, ls=(0, (1, 0)), zorder=0.4, alpha=0.35)
        ax.plot([sx * h, sx * h], [-LEG, 0], color=DARK, lw=1.6)
        ax.add_patch(Circle((sx * h, -LEG), 6, fc="#e67e22", ec="k"))
    ax.plot([-h, h], [0, 0], color=DARK, lw=1.2, ls=(0, (6, 3)))
    ax.plot([-h, h], [L, L], color=DARK, lw=1.2, ls=(0, (6, 3)))
    z_nose = [L - TOP_MARGIN - i * S for i in range(N + 1)]
    for i, zn in enumerate(z_nose):
        col = RED if i % 2 == 0 else BLUE
        for r in range(ROWS):
            y = Y0 + r * PITCH
            ax.plot([y, y], [zn, zn - CHORD], color=col, lw=1.1)
            if i < N:
                ax.plot([y], [zn - CHORD], marker="o", ms=1.6, color="k")
        sign = "+" if i % 2 == 0 else "−"
        ax.text(-h - 18, zn - CHORD / 2, f"G{i+1}  {sign}{VHALF:.0f} kV", color=col, fontsize=9, va="center", ha="right", weight="bold")
    for y in (-490, -245, 245, 490):
        ax.plot([y, y], [8, z_nose[0] + 6], color=GREY, lw=0.5)
    ax.plot([0, 0], [-45, L + 50], color=DARK, lw=1.8)
    ax.add_patch(Rectangle((-45, L - 15), 90, 60, fc="#2b2b30", ec="k", zorder=5))
    ax.text(0, L + 15, "POD", color="w", ha="center", va="center", fontsize=8, zorder=6)
    ax.add_patch(Rectangle((-20, -50), 40, 30, fc="#27ae60", ec="k", zorder=5))
    note(ax, (20, -35), (130, -75), f"payload, {REP['assumptions']['m_payload']*1e3:.0f} g", fs=8)
    note(ax, (45, L + 30), (150, L + 105), "pod: battery, 4 HV modules, flight controller\n(45 mm clear of the +10 kV strips, film shield under)", fs=8)
    # flow arrows
    for y in (-380, -190, 190, 380):
        ax.add_patch(FancyArrowPatch((y, L + 95), (y, L + 20), arrowstyle="-|>", mutation_scale=12, color="#4a90c2", lw=1.2))
        ax.add_patch(FancyArrowPatch((y, -15), (y, -95), arrowstyle="-|>", mutation_scale=12, color="#4a90c2", lw=1.2))
    for sx in (-1, 1):
        ax.add_patch(FancyArrowPatch((sx * (h + 150), L + 60), (sx * (h - 45), L - 20), connectionstyle=f"arc3,rad={-0.45*sx}",
                                     arrowstyle="-|>", mutation_scale=12, color="#4a90c2", lw=1.2))
    ax.text(-h - 150, L + 120, "air in (drawn round the lip)", color="#2f6f9f", fontsize=8.5, ha="left")
    ax.text(-380, -112, f"air out, about {B['v_duct']:.1f} m/s", color="#2f6f9f", fontsize=8.5, ha="left")
    ax.add_patch(FancyArrowPatch((-h - 215, 60), (-h - 215, 230), arrowstyle="-|>", mutation_scale=18, color="#c0392b", lw=2.2))
    ax.text(-h - 228, 145, f"thrust\n{B['T']:.1f} N max", color=RED, fontsize=9, ha="right", va="center", weight="bold")
    dim(ax, (-h, -LEG), (h, -LEG), f"{W:.0f}", off=-40)
    dim(ax, (h, 0), (h, L), f"{L}", off=190, horiz=False)
    dim(ax, (h, -LEG), (h, 0), f"{LEG:.0f}", off=190, horiz=False)
    dim(ax, (h, z_nose[1]), (h, z_nose[0]), f"stage pitch {S:.0f}", off=95, horiz=False, fs=7.5)
    ax.text(0, L + 160, "SECTION  (cut across the strips; strips seen end-on)", ha="center", fontsize=11, weight="bold")
    ax.text(-h, -200, f"All-up mass {B['m_total']*1e3:.0f} g  ·  max thrust {B['T']/9.80665*1e3:.0f} gf  ·  "
            f"hover {REP['hover']['V']/1e3:.1f} kV, {REP['hover']['P_batt']:.0f} W from a {REP['battery']['name']}",
            fontsize=9)
    ax.set_xlim(-860, 760)
    ax.set_ylim(-230, L + 190)
    ax.set_aspect("equal")
    ax.axis("off")
    title_block(fig, "Ion Lantern — general arrangement", 1)
    save(fig, "01_general_arrangement")


# ------------------------------------------------------------------ 2. electrode detail
def strip(ax, y, zn, col, wire=True, scale=1.0, lw=1.0, labels=False):
    """One integrated strip: foam nose stick, metallised film collector, bare tail, emitter wire."""
    t = TC
    ax.add_patch(Rectangle((y - t / 2, zn - 3.0), t, 3.0, fc=FOAM, ec="k", lw=0.4))
    ax.add_patch(Circle((y, zn - t / 2), t / 2, fc=FOAM, ec="k", lw=0.4))
    # metallised film wrapped over the nose (both faces) down to CC
    ax.plot([y - t / 2 - 0.15, y - t / 2 - 0.15, y - 0.05], [zn - t / 2, zn - 5, zn - CC], color=col, lw=lw * 1.6, solid_capstyle="round")
    ax.plot([y + t / 2 + 0.15, y + t / 2 + 0.15, y + 0.05], [zn - t / 2, zn - 5, zn - CC], color=col, lw=lw * 1.6, solid_capstyle="round")
    th = np.linspace(0, np.pi, 20)
    ax.plot(y + (t / 2 + 0.15) * np.cos(th), zn - t / 2 + (t / 2 + 0.15) * np.sin(th), color=col, lw=lw * 1.6)
    ax.plot([y, y], [zn - CC, zn - CHORD], color="#7fb2d9", lw=lw * 1.1)
    if wire:
        ax.add_patch(Circle((y, zn - CHORD), 0.35 * scale, fc="k", ec="k"))


def electrode_detail():
    fig = plt.figure(figsize=(16, 10))
    ax = fig.add_axes([0.02, 0.04, 0.50, 0.89])
    z_nose = [L - TOP_MARGIN - i * S for i in range(N + 1)]
    ys = (-PITCH, 0.0, PITCH)
    for i, zn in enumerate(z_nose):
        col = RED if i % 2 == 0 else BLUE
        for y in ys:
            strip(ax, y, zn, col, wire=i < N)
        if i < N:
            pos = i % 2 == 0
            c = "#f5c6c0" if pos else "#c5d8f3"
            for y in ys:  # ion plume from the wire to the next strip
                ax.add_patch(Polygon([(y, zn - CHORD), (y - PITCH * 0.47, zn - CHORD - D * 0.9), (y - 2, zn - S - CC),
                                      (y + 2, zn - S - CC), (y + PITCH * 0.47, zn - CHORD - D * 0.9)], fc=c, ec="none", alpha=0.55, zorder=0))
            ax.text(PITCH * 1.75, zn - CHORD - D / 2, ("positive corona\n+ ions drift down" if pos else "negative corona\n− ions drift down"),
                    fontsize=9, color=RED if pos else BLUE, va="center")
            ax.annotate("", (-PITCH * 1.72, zn - S + 1), (-PITCH * 1.72, zn - CHORD - 1),
                        arrowprops=dict(arrowstyle="-|>" if pos else "<|-", color="k", lw=1.0))
            ax.text(-PITCH * 1.8, zn - CHORD - D / 2, "E", fontsize=10, style="italic", ha="right", va="center")
        sign = "+" if i % 2 == 0 else "−"
        ax.text(PITCH * 1.75, zn - 8, f"grille G{i+1}:  {sign}{VHALF:.0f} kV rail", fontsize=9.5, color=col, weight="bold", va="center")
    for y in (-PITCH / 2, PITCH / 2):
        ax.add_patch(FancyArrowPatch((y, z_nose[0] + 34), (y, z_nose[0] + 6), arrowstyle="-|>", mutation_scale=11, color="#4a90c2", lw=1.1))
        ax.add_patch(FancyArrowPatch((y, z_nose[-1] - CHORD - 6), (y, z_nose[-1] - CHORD - 34), arrowstyle="-|>", mutation_scale=11, color="#4a90c2", lw=1.1))
    ax.text(0, z_nose[0] + 40, "air flow", color="#2f6f9f", ha="center", fontsize=9)
    dim(ax, (-PITCH, z_nose[0] + 3), (0, z_nose[0] + 3), f"{PITCH:.0f}", off=10, tick=3)
    dim(ax, (-PITCH, z_nose[0]), (-PITCH, z_nose[0] - CHORD), f"{CHORD:.0f}", off=-14, horiz=False, tick=3)
    dim(ax, (-PITCH, z_nose[0] - CHORD), (-PITCH, z_nose[1]), f"gap {D:.0f}", off=-14, horiz=False, tick=3)
    dim(ax, (-PITCH, z_nose[1]), (-PITCH, z_nose[2]), f"stage pitch {S:.0f}", off=-30, horiz=False, tick=3)
    ax.text(-PITCH * 2.6, z_nose[-1] - CHORD - 42,
            "Each strip hangs from the one above on two fine aramid threads per support rail.\n"
            "Neighbouring grilles sit on opposite rails, so every stage sees the full rail-to-rail\n"
            "voltage, and each emitter wire is at the potential of the strip it is bonded to —\n"
            "there is nothing upstream for it to discharge back to.",
            fontsize=8.5, va="top")
    ax.set_xlim(-PITCH * 2.7, PITCH * 4.6)
    ax.set_ylim(z_nose[-1] - CHORD - 88, z_nose[0] + 50)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("STACK  (3 strip columns of 32 shown, true proportions)", fontsize=11, weight="bold")

    # enlarged single strip
    ax = fig.add_axes([0.55, 0.30, 0.20, 0.63])
    strip(ax, 0, 0, RED, wire=True, scale=0.25, lw=1.8)
    note(ax, (0, -TC / 2), (6, 3), f"nose: foam stick {TC:.1f} × 3,\nedge sanded round (R ≥ {TC/2:.2f})", fs=8.5)
    note(ax, (TC / 2 + 0.2, -9), (6, -7), f"collector: 6 µm metallised film,\nmetal side out, wrapped over the\nnose, {CC:.0f} long on both faces", fs=8.5)
    note(ax, (0.1, -CC - CD / 2), (6, -22), f"tail: clear 6 µm film (no metal), bonded\nbetween the two collector faces, {CD:.0f} long", fs=8.5)
    note(ax, (0.2, -CHORD), (6, -CHORD - 1), f"emitter: Ø{2*CELL['r_e']*1e6:.0f} µm tungsten wire\nbonded along the trailing edge,\njumpered to this strip's collector", fs=8.5)
    dim(ax, (-TC / 2 - 0.2, 0), (-TC / 2 - 0.2, -CC), f"{CC:.0f}", off=-4, horiz=False, tick=1.0)
    dim(ax, (-TC / 2 - 0.2, -CC), (-TC / 2 - 0.2, -CHORD), f"{CD:.0f}", off=-4, horiz=False, tick=1.0)
    ax.set_xlim(-9, 22)
    ax.set_ylim(-CHORD - 5, 6)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("ONE STRIP  (section, ×5)", fontsize=11, weight="bold")

    ax = fig.add_axes([0.77, 0.04, 0.22, 0.89])
    ax.axis("off")
    b = B
    rows = [
        ("Strip rows per grille", f"{ROWS}  (×2 half-strips)"),
        ("Half-strip length", f"{HALF:.0f} mm"),
        ("Grilles / stages", f"{N+1} / {N}"),
        ("Strip length, total", f"{(N+1)*ROWS*2*HALF/1e3:.0f} m"),
        ("Emitter wire, total", f"{N*ROWS*2*HALF/1e3:.0f} m"),
        ("Strip mass", f"{b['m_strips']*1e3/((N+1)*ROWS*2*HALF/1e3):.2f} g/m"),
        ("Stage voltage, max", f"{b['V']/1e3:.0f} kV"),
        ("Mean field in gap", f"{b['E_avg']/1e5:.1f} kV/cm"),
        ("Corona onset", f"{b['V_onset']/1e3:.1f} kV"),
        ("Current per stage", f"{b['I_stage']*1e3:.1f} mA (+), ×1.3 (−)"),
        ("Peak field at nose", f"{b['E_coll']/1e6:.1f} MV/m"),
        ("Pressure rise / stage", f"{b['dp_stage']:.2f} Pa"),
        ("Duct velocity", f"{b['v_duct']:.2f} m/s"),
        ("Duct gain on body force", f"×{b['aug']:.2f}"),
    ]
    ax.text(0, 0.985, "Electrode data at full throttle", fontsize=11, weight="bold", va="top")
    for k, (a_, v) in enumerate(rows):
        y = 0.94 - 0.036 * k
        ax.text(0.0, y, a_, fontsize=9, va="top")
        ax.text(1.0, y, v, fontsize=9, va="top", ha="right")
        ax.plot([0, 1], [y - 0.027, y - 0.027], color="#dddddd", lw=0.6)
    ax.text(0, 0.40, "Build notes", fontsize=11, weight="bold", va="top")
    ax.text(0, 0.365,
            "• Nose radius matters: a sharp or frayed nose starts a\n  reverse discharge. Keep it smooth; no film edges or\n  wire ends pointing up-stream.\n"
            "• The wire must lie exactly on the trailing edge, free of\n  glue on its lower half. Use dots of CA every 40 mm.\n"
            "• Strips on the +rail and −rail are identical; polarity is\n  set only by which bus a grille is wired to.\n"
            "• G1 needs no metallised nose (nothing above it) and\n  G4 needs no wire; building them identical is fine.\n"
            "• Keep 20 mm between strip ends and the film wall,\n  and 20 mm between the two half-strips of a row.",
            fontsize=8.5, va="top", linespacing=1.35)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    title_block(fig, "Ion Lantern — electrode stack and strip construction", 2)
    save(fig, "02_electrode_detail")


# ------------------------------------------------------------------ 3. HV module schematic
def hv_schematic():
    import schemdraw
    import schemdraw.elements as elm

    schemdraw.use("matplotlib")
    d = schemdraw.Drawing(unit=2.2, fontsize=10)
    d.config(lw=1.3)
    yt, yb, ym = 8.0, 0.0, 4.0
    # --- input
    d += elm.Battery().endpoints((0, yb), (0, yt)).label("BT1\n6S LiPo\n18–25.2 V", loc="left").reverse()
    d += elm.Fuse().endpoints((0, yt), (2.4, yt)).label("F1  15 A")
    d += elm.Line().endpoints((2.4, yt), (9.2, yt))
    d += elm.Line().endpoints((0, yb), (9.2, yb))
    d += elm.Capacitor().endpoints((3.2, yt), (3.2, yb)).label("C1\n6×22 µF\n35 V X7R", loc="bottom")
    d += elm.Dot().at((3.2, yt))
    d += elm.Dot().at((3.2, yb))
    # --- full bridge
    for x, names in ((5.6, ("Q1", "Q2")), (9.2, ("Q3", "Q4"))):
        hi = elm.NFet(bulk=False).at((x, ym + 2.0)).anchor("source").theta(0).reverse()
        d += hi
        d += elm.Line().endpoints(hi.drain, (x, yt))
        lo = elm.NFet(bulk=False).at((x, ym)).anchor("drain").theta(0).reverse()
        d += lo
        d += elm.Line().endpoints(lo.source, (x, yb))
        d += elm.Line().endpoints(hi.source, lo.drain)
        d += elm.Dot().at((x, ym + 1.0))
        if x < 9:
            d += elm.Dot().at((x, yt))
            d += elm.Dot().at((x, yb))
        d += elm.Label().at(hi.gate).label(names[0], loc="left", ofst=0.1)
        d += elm.Label().at(lo.gate).label(names[1], loc="left", ofst=0.1)
    d += elm.Label().at((8.6, -0.75)).label("Q1–Q4: 40 V, ≤ 5 mΩ MOSFET (3 × 3 mm), two half-bridge gate drivers", fontsize=9)
    # --- resonant tank and transformer
    sa, sb = (5.6, ym + 1.0), (9.2, ym + 1.0)
    d += elm.Line().endpoints(sa, (6.6, ym + 1.0))
    d += elm.Line().endpoints((6.6, ym + 1.0), (6.6, 10.4))
    d += elm.Inductor2().endpoints((6.6, 10.4), (9.2, 10.4)).label("Lr  ≈ 1.5 µH\n(leakage + trim)")
    d += elm.Capacitor().endpoints((9.2, 10.4), (11.4, 10.4)).label("Cr  220 nF\nC0G / film")
    d += elm.Line().endpoints((11.4, 10.4), (12.6, 10.4))
    d += elm.Line().endpoints(sb, (10.6, ym + 1.0))
    d += elm.Line().endpoints((10.6, ym + 1.0), (10.6, 9.0))
    d += elm.Line().endpoints((10.6, 9.0), (12.6, 9.0))
    xfp, xfs = 12.6, 14.2
    d += elm.Inductor2(loops=3).endpoints((xfp, 10.4), (xfp, 9.0)).flip()
    d += elm.Inductor2(loops=5).endpoints((xfs, 10.9), (xfs, 8.5))
    d += elm.Line().endpoints((13.3, 10.9), (13.3, 8.5)).linewidth(0.8)
    d += elm.Line().endpoints((13.5, 10.9), (13.5, 8.5)).linewidth(0.8)
    d += elm.Label().at((11.9, 12.0)).label("T1  RM8 ferrite (N49 / 3F36)\n4 : 500, sectioned HV bobbin\n≈ 2.6 kV peak at the secondary", fontsize=9)
    # --- Cockcroft-Walton, positive ladder up, negative ladder down
    x0 = 16.4
    d += elm.Line().endpoints((xfs, 8.5), (x0, 8.5))
    d += elm.Dot().at((x0, 8.5))
    d += elm.Label().at((15.25, 7.7)).label("HV return\n(0 V)", fontsize=9)
    d += elm.Dot().at((xfs, 10.9))
    d += elm.Line().endpoints((xfs, 10.9), (xfs, 12.6))
    d += elm.Line().endpoints((xfs, 12.6), (x0 - 0.4, 12.6))
    d += elm.Line().endpoints((xfs, 10.9), (xfs + 0.7, 10.9))
    d += elm.Line().endpoints((xfs + 0.7, 10.9), (xfs + 0.7, 4.4))
    d += elm.Line().endpoints((xfs + 0.7, 4.4), (x0 - 0.4, 4.4))

    def ladder(y_ac, y_dc, sgn, tag):
        n1, n3 = (x0 + 2.2, y_ac), (x0 + 6.6, y_ac)
        g0, n2, n4 = (x0, y_dc), (x0 + 4.4, y_dc), (x0 + 8.8, y_dc)
        up, dn = ("top", "bottom") if sgn > 0 else ("bottom", "top")
        d.add(elm.Capacitor().endpoints((x0 - 0.4, y_ac), n1).label(f"C{tag}1", loc=up))
        d.add(elm.Capacitor().endpoints(n1, n3).label(f"C{tag}3", loc=up))
        d.add(elm.Capacitor().endpoints(g0, n2).label(f"C{tag}2", loc=dn))
        d.add(elm.Capacitor().endpoints(n2, n4).label(f"C{tag}4", loc=dn))
        for k, (p, q) in enumerate([(g0, n1), (n1, n2), (n2, n3), (n3, n4)]):
            e = elm.Diode().endpoints(p, q) if sgn > 0 else elm.Diode().endpoints(q, p)
            d.add(e.label(f"D{tag}{k+1}", loc="top", ofst=0.05, fontsize=8))
        for p in (n1, n3, n2):
            d.add(elm.Dot().at(p))
        return n4

    d += elm.Line().endpoints((x0, 8.5), (x0, 9.9))
    d += elm.Line().endpoints((x0, 8.5), (x0, 7.1))
    pos = ladder(12.6, 9.9, +1, "p")
    neg = ladder(4.4, 7.1, -1, "n")
    for node, y_out, nm, col in ((pos, 9.9, f"+HV  0 … +{VHALF:.0f} kV", RED), (neg, 7.1, f"−HV  0 … −{VHALF:.0f} kV", BLUE)):
        d += elm.Resistor().endpoints(node, (node[0] + 2.4, y_out)).label("Ro 100 kΩ\n3 W, HV type", loc="top" if y_out > 8 else "bottom", fontsize=8)
        d += elm.Dot(open=True).at((node[0] + 2.4, y_out)).label(nm, loc="right", color=col)
    d += elm.Resistor().endpoints((pos[0], 9.9), (pos[0], 8.5)).label("Rb", loc="bottom", fontsize=9)
    d += elm.Resistor().endpoints((neg[0], 7.1), (neg[0], 8.5)).label("Rb", loc="top", fontsize=9)
    d += elm.Dot().at((pos[0], 8.5))
    d += elm.Line().endpoints((x0, 8.5), (pos[0], 8.5)).linestyle(":")
    d += elm.Dot().at(pos)
    d += elm.Dot().at(neg)
    d += elm.Label().at((x0 + 4.4, 14.2)).label(
        "Bipolar 2-stage Cockcroft–Walton multiplier.   C: 1 nF, ≥ 9 kV (3 × 3 kV C0G 1812 in series).\n"
        "D: ≥ 8 kV, ≥ 30 mA, trr ≤ 100 ns (HV fast-recovery diode, or 8 × US1M in series with 10 MΩ grading resistors).\n"
        "Rb: 100 MΩ HV resistor per rail to return, with a 100 kΩ low-side tap for voltage sense (1000 : 1); bleeds the rails in < 1 s.",
        fontsize=9)
    # --- controller block
    bx0, bx1, by0, by1 = -0.6, 14.2, -4.6, -1.7
    for p, q in (((bx0, by1), (bx1, by1)), ((bx1, by1), (bx1, by0)), ((bx1, by0), (bx0, by0)), ((bx0, by0), (bx0, by1))):
        d += elm.Line().endpoints(p, q)
    d += elm.Label().at(((bx0 + bx1) / 2, (by0 + by1) / 2)).label(
        "U1  module controller (STM32G0-class MCU + two half-bridge gate drivers)\n"
        "• reads the thrust command (PWM 1–2 ms) from the flight controller and maps it to a rail-voltage set-point\n"
        "• closes the loop on rail voltage via the Rb taps; sets bridge frequency / phase shift (≈ 250 kHz, zero-voltage switching)\n"
        "• arc / over-current trip from a shunt in the HV return: gates off within 20 µs, automatic retry after 50 ms, FAULT line to the FC", fontsize=9)
    d += elm.Line().endpoints((2.0, by1), (2.0, -1.0)).linestyle("--")
    d += elm.Label().at((2.0, -0.75)).label("gate drive to Q1–Q4", fontsize=8)
    d += elm.Line().endpoints((bx1, -3.15), (pos[0], -3.15)).linestyle("--")
    d += elm.Line().endpoints((pos[0], -3.15), (pos[0], 6.4)).linestyle("--")
    d += elm.Label().at((pos[0] - 3.0, -2.8)).label("rail-voltage sense, return-current sense", fontsize=8)
    d += elm.Line().endpoints((-2.6, -3.15), (bx0, -3.15)).linestyle("--")
    d += elm.Label().at((-1.7, -2.8)).label("PWM in", fontsize=8)
    d += elm.Label().at((14.0, 16.3)).label(
        "Ion Lantern — HV module (one per zone, four identical): 150 W, 0 … ±10 kV — reference schematic", fontsize=15)
    d += elm.Label().at((14.0, -5.4)).label(
        "Reference design, not a tested board: verify every voltage, current and creepage rating before building.  "
        "Short-duration rating (90 s on / 10 min off).", fontsize=8.5, color=GREY)
    fig = d.draw(show=False).fig
    fig.set_size_inches(19, 11.5)
    save(fig, "03_hv_module_schematic")


# ------------------------------------------------------------------ 4. system wiring
def box(ax, xy, wh, text, fc="#f4f4f4", ec="k", fs=9, tc="k"):
    x, y = xy
    w, h = wh
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec=ec, lw=1.1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc)
    return (x, y, w, h)


def wire(ax, pts, color="k", lw=1.4, ls="-", label=None, lpos=0.5, fs=7.5, dy=0.08):
    xs, ys = zip(*pts)
    ax.plot(xs, ys, color=color, lw=lw, ls=ls, solid_capstyle="round")
    if label:
        i = max(0, min(len(pts) - 2, int(lpos * (len(pts) - 1))))
        ax.text((xs[i] + xs[i + 1]) / 2, (ys[i] + ys[i + 1]) / 2 + dy, label, fontsize=fs, ha="center", va="bottom", color=color)


def system_wiring():
    fig = plt.figure(figsize=(17, 9.5))
    ax = fig.add_axes([0.01, 0.03, 0.98, 0.90])
    ax.set_xlim(0, 17)
    ax.set_ylim(0, 9.4)
    ax.axis("off")
    bat = box(ax, (0.3, 6.5), (2.2, 1.3), f"BATTERY\n6S 650 mAh 95C LiPo\n14.4 Wh, {REP['battery']['mass']*1e3:.0f} g, XT30", fc="#fde9d9", fs=8.5)
    arm = box(ax, (3.1, 6.7), (1.7, 0.9), "ARMING PLUG\n+ 40 A fuse", fc="#fff3c4")
    bus = box(ax, (5.4, 6.5), (1.3, 1.3), "POWER\nBUS\n22 V", fc="#eeeeee")
    wire(ax, [(2.5, 7.15), (3.1, 7.15)], RED, 2.4)
    wire(ax, [(4.8, 7.15), (5.4, 7.15)], RED, 2.4)
    bec = box(ax, (5.4, 4.7), (1.3, 0.9), "5 V BEC\n1 A", fc="#eeeeee")
    wire(ax, [(6.05, 6.5), (6.05, 5.6)], RED, 1.6)
    fc_ = box(ax, (3.0, 2.2), (2.6, 1.7), "FLIGHT CONTROLLER\nINAV-class board, IMU + baro\nmixer: quad-X, 4 PWM outputs", fc="#e4e0f5")
    rx = box(ax, (0.3, 3.2), (1.9, 0.8), "RC RECEIVER\n(ELRS nano)", fc="#e4e0f5")
    of = box(ax, (0.3, 1.9), (1.9, 0.9), "OPTICAL FLOW +\nToF RANGE SENSOR\n(looks down)", fc="#e4e0f5", fs=8)
    wire(ax, [(2.2, 3.6), (3.0, 3.4)], "#5b4ea8", 1.2, label="CRSF", fs=7)
    wire(ax, [(2.2, 2.35), (3.0, 2.6)], "#5b4ea8", 1.2, label="UART", fs=7)
    wire(ax, [(6.05, 4.7), (6.05, 4.3), (4.3, 4.3), (4.3, 3.9)], RED, 1.2, label="5 V", lpos=0.5, fs=7)
    # HV modules and zones
    zones = [("A", 7.6), ("B", 5.6), ("C", 3.6), ("D", 1.6)]
    for k, (nm, y) in enumerate(zones):
        m = box(ax, (8.3, y - 0.55), (2.4, 1.1), f"HV MODULE {nm}\n22 V → 0 … ±{VHALF:.0f} kV\n150 W, ≤ 54 g", fc="#dbe9f6")
        wire(ax, [(6.7, 7.15 - 0.08 * k), (7.3 + 0.12 * k, 7.15 - 0.08 * k), (7.3 + 0.12 * k, y + 0.25), (8.3, y + 0.25)], RED, 1.8)
        xc = 6.85 + 0.1 * (3 - k)
        wire(ax, [(5.6, 3.5 - 0.3 * k), (xc, 3.5 - 0.3 * k), (xc, y - 0.25), (8.3, y - 0.25)], "#5b4ea8", 1.1, ls="--")
        z = box(ax, (12.9, y - 0.62), (3.7, 1.24), "", fc="white")
        ax.text(12.98, y + 0.47, f"ZONE {nm}  (16 rows × {HALF:.0f} mm)", fontsize=8, va="center", weight="bold")
        for g in range(N + 1):
            col = RED if g % 2 == 0 else BLUE
            yy = y + 0.25 - 0.22 * g
            ax.plot([13.9, 16.4], [yy, yy], color=col, lw=2.0)
            ax.text(16.45, yy, f"G{g+1}", fontsize=7, va="center", color=col)
        # +HV to G1,G3 ; -HV to G2,G4
        wire(ax, [(10.7, y + 0.3), (11.5, y + 0.3)], RED, 1.6, label=f"+{VHALF:.0f} kV", fs=7, dy=0.04)
        wire(ax, [(10.7, y - 0.3), (11.5, y - 0.3)], BLUE, 1.6, label=f"−{VHALF:.0f} kV", fs=7, dy=-0.3)
        wire(ax, [(11.5, y + 0.3), (13.3, y + 0.3), (13.3, y + 0.25), (13.9, y + 0.25)], RED, 1.6)
        wire(ax, [(13.3, y + 0.25), (13.3, y - 0.19), (13.9, y - 0.19)], RED, 1.6)
        wire(ax, [(11.5, y - 0.3), (13.6, y - 0.3), (13.6, y + 0.03), (13.9, y + 0.03)], BLUE, 1.6)
        wire(ax, [(13.6, y - 0.3), (13.6, y - 0.41), (13.9, y - 0.41)], BLUE, 1.6)
    ax.text(6.6, 4.05, "4 × PWM thrust\ncommand (dashed)", fontsize=7.5, color="#5b4ea8", ha="center")
    ax.text(11.8, 8.75, "HV leads: silicone wire rated ≥ 15 kV DC, as short as possible,\n+ and − of a zone routed ≥ 40 mm apart, no sharp solder points",
            fontsize=8, ha="center")
    ax.text(0.3, 0.05,
            "Mixer (same as a quad-X drone): zone thrust = collective ± roll ± pitch ± yaw.  Roll and pitch come from the thrust difference across the duct;\n"
            "yaw from small canted exit tabs under each zone (A, C one way; B, D the other).  Each module also returns a FAULT line.\n"
            "Failsafe: loss of RC link or any module FAULT → all PWM to zero.\n"
            "The pod electronics are the 0 V reference of all four HV modules; the airframe carries no ground conductor.",
            fontsize=8.5, va="bottom")
    ax.text(0.3, 8.75, "Low-voltage side: 16 AWG silicone, XT30.  Total draw in hover ≈ "
            f"{REP['hover']['P_batt']/21.5:.0f} A ({REP['hover']['P_batt']:.0f} W), {B['P_batt']/21.5:.0f} A at full throttle.", fontsize=8.5)
    title_block(fig, "Ion Lantern — system wiring", 4)
    save(fig, "04_system_wiring")


if __name__ == "__main__":
    general_arrangement()
    electrode_detail()
    system_wiring()
    hv_schematic()
    print("drawings written to", HERE)
