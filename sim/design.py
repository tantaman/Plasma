"""
Frozen design point of the vehicle ("Ion Lantern").  Everything downstream -- the
flight simulation, the drawings, the bill of materials -- reads its numbers from here
and from results/design.json, which sizing.py computes from these inputs.
"""

# --- thruster ------------------------------------------------------------------------
D_GAP = 0.030      # m, emitter wire -> nose of the next collector strip (= strip pitch)
R_WIRE = 25e-6     # m, emitter radius (50 um tungsten, 44 AWG -- the wire validated against MIT data)
N_STAGES = 3       # acceleration stages (4 grilles of strips)
DUCT_W = 1.00      # m, inner side of the square duct
V_MAX = 20.0e3     # V across every stage at full throttle (rails at +10 kV and -10 kV)
N_QUADS = 4        # independently throttled HV zones (2 x 2)

# --- battery: Tattu R-Line 6S1P 650 mAh 95C (retail listings: 115-118 g) ---------------
BATTERY = dict(
    name="Tattu R-Line 6S 650 mAh 95C (XT30)",
    cells=6,
    capacity_ah=0.65,
    mass=0.116,          # kg
    r_cell=0.0092,       # ohm per cell (6 mOhm*Ah / 0.65 Ah), typical of race packs when warm
    r_leads=0.006,       # ohm, leads + connector
)

# --- fixed equipment -------------------------------------------------------------------
M_AVIONICS = 0.012   # kg: flight controller + receiver + ToF/optical-flow sensor + wiring
M_PAYLOAD = 0.020    # kg: design payload


def cell():
    """Electrode unit-cell geometry (kwargs for corona.mesh_cell)."""
    from sizing import cell_geometry

    return cell_geometry(D_GAP, R_WIRE)


# kept for validate_corona.py
STAGE_V = V_MAX
