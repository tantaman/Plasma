# Build and test guide

This is a paper design. Nothing here has been built or flown, and the predicted thrust
margin is small (see the README). The guide is therefore ordered so that the cheap,
decisive measurements come first: each test gate has a number to hit, and missing it
means stopping to fix the design, not carrying on.

Drawings: `hardware/01_general_arrangement`, `02_electrode_detail`,
`03_hv_module_schematic`, `04_system_wiring` (SVG and PNG). Parts and masses:
`hardware/BOM.md`.

## Safety — read first

- **High voltage.** Four supplies of ±10 kV, together able to deliver about 30 mA
  continuously. That is enough to kill. Never touch, adjust or approach the craft while
  the battery is connected. Arm and disarm only through the arming plug on a 1 m
  insulated lanyard or by radio. After disarming, wait 5 s, then touch every grille bus
  with an earthed discharge stick before handling.
- **Work alone? Don't.** Have a second person present who knows where the arming plug is.
- **Ozone and nitrogen oxides.** A 470 W corona, half of it negative polarity, produces
  ozone at a rate that will exceed workplace limits in an ordinary room within one
  flight. Run only outdoors in still air or in a large, well-ventilated hall; keep runs
  to about a minute; leave and ventilate afterwards. Stop at once if you smell ozone
  strongly or feel throat or chest irritation.
- **Fire.** Arcs can ignite 6 µm film and foam. Keep a CO₂ extinguisher at hand, fly
  over a non-flammable floor, and do not defeat the arc trip in the HV modules.
- **Battery.** The pack is discharged at 40–50C. Use a pack in good condition, charge
  it in a LiPo bag, let it cool between flights, and retire it if it puffs.
- **Bystanders and electronics.** Keep people 3 m away. The corona radiates broadband
  radio noise; expect interference with nearby receivers.

## Test gates at a glance

| Gate | What you measure | Pass |
|---|---|---|
| G-A | One zone-sized single stage on a scale | Current within 15 % of the table in step 3; thrust at 20 kV ≥ 0.43 N |
| G-B | One HV module on a resistive load | 145 W at ±10 kV for 90 s, ≥ 88 % efficient, ≤ 54 g |
| G-C | One full 3-stage zone with its module | ≥ 1.35 N at 20 kV, no sparking |
| G-D | Complete craft weighed | ≤ 476 g without payload |
| G-E | Complete craft on a scale, all zones | ≥ 5.4 N at 20 kV |
| G-F | Tethered hover | Stable for 30 s with slack tethers |

## 1. Strips

1. Cut foam sticks 1.5 × 3 × 470 mm. Sand one 1.5 mm edge to a full half-round along
   the whole length; this is the collector nose. Reject any stick with a nick.
2. Cut metallised film strips 32 × 470 mm (the collector) and clear, unmetallised film
   strips 20 × 470 mm (the tail).
3. Wrap the metallised strip, metal outwards, over the rounded edge of the stick so that
   15 mm hangs down each face; bond lightly to the foam. The film must be smooth round
   the nose, with no wrinkle or lifted edge facing upwards.
4. Slip the top 2 mm of the clear strip between the two metallised faces just below the
   stick and bond the three layers together. The clear tail now hangs 18 mm below the
   end of the metal.
5. Lay the Ø 50 µm tungsten wire along the bottom edge of the tail under light tension
   and fix it with a dot of thin CA every 40 mm. The lower half of the wire must stay
   clean. Leave 15 mm of wire free at one end and jumper it to the metallised film of
   the same strip with a dab of conductive paint.
6. Weigh: a finished half-strip should be 0.24 g. Make 256 plus spares.

## 2. One-stage test rig (before building anything else)

Build a 470 × 480 mm open frame holding two grilles of 16 half-strips at 30 mm pitch,
the upper grille's wires 30 mm above the lower grille's noses, inside a 330 mm tall
film-walled duct with the R35 lip. This is one zone with one stage.

## 3. Gate G-A — does the physics match?

Stand the rig on an insulating column on a 0.01 g scale, exit at least 300 mm above the
table. Feed it from a bench HV supply (+rail to the upper grille, −rail to the lower)
through a micro-ammeter. The model predicts, for this rig (one positive stage, 7.5 m of
strip):

| Stage voltage | Current | Thrust |
|---:|---:|---:|
| 12 kV | 0.48 mA | 0.11 N (11 gf) |
| 16 kV | 1.21 mA | 0.28 N (29 gf) |
| 20 kV | 2.20 mA | 0.51 N (52 gf) |

A single stage has a higher duct gain than three in series, so these thrusts are per
stage about 12 % above what each stage contributes in the full stack. If current is
right but thrust is low, the duct is under-performing: check the lip, leaks, and that
the exit is unobstructed. If both are low, check the gap and wire condition. Repeat
with the polarity reversed; the negative stage should draw about 30 % more current for
about 10 % more thrust. Raise the voltage slowly to find the sparking limit; it must be
above 22 kV.

**Stop here if thrust at 20 kV is below 0.43 N** (85 % of the prediction). The sizing
model says the full craft then no longer lifts its own weight.

## 4. HV modules and gate G-B

Build one module to `03_hv_module_schematic`. Wind the secondary in at least three
sections with interlayer tape and pot it under vacuum. Space the multiplier parts for
10 kV creepage (≥ 15 mm between rails on the board, slots milled between stages) and
round every solder joint. Test on a 2 × 1.3 MΩ HV resistor load (±10 kV, 7.5 mA) for
90 s while logging input power and case temperature, then test the arc trip by drawing
an arc to a grounded ball through 1 MΩ. Weigh it. Build the other three only when this
one passes.

## 5. Frame, skin and lip

1. On the building board, lash and glue the two 1.0 × 1.0 m perimeter frames from
   1.5 mm tube, each with three interior rails at 245 mm spacing and diagonal stays to
   a centre hub.
2. Join them with the four corner longerons and the centre mast at 306 mm spacing; the
   longerons run 120 mm below the lower frame as legs.
3. Check squareness by measuring diagonals, then fix with thread bracing across each
   side face.
4. Glue the foam lip formers round the top frame and skin the lip, then the four walls,
   with 6 µm PET film. Seal all seams; the walls carry a few pascals of suction.

## 6. Hanging the grilles

1. Each strip hangs by two aramid threads from every rail it crosses (edge rail, one
   interior rail, centre rail), tied round the foam stick. Hang grille G1 with its
   noses 60 mm below the top frame.
2. Hang G2 from G1, G3 from G2, G4 from G3 with thread drops that leave exactly 30 mm
   from each wire to the nose below (63 mm nose to nose). Use a 30 mm gauge block.
3. Tie the G4 sticks down to the bottom rails so that thrust cannot lift the stack.
4. Run one Ø 0.1 mm bus wire along the wall end of every grille half and bond each
   strip's metallised film to it with conductive paint. G1 and G3 of a zone go to that
   zone's +rail, G2 and G4 to its −rail.
5. Sight along every row: wires and noses must line up within ±2 mm.

## 7. Pod and wiring

Mount battery, the four modules, flight controller, receiver and regulator on the mast
so the pod's underside is 45 mm above the G1 noses, with a polyimide sheet beneath it.
Route each zone's two HV leads down opposite corners of that zone, at least 40 mm apart
and 20 mm from any strip of the other polarity. Fit the optical-flow sensor and the
payload under the mast below the exit. Fix a small film tab under each zone's exit,
canted 10°, zones A and C one way and B and D the other, for yaw control.

## 8. Gates G-C to G-F

- **G-C.** Power one zone alone from its module on the scale (other zones unplugged).
- **G-D.** Weigh the complete craft.
- **G-E.** Whole craft on the insulating column on the scale; ramp all zones together
  to 20 kV. Log thrust, pack current and voltage. Hover needs 4.9 N; the model says
  5.7 N. Below 5.4 N there is not enough margin to control it.
- **G-F.** Configure the flight controller as a quad-X with the four PWM outputs mapped
  to zones A–D, low rates and a 7° tilt limit. Fly from an open grating 0.3 m above the
  floor so the exit is not blocked, with four slack threads from the leg tips to floor
  anchors limiting travel to 0.3 m. Only then try free hover, indoors, in still air:
  the simulation holds position against a 0.18 m/s draught only with about 0.4 m of
  drift.

## Expected flight envelope (from the simulation)

Hover at 18.9 kV drawing 24.7 A (540 W) from the pack; about 75 s of hover before the
22 % reserve; climb limited to roughly 0.2 m/s by the 17 % thrust margin; translation
about 0.1 m/s; no tolerance for wind.
