# Bill of materials and mass budget

Masses are the ones the sizing model uses (`results/design.json`); the flight simulation
is built from the same numbers. Anything heavier than its line here comes straight out of
the 83 g of thrust margin, so weigh every sub-assembly as you build it.

| Group | Mass | Share |
|---|---:|---:|
| Electrode strips (4 grilles) | 61 g | 12 % |
| Frame | 41 g | 8 % |
| Duct skin and inlet lip | 20 g | 4 % |
| Threads, HV leads, joints, feet | 12 g | 2 % |
| HV modules (4) | 214 g | 43 % |
| Battery | 116 g | 23 % |
| Avionics | 12 g | 2 % |
| Payload | 20 g | 4 % |
| **All-up** | **496 g** | |
| Maximum thrust (model) | 579 gf | T/W = 1.17 |

## Electrodes — 61 g

| Item | Specification | Quantity | Mass |
|---|---|---:|---:|
| Collector and tail film | 6 µm polypropylene capacitor film, ≈ 5.4 g/m². Collector: metallised one side, 32 mm wide, wrapped over the nose. Tail: clear (no metal), 20 mm wide. | 120 m of strip, 5.8 m² in all | 31 g |
| Nose stick | XPS/Depron foam, ≈ 35 kg/m³, 1.5 × 3 mm, top edge sanded to a full round | 120 m | 19 g |
| Emitter wire | Tungsten, Ø 50 µm (44 AWG), plain or gold-plated | 90 m (+ spare) | 5 g |
| Adhesive | Thin CA in dots; 3M 77-type spray for film to foam, used sparingly | — | 6 g |

256 half-strips in total: 32 rows × 2 halves × 4 grilles, each 470 mm long.

## Frame — 41 g

| Item | Specification | Quantity | Mass |
|---|---|---:|---:|
| Perimeter frames, top and bottom | Pultruded carbon tube 1.5 × 0.7 mm, ≈ 2.2 g/m | 8 × 1.0 m | 17.6 g |
| Support rails and hub struts | Same tube; per level 3 interior rails + 1 m of diagonal stays to the mast | 8 × 1.0 m | 17.6 g |
| Corner longerons and centre mast | Pultruded carbon tube 2 × 1 mm, ≈ 3.8 g/m | 5 × 0.306 m | 5.8 g |

## Duct skin and lip — 20 g

| Item | Specification | Quantity | Mass |
|---|---|---:|---:|
| Wall film | 6 µm PET (Mylar), ≈ 8.4 g/m², clear, **not** metallised | 4 × 1.0 × 0.306 m | 10.3 g |
| Lip film | Same film over the formers, half-round R35 | 0.44 m² | 3.7 g |
| Seams | 6 mm PET tape | ≈ 6 m | 2 g |
| Lip formers | 1 mm foam half-discs R35 every 100 mm | 40 | 4 g |

## Threads, leads, joints, feet — 12 g

| Item | Specification | Mass |
|---|---|---:|
| Hanging threads | Aramid (Kevlar) thread, 20–40 dtex | 1 g |
| HV leads | Silicone-insulated wire rated ≥ 15 kV DC, 28 AWG, 8 runs × 0.3 m | 3 g |
| Grille bus wires | Enamelled copper Ø 0.1 mm along one end of every grille half | 0.5 g |
| Joints | Thread-and-CA lashings at 28 tube joints | 3 g |
| Leg extensions and feet | 4 × 120 mm of longeron tube, foam foot pads | 2.5 g |
| Pod mount and shield | 25 µm polyimide sheet under the pod, foam cradle | 2 g |

## HV modules — 4 × 54 g = 214 g

One per zone, identical. See `03_hv_module_schematic`. The budget is the modelling
assumption of 2.7 kW/kg at 145 W per module; the bottom-up estimate below leaves about a
quarter of that as contingency.

| Item | Specification | Mass each |
|---|---|---:|
| Transformer T1 | RM8 ferrite core set (N49 or 3F36), 4-turn litz primary, 500-turn secondary of Ø 0.1 mm wire in ≥ 3 sections, vacuum-potted | 20 g |
| Multiplier board | 8 capacitor positions (3 × 1 nF 3 kV C0G 1812 in series each), 8 diode positions (≥ 8 kV, ≥ 30 mA, trr ≤ 100 ns), 0.8 mm FR-4, conformal coat | 10 g |
| Inverter and controller | 4 × 40 V ≤ 5 mΩ MOSFET, 2 half-bridge drivers, small MCU, 6 × 22 µF 35 V input capacitors, resonant capacitor, 0.8 mm FR-4 | 8 g |
| Output and bleed resistors, leads, potting | 2 × 100 kΩ HV, 2 × 100 MΩ HV, silicone leads | 5 g |
| Bottom-up subtotal | | 43 g |
| Budget (model) | | 54 g |

Required rating: 18–25.2 V in; 0 to ±10 kV out (20 kV rail-to-rail) at 7.5 mA; 150 W
for 90 s on / 10 min off; efficiency ≥ 88 %; regulation bandwidth ≥ 30 Hz.

## Battery — 116 g

Tattu R-Line 6S1P 650 mAh 95C, XT30 (retail listings give 115–118 g). Any 6S 650 mAh
pack of ≤ 118 g rated ≥ 75C continuous is equivalent. The pack runs at about 38C in hover
and 46C at full throttle.

## Avionics — 12 g

| Item | Specification | Mass |
|---|---|---:|
| Flight controller | Small INAV-capable board with IMU and barometer, 4 PWM outputs | 4.5 g |
| Receiver | ELRS nano | 0.6 g |
| Optical-flow and ToF range sensor | Downward-looking module (PMW3901 + VL53L1X class) | 2 g |
| 5 V regulator | 1 A switching BEC | 1 g |
| Power wiring | XT30, arming plug with 40 A fuse, 16 AWG silicone | 4 g |

## Payload — 20 g

Hung under the centre mast, below the exit plane. The model predicts lift-off with up to
about 50 g (T/W 1.10), but with very little control margin left.

## Tools and consumables

Flat building board ≥ 1.1 × 1.1 m; strip-assembly jig (two parallel fences 470 mm apart
with a 1.5 mm slot for the foam stick); foam cutter or steel rule and fresh blades; 0.01 g
scale; kitchen scale reading 0.1 g up to 1 kg; HV probe (1000:1, ≥ 30 kV) and a
micro-ammeter; earthed discharge stick; insulating stand; CO₂ extinguisher; LiPo bag.
