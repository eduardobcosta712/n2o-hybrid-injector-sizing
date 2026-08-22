# Example 1 — Sizing Mode: Predict Mass Flow for a Known Injector Geometry

This example shows how to use **Sizing mode** to predict the real oxidiser mass flow delivered by an existing injector design, accounting for two-phase effects in the feed path.

---

## Scenario

A hybrid motor test bench has the following configuration:

- N₂O tank at **20 °C**, pressurised to **58 bar**
- Feed line: 800 mm of 10 mm ID tubing → ball valve → 500 mm of 10 mm ID tubing → 90° elbow
- Injector: **6 orifices**, each **1.5 mm diameter** (Cd = 0.65)
- Chamber pressure: **22 bar**

The team sized the injector using the SPI model and expects approximately 514 g/s. The question is: what does the system actually deliver?

---

## Inputs

| Parameter | Value |
|---|---|
| Tank temperature | 20 °C |
| Tank pressure | 58 bar |
| P_sat at 20 °C | 51.37 bar |
| Subcooling margin | **6.63 bar** |
| Chamber pressure | 22 bar |
| Discharge coefficient Cd | 0.65 |
| Number of orifices | 6 |
| Orifice diameter | 1.5 mm |
| Total orifice area | 10.603 mm² |

**Feed line segments:**

| # | Type | Length | ID | K |
|---|---|---|---|---|
| 1 | Pipe | 800 mm | 10 mm | — |
| 2 | Fitting | — | 10 mm | 0.05 (ball valve, fully open) |
| 3 | Pipe | 500 mm | 10 mm | — |
| 4 | Fitting | — | 10 mm | 0.90 (elbow 90°, standard) |

---

## Results

**Feed line evaluation:**
- No flashing detected — pressure at the injector inlet is **57.51 bar**, still well above P_sat(20 °C) = 51.37 bar.
- The feed line is short and uses low-K fittings, so it preserves most of the tank subcooling margin.

**Injector model selection:**
- P_downstream = 22 bar < P_sat(20 °C) = 51.37 bar → the pressure drop across the orifice crosses the saturation curve inside the orifice.
- SPI alone is **not sufficient** → **Dyer model** applied.

| Model | Predicted mass flow |
|---|---|
| SPI (single-phase, reference) | 514.5 g/s |
| HEM (full equilibrium limit) | — |
| **Dyer (adopted)** | **368.1 g/s** |
| SPI over-prediction | **28.5%** |

**Dyer model parameters:**
- Non-equilibrium parameter κ = 1.100
- Exit vapour quality x = 0.299 (≈ 30% of the mass exits as vapour)

---

## Interpretation

The SPI model over-predicts the mass flow by 28.5%. If the motor was designed around 514 g/s, the real O/F ratio will be lower than expected, reducing specific impulse and potentially compromising combustion stability.

**Combustion stability check:**
- Injector pressure drop (57.51 − 22 = 35.51 bar) / P_chamber (22 bar) = **161%** — well above the 15% minimum. Stable.

**Corrective options if 514 g/s is required:**
1. Increase total orifice area — size with Dyer instead of SPI (see Example 2).
2. Reduce chamber pressure to lower the required mass flow at the design O/F.

---

## How to reproduce in the tool

1. Open **Sizing mode**.
2. Set sidebar: T = 20 °C, P = 58 bar, Cd = 0.65, P_chamber = 22 bar.
3. Set orifice: 6 holes, 1.5 mm diameter.
4. Add feed line segments as listed above.
5. Read the results — Dyer prediction should show **368.1 g/s** and SPI over-prediction **28.5%**.
