# Example 1 — Sizing Mode: Predict Mass Flow for a Known Injector Geometry

> **Property-backend status: CONFIRMED (September 2026).** All figures below
> were regenerated after installing the real CoolProp package (8.0.0) and
> confirming the full test suite (262 tests) and `validation/waxman_2013_validation.py`
> against it. They replace the earlier, Perry/McGill-backed figures (381.8 g/s,
> 25.9% over-prediction), which are now superseded — the differences (about
> 1 percentage point) are consistent with the ~2–5% error the old correlations
> were known to carry.

This example shows how to use **Sizing mode** to predict the real oxidiser mass flow delivered by an existing injector design, accounting for two-phase effects in the feed path.

---

## Scenario

A hybrid motor test bench has the following configuration:

- N₂O tank at **20 °C**, pressurised to **58 bar**
- Feed line: 800 mm of 10 mm ID tubing → ball valve → 500 mm of 10 mm ID tubing → 90° elbow
- Injector: **6 orifices**, each **1.5 mm diameter** (Cd = 0.65)
- Chamber pressure: **22 bar**

The team sized the injector using the SPI model and expects approximately 515 g/s. The question is: what does the system actually deliver?

---

## Inputs

| Parameter | Value |
|---|---|
| Tank temperature | 20 °C |
| Tank pressure | 58 bar |
| P_sat at 20 °C | 50.53 bar |
| Subcooling margin | **7.47 bar** |
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

**Feed line evaluation (at the converged flow):**
- No flashing detected — pressure at the injector inlet is **57.56 bar**, still well above P_sat(20 °C) = 50.53 bar.
- The feed line is short and uses low-K fittings, so it preserves most of the tank subcooling margin.

**Injector model selection:**
- P_downstream = 22 bar < P_sat(20 °C) = 50.53 bar → the pressure drop across the orifice crosses the saturation curve inside the orifice.
- SPI alone is **not sufficient** → **Dyer model** applied.

| Model | Predicted mass flow |
|---|---|
| SPI (single-phase, reference) | 515.0 g/s |
| HEM (full equilibrium limit) | 242.0 g/s |
| **Dyer (adopted)** | **386.1 g/s** |
| SPI over-prediction | **25.0%** |

**Dyer model parameters:**
- Non-equilibrium parameter κ = 1.117
- Exit vapour quality x = 0.291 (≈ 29% of the mass exits as vapour)

---

## Non-equilibrium choking diagnostic

Running this exact case in the tool also shows a **choking warning**: the Dyer prediction exceeds the Henry-Fauske (1971) non-equilibrium choking ceiling, which at these tank/chamber conditions is **322.1 g/s** — 17% below the 386.1 g/s Dyer prediction.

This does **not** change the 386.1 g/s figure — see `docs/future_work.md`, Priority 1, for why the ceiling is surfaced as a diagnostic warning rather than applied automatically: it is theoretically sound (a primary-source non-equilibrium critical flow model) but has only been confirmed *not* to interfere with the Waxman-validated cases (8–14 bar pressure drop); it has not itself been experimentally validated at this example's conditions (here, roughly 35.6 bar pressure drop, well outside that validated band).

**How to read this result honestly:** the Dyer model itself is validated only at 8–14 bar. At 35 bar the prediction is a model estimate, and the choking diagnostic says the true flow may be as low as ≈322 g/s. A defensible planning range is therefore roughly **322–386 g/s**, to be narrowed by a cold-flow or hot-fire measurement — but in either case far below the SPI figure the motor was designed around.

---

## Interpretation

The SPI model over-predicts the mass flow by 25.0% against the Dyer estimate (and by about 37% if the choking ceiling turns out to be the binding one). If the motor was designed around 515 g/s, the real O/F ratio will be lower than expected, reducing specific impulse and potentially compromising combustion stability.

**Combustion stability check:**
- Injector pressure drop (57.56 − 22 = 35.56 bar) / P_chamber (22 bar) = **161.7%** — well above the 15% minimum. Stable.

**Corrective options if 515 g/s is required:**
1. Increase total orifice area — size with Dyer instead of SPI (see Example 2).
2. Reduce chamber pressure to lower the required mass flow at the design O/F.

---

## How to reproduce in the tool

1. Open **Sizing mode**.
2. Set sidebar: T = 20 °C, P = 58 bar, Cd = 0.65, P_chamber = 22 bar.
3. Set orifice: 6 holes, 1.5 mm diameter.
4. Add feed line segments as listed above.
5. Read the results — Dyer prediction should show **386.1 g/s** and SPI over-prediction **25.0%**. An amber choking-ceiling warning (322.1 g/s) should also appear underneath the Dyer caption.
