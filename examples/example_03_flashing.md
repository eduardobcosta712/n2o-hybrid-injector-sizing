# Example 3 — Flashing Case: Diagnosis and Design Correction

> **Property-backend note (September 2026).** All figures below were computed with the Perry/McGill + NIST property set used before the CoolProp integration (`docs/future_work.md`, Priority 4). They have not yet been regenerated with the real CoolProp package (unavailable in the development sandbox); expect small shifts (a few per cent) once that is done. Re-run this scenario after `pip install CoolProp` to confirm.

This example shows how to use the tool to **diagnose a flashing feed line**, understand the cause, and apply corrections until the system operates safely. It illustrates the physical mechanism that caused unexpected results in several documented university hybrid motor tests.

> **Regenerated in the September 2026 audit.** The previous version of this example was written for an earlier model (single-phase line, uncorrected Dyer weights, no coupled solver) and contained statements that the current model does not support: that raising the tank pressure alone would still leave the line flashing, that a tank at 45 bar could have "5 bar of margin" at 22 °C (impossible: P_sat(22 °C) = 53.68 bar), and that pressure drop scales as 1/D⁴ so that doubling the diameter cuts friction 16×. All figures below were produced by running the current code.

---

## Scenario

A team designs a motor with N₂O at **22 °C**, tank pressure **54 bar**, and the following feed line:

- 2500 mm of 8 mm ID tubing
- Needle valve (fully open, K = 2.0)
- 90° elbow (standard, K = 0.9)
- Injector: 4 holes, 1.8 mm diameter (Cd = 0.65), total area 10.18 mm²
- Chamber pressure: 18 bar

The SPI model was used for the initial design (a naive SPI calculation with the full 36 bar tank-to-chamber drop gives **492 g/s**). On test, the motor produces less thrust than expected and combustion is unstable.

---

## Initial configuration — what goes wrong

| Parameter | Value |
|---|---|
| Tank temperature | 22 °C |
| Tank pressure | 54 bar |
| P_sat at 22 °C | **53.68 bar** |
| Subcooling margin | **only 0.32 bar** |

With only 0.32 bar of margin, even small friction losses are enough to push the fluid below its saturation pressure.

**Feed line trace (at the converged operating point, 189 g/s):**

| Segment | Type | P after (bar) | ΔT_sub (K) | x at segment start |
|---|---|---|---|---|
| Tank exit | — | 54.00 | +0.27 | — |
| After pipe (2.5 m) | Pipe | 53.55 | −0.11 | 0 |
| After needle valve | Fitting (K=2.0) | 53.37 | −0.27 | 0.0021 |
| After elbow | Fitting (K=0.9) | 53.28 | −0.34 | 0.0048 |

**Result:** flashing is detected already after the first pipe segment (a friction loss of only 0.45 bar exhausts the 0.32 bar margin). The fluid reaches the injector with vapour quality x = 0.006. The tool applies the **HEM two-phase-inlet** model and predicts **189 g/s** — about **2.6× below** the naive SPI figure of 492 g/s. (This two-phase-inlet result is implemented and unit tested but has not been validated against experimental data; see "Caveats" below.)

---

## Diagnosis

The tool's diagnostic panel identifies four contributing factors:

1. **Subcooling margin below 5 bar** — only 0.32 bar. Any friction loss causes flashing.
2. **Long feed line** — 2.5 m total; the panel targets under 1.5 m.
3. **Small pipe diameter** — 8 mm ID; the panel suggests 10 mm or more.
4. **High-K fitting** — the needle valve (K = 2.0). A fully open ball valve (K = 0.05) is 40 times lower.

---

## Correction — step by step

The three changes below are cumulative.

| Step | Change | Flashing | Injector inlet | Mass flow | Model | Line loss |
|---|---|---|---|---|---|---|
| 0 | Original | **YES** | 53.28 bar | 189.1 g/s | HEM two-phase inlet | 0.72 bar |
| 1 | Needle valve → ball valve (K = 0.05) | **YES** | 53.46 bar | 190.0 g/s | HEM two-phase inlet | 0.54 bar |
| 2 | Raise tank pressure to **60 bar** | NO | 58.09 bar | 365.6 g/s | Dyer | 1.91 bar |
| 3 | Pipe ID 8 mm → **12 mm** | NO | 59.72 bar | 374.6 g/s | Dyer | 0.28 bar |

**Step 1** helps only marginally: the 2.5 m pipe alone (0.45 bar) already consumes the 0.32 bar margin.

**Step 2 is the decisive one.** Raising the regulator to 60 bar increases the subcooling margin from 0.32 bar to 6.32 bar, which the line losses (1.91 bar at 365 g/s) no longer exhaust: the fluid arrives as subcooled liquid (ΔT_sub = +3.64 K at the injector inlet) and the flow more than **doubles** (190 → 366 g/s). (Keeping the needle valve at 60 bar would still avoid flashing: 362 g/s, 2.53 bar line loss.)

**Step 3** is not needed to avoid flashing at 60 bar, but it recovers pressure and margin: line loss falls from 1.91 to 0.28 bar and ΔT_sub at the inlet rises to +4.94 K, giving robustness against a warmer day or a slightly lower tank pressure. Friction in a pipe scales roughly as 1/D⁵ at fixed flow (fitting losses as 1/D⁴), so widening from 8 to 12 mm (×1.5) cuts the pipe loss by a factor of about 7.

**Final design:** 60 bar, 12 mm line, ball valve — no flashing, comfortable margin, **374.6 g/s** (Dyer).

---

## Caveats that the tool also reports

1. **Choking warning still appears on the corrected design.** At 60 bar → 18 bar the Henry-Fauske ceiling is 298.7 g/s (286.5 g/s after Step 2), *below* the Dyer prediction of 374.6 g/s (365.6 g/s). As in Examples 1 and 2, this is a diagnostic, not a correction; the pressure drop (about 42 bar) is far outside the 8–14 bar band where the Dyer model has been validated. Read the result as "roughly 300–375 g/s", to be confirmed by measurement.
2. **The jump from 190 to 366 g/s between Steps 1 and 2 is partly a model artefact.** Dyer (liquid inlet) and HEM (two-phase inlet) are different models, and HEM at vanishing vapour quality predicts about half of what Dyer predicts at the same conditions. The qualitative message — crossing into flashing in the line is very costly — is robust; the exact size of the drop is not validated.

---

## Key lessons

1. **Subcooling margin matters more than absolute pressure.** At a given temperature, the tank pressure fixes the margin: 54 bar at 22 °C is only 0.32 bar above P_sat, while 60 bar at the same temperature gives 6.32 bar. (Equivalently, the same pressure at a lower temperature gives more margin.)
2. **The needle valve is the worst fitting for N₂O systems.** Its high K value generates large local pressure losses. Use ball valves for isolation and flow control.
3. **Pipe diameter has a disproportionate effect.** At fixed flow, pipe friction scales as roughly 1/D⁵ (fittings 1/D⁴): going from 8 to 12 mm cuts pipe friction by about 7×, and doubling the diameter would cut it by about 30×. Undersized tubing is a common and avoidable cause of feed line flashing.
4. **The SPI model cannot predict this.** It has no mechanism to detect or account for flashing — it always produces a number (492 g/s here), even when the real flow is a fraction of that number.

---

## How to reproduce in the tool

**Original (flashing) case:**
1. Open **Sizing mode**.
2. T = 22 °C, P = 54 bar, Cd = 0.65, P_chamber = 18 bar. Orifices: 4 × 1.8 mm.
3. Segments: pipe 2500 mm / 8 mm → needle valve / 8 mm → elbow 90° / 8 mm.
4. Observe: Flashing YES, x_in = 0.006, HEM two-phase-inlet mass flow **189.1 g/s**, and the diagnostic panel with the four suggestions above.

**Corrected case:**
1. Change T = 22 °C, P = **60 bar**.
2. Change pipe segments to **12 mm ID**.
3. Change needle valve to **Ball valve (fully open)**.
4. Observe: Flashing NO, Dyer model applied, mass flow **374.6 g/s**, plus the amber choking-ceiling warning (298.7 g/s).
