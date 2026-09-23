# Example 3 -- Flashing Case: Diagnosis and Design Correction

> **Property-backend status.** The two headline mass-flow figures in this example (193.9 g/s and 344.4 g/s) are the confirmed CoolProp results for this exact geometry, as reported in `docs/future_work.md`, Priority 2. The starting tank pressure used here (53.5 bar, instead of the original 54 bar) was chosen specifically to sit outside a narrow non-convergence band -- roughly 54.0-55.5 bar for this geometry -- that the coupled solver runs into right at the flashing threshold once the real (lower) CoolProp saturation pressure is used; see "A solver note" below. The supporting detail in this example (the segment-by-segment trace, kappa, exit quality, and the Henry-Fauske ceiling at the corrected point) has not yet been regenerated with CoolProp and is shown as a legacy-backend (Perry/McGill + NIST) cross-check, clearly marked wherever it appears -- treat those specific numbers as indicative pending a full CoolProp re-run.

This example shows how to use the tool to **diagnose a flashing feed line**, understand the cause, and apply corrections until the system operates safely. It illustrates the physical mechanism that caused unexpected results in several documented university hybrid motor tests.

---

## Scenario

A team designs a motor with N2O at **22 degC**, tank pressure **53.5 bar**, and the following feed line:

- 2500 mm of 8 mm ID tubing
- Needle valve (fully open, K = 2.0)
- 90 degree elbow (standard, K = 0.9)
- Injector: 4 holes, 1.8 mm diameter (Cd = 0.65), total area 10.18 mm$^2$
- Chamber pressure: 18 bar

The SPI model was used for the initial design. A naive SPI calculation with the full tank-to-chamber drop, at the legacy-backend liquid density for 22 degC, gives roughly 489 g/s -- this figure moves slightly under CoolProp's slightly different liquid density, but not enough to change the conclusion below. On test, the motor produces less thrust than expected and combustion is unstable.

---

## Initial configuration -- what goes wrong

| Parameter | Value |
|---|---|
| Tank temperature | 22 degC |
| Tank pressure | 53.5 bar |
| Subcooling margin at the tank exit | small -- on the order of a few tenths of a bar, close enough to zero that line friction alone is enough to cross saturation |

With so little margin, even the first metre of tubing is enough to push the fluid below its saturation pressure.

**Feed line trace:**

| Segment | Type | P after (bar) | Vapour quality x at segment start |
|---|---|---|---|
| Tank exit | -- | 53.50 | -- |
| After pipe (2.5 m) | Pipe | 53.06 | 0.003 |
| After needle valve | Fitting (K=2.0) | 52.87 | 0.009 |
| After elbow | Fitting (K=0.9) | 52.79 | 0.012 |

**Result -- CONFIRMED (CoolProp):** flashing is detected in the line, the tool applies the **HEM two-phase-inlet** model, and the converged real mass flow is **193.9 g/s** -- roughly 2.5x below the naive SPI figure. (This two-phase-inlet result is implemented and unit tested but has not been validated against experimental data; see "Caveats" below.)

---

## Diagnosis

The tool's diagnostic panel identifies the contributing factors:

1. **Subcooling margin too small at the tank exit** -- any friction loss at all is enough to cause flashing.
2. **Long feed line** -- 2.5 m total; the panel targets under 1.5 m.
3. **Small pipe diameter** -- 8 mm ID; the panel suggests 10 mm or more.
4. **High-K fitting** -- the needle valve (K = 2.0). A fully open ball valve (K = 0.05) is 40 times lower.

---

## Correction -- step by step

The changes below are cumulative. The mass-flow figures for Steps 0 and 2 are the CoolProp-confirmed headline numbers; the other columns and Steps 1 and 3 are the legacy-backend cross-check, kept for the qualitative trend until they too are regenerated.

| Step | Change | Flashing | Mass flow | Model | Status |
|---|---|---|---|---|---|
| 0 | Original (53.5 bar, needle valve, 8 mm) | YES | **193.9 g/s** | HEM two-phase inlet | CoolProp-confirmed |
| 1 | Needle valve to ball valve (K = 0.05) | YES | roughly 188 g/s | HEM two-phase inlet | legacy-backend cross-check |
| 2 | Raise tank pressure to **56 bar** | NO | **344.4 g/s** | Dyer | CoolProp-confirmed |
| 3 | Pipe ID 8 mm to **12 mm** (still at 56 bar) | NO | roughly 352 g/s | Dyer | legacy-backend cross-check |

**Step 1** helps only marginally: the 2.5 m pipe alone already consumes almost all of the available margin, so swapping the valve barely moves the needle.

**Step 2 is the decisive one.** Raising the regulator to 56 bar restores a positive subcooling margin that the line losses no longer exhaust: the fluid arrives at the injector still liquid, and the flow **nearly doubles** relative to Step 0 (194 to 344 g/s).

**Step 3** is not needed to avoid flashing at 56 bar, but it recovers pressure margin and gives robustness against a warmer day or a slightly lower tank pressure, at the cost of a small further increase in flow that then needs to be accounted for in the injector sizing.

**Final design:** 56 bar, needle valve replaced by a ball valve, no flashing, **344.4 g/s** (Dyer, CoolProp-confirmed) as the number to design the injector around.

---

## Caveats that the tool also reports

1. **Choking warning still appears on the corrected design.** At the legacy-backend cross-check, the Henry-Fauske ceiling at 56 bar to 18 bar comes out below the Dyer prediction (`choked = True`), consistent with the same pattern seen in Examples 1 and 2 at comparably large pressure drops. As in those examples, this is a diagnostic, not a correction: the pressure drop here (about 38 bar) is far outside the 8-14 bar band where the Dyer model has been validated against experiment. Treat the 344.4 g/s figure as a model estimate pending a cold-flow or hot-fire measurement, and expect the exact ceiling value to move once it too is regenerated with CoolProp.
2. **The jump between Steps 0 and 2 is partly a model artefact, not purely a physical effect.** Dyer (liquid inlet) and HEM (two-phase inlet) are different models, and HEM at vanishing vapour quality predicts noticeably less than Dyer at the same conditions -- the already-documented Dyer-to-HEM discontinuity (`docs/future_work.md`, Priority 2). The qualitative message -- crossing into flashing in the line is very costly -- is robust; the exact size of the jump is not fully pinned down.

---

## A solver note

The original version of this example used a 54 bar starting point. With the real CoolProp saturation pressure (lower, at 22 degC, than the older Perry/McGill correlation predicted), that specific starting point falls inside a narrow tank-pressure band -- roughly 54.0-55.5 bar for this exact geometry -- where the coupled solver in `full_system.py` alternates between the Dyer and HEM branches instead of settling, and correctly raises `RuntimeError` rather than returning an unconverged answer. The two operating points used in this example, 53.5 bar and 56 bar, sit cleanly on either side of that band, which is why they were chosen for the regenerated numbers above (see `docs/future_work.md`, Priority 2, and `docs/04_implementation.md`, Section 4.5, for the full account). The supporting detail (Step 1, Step 3, the exact trace, kappa, exit quality, and the choking ceiling) still needs a direct CoolProp run to replace the legacy-backend cross-check values shown here.

---

## Key lessons

1. **Subcooling margin matters more than absolute pressure.** At a given temperature, the tank pressure fixes the margin: 53.5 bar at 22 degC leaves almost none, while 56 bar at the same temperature restores a workable margin before line losses.
2. **The needle valve is the worst fitting for N2O systems.** Its high K value generates large local pressure losses. Use ball valves for isolation and flow control.
3. **Pipe diameter has a disproportionate effect.** At fixed flow, pipe friction scales as roughly $1/D^5$ (fittings as $1/D^4$): widening the line is one of the most effective single changes available.
4. **The SPI model cannot predict this.** It has no mechanism to detect or account for flashing -- it always produces a number, even when the real flow is a fraction of that number.

---

## How to reproduce in the tool

**Original (flashing) case:**
1. Open **Sizing mode**.
2. T = 22 degC, P = 53.5 bar, Cd = 0.65, P_chamber = 18 bar. Orifices: 4 by 1.8 mm.
3. Segments: pipe 2500 mm / 8 mm, needle valve / 8 mm, elbow 90 degrees / 8 mm.
4. Observe: Flashing YES, HEM two-phase-inlet mass flow **193.9 g/s**, and the diagnostic panel with the four suggestions above.

**Corrected case:**
1. Change T = 22 degC, P = **56 bar**.
2. Change the needle valve to **Ball valve (fully open)**.
3. Observe: Flashing NO, Dyer model applied, mass flow **344.4 g/s**, plus the amber choking-ceiling warning.
