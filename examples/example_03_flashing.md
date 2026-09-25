# Example 3 — Flashing Case: Diagnosis and Design Correction

> **Property-backend status: CONFIRMED (September 2026).** Every figure below
> was regenerated with the real CoolProp package installed, by running
> `tests/generate_example_03.py` (same geometry, same solver settings as the
> tool). It replaces the earlier version of this example, which mixed two
> CoolProp headline numbers with legacy Perry/McGill cross-check values. Two
> statements of that earlier version turned out to be wrong and were corrected
> here: (i) Step 1 (ball valve at 53.5 bar) does **not** give "roughly 188 g/s",
> the solver does not converge there; (ii) the 344.4 g/s figure belongs to the
> *needle-valve* geometry at 56 bar, not to the ball-valve design.

This example shows how to use the tool to **diagnose a flashing feed line**, understand the cause, and apply corrections until the system operates safely. It illustrates the physical mechanism behind unexpected results in several documented university hybrid motor tests.

---

## Scenario

A team designs a motor with N₂O at **22 °C**, tank pressure **53.5 bar**, and the following feed line:

- 2500 mm of 8 mm ID tubing
- Needle valve (fully open, K = 2.0)
- 90° elbow (standard, K = 0.9)
- Injector: 4 holes, 1.8 mm diameter (Cd = 0.65), total area 10.18 mm²
- Chamber pressure: 18 bar

The SPI model was used for the initial design: with the full tank-to-chamber pressure drop it predicts **489.0 g/s**. On test, the motor produces less thrust than expected and combustion is unstable.

---

## Initial configuration — what goes wrong

| Parameter | Value |
|---|---|
| Tank temperature | 22 °C |
| Tank pressure | 53.5 bar |
| P_sat at 22 °C | 52.859 bar |
| Subcooling margin at the tank exit | **0.641 bar** |

The margin is under 1 bar, so it does not survive the line.

**Feed line trace (at the converged flow):**

| Segment | Type | ΔP (bar) | P after (bar) | x at segment start | ρ_eff (kg/m³) |
|---|---|---|---|---|---|
| Tank exit | — | — | 53.500 | — | — |
| Pipe, 2.5 m | Pipe | 0.468 | 53.032 | 0.0000 | 769.2 |
| Needle valve | Fitting (K = 2.0) | 0.193 | 52.839 | 0.0000 | 769.2 |
| Elbow | Fitting (K = 0.9) | 0.087 | 52.752 | 0.0003 | 768.4 |

The pipe alone uses 0.468 bar of the 0.641 bar margin and leaves 0.173 bar. The needle valve then drops 0.193 bar, which is more than what is left, so the pressure falls below P_sat(22 °C) = 52.859 bar **at the needle valve**, and the elbow already runs (barely) two-phase. The vapour quality at the injector inlet is x_in = 0.0016, small in the line but with a large effect at the orifice.

**Result:** flashing is detected, the tool applies the **HEM two-phase-inlet** model, and the solver converges (13 iterations) to a real mass flow of **193.9 g/s**, about 2.5× below the naive SPI figure (489.0 g/s). The vapour quality at the orifice exit is x = 0.348. This two-phase-inlet result is implemented and unit tested but has not been validated against experimental data (see "Caveats").

---

## Diagnosis

The tool's diagnostic panel identifies the contributing factors:

1. **Subcooling margin too small at the tank exit** — 0.64 bar; any friction loss at all is nearly enough to cause flashing.
2. **Long feed line** — 2.5 m total; the panel targets under 1.5 m.
3. **Small pipe diameter** — 8 mm ID; the panel suggests 10 mm or more.
4. **High-K fitting** — the needle valve (K = 2.0) is the segment that actually tips the fluid across saturation. A fully open ball valve (K = 0.05) is 40 times lower.

---

## Correction — step by step

The changes are cumulative. All rows are CoolProp results.

| Step | Change | Flashing | Mass flow | Model | Notes |
|---|---|---|---|---|---|
| 0 | Original (53.5 bar, needle valve, 8 mm) | YES | **193.9 g/s** | HEM two-phase inlet | converged, 13 iterations |
| 1 | Needle valve → ball valve (K = 0.05) | — | **no convergence** | — | solver oscillates, see below |
| 2 | Raise tank pressure to **56 bar** | NO | **347.7 g/s** | Dyer | 6 iterations; κ = 1.020, x = 0.347 |
| 3 | Line ID 8 mm → **12 mm** (still 56 bar) | NO | **356.2 g/s** | Dyer | 8 iterations; κ = 1.041, x = 0.347 |

**Step 1 is the instructive failure.** Swapping the valve removes most of the local loss, and that is exactly what puts the operating point *on* the flashing threshold. The coupled solver then cycles instead of settling: its iterates repeat the sequence 271.55, 231.90, 212.77, 203.50 g/s indefinitely. Working back through the damped update (α = 0.5), the two targets it alternates between are about 340 g/s (Dyer branch, no flashing) and about 193 g/s (HEM branch, flashing). At the higher flow the line loses enough pressure to flash, which sends the model to the low-flow HEM branch, where the line loses less pressure and no longer flashes, which sends it back. Neither branch is self-consistent, so the tool raises a `RuntimeError` ("Coupled solver did not converge…", shown as "Model error" in the interface) instead of returning a number. That is the correct behaviour, and it is a signal that the point is too close to the threshold to be predicted by this model.

**Step 2 is the decisive one.** At 56 bar the tank-exit margin is 3.141 bar. The line now loses 1.733 bar at the converged flow, the fluid reaches the injector still liquid at 54.267 bar (1.41 bar above P_sat), and the flow rises from 193.9 to 347.7 g/s.

**Step 3** is not needed to avoid flashing at 56 bar, but it cuts the line loss from 1.733 bar to 0.257 bar and raises the margin at the injector inlet from 1.41 bar to 2.88 bar (inlet at 55.743 bar). That gives robustness against a warmer day or a slightly lower tank pressure. The flow rises a further 2.4 %, which must be accounted for in the injector sizing.

**Final design:** 56 bar, ball valve, 12 mm line, no flashing, **356.2 g/s** (Dyer). The minimum change that removes flashing is Step 2 (347.7 g/s).

---

## Where the solver does and does not converge

Because the switch from Dyer to HEM at the flashing threshold is discontinuous, there is a band of tank pressures where the coupled solver cannot settle. Scan for the **original geometry** (needle valve, 8 mm line), same chamber and injector:

| Tank pressure (bar) | Result |
|---|---|
| 52.00 – 53.50 | converges, HEM two-phase inlet: 186.7 → 193.9 g/s |
| 53.75 – 55.00 | **no convergence** |
| 55.25 – 58.00 | converges, Dyer: 340.3 → 355.1 g/s |

The band lies between 53.50 bar (converges, HEM) and 55.25 bar (converges, Dyer). Its position depends on the line geometry, as Step 1 shows: the same 53.5 bar that converges with the needle valve does not converge with the ball valve. The starting point of this example (53.5 bar) sits right at the lower edge of the band for the original geometry, so a slightly higher tank pressure would have made the tool refuse to answer. The scan also gives a consistency check on the Step 2 value: at 56 bar with the needle valve the flow is 344.4 g/s, slightly below the 347.7 g/s obtained with the ball valve, as expected from the lower line loss.

Across the band the flow changes from 193.9 g/s (53.5 bar) to 340.3 g/s (55.25 bar), a factor of 1.75 for 1.75 bar of tank pressure.

---

## Caveats that the tool also reports

1. **The choking warning appears on the corrected design.** The Henry-Fauske (1971) non-equilibrium ceiling is **267.8 g/s** at Step 2 (Dyer 347.7 g/s, `choked = True`, ceiling 23 % below Dyer) and **278.1 g/s** at Step 3 (Dyer 356.2 g/s, ceiling 22 % below). This is the same pattern as in Examples 1 and 2 at comparably large pressure drops (36–38 bar here). It is a diagnostic, not a correction: the Dyer model is validated only at 8–14 bar, and the ceiling itself has no experimental confirmation where it binds. A defensible planning range for the final design is roughly **278–356 g/s**, to be narrowed by a cold-flow or hot-fire measurement. In either case it is far below the ≈ 506 g/s that SPI predicts at 56 bar.
2. **The jump between Steps 0 and 2 is partly a model artefact.** Dyer (liquid inlet) and HEM (two-phase inlet) are different models, and HEM at vanishing vapour quality predicts noticeably less than Dyer at the same conditions (`docs/future_work.md`, Priority 2). The qualitative message, that crossing into flashing in the line is very costly, is robust. The exact size of the jump, and the behaviour inside the non-convergent band, are not pinned down by any experimental data.
3. **The tool's own margin guideline is 5 bar before line losses.** 56 bar gives 3.14 bar at 22 °C. If more margin is wanted, the scan above shows a converging Dyer solution at 58 bar (5.14 bar margin) for the original needle-valve geometry, 355.1 g/s; re-run the final geometry in the tool to get its own figure.

---

## Key lessons

1. **Subcooling margin matters more than absolute pressure.** At a given temperature the tank pressure fixes the margin: 53.5 bar at 22 °C leaves 0.64 bar (used up by the first pipe and the needle valve), while 56 bar restores 3.14 bar before line losses.
2. **The needle valve is the worst fitting for N₂O systems.** In the original line it is the element that tips the fluid across saturation (0.193 bar against 0.173 bar of margin left). Use ball valves for isolation and flow control.
3. **Fixing one thing near the threshold can make the model unusable, not better.** Step 1 removes a loss and lands on the flashing threshold, where the model has no consistent answer. Move away from the threshold (more tank pressure, wider line) rather than sitting just above it.
4. **Pipe diameter has a disproportionate effect.** At fixed flow, pipe friction scales as roughly $1/D^5$ and fittings as $1/D^4$. Going from 8 to 12 mm cut the pipe loss from 1.437 bar to 0.195 bar in Step 3.
5. **The SPI model cannot predict any of this.** It has no mechanism to detect or account for flashing: it always produces a number (489.0 g/s here), even when the real flow is 194 g/s.

---

## How to reproduce in the tool

**Original (flashing) case:**
1. Open **Sizing mode**.
2. T = 22 °C, P = 53.5 bar, Cd = 0.65, P_chamber = 18 bar. Orifices: 4 × 1.8 mm.
3. Segments: pipe 2500 mm / 8 mm, needle valve / 8 mm, elbow 90° / 8 mm.
4. Observe: Flashing YES, x_in = 0.0016, HEM two-phase-inlet mass flow **193.9 g/s**, and the diagnostic panel with the four suggestions above.

**Step 1 (no convergence):** change the needle valve to **Ball valve (fully open)**, keep 53.5 bar. The tool shows a "Model error: Coupled solver did not converge" message.

**Corrected case:**
1. Set P = **56 bar**, ball valve, 8 mm line: Flashing NO, Dyer, **347.7 g/s**, plus the amber choking-ceiling warning (267.8 g/s).
2. Change the line ID to **12 mm**: **356.2 g/s**, ceiling 278.1 g/s.

The scripted version of all of the above is `tests/generate_example_03.py`.
