# User Manual — N₂O Injector Sizing Tool

This document walks through the interactive Streamlit tool step by step. For the underlying physics, see the theory documents (`01_n2o_thermodynamics.md` through `03b_grain_sizing.md`).

---

## Installation

From the repository root, install the required packages:

```bash
pip install streamlit plotly matplotlib numpy reportlab CoolProp
```

`CoolProp` provides the N₂O equation of state the model itself needs (not just the interface); the tool will not start without it. See the note below if it fails to install.

Then launch the tool:

```bash
streamlit run src/interface/app.py
```

The tool opens automatically in your browser at `http://localhost:8501`. It runs entirely locally — no data is sent anywhere.

---

## A note on the N₂O properties (CoolProp)

As of September 2026 the tool's saturation properties (pressure, densities, enthalpy, entropy) come from **CoolProp**, an accurate, open-source equation-of-state library, replacing an earlier set of approximate correlations. `pip install CoolProp` should work on a normal machine with internet access (pre-built wheels exist for Windows, macOS and Linux). If it fails on your machine, see `docs/references.md`, "CoolProp and alternatives", for what to try instead. Every number quoted elsewhere in this manual, in `validation/`, and in `examples/` predates this change and has not yet been re-confirmed with CoolProp installed — treat them as provisional (see `docs/future_work.md`, Priority 4).

---

## Landing page

The landing page presents two modes and two reference sections:

- **Sizing mode** — you have an existing or proposed orifice geometry and want to predict the real mass flow the system will deliver.
- **Design mode** — you have a target mass flow and want the orifice area required to deliver it, with the Dyer correction applied.
- **Model assumptions** — lists what the model does and does not account for. Read this before interpreting results.
- **Disclaimer** — the tool is for predictive purposes only. See the disclaimer for full details.

---

## What the results are worth

Before using any number from the tool, keep the validation status in mind:

- The Dyer injector model is validated against experiment (Waxman) for injector pressure drops of **8–14 bar**. Typical motor designs use larger drops (20–50 bar), where the tool's output is an **unvalidated model estimate**.
- When the feed line flashes, the injector is evaluated with a HEM two-phase-inlet model that is implemented and unit tested but **not validated** against experimental data.
- The amber **non-equilibrium choking warning** (below) marks operating points where a second, independent model predicts a lower flow than Dyer. Treat the pair of numbers as a plausible range.

A cold-flow or hot-fire measurement should confirm any hardware decision.

---

## Sidebar — common inputs

Both modes share the same sidebar inputs.

### Quick start — presets

The **Load preset configuration** dropdown provides four generic starting-point configurations:

| Preset | Typical use |
|---|---|
| Compact lab motor (low pressure) | Small motors, P_chamber ≈ 12 bar |
| Mid-scale hybrid (moderate pressure) | General university motor, P_chamber ≈ 20 bar |
| High-pressure research motor | Larger motors, P_chamber ≈ 30 bar |
| Near-critical conditions (risk case) | Demonstrates the sensitivity of N₂O near its critical point |

Loading a preset fills all sidebar fields and the feed line segment list automatically. Some presets deliberately start with a tank at or below saturation (for example the compact lab motor: 45 bar at 15 °C, where P_sat = 45.9 bar) to show the flashing behaviour.

### Tank conditions

**Tank temperature (deg C)** — N₂O temperature at the tank outlet (−10 to 35 °C). The tool displays a live badge showing the current subcooling margin (P_tank − P_sat(T_tank)). A green badge means the fluid is safely subcooled; a red badge means the tank is already at or below saturation, and flashing will occur regardless of feed line geometry.

> Practical note: N₂O's saturation pressure at 20 °C is approximately 51.4 bar (correlation used by the tool). A tank at 60 bar and 20 °C has about 8.6 bar of subcooling margin. Between 20 °C and 35 °C the saturation pressure rises to about 70 bar, so keeping a positive margin at warm temperatures requires a supercharged tank.

**Tank pressure (bar)** — 5 to 71 bar; must exceed P_sat(T_tank) for the fluid to remain liquid at the outlet.

### Injector

**Discharge coefficient Cd** — empirical efficiency of the orifice. Typical values:
- Sharp-edged orifice: 0.61
- Well-rounded orifice: 0.82
- Typical machined injector: 0.65–0.70

**Chamber pressure (bar)** — combustion chamber pressure. Sets the pressure drop across the injector. A minimum injector pressure drop of about 15% of chamber pressure is recommended for combustion stability — the tool checks this automatically and flags a warning if the margin is insufficient.

### Combustion stability indicator

Below the chamber pressure field, a badge shows the injector pressure drop (tank pressure minus chamber pressure, before line losses) as a percentage of chamber pressure:

- Green (≥ 20%): stable
- Blue (15–20%): marginal
- Red (< 15%): instability risk — consider raising tank pressure or reducing orifice area

After the calculation, the results section repeats the check with the *actual* injector inlet pressure.

### Advanced — pipe roughness

Absolute roughness of the pipe wall in micrometres. Default 1.5 µm is appropriate for smooth stainless steel tubing. Commercial steel pipe: approximately 46 µm.

---

## Feed line geometry

Both modes include an editable feed line segment list, displayed below the sidebar inputs.

### Adding segments

Use the **+ Pipe** and **+ Fitting** buttons. Insert buttons appear between every pair of existing segments, so you can place a new segment at any position in the line without removing and re-adding components.

### Segment types

**Pipe segment:**
- Length (mm) — physical length of the straight section
- Inner diameter (mm) — internal diameter of the pipe

**Fitting:**
- Type — choose from a dropdown of common fittings with pre-filled K values:

| Fitting | K value |
|---|---|
| Ball valve (fully open) | 0.05 |
| Ball valve (½ open) | 5.5 |
| Needle valve (fully open) | 2.0 |
| Globe valve (fully open) | 10.0 |
| Check valve | 2.5 |
| Elbow 90° (standard) | 0.9 |
| Elbow 90° (long radius) | 0.4 |
| Elbow 45° | 0.4 |
| Tee (flow-through) | 0.4 |
| Tee (branch) | 1.5 |
| Union / coupling | 0.04 |
| Custom | enter manually |

- Inner diameter (mm) — diameter of the adjacent pipe at that fitting

### Line schematic

Below the segment list, a live schematic shows the feed line layout: pipes as horizontal rectangles (proportional to length), fittings as diamond symbols, with tank and injector at the ends. This updates immediately as you add or modify segments.

---

## Sizing mode

### Orifice geometry

- **Number of orifices** — total number of identical holes in the injector plate.
- **Orifice diameter (mm)** — diameter of each individual hole. The total orifice area (N × π(d/2)²) is shown as a badge and updates live.

### Initial mass-flow guess for the coupled solver

This field is only the **starting point** of the iteration. The model iterates automatically until the feed-line losses and injector flow are mutually consistent, so the converged result does not depend on it (it only affects convergence speed). Set it to your design target mass flow as a reasonable starting value.

### Results

After entering all inputs, results appear immediately below the feed line:

| Indicator | Meaning |
|---|---|
| **Flashing in line — NO** | Feed line stays subcooled liquid throughout; SPI or Dyer is used at the injector |
| **Flashing in line — YES** | Pressure dropped below P_sat somewhere in the line; the injector is evaluated with the HEM two-phase-inlet model |
| **Inlet vapour quality** | Vapour fraction at the injector inlet (shown when flashing occurs) |
| **Exit vapour quality** | Vapour fraction at the orifice exit (flashing case) |
| **SPI sufficient — YES** | Single-phase model valid; no two-phase correction needed |
| **SPI sufficient — NO** | Dyer model used; two-phase correction applied |
| **Real mass flow** | Model-selected prediction in g/s |
| **SPI over-prediction** | How much the SPI model would have over-predicted (%) |

When Dyer is used, a caption shows kappa (non-equilibrium parameter), exit vapour quality x, and the HEM prediction for reference.

### Non-equilibrium choking warning

Directly beneath that caption, an amber **"Non-equilibrium choking ceiling exceeded"** badge and warning box appear whenever the Dyer prediction exceeds the Henry-Fauske (1971) non-equilibrium choking ceiling for the current tank/chamber conditions. This is a **diagnostic warning only** — it does not change the displayed "Real mass flow" figure. The ceiling is theoretically sound (a primary-source non-equilibrium critical-flow model) but has only been confirmed not to interfere with the tool's validated reference case (Waxman 2013/2014, 8–14 bar pressure drop); outside that band it is an unconfirmed, conservative alternative estimate, not a certainty. If you see this warning, treat it as a prompt to think carefully about the operating point (and, ideally, to seek experimental confirmation) rather than as a correction to apply by hand. See `docs/future_work.md`, Priority 1, for the full reasoning.

If the ceiling could not be computed at all (tank temperature above ≈307 K, near the critical point, where the required entropy data is unavailable), a grey caption states this instead of the warning.

### Diagrams

Four interactive Plotly charts appear:

1. **Pressure along the feed line** — pressure at each segment, with the P_sat(T_tank) threshold as a dashed red line. Use scroll to zoom, click-drag to pan.
2. **P-T diagram** — saturation curve with liquid/vapour regions shaded, and the three operating points (tank, injector inlet, chamber) marked. The chamber point is plotted at T_sat(P_chamber), the saturation temperature at chamber pressure.
3. **Subcooling margin along the line** — ΔT_sub at each segment. Green fill = margin remaining; red fill = flashing zone.
4. **Sensitivity analysis** — how the real mass flow varies as tank temperature changes ±5 °C from the design point. The design point is marked with a star. Any temperatures where flashing is predicted are highlighted as a red zone.

A **sensitivity tornado** expander is available below the charts: each input is varied one-at-a-time (±5 K temperature, ±5% pressure, ±0.05 Cd, ±0.05 mm orifice diameter, ±20% line length) to show which input most affects the result.

A fifth chart is available in the **Pressure drop by segment** expander: a horizontal bar chart showing which segments consume the most pressure.

### PDF export

The **Download PDF report** button generates a one-page A4 summary including all inputs, feed line geometry, results, and the two main charts. When the non-equilibrium choking ceiling is exceeded, the results table and a short note both record this, alongside the ceiling value. Suitable for attaching to a project report. (The grain sizing panel is not part of the PDF.)

---

## Design mode

### Target inputs

- **Target oxidiser mass flow (g/s)** — the mass flow the system must deliver at the design O/F ratio and chamber pressure.
- **Number of orifices** — the total number of holes. The tool computes the required diameter per hole for both the SPI and Dyer predictions.

The feed line is evaluated at the target flow (the flow the designed injector is meant to deliver).

### Results

| Output | Description |
|---|---|
| **SPI area** | Total orifice area computed from the SPI model (under-predicts required area) |
| **SPI hole diameter** | Diameter per hole for the SPI area |
| **Dyer area** | Total orifice area with Dyer correction (larger than SPI) |
| **Dyer hole diameter** | Diameter per hole for the Dyer area — this is the value to manufacture |

A caption below the cards shows the percentage by which the Dyer area exceeds the SPI area, explaining why the correction is necessary.

> **Why Dyer area > SPI area?** The Dyer model predicts a lower mass flow than SPI for the same geometry (two-phase effects reduce effective density). To hit the target mass flow, the orifice area must be increased relative to what SPI would specify. If you size with SPI only, the system will under-deliver.

**When the chamber pressure is at or above P_sat(T_tank)**, the flow stays single-phase through the orifice: the cards are labelled "Recommended (SPI)", the recommended area equals the SPI area, and an information box explains why no Dyer correction is applied. The model-comparison chart is replaced by the pressure-along-line chart.

### Non-equilibrium choking warning

Immediately below that caption, the same amber choking warning described for Sizing mode can appear here too — checked at the recommended area and the tank/chamber conditions, before the combustion-stability check. **Important:** this condition is independent of the computed orifice area (both the Dyer prediction and the choking ceiling scale linearly with area), so if it appears, resizing the orifice will not resolve it — the message suggests adjusting tank or chamber pressure instead. As in Sizing mode, this is a diagnostic warning, not a correction applied to the area shown above.

### Diagrams

1. **Model comparison bar chart** — SPI, HEM, Dyer, and target mass flow side by side for the recommended area. When available, the Henry-Fauske ceiling is drawn as an additional dash-dot reference line (orange normally, red if exceeded).
2. **P-T diagram** — same as Sizing mode.
3. **Sensitivity analysis** — for the recommended area, how mass flow varies with tank temperature.
4. **Subcooling margin** — along the feed line.

---

## Grain sizing (fuel side)

Both modes end with a collapsed **Grain sizing (fuel side)** panel. It sizes the *initial* fuel grain for a target O/F ratio from the oxidiser mass flow computed above (the real mass flow in Sizing mode, the target in Design mode). It does **not** simulate the transient burn — port radius, O/F and thrust all drift as the grain regresses.

Inputs:

| Input | Notes |
|---|---|
| **Target O/F ratio** | From your thermochemical sizing (CEA, RPA, …); not computed by this tool |
| **Fuel** | Paraffin wax, HTPB, ABS or PMMA — fills in the fuel **density** only, which you may override |
| **Grain length (mm)** | A direct input, usually fixed by your motor case |
| **Number of ports** | Circular ports only; more ports add burning perimeter (∝ √N at fixed total port area) |
| **Burn duration (s)** | Optional; adds a first-order, conservative burnback estimate |
| **Regression-rate data point** | **Required** — three numbers: a regression rate in mm/s, the oxidiser mass flux G_o (kg/(m²·s)) it was measured at, and the exponent n |

**Why a data point and not the coefficient a?** The Marxman coefficient's units depend on n, and literature values mix unit systems, so entering `a` by hand is the easiest way to get a nonsense result. Instead, enter one regression-rate point exactly as it appears in your test data or in a paper's plot/table, and the tool converts it. The fuel's literature note is shown as orientation only: published rates for the same fuel/oxidiser pair differ by 2–3× between studies, so use your own data or a source you have vetted for your scale and injector.

Results: total fuel mass flow, initial port radius (per port), initial G_o and initial regression rate, plus the burnback estimate when a burn time is given. If the solved port radius is outside 5–300 mm the tool shows an error instead of a result, with the regression rate your inputs imply; this almost always means a unit problem in the data point or a target O/F, length or port count that does not fit.

---

## When flashing is detected in the feed line

If the feed line pressure drops below P_sat(T_tank) before reaching the injector, both modes show a **diagnostic panel** identifying the likely causes and suggesting specific corrective actions:

- **Subcooling margin below 5 bar** (or tank pressure at/below P_sat) — raise tank pressure or lower the tank temperature.
- **Line longer than 1.5 m** — shorten it.
- **Pipe inner diameter below 10 mm** — widen it (friction scales roughly as 1/D⁵ in pipe, 1/D⁴ in fittings).
- **High-K fittings** (K > 1) — replace with fully open ball valves.
- **Tank above 25 °C** — pre-cool the oxidiser.

What else happens depends on the mode:

- **Sizing mode** — the injector is evaluated with the **HEM two-phase-inlet** model, and the result cards show the inlet and exit vapour quality and the resulting mass flow (labelled "HEM two-phase inlet"). This estimate is not validated against experimental data, and the model switches abruptly between Dyer (no flashing) and HEM (flashing), so the size of the drop at the flashing threshold is uncertain. A naive SPI reference value is shown for comparison.
- **Design mode** — area sizing is **not** offered; fix the line first.

The most common fix is to increase tank pressure so that P_tank − P_sat(T_tank) > 5 bar before line losses, providing a safety margin. See `examples/example_03_flashing.md`.

---

## Worked examples

See the `examples/` folder for three complete scenarios, including all inputs and interpretation of results:

- [`example_01_sizing.md`](../examples/example_01_sizing.md) — Predict mass flow for a known injector geometry
- [`example_02_design.md`](../examples/example_02_design.md) — Design an injector for a target mass flow
- [`example_03_flashing.md`](../examples/example_03_flashing.md) — Diagnose and correct a flashing feed line
