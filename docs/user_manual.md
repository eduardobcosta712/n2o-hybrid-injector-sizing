# User Manual — N₂O Injector Sizing Tool

This document walks through the interactive Streamlit tool step by step. For the underlying physics, see the theory documents (`01_n2o_thermodynamics.md` through `03_two_phase_flow.md`).

---

## Installation

From the repository root, install the required packages:

```bash
pip install streamlit plotly matplotlib numpy reportlab
```

Then launch the tool:

```bash
streamlit run src/interface/app.py
```

The tool opens automatically in your browser at `http://localhost:8501`. It runs entirely locally — no data is sent anywhere.

---

## Landing page

The landing page presents two modes and three reference sections:

- **Sizing mode** — you have an existing or proposed orifice geometry and want to predict the real mass flow the system will deliver.
- **Design mode** — you have a target mass flow and want the orifice area required to deliver it, with the Dyer correction applied.
- **Model assumptions** — lists what the model does and does not account for. Read this before interpreting results.
- **Disclaimer** — the tool is for predictive purposes only. See the disclaimer for full details.

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

Loading a preset fills all sidebar fields and the feed line segment list automatically.

### Tank conditions

**Tank temperature (deg C)** — N₂O temperature at the tank outlet. The tool displays a live badge showing the current subcooling margin (P_tank − P_sat(T_tank)). A green badge means the fluid is safely subcooled; a red badge means the tank is already at or below saturation, and flashing will occur regardless of feed line geometry.

> Practical note: N₂O's saturation pressure at 20 °C is approximately 51.4 bar. A tank at 60 bar and 20 °C has about 8.6 bar of subcooling margin.

**Tank pressure (bar)** — must exceed P_sat(T_tank) for the fluid to remain liquid at the tank outlet.

### Injector

**Discharge coefficient Cd** — empirical efficiency of the orifice. Typical values:
- Sharp-edged orifice: 0.61
- Well-rounded orifice: 0.82
- Typical machined injector: 0.65–0.70

**Chamber pressure (bar)** — combustion chamber pressure. Sets the pressure drop across the injector. A minimum injector pressure drop of 15% of chamber pressure is recommended for combustion stability — the tool checks this automatically and flags a warning if the margin is insufficient.

### Combustion stability indicator

Below the chamber pressure field, a badge shows the current injector pressure drop as a percentage of chamber pressure:

- Green (≥ 20%): stable
- Blue (15–20%): marginal
- Red (< 15%): instability risk — consider raising tank pressure or reducing orifice area

### Advanced — pipe roughness

Absolute roughness of the pipe wall in micrometres. Default 1.5 µm is appropriate for smooth stainless steel tubing. Commercial steel pipe: approximately 46 µm.

---

## Feed line geometry

Both modes include an editable feed line segment list, displayed below the sidebar inputs.

### Adding segments

Use the **Add pipe segment** and **Add fitting** buttons. Insert buttons appear between every pair of existing segments, so you can place a new segment at any position in the line without removing and re-adding components.

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

### Mass flow for line evaluation

This field sets the mass flow used to compute the feed line velocity and friction losses. It does not directly affect the injector result — the injector model uses the orifice geometry and the pressure conditions at the injector inlet, not this value. Set it to your design target mass flow as a reasonable approximation.

### Results

After entering all inputs, results appear immediately below the feed line:

| Indicator | Meaning |
|---|---|
| **Flashing in line — NO** | Feed line stays subcooled liquid throughout; injector models are valid |
| **Flashing in line — YES** | Pressure dropped below P_sat somewhere in the line; injector models not evaluated |
| **SPI sufficient — YES** | Single-phase model valid; no two-phase correction needed |
| **SPI sufficient — NO** | Dyer model used; two-phase correction applied |
| **Real mass flow** | Model-selected prediction in g/s |
| **SPI over-prediction** | How much the SPI model would have over-predicted (%) |

When Dyer is used, a caption shows kappa (non-equilibrium parameter), exit vapour quality x, and the HEM prediction for reference.

### Diagrams

Four interactive Plotly charts appear:

1. **Pressure along the feed line** — pressure at each segment, with the P_sat(T_tank) threshold as a dashed red line. The green shading shows the subcooling margin. Use scroll to zoom, click-drag to pan.
2. **P-T diagram** — saturation curve with liquid/vapour regions shaded, and the three operating points (tank, injector inlet, chamber) marked. The chamber point is plotted at T_sat(P_chamber), the saturation temperature at chamber pressure.
3. **Subcooling margin along the line** — delta T_sub at each segment. Green fill = margin remaining; red fill = flashing zone.
4. **Sensitivity analysis** — how the real mass flow varies as tank temperature changes ±5 °C from the design point. The design point is marked with a star. Any temperatures where flashing is predicted are highlighted as a red zone.

A fifth chart is available in the **Pressure drop by segment** expander: a horizontal bar chart showing which segments consume the most pressure.

### PDF export

The **Download PDF report** button generates a one-page A4 summary including all inputs, feed line geometry, results, and the two main charts. Suitable for attaching to a project report.

---

## Design mode

### Target inputs

- **Target oxidiser mass flow (g/s)** — the mass flow the system must deliver at the design O/F ratio and chamber pressure.
- **Number of orifices** — the total number of holes. The tool computes the required diameter per hole for both the SPI and Dyer predictions.

### Results

| Output | Description |
|---|---|
| **SPI area** | Total orifice area computed from the SPI model (under-predicts required area) |
| **SPI hole diameter** | Diameter per hole for the SPI area |
| **Dyer area** | Total orifice area with Dyer correction (larger than SPI) |
| **Dyer hole diameter** | Diameter per hole for the Dyer area — this is the value to manufacture |

A caption below the cards shows the percentage by which the Dyer area exceeds the SPI area, explaining why the correction is necessary.

> **Why Dyer area > SPI area?** The Dyer model predicts a lower mass flow than SPI for the same geometry (two-phase effects reduce effective density). To hit the target mass flow, the orifice area must be increased relative to what SPI would specify. If you size with SPI only, the system will under-deliver.

### Diagrams

1. **Model comparison bar chart** — SPI, HEM, Dyer, and target mass flow side by side for the Dyer area, showing the spread between model predictions.
2. **P-T diagram** — same as Sizing mode.
3. **Sensitivity analysis** — for the Dyer area, how mass flow varies with tank temperature.
4. **Subcooling margin** — along the feed line.

---

## When flashing is detected in the feed line

If the feed line pressure drops below P_sat(T_tank) before reaching the injector, the injector models are not evaluated. Instead, the tool shows:

- A **diagnostic panel** identifying the likely cause (tank pressure too low, line too long, pipe diameter too small, high-K fittings) and suggesting specific corrective actions.
- A **reference SPI value** — what the SPI model would predict if the fluid were single-phase throughout. This is an upper bound; the real mass flow will be substantially lower due to two-phase flow in the line.

The most common fix is to increase tank pressure so that P_tank − P_sat(T_tank) > 5 bar before line losses, providing a safety margin.

---

## Worked examples

See the `examples/` folder for three complete sizing scenarios, including all inputs and interpretation of results:

- [`example_01_sizing.md`](../examples/example_01_sizing.md) — Predict mass flow for a known injector geometry
- [`example_02_design.md`](../examples/example_02_design.md) — Design an injector for a target mass flow
- [`example_03_flashing.md`](../examples/example_03_flashing.md) — Diagnose and correct a flashing feed line
