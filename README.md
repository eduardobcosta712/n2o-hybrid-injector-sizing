# n2o-hybrid-injector-sizing

A physics-based tool and technical reference for predicting premature vaporisation (*flashing*) of N₂O in the feed systems of paraffin/N₂O hybrid rocket motors, and for sizing injector orifices that account for two-phase flow effects.

> **Disclaimer:** This tool is provided for predictive and educational purposes only. It is an academic project under active development and has not been independently validated against a comprehensive experimental dataset. Results must not be used as the sole basis for engineering decisions or hardware fabrication. The author accepts no responsibility for any damages arising from the use of this tool.

---

## The problem

N₂O is the most common oxidiser in university-built hybrid motors — self-pressurising, relatively safe to handle, and commercially available. It has, however, a critical point at only ≈36.4 °C, dangerously close to typical ambient launch temperatures. This makes it exceptionally prone to partial vaporisation before it reaches the combustion chamber, whenever pressure in the feed system drops due to line friction, fittings, or acceleration through the injector orifice itself.

When this happens, the model most commonly used by university teams — SPI (*Single Phase Incompressible*, Bernoulli's equation with an empirical discharge coefficient) — is no longer valid: it implicitly assumes the fluid stays liquid from tank to chamber. Documented results in the literature show real mass flow far below the SPI prediction, by a factor of several times, shifting the motor's real O/F far from the design point and directly affecting thrust and combustion stability.

This project gives a hybrid propulsion team a way to predict, before testing, whether their feed system is at risk — and to size their injector orifices so the system actually delivers the intended mass flow.

---

## What this project delivers

| Component | Description |
|---|---|
| **Theory docs** (`docs/`) | N₂O thermodynamics, SPI model derivation, HEM and Dyer two-phase models — written from first principles, no prior two-phase flow knowledge required |
| **Calculation model** (`src/model/`) | Python implementation of the full tank → feed line → injector path: Darcy-Weisbach friction losses, SPI/HEM/Dyer injector models, validated against published data |
| **Interactive tool** (`src/interface/`) | Streamlit web app — two modes (Sizing and Design), live diagrams, combustion stability check, sensitivity analysis, PDF report export |
| **Practical examples** (`examples/`) | Worked cases showing how to use the tool for real sizing scenarios |

---

## Interactive tool — quick start

### Requirements

```bash
pip install streamlit plotly matplotlib numpy reportlab
```

### Running

From the repository root:

```bash
streamlit run src/interface/app.py
```

Opens automatically at `http://localhost:8501`.

### Two modes

**Sizing mode** — you have a known orifice geometry and want to know the real mass flow the system will deliver, accounting for two-phase effects along the feed path.

**Design mode** — you have a target mass flow and want the orifice area (and diameter per hole) required to deliver it, with the Dyer correction applied so SPI under-sizing is avoided.

See [`docs/user_manual.md`](docs/user_manual.md) for a full walkthrough.

---

## Model overview

The tool chains three physics blocks:

```
Tank
  └─→ Feed line (Darcy-Weisbach friction + fitting losses)
        └─→ Injector sufficiency check (SPI criterion)
              ├─→ [SPI valid] Single-phase Bernoulli prediction
              └─→ [Two-phase] Dyer/NHNE model (blend of SPI and HEM limits)
```

The Dyer model is a weighted combination of the SPI limit ("no time to vaporise") and the HEM limit ("full thermodynamic equilibrium"), with the weight determined by how close the upstream pressure already is to saturation. It is the standard two-phase injector model in the university rocketry literature for this class of problem.

---

## Scope and known limitations

- Tank temperature is a **direct user input** — thermal balance with the environment (solar radiation, convection) is not modelled. See [`docs/future_work.md`](docs/future_work.md).
- Feed line assumed **adiabatic** and **steady-state** — no transient start-up effects.
- When flashing is detected in the feed line, the injector models are not evaluated (fluid arrives two-phase at the orifice inlet, outside their domain). A qualitative diagnostic and the SPI reference value are shown instead.
- Discharge coefficients and Dyer model parameters use **literature reference values**, not team-calibrated data.

---

## Repository structure

```
n2o-hybrid-injector-sizing/
├── README.md
├── docs/
│   ├── 01_n2o_thermodynamics.md     — N₂O saturation properties and flashing
│   ├── 02_spi_model.md              — SPI model derivation and assumptions
│   ├── 03_two_phase_flow.md         — HEM and Dyer models
│   ├── 04_implementation.md         — Module-by-module implementation notes
│   ├── future_work.md               — Identified extensions out of current scope
│   ├── user_manual.md               — Step-by-step guide to the interactive tool
│   └── references.md                — Sources and citations
├── src/
│   ├── model/
│   │   ├── n2o_properties.py        — Saturated N₂O thermophysical properties
│   │   ├── n2o_saturation_table.csv — Tabulated data (Table A.1, arXiv:2302.06725)
│   │   ├── feed_line.py             — Feed line pressure drop (Darcy-Weisbach)
│   │   ├── injector_spi.py          — SPI injector model
│   │   ├── injector_two_phase.py    — HEM and Dyer two-phase models
│   │   └── full_system.py           — Full path orchestration
│   └── interface/
│       ├── app.py                   — Streamlit interactive tool
│       ├── plotting.py              — Plotly interactive charts
│       └── export.py                — PDF report generation
└── examples/
    ├── example_01_sizing.md         — Sizing mode: predict mass flow for known geometry
    ├── example_02_design.md         — Design mode: find orifice area for target flow
    └── example_03_flashing.md       — Flashing case: diagnosis and design correction
```

---

## Author

Eduardo Costa — Aerospace Engineering, Instituto Superior Técnico (Técnico Lisboa).
Developed independently as a personal project in hybrid rocket propulsion.
