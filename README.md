# n2o-hybrid-injector-sizing

![Tests](https://github.com/eduardobcosta712/n2o-hybrid-injector-sizing/actions/workflows/tests.yml/badge.svg)

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
| **Calculation model** (`src/model/`) | Coupled iterative solver for the full tank → feed line → injector path: Darcy-Weisbach friction losses, SPI/HEM/Dyer injector models, self-consistent operating point, non-equilibrium choking diagnostic, validated against published data |
| **Interactive tool** (`src/interface/`) | Streamlit web app — two modes (Sizing and Design), live diagrams, combustion stability check, non-equilibrium choking warning, sensitivity analysis, PDF report export |
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

The tool uses a **coupled iterative solver** that finds the self-consistent operating point where feed-line losses and injector flow are mutually consistent:

```
Initial mass-flow guess
         ↓
   Feed-line model  ──→  Injector inlet pressure
         ↑                        ↓
         └──────  Injector model (SPI / Dyer / HEM)
                          ↓
                     Converged?  →  Result
```

Three injector regimes are handled automatically:
- **SPI** — single-phase throughout the orifice
- **Dyer/NHNE** — liquid at the inlet, two-phase inside the orifice
- **HEM two-phase inlet** — fluid arrives partially vaporised (flashing in the feed line)

The Dyer model is a weighted combination of the SPI limit ("no time to vaporise") and the HEM limit ("full thermodynamic equilibrium"), with the weight determined by how close the upstream pressure already is to saturation. The formula uses the corrected weight convention of Solomon (2011) and Waxman (2013, Eq. 9) — large κ weights toward SPI (less equilibrium), small κ toward HEM.

**Non-equilibrium choking diagnostic.** Every Dyer evaluation is also checked against a Henry-Fauske (1971) non-equilibrium critical-flow ceiling — the physically correct bound for a non-equilibrium prediction, as opposed to the *equilibrium* HEM ceiling (which Dyer legitimately and correctly exceeds by design). This ceiling is surfaced as a warning (`choked` flag, both in the interface and the PDF export) when exceeded, rather than silently applied to the reported mass flow — see "Validation" and "Scope and known limitations" below for why.

---

## Validation

The full coupled model has been validated against experimental data from Waxman (2013/2014) for supercharged N₂O injectors — the correct domain for the model (subcooled liquid at the injector inlet, QF_upstream = 0). Four operating points at moderate pressure drops (8–14 bar), representative of real motor design conditions, give a **MAPE of 3.51%** (mean error −1.9%), with all points within ±5%. This compares favourably with the Niño & Razavi (2019) reference result of 3.91% for the same dataset. Full report in [`validation/waxman_2013_results.md`](validation/waxman_2013_results.md).

The `hem_critical_flow()` and `hem_critical_flow_isentropic()` functions provide the *equilibrium* two-phase choking ceiling as standalone diagnostics (Waxman 2013 Eq. 5, isenthalpic and isentropic paths respectively — they agree to within 1.4%). At Waxman conditions they give ≈41 g/s; the Dyer predictions (42–50 g/s) sit above this, consistent with the non-equilibrium correction accounting for partial vaporisation inside the orifice.

The `henry_fauske_critical_flow()` function provides the physically appropriate **non-equilibrium** choking ceiling (Henry & Fauske, 1971). At Waxman conditions it gives 50.67 g/s — above all 4 validated Dyer predictions, so applying it as a diagnostic does not change the validated MAPE = 3.51% result. At other operating points further from the Waxman geometry, this ceiling *can* bind (roughly 13–17% below the uncapped Dyer prediction in one tested case) — there is currently no experimental confirmation in this project's validation set at conditions where it actually changes the answer, which is why it is surfaced as a warning rather than applied automatically. See `docs/future_work.md`, Priority 1, for the full reasoning.

**On validation coverage.** The Waxman dataset is the only open-access experimental dataset in the correct domain (supercharged N₂O, tabulated operating points at design-relevant pressure drops). The two-phase inlet path (feed-line flashing) is physically implemented and tested but has not been validated against published data — no suitable open-access dataset was found. Neither has the Henry-Fauske choking diagnostic, outside the regime where it happens not to bind. These are noted as known limitations.

---

> **182 automated tests** pass locally (137 + 12 new in `test_n2o_properties.py` + 33 new in `test_injector_two_phase.py`, added September 2026) on Python 3.10 and 3.12 via GitHub Actions CI (pytest, 5 test modules covering all model components).

## Scope and known limitations

- Tank temperature is a **direct user input** — thermal balance with the environment (solar radiation, convection) is not modelled. See [`docs/future_work.md`](docs/future_work.md).
- Feed line assumed **adiabatic** and **steady-state** — no transient start-up effects.
- When flashing is detected in the feed line, the model estimates the vapour quality at the injector inlet via isenthalpic flash and applies HEM with a two-phase inlet enthalpy. The Dyer blend is not used in this regime (it collapses to SPI when P_upstream ≈ P_sat, which is physically incorrect — HEM is the appropriate limit).
- Discharge coefficients use **literature reference values**, not team-calibrated data.
- The Dyer model is validated for moderate pressure drops (design regime, ΔP = 20–50 bar in general use, 8–14 bar at the specific Waxman validation points). A **non-equilibrium choking ceiling** (Henry-Fauske, 1971) is checked automatically and surfaced as a warning (not an automatic cap) when the Dyer prediction exceeds it — see "Validation" above. The equilibrium HEM-based ceilings (`hem_critical_flow`, `hem_critical_flow_isentropic`) remain available as standalone diagnostics but are *not* the correct bound for a non-equilibrium prediction like Dyer, since real non-equilibrium two-phase flow chokes at a higher mass flux than full equilibrium allows — using them as a cap was tried, found to contradict the validated Waxman results, and reverted (see `docs/future_work.md`, Priority 1, for the full history).
- N₂O thermophysical properties sourced from McGill/Perry (Tables A.1–A.2) and NIST WebBook/Lemmon & Span 2006 (Tables A.3–A.4: μ_v, C_pl, μ_l, entropy). Functions include `mu_liquid_sat(T)`, `cp_liquid_sat(T)`, `mu_vapor_sat(T)`, `s_liquid_sat(T)`, `s_vapor_sat(T)`. Accuracy ~1–2% in design range; higher uncertainty near the critical point. Table A.4 (and therefore entropy-dependent functions, including the choking diagnostics) covers a narrower temperature range (182.33–307.33 K) than the main correlations (up to 309.52 K) — calls above 307.33 K raise a clear `ValueError` rather than silently failing.

---

## Repository structure

```
n2o-hybrid-injector-sizing/
├── README.md
├── docs/
│   ├── 01_n2o_thermodynamics.md
│   ├── 02_spi_model.md
│   ├── 03_two_phase_flow.md
│   ├── 04_implementation.md
│   ├── future_work.md
│   ├── user_manual.md
│   └── references.md
├── src/
│   ├── model/
│   │   ├── n2o_properties.py
│   │   ├── n2o_saturation_table.csv
│   │   ├── feed_line.py
│   │   ├── injector_spi.py
│   │   ├── injector_two_phase.py
│   │   └── full_system.py
│   └── interface/
│       ├── app.py
│       ├── plotting.py
│       └── export.py
├── validation/
│   ├── waxman_2013_results.md
│   ├── waxman_2013_validation.py
│   └── waxman_2013_experimental_data.csv
├── tests/
│   ├── conftest.py
│   ├── test_n2o_properties.py
│   ├── test_feed_line.py
│   ├── test_injector_spi.py
│   ├── test_injector_two_phase.py
│   └── test_full_system.py
└── examples/
    ├── example_01_sizing.md
    ├── example_02_design.md
    └── example_03_flashing.md
```

---

## Author

Eduardo Costa — Aerospace Engineering, Instituto Superior Técnico (Técnico Lisboa).
Developed independently as a personal project in hybrid rocket propulsion.
