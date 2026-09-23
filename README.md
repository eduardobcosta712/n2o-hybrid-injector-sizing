# n2o-hybrid-injector-sizing

![Tests](https://github.com/eduardobcosta712/n2o-hybrid-injector-sizing/actions/workflows/tests.yml/badge.svg)

A physics-based tool and technical reference for predicting premature vaporisation (*flashing*) of N₂O in the feed systems of paraffin/N₂O hybrid rocket motors, and for sizing injector orifices that account for two-phase flow effects.

> **Disclaimer:** This tool is provided for predictive and educational purposes only. It is an academic project under active development. Its injector model has been validated against experiment only for injector pressure drops of 8–14 bar (see "Validation" below); results at other conditions are unvalidated model estimates. Results must not be used as the sole basis for engineering decisions or hardware fabrication. The author accepts no responsibility for any damages arising from the use of this tool.

---

## The problem

N₂O is the most common oxidiser in university-built hybrid motors — self-pressurising, relatively safe to handle, and commercially available. It has, however, a critical point at only ≈36.4 °C, dangerously close to typical ambient launch temperatures. This makes it exceptionally prone to partial vaporisation before it reaches the combustion chamber, whenever pressure in the feed system drops due to line friction, fittings, or acceleration through the injector orifice itself.

When this happens, the model most commonly used by university teams — SPI (*Single Phase Incompressible*, Bernoulli's equation with an empirical discharge coefficient) — is no longer valid: it implicitly assumes the fluid stays liquid from tank to chamber. Documented results in the literature show real mass flow far below the SPI prediction, by a factor of several times, shifting the motor's real O/F far from the design point and directly affecting thrust and combustion stability.

This project gives a hybrid propulsion team a way to predict, before testing, whether their feed system is at risk — and to size their injector orifices so the system actually delivers the intended mass flow.

---

## What this project delivers

| Component | Description |
|---|---|
| **Theory docs** (`docs/`) | N₂O thermodynamics, SPI model derivation, HEM and Dyer two-phase models, fuel grain sizing — written from first principles, no prior two-phase flow knowledge required |
| **Calculation model** (`src/model/`) | Coupled iterative solver for the full tank → feed line → injector path: Darcy-Weisbach friction losses, SPI/HEM/Dyer injector models, self-consistent operating point, non-equilibrium choking diagnostic, fuel grain sizing (Marxman), validated against published data in the band described below |
| **Interactive tool** (`src/interface/`) | Streamlit web app — two modes (Sizing and Design), live diagrams, combustion stability check, flashing diagnostics, non-equilibrium choking warning, fuel grain sizing panel, sensitivity analysis, PDF report export |
| **Practical examples** (`examples/`) | Worked cases showing how to use the tool for real sizing scenarios |

---

## Interactive tool — quick start

### Requirements

```bash
pip install streamlit plotly matplotlib numpy reportlab CoolProp
```

### Running

From the repository root:

```bash
streamlit run src/interface/app.py
```

Opens automatically at `http://localhost:8501`.

### Two modes

**Sizing mode** — you have a known orifice geometry and want to know the real mass flow the system will deliver, accounting for two-phase effects along the feed path.

**Design mode** — you have a target mass flow and want the orifice area (and diameter per hole) required to deliver it, with the Dyer correction applied so SPI under-sizing is avoided (or the SPI area itself when the flow stays single-phase through the orifice).

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

The Dyer model is a weighted combination of the SPI limit ("no time to vaporise") and the HEM limit ("full thermodynamic equilibrium"), with the weight determined by the non-equilibrium parameter κ = √[(P_up − P_down)/(P_sat − P_down)]. The formula uses the corrected weight convention of Solomon (2011) and Waxman (2013, Eq. 9) — large κ weights toward SPI (less equilibrium), small κ toward HEM.

**Non-equilibrium choking diagnostic.** Every Dyer evaluation is also checked against a Henry-Fauske (1971) non-equilibrium critical-flow ceiling — the physically correct bound for a non-equilibrium prediction, as opposed to the *equilibrium* HEM ceiling (which Dyer legitimately and correctly exceeds by design). This ceiling is surfaced as a warning (`choked` flag, both in the interface and the PDF export) when exceeded, rather than silently applied to the reported mass flow — see "Validation" and "Scope and known limitations" below for why.

---

## N₂O properties (CoolProp)

Saturation thermodynamics (P_sat, T_sat, densities, enthalpy, entropy) are computed by **CoolProp** (Lemmon & Span 2006 equation of state) instead of the earlier Perry/McGill correlations, which is more accurate (the earlier correlations carried a known ~2–5 % error in P_sat, largest away from the critical point — see `docs/references.md`). This has now been confirmed with the real CoolProp package installed (September 2026): `pip install CoolProp` succeeded, the full 262-test suite passes, and `validation/waxman_2013_results.md` was regenerated end to end. Viscosity (used in the feed line) still comes from NIST WebBook tables, unaffected by this change — CoolProp has no N₂O viscosity model.

## Validation

The full coupled model has been validated against experimental data from Waxman (2013/2014) for supercharged N₂O injectors — the domain for which the Dyer model is intended (subcooled liquid at the injector inlet, QF_upstream = 0). Four operating points at pressure drops of **8–14 bar** give a **MAPE of 2.76%** (mean error −0.3%), with all points within ±5%. This compares favourably with the Niño & Razavi (2019) reference result of 3.91% for the same dataset. The experimental values are read off a graph (±3% read-off uncertainty). Full report in [`validation/waxman_2013_results.md`](validation/waxman_2013_results.md); the script `validation/waxman_2013_validation.py` reproduces it.

The `hem_critical_flow()` and `hem_critical_flow_isentropic()` functions provide the *equilibrium* two-phase choking ceiling as standalone diagnostics (Waxman 2013 Eq. 5, isenthalpic and isentropic paths respectively — they agree to within 1.3%). At Waxman conditions they give ≈42–43 g/s; the Dyer predictions (43–50 g/s) sit above this, consistent with the non-equilibrium correction accounting for partial vaporisation inside the orifice.

The `henry_fauske_critical_flow()` function provides the physically appropriate **non-equilibrium** choking ceiling (Henry & Fauske, 1971). At Waxman conditions it gives 52.15 g/s — above all 4 validated Dyer predictions, so applying it as a diagnostic does not change the validated MAPE = 2.76% result. At other operating points further from the Waxman geometry this ceiling *does* bind (17% below the Dyer prediction in Example 1, see `examples/`); there is currently no experimental confirmation in this project's validation set at conditions where it actually changes the answer, which is why it is surfaced as a warning rather than applied automatically. See `docs/future_work.md`, Priority 1, for the full reasoning.

**What is *not* validated.** The Waxman dataset is the only open-access experimental dataset found in the correct domain. Consequently:
- **Pressure drops of 20–50 bar**, typical of motor designs and used in every worked example, are **outside the validated band**; there the model gives estimates with unquantified error (and the Henry-Fauske diagnostic flags them).
- The **two-phase feed-line model** and the **HEM two-phase-inlet injector path** (flashing in the line) are implemented and unit tested but have not been validated against any published data. Moreover, the switch from Dyer to HEM at the flashing threshold is discontinuous (HEM at vanishing vapour quality predicts roughly half the Dyer flow) — and, confirmed during the September 2026 CoolProp validation run, this discontinuity can also make the coupled fixed-point solver fail to converge for operating points that sit close enough to the flashing threshold, oscillating between the two branches instead of settling (see `docs/future_work.md`, Priority 2 and the "Audit follow-ups" section).
- The **coupled solver** has not been checked against a measured tank-to-chamber flow with a significant line.

---

> **262 automated tests** on Python 3.14 via `pytest tests/ -v` (6 test modules): `test_n2o_properties` 70, `test_feed_line` 32, `test_injector_spi` 13, `test_injector_two_phase` 55, `test_full_system` 45, `test_grain_sizing` 47. Confirmed passing with the real CoolProp package installed (September 2026).

## Scope and known limitations

- Tank temperature is a **direct user input** — thermal balance with the environment (solar radiation, convection) is not modelled. See [`docs/future_work.md`](docs/future_work.md).
- Feed line assumed **adiabatic** and **steady-state** — no transient start-up effects.
- When flashing is detected in the feed line, the model estimates the vapour quality at the injector inlet via isenthalpic flash and applies HEM with a two-phase inlet enthalpy. The Dyer blend is not used in this regime: its SPI branch represents delayed nucleation in a *liquid*, and its inlet state is undefined once vapour is present (κ itself is not the problem — it equals 1 at a saturated inlet). The switch is discontinuous and the path is unvalidated (see above) — and, as of the September 2026 CoolProp confirmation, can cause the coupled solver to fail to converge right at the threshold.
- Discharge coefficients use **literature reference values**, not team-calibrated data.
- The Dyer model is validated for pressure drops of 8–14 bar only (Waxman). A **non-equilibrium choking ceiling** (Henry-Fauske, 1971) is checked automatically and surfaced as a warning (not an automatic cap) when the Dyer prediction exceeds it — see "Validation" above.
- N₂O thermophysical properties: saturation thermodynamics from CoolProp (Lemmon & Span 2006 equation of state, confirmed installed and passing); viscosity (μ_v, μ_l) and the legacy Perry/McGill Table A.1 (kept only as test cross-checks) still come from `n2o_saturation_table.csv`. The literature attribution of the NIST *viscosity* correlations is unverified (see `docs/references.md`, "Open bibliographic points").
- Fuel grain sizing sizes the **initial** ($t=0$) circular-port geometry only; the regression-rate data point is a required user input (see `docs/03b_grain_sizing.md`).

---

## Repository structure

```
n2o-hybrid-injector-sizing/
├── README.md
├── docs/
│   ├── 01_n2o_thermodynamics.md
│   ├── 02_spi_model.md
│   ├── 03_two_phase_flow.md
│   ├── 03b_grain_sizing.md
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
│   │   ├── full_system.py
│   │   └── grain_sizing.py
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
│   ├── test_full_system.py
│   └── test_grain_sizing.py
└── examples/
    ├── example_01_sizing.md
    ├── example_02_design.md
    └── example_03_flashing.md
```

---

## Author

Eduardo Costa — Aerospace Engineering, Instituto Superior Técnico (Técnico Lisboa).
Developed independently as a personal project in hybrid rocket propulsion.
