# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project aims to follow [Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-09-30

First stable release. The open items of the release candidate (bibliographic
checks, clean-clone check, fate of `apply_choking_limit()`) are all closed.

### Added
- **Coupled solver** (`full_system.py`): damped fixed-point iteration for the
  tank -> feed line -> injector operating point, with automatic SPI / Dyer
  (NHNE) / HEM two-phase-inlet regime selection; non-convergence raises a
  `RuntimeError` with the iteration history.
- **Design-mode sizing function** `design_injector_area()` (unit tested, no
  Streamlit dependency).
- **Feed line model** (`feed_line.py`): Darcy-Weisbach with Swamee-Jain,
  fitting losses, temperature-dependent liquid viscosity, and a HEM
  two-phase pressure-drop model after flashing onset.
- **Injector models** (`injector_spi.py`, `injector_two_phase.py`): SPI,
  HEM, Dyer with the corrected Solomon (2011) / Waxman (2013, Eq. 9)
  weights, HEM with two-phase inlet, isenthalpic and isentropic HEM critical
  flow.
- **Non-equilibrium choking diagnostic** (Henry-Fauske, 1971): returned as
  `m_dot_crit_HF` / `choked` next to the Dyer result and shown as a warning
  in the interface and the PDF; deliberately not applied as a cap.
- **Fuel grain sizing** (`grain_sizing.py`): Marxman regression rate, initial
  circular-port radius, multi-port, conservative burnback estimate;
  `a_from_reference_rate()` and a plausibility guard (5-300 mm).
- **Streamlit interface**: Sizing and Design modes, presets, Plotly diagrams,
  tornado sensitivity, combustion-stability check, low-supercharge notice,
  grain sizing panel, PDF export.
- **Validation** against Waxman et al. (2013/2014): original 4 points plus
  104 digitised points (Figs. 11-16, injector 3, dP up to 46 bar).
- **Exploratory**: `dyer_mass_flow_corrected()` (supercharge-gated kappa
  correction) and `validation/explore_supercharge_correction.py`; not called
  by `full_system.py` or the interface.
- Theory documents (`docs/01`-`04`), three worked examples, user manual,
  262 automated tests, GitHub Actions CI (Python 3.10 and 3.12).

### Changed
- N2O saturation properties now come from **CoolProp** (Lemmon & Span 2006)
  instead of the Perry / McGill correlations (Perry P_sat differed by up to
  4.81 % at 230 K). Validation report and all three examples regenerated.
  Part A MAPE: 3.51 % (old backend) -> **2.76 %**.
- Valid thermodynamic range is now [T_MIN, T_MAX] (~182.3-309.5 K); the old
  307.33 K limit on the entropy functions no longer applies (viscosity tables
  still stop at 307.33 K).

### Fixed
- Design mode crashed with a "math domain error" when the flow stays
  single-phase through the orifice; SPI is now selected in that regime.
- `dyer_non_equilibrium_parameter()` now raises a clear `ValueError` when
  P_downstream >= P_sat(T_upstream).
- Dyer weighting sign error of the original paper (Dyer et al. 2007)
  corrected.
- Viscosity tables (A.3/A.4) in `n2o_saturation_table.csv` no longer cite
  papers that do not support them (Millat et al. 1991 covers only the
  zero-density limit; Laesecke & Hafer 1998 is about fluorinated propanes);
  the values are documented as NIST WebBook data. No numerical change.
- Bibliography re-checked: Dyer et al. (2007) author list and Niño & Razavi
  (2019) citations verified; `docs/references.md` updated.

### Deprecated
- `apply_choking_limit()` (equilibrium-ceiling cap, retracted plan); emits a
  `DeprecationWarning` and is not used by the model. Kept on purpose, for
  reference: it documents the retracted approach and why the equilibrium
  ceiling is not a valid cap for Dyer (see `docs/future_work.md`, Priority 1).

### Known limitations
- Two-phase feed-line model and the HEM two-phase-inlet path are **not
  validated** against experimental data. The Dyer -> HEM switch at the
  flashing threshold is discontinuous, and the coupled solver can fail to
  converge in a narrow band of operating points close to it.
- Dyer accuracy depends on tank supercharge: MAPE ~2 % above ~14 bar
  (200 psi), up to ~12-18 % below it (always over-predicting). Validation
  covers one injector geometry, plus 4 points of a second.
- Tank temperature is a user input (no thermal model); steady-state only.
- Grain sizing is initial (t = 0), circular ports only.

[1.0.0]: https://github.com/eduardobcosta712/n2o-hybrid-injector-sizing/releases/tag/v1.0.0
