# 4. Implementation

This document tracks the computational implementation of the model described in Sections 1-3, module by module, in the order they were built. Each section documents: the theory being translated into code, the implementation itself, its validation, and any issues encountered along the way -- including issues that turned out to be informative rather than mere bugs.

Numbers quoted in the "Validation" paragraphs were regenerated with the confirmed CoolProp backend (September 2026); see each module's section below for what changed and why.

## 4.1 `n2o_properties.py` -- Saturated N2O properties

### Purpose

Every other module in this project depends on this one: the subcooling margin (Section 1.4), the flashing criterion (Section 1.6), and the SPI/HEM/Dyer models (Sections 2-3) all require, at minimum, $`P_{sat}(T)`$ and its inverse $`T_{sat}(P)`$. This module implements those, plus saturated liquid/vapour density, enthalpy, entropy, and the molar mass constant `M_N2O` (defined once here and imported by every other module).

### CoolProp integration (September 2026, Priority 4) -- confirmed

**What changed.** The closed-form Perry/McGill correlations and the McGill Table A.1 (nu_v, h_l, h_v) that the thermodynamic functions used before September 2026 have been replaced by **CoolProp** (Bell et al., 2014), which implements the **Lemmon & Span (2006)** fundamental equation of state for nitrous oxide directly -- the same equation of state the NIST WebBook itself uses. `P_sat`, `T_sat`, `rho_liquid_sat`, `nu_vapor_sat`, `h_liquid_sat`/`h_vapor_sat`/`h_fg`, `cp_liquid_sat`, and `s_liquid_sat`/`s_vapor_sat`/`s_fg` all now call `CoolProp.CoolProp.PropsSI` for a saturated state (`Q=0` or `Q=1`) at a given temperature or pressure.

**Why.** The Perry/McGill correlations carried a documented P_sat error that, once checked against the real CoolProp output, turned out to run largest at low temperature (up to 4.81 % at T = 230 K) and smallest near the critical point (about 0.9 % near T_MAX) -- the opposite trend from what this section originally assumed. Both P_sat and h_fg feed directly into the Dyer $`\kappa`$ weighting and the vapour-quality calculations that are the whole point of this project. CoolProp removes that error source and, since entropy now comes from the same equation of state as everything else, also removes the old 307.33 K ceiling on the isentropic choking models (Section 4.4): the valid range is now `[T_MIN, T_MAX]`, about 0.01 K above the triple point to 0.02 K below the critical point (evaluating a saturated state exactly at either endpoint is ill-conditioned for most equation-of-state implementations, hence the small margins).

**What did NOT change.** Viscosity ($`\mu_l`$, $`\mu_v`$) still comes from the NIST WebBook tables (A.3, A.4), interpolated exactly as before: CoolProp's fluid page for nitrous oxide lists only the equation of state and a surface-tension correlation, not a viscosity model, so there is nothing to switch to there. These tables still stop at 307.33 K. The molar-quantity API (`nu_vapor_sat` in m3/kmol, `h_liquid_sat` in kJ/kmol, ...) is unchanged, so `feed_line.py`, `injector_spi.py`, `injector_two_phase.py` and `full_system.py` needed no change to their own formulas -- only `injector_two_phase.py`'s domain-check messages, which used to name "Table A.4" for the entropy range, were updated to name the new, wider range.

**Requirement.** `pip install CoolProp`. Importing `n2o_properties.py` without it raises a clear `ImportError` naming the install command, rather than failing deeper in the module with an opaque `NameError`.

**Reference-state note.** CoolProp's zero of enthalpy and entropy is not the one the NIST tables in the CSV use. Every physical quantity in this project is a *difference* of enthalpy or entropy taken from the *same* source (e.g. $`x = (h_{up} - h_l)/h_{fg}`$), so the arbitrary reference cancels and this has no effect on any result.

**Verification status.** `pip install CoolProp` (8.0.0) succeeds on a normal machine with internet access. With the real package installed:

```bash
python src/model/n2o_properties.py        # self-check against literature references
pytest tests/ -v                          # full 262-test suite, including the legacy cross-checks
python validation/waxman_2013_validation.py
```

all pass, and the module's own self-check reproduces the literature reference values it targets. One test needed a genuine fix rather than a re-run: `test_agrees_with_perry_correlation` failed at the original `rel_tol=0.04` because the real CoolProp P_sat differs from the Perry correlation by up to 4.81 % at T = 230 K, well outside that tolerance; the test's tolerance was widened to `rel_tol=0.06`, with a comment recording the finding, rather than silently loosened. Every figure downstream of P_sat and h_fg moved by roughly 1-5 % relative to the earlier Perry/McGill-backed numbers, in the direction that error was known to point. `validation/waxman_2013_results.md`, `examples/example_01_sizing.md`, and `examples/example_02_design.md` have all been regenerated against the confirmed backend. `examples/example_03_flashing.md` needed a different starting operating point before it could be regenerated -- see that file and Section 4.2's "Coupled-solver note" below for why.

The Perry/McGill correlations and the McGill Table A.1 remain in the codebase as private helper functions inside `tests/test_n2o_properties.py` and in `n2o_saturation_table.csv`, used only as independent cross-checks in the test suite (e.g. `TestPsat.test_agrees_with_perry_correlation`), never by the model itself.

**If CoolProp cannot be installed** (e.g. a fully offline machine): see `docs/references.md`, "CoolProp and alternatives", for a short comparison of the fallback options (REFPROP, reverting to the Perry/McGill correlations, or implementing the Lemmon & Span equation of state directly).

### Tabulated viscosity: $`\mu_v(T)`$, $`\mu_l(T)`$

**Table A.3** (NIST WebBook): saturated vapour dynamic viscosity $`\mu_v(T)`$ at 26 points (182-307 K), linearly interpolated in `mu_vapor_sat(T)`. **Table A.4** (NIST WebBook): saturated liquid dynamic viscosity $`\mu_l(T)`$, used the same way in `mu_liquid_sat(T)`. Both raise a named `ValueError` outside their range (182.33-307.33 K), narrower than the thermodynamic range now available from CoolProp; `feed_line.py` clamps the local saturation temperature to 307.33 K before calling `mu_vapor_sat` when it is exceeded (Section 4.2). The literature attribution of the underlying viscosity correlations is unverified -- see `docs/references.md`, "Open bibliographic points".

Table A.4 also still carries $`c_{pl}(T)`$ and $`s_l(T)`$/s_v(T) columns (NIST); these are not read by the model any more (CoolProp supplies $`c_{pl}`$ and entropy) and are used in the tests purely as an independent cross-check, since they come from the same equation of state as CoolProp and should therefore agree tightly.

### Validation

`python src/model/n2o_properties.py` reproduces the literature reference values it targets (0 degC, 20 degC, and near the critical point), and every reference-value test in `tests/test_n2o_properties.py` passes against the real CoolProp backend (see "Verification status" above for the one tolerance that needed adjusting and why).

### File location

`src/model/n2o_properties.py`, reading `src/model/n2o_saturation_table.csv` (viscosity and legacy cross-check data only)

---
*Next section: 4.2 `feed_line.py` -- pressure drop and subcooling margin along the feed line (Darcy-Weisbach friction losses, fitting losses).*

## 4.2 `feed_line.py` -- Feed line pressure drop and subcooling margin

### Purpose

Implements Module 1 from the implementation plan (Section 3.5 synthesis): tracks pressure and subcooling margin from the tank exit to the injector inlet, through an arbitrary sequence of straight pipe segments and fittings (valves, elbows), flagging the point (if any) where the fluid crosses the saturation curve.

### Theory implemented

- **Reynolds number**, $`Re = \rho v D / \mu`$, to classify the flow regime.
- **Darcy friction factor** $`f`$: exact laminar solution $`f = 64/Re`$ for $`Re < 2300`$; the explicit **Swamee-Jain approximation** to the (implicit) Colebrook equation for $`Re \geq 2300`$,

```math
f = \frac{0.25}{\left[\log_{10}\left(\dfrac{\varepsilon/D}{3.7} + \dfrac{5.74}{Re^{0.9}}\right)\right]^2}
```

chosen over solving Colebrook directly because it is explicit (no iteration needed) while remaining within about 1 percent of the implicit solution for the range of interest. The 2300-4000 transitional regime is conservatively treated with the turbulent formula, documented as a deliberate simplification in the function's docstring.
- **Darcy-Weisbach friction loss**, $`\Delta P = f (L/D)(\rho v^2/2)`$, and **fitting (minor) losses**, $`\Delta P = K(\rho v^2/2)`$.
- The line is assumed **adiabatic** (constant temperature, per the Section 1.6/2.1 justification): only pressure is tracked segment by segment; temperature stays fixed at the tank value throughout.

**Scaling with diameter.** At fixed mass flow the velocity scales as $`1/D^2`$, so a fitting loss scales as $`1/D^4`$ and a pipe friction loss as roughly $`1/D^5`$ (slightly less, since $`f`$ rises weakly as $`Re`$ falls). Widening a line is therefore an even more effective remedy than the $`1/D^4`$ rule of thumb suggests.

### Liquid viscosity

The liquid viscosity used for $`Re`$ and for the liquid part of $`\mu_{mix}`$ is, by default, the temperature-dependent saturated-liquid value $`\mu_l(T_{tank})`$ = `mu_liquid_sat(T_tank)` (Table A.4). A constant `MU_LIQUID_N2O = 1.5e-4` Pa s is kept only as a fallback for $`T_{tank} > 307.33`$ K (outside Table A.4) or when the caller passes `mu=` explicitly.

### Two-phase HEM model

The feed line model handles two distinct flow regimes:

**Single-phase region** ($`P > P_{sat}(T_{tank})`$): Darcy-Weisbach with pure liquid properties -- $`\rho_l(T_{tank})`$ and $`\mu_l(T_{tank})`$.

**Two-phase region** ($`P \leq P_{sat}(T_{tank})`$): once flashing is detected, all subsequent segments use HEM mixture properties updated at each segment's local pressure:

```math
x(s) = \frac{h_l(T_{tank}) - h_l(T_{sat}(P(s)))}{h_{fg}(T_{sat}(P(s)))}
```

```math
\rho_{mix} = \frac{1}{\dfrac{1-x}{\rho_l} + \dfrac{x}{\rho_v}}
```

```math
\mu_{mix} = (1-x)\,\mu_l + x\,\mu_v
```

where $`\mu_v`$ comes from `mu_vapor_sat(T)` (NIST WebBook, Table A.3). For a near-critical local saturation temperature above Table A.3's 307.33 K limit, $`\mu_v`$ is evaluated at the table's last point (a small, documented approximation). The per-segment trace includes `x_quality`, `rho_eff_kg_m3`, and `mu_eff_Pa_s` for inspection and plotting.

**Physical significance.** With typical conditions (N2O at 20 degC, $`x \approx 0.05`$), $`\rho_{mix}`$ falls to roughly 60-100 kg/m$`^3`$ against $`\rho_l \approx 785`$ kg/m$`^3`$ -- a factor of 8-13x reduction. Since $`\Delta P \propto \dot{m}^2 / (\rho_{mix} A^2)`$, the two-phase friction losses in this region are correspondingly 8-13x larger than the single-phase model would predict. This is the dominant effect; the viscosity correction ($`\mu_{mix}`$ vs $`\mu_l`$) is secondary.

**`x_inlet` output.** The final `x_inlet` is computed via isenthalpic flash at $`P_{final}`$, which feeds directly into `hem_mass_flow_two_phase_inlet` in `full_system.py` when `flashing_detected = True`.

### Validation

Two cases were run with an identical line geometry (2 m of 8 mm ID tubing, one ball valve, one 90 degree elbow, 0.5 kg/s), differing only in initial tank subcooling:

- **Case A -- tank exactly at saturation** ($`\Delta T_{sub} = 0`$ initially): the very first segment already pushes the margin negative, and the flag `flashing_detected` correctly triggers. This is expected, not a bug: with zero initial margin, any pressure drop, however small, crosses the saturation curve -- the code is correctly enforcing the definition from Section 1.4 in the least forgiving case.
- **Case B -- tank with 5 bar of initial subcooling**: the same line loses about 2.5 bar but the final margin stays positive, and no flashing is flagged. Same hardware, different outcome -- driven entirely by how much margin the tank started with.

**Known limitations.** (i) `rho_liquid_sat(T_tank)` depends only on temperature, not pressure -- consistent with the incompressibility assumption underlying the SPI model (Section 2.4), but for lines with much larger pressure excursions real liquid density would vary somewhat with pressure too. (ii) The two-phase line model has been unit tested but not validated against experimental data. (iii) If the pressure reaches zero or below inside the line (an impossibly restrictive line at the requested flow), the model does not raise an error at that point; the downstream injector evaluation then fails with a `ValueError`.

**Coupled-solver note (added with the confirmed CoolProp backend).** The Dyer-to-HEM switch at the flashing threshold (Section 4.4) is discontinuous, and with the real CoolProp saturation pressure this discontinuity was found to create a narrow tank-pressure band, for a specific line geometry, where the coupled fixed-point solver in `full_system.py` oscillates between the two branches and fails to converge rather than settling on either one -- see Section 4.5 and `docs/future_work.md`, Priority 2, for the full account and for how `examples/example_03_flashing.md` was adjusted to sit clearly outside that band.

### File location

`src/model/feed_line.py`

---
*Next section: 4.3 `injector_spi.py` -- the SPI injector model, and the criterion for when it is sufficient on its own.*

## 4.3 `injector_spi.py` -- SPI injector model and sufficiency criterion

### Purpose

Implements the SPI model derived in Section 2 -- both directions (mass flow from a known orifice area, and required orifice area from a target mass flow) -- plus the decision criterion (planned in the pre-implementation discussion, "Module 2" logic) that determines whether SPI alone is valid at a given operating point, or whether two-phase effects must be considered.

### Theory implemented

```math
\dot{m}_{SPI} = C_d \, A \sqrt{2 \rho \Delta P}
```

and its algebraic inverse, $`A = \dot m_{target} / (C_d \sqrt{2\rho\Delta P})`$ -- the direction more commonly needed in practice, since the target mass flow is usually fixed by the motor's thermochemical sizing (O/F ratio, chamber pressure) and the orifice area is the unknown being solved for.

**Sufficiency criterion.** Per Section 3.1 (pressure is minimum, and therefore flashing risk is greatest, at the orifice's narrowest cross-section), SPI alone is valid only if even the downstream pressure stays at or above the saturation pressure evaluated at the upstream temperature:

```math
\text{SPI sufficient} \iff P_{downstream} \geq P_{sat}(T_{upstream})
```

If this fails, the flow crosses the saturation curve somewhere inside the orifice, and the two-phase model (Section 4.4, Dyer) must be used instead.

### Validation

A deliberately demanding illustrative case: N2O at 20 degC, $`P_{upstream} = 50`$ bar (already slightly below $`P_{sat}(20\ \text{degC})`$, i.e. entering the orifice already at the edge of saturation), $`P_{downstream} = 20`$ bar, $`A = 3.79`$ mm$`^2`$. The `spi_sufficient` criterion correctly returns `False`, so a single-phase SPI evaluation at this point would be reporting a number that is not the physically correct flow rate for this operating point, but rather the value SPI would (incorrectly) predict by assuming single-phase liquid throughout -- this is exactly the comparison baseline the Dyer model in Section 4.4 is validated against.

### File location

`src/model/injector_spi.py`

---
*Next section: 4.4 `injector_two_phase.py` -- HEM and Dyer models.*

## 4.4 `injector_two_phase.py` -- HEM and Dyer two-phase models

### Purpose

Implements the two-phase injector models from Section 3.4: HEM (full thermodynamic equilibrium) and Dyer (weighted blend of SPI and HEM), used when `spi_sufficient` returns `False`, plus the choking diagnostics.

### Theory implemented

Vapor quality at the orifice exit is obtained assuming an isenthalpic process and full equilibrium at the exit (exit sits on the saturation curve at $`T_{down} = T_{sat}(P_{down})`$):

```math
x = \frac{h_{upstream} - h_l(T_{down})}{h_{fg}(T_{down})}
```

The HEM mixture density follows from the mass-weighted average of the two phases' specific volumes, and $`\dot m_{HEM}`$ from the same orifice equation used throughout the project, with $`\rho_{HEM}`$ in place of the pure-liquid density. Dyer blends $`\dot m_{SPI}`$ and $`\dot m_{HEM}`$ via the non-equilibrium parameter $`\kappa`$:

```math
\kappa = \sqrt{\frac{P_{upstream} - P_{downstream}}{P_{sat}(T_{upstream}) - P_{downstream}}}
```

```math
\dot m_{Dyer} = \frac{\kappa}{1+\kappa}\,\dot m_{SPI} + \frac{1}{1+\kappa}\,\dot m_{HEM}
```

This uses the corrected weight convention of Solomon (2011) and Waxman (2013, Eq. 9): large $`\kappa`$ weights toward SPI (less time to reach equilibrium), small $`\kappa`$ toward HEM.

**Domain restrictions**, both enforced with an explicit `ValueError`:

1. $`P_{upstream} > P_{sat}(T_{upstream})`$ -- the fluid must still be liquid at the orifice inlet, consistent with this model addressing vaporization *inside* the orifice (Section 3.1), not an already-two-phase feed line (that case is handled by `hem_mass_flow_two_phase_inlet`, below).
2. $`P_{downstream} < P_{sat}(T_{upstream})`$ -- the pressure must actually cross saturation inside the orifice; otherwise $`\kappa`$ is undefined and SPI applies. Design mode routes through `full_system.design_injector_area`, which selects SPI automatically in that regime.

### Validation

Against the Waxman (2013/2014) dataset, with the confirmed CoolProp backend (`validation/waxman_2013_results.md`): four operating points at injector pressure drops of 8-14 bar give Dyer predictions of 42.91-50.37 g/s against experimental values of 44.0-48.0 g/s, for a MAPE of 2.76 percent (mean error -0.3 percent). No experimental validation exists outside that pressure-drop band.

### HEM with a two-phase inlet -- `hem_mass_flow_two_phase_inlet()`

When the feed line flashes, the fluid reaches the orifice with quality $`x_{inlet} > 0`$; the upstream enthalpy becomes $`h_{up} = h_l(T_{tank}) + x_{inlet}\,h_{fg}(T_{tank})`$ and HEM is applied from there. The Dyer blend is deliberately not used: its SPI branch encodes delayed nucleation in a *liquid*, which no longer applies once vapour is present, and its inlet state is undefined ($`\kappa`$ equals 1 at a saturated inlet, not infinity -- it only diverges as $`P_{down} \to P_{sat}`$).

**Known limitation.** The switch from Dyer to this HEM model at the flashing threshold is discontinuous: HEM at $`x_{inlet} \to 0`$ gives roughly half the Dyer flow at the same conditions. The model is unit tested but has no experimental validation, and (Section 4.2, "Coupled-solver note") the discontinuity can also make the coupled solver fail to converge for an operating point sitting close enough to the flashing threshold.

### HEM critical flow -- `hem_critical_flow()`

Available as a standalone function. Implements Waxman (2013) Eq. (5): scans $`P_2`$ from $`P_{sat}(T_{upstream})`$ downward and finds the maximum of the isenthalpic HEM mass-flow curve:

```math
\dot{m}_{crit} = \max_{P_2 < P_{sat}} \left[ C_d A \sqrt{2\,\rho_{mix}(P_2)\,\Delta P} \right]
```

This maximum is the choking limit of the equilibrium model. At Waxman conditions ($`T_1 = 280`$ K, $`P_1 = 4.36`$ MPa): 42.48 g/s at $`P_{2,crit} = 30.5`$ bar, $`x_{crit} = 0.0744`$. The function is not applied automatically to Dyer, because the Dyer non-equilibrium correction legitimately predicts above the HEM-only ceiling, confirmed by the Waxman experimental data (44-48 g/s vs. the equilibrium ceiling of 42.48 g/s).

### Isentropic choking scan -- `hem_critical_flow_isentropic()`

`hem_critical_flow()` above uses the isenthalpic path ($`h =`$ const), which correctly describes the real thermodynamic state of the fluid at the orifice exit (an orifice is adiabatic, so the first law gives $`h_{up}=h_{down}`$ regardless of internal irreversibility), but is only an approximation to the true choking condition. Choking is set by the two-phase speed of sound, $`c^2 = (\partial P/\partial\rho)_s`$ -- a derivative taken at constant entropy, because an acoustic disturbance is a small, fast, essentially reversible perturbation on top of the (possibly irreversible) mean flow.

`hem_critical_flow_isentropic()` implements this more rigorous scan, mirroring `hem_critical_flow()` exactly except along $`s =`$ const:

```math
x_{is}(P_2) = \frac{s_{up} - s_l(T_{sat}(P_2))}{s_v(T_{sat}(P_2)) - s_l(T_{sat}(P_2))}
```

```math
\dot m_{crit,\,is} = \max_{P_2 < P_{sat}} \left[C_d A \sqrt{2\,\rho_{HEM}(x_{is})\,(P_{up}-P_2)}\right]
```

The entropy functions `s_liquid_sat(T)`, `s_vapor_sat(T)` and `s_fg(T) = s_v(T) - s_l(T)` now come from CoolProp across the full thermodynamic range, so the scan is no longer limited to 307.33 K as it was under the older NIST-table-only entropy source.

**Validation.** At Waxman conditions ($`T_1=280`$ K, $`P_1=4.36`$ MPa, $`D=1.5`$ mm, $`C_d=0.65`$), with the confirmed CoolProp backend:

| Path | $`\dot m_{crit}`$ | $`P_{2,crit}`$ | $`x_{crit}`$ |
|---|---|---|---|
| Isenthalpic (`hem_critical_flow`) | 42.48 g/s | 30.5 bar | 0.0744 |
| Isentropic (`hem_critical_flow_isentropic`) | 43.02 g/s | 29.8 bar | 0.0774 |

+1.3 percent difference -- both remain below the experimental range (44.0-48.0 g/s), consistent with the interpretation that Dyer's non-equilibrium correction legitimately predicts above either HEM-only ceiling.

### Why the equilibrium ceiling is not a cap for Dyer

The original plan -- applying an equilibrium HEM ceiling as an automatic minimum cap on `dyer_mass_flow()` -- was tested numerically and found to be wrong:

- All 4 validated Waxman operating points have `m_dot_Dyer` sitting above the equilibrium ceiling (42.48 g/s), matching experiment within MAPE = 2.76 percent. Capping there would break this validated behaviour.
- Capping only the HEM term before blending does not help either: `hem_mass_flow()` evaluated at the actual `P_downstream` already falls below the critical value once past the choke point (it follows the descending branch of the curve), so a minimum cap there is a no-op exactly where it would be needed.

The function `apply_choking_limit()` implementing that retracted plan is kept only for backward compatibility; it is marked deprecated, emits a `DeprecationWarning`, and nothing in the repository calls it.

### Henry-Fauske non-equilibrium critical flow -- `henry_fauske_critical_flow()`

**Source.** Henry, R.E. & Fauske, H.K. (1971), "The Two-Phase Critical Flow of One-Component Mixtures in Nozzles, Orifices, and Short Tubes," ASME J. Heat Transfer, 93(2), 179-187. Equations transcribed from the simplified form in Simoneau, R.J., Henry, R.E., Hendricks, R.C. & Watterson, R. (1971), "Two-Phase Critical Discharge of High Pressure Liquid Nitrogen," NASA TM X-67863, Eqs. (2)-(5).

**Theory.** For saturated/subcooled liquid at the nozzle inlet ($`P_0/P_c > 0.05`$, comfortably true here), five assumptions apply: negligible vapour before the throat (so the inlet-to-throat momentum balance is single-phase Bernoulli), incompressible liquid, equilibrium vapour formation at the throat, equal liquid/vapour velocity at the throat, and -- the key non-equilibrium closure -- a fractional mass-transfer rate:

```math
\eta = \frac{P_t}{P_0} = 1 - \frac{v_{l0}\,G_c^2}{2P_0} \qquad \text{(momentum, Eq. 2)}
```

```math
N = \min\left(1,\ \frac{x_E}{0.14}\right) \qquad \text{(Henry 1970 fit to Starkman et al. steam-water data)}
```

```math
G_c^2 = \left[\frac{N\,(v_{gE}-v_{l0})}{s_{gE}-s_{lE}}\,\frac{ds_{lE}}{dP}\right]^{-1} \qquad \text{(mass-transfer closure, Eq. 5)}
```

where $`x_E`$ is the equilibrium quality at the throat, computed by `vapor_quality_isentropic()`. $`ds_{lE}/dP`$ comes from the chain rule $`ds_{lE}/dP = (ds_l/dT) / (dP_{sat}/dT)`$, with $`ds_l/dT`$ from a small central finite difference on `s_liquid_sat(T)` and $`dP_{sat}/dT`$ from the analytical `dP_sat_dT`. All volumes and entropies are converted to specific (per unit mass) quantities for Eq. 5's units to work out to a mass flux.

**Solved by bisection, not fixed-point iteration.** Eqs. (2) and (5) are coupled ($`G_c`$ depends on properties at the unknown throat pressure $`P_t`$, which itself depends on $`G_c`$ via Eq. 2). A first, naive fixed-point implementation diverged: $`G_c`$ from Eq. 5 blows up as $`P_t \to P_{sat}(T_0)`$ from below (since $`N\to0`$ there). Framing it as a residual $`f(P_t) = P_{t,\,momentum}(G_c(P_t)) - P_t`$ and bisecting is robust.

**Validation.** At Waxman conditions, with the confirmed CoolProp backend: $`\dot m_{crit} = 52.15`$ g/s, at $`P_{2,crit} = 31.76`$ bar, $`x_E = 0.0576`$, $`N = 0.412`$ -- above the equilibrium ceilings (correct direction), and above all 4 validated Dyer predictions (42.91-50.37 g/s). Re-running the full Waxman validation through the coupled solver confirms MAPE = 2.76 percent, `choked = False` at all 4 points.

**Important caveat -- surfaced, not hidden.** At operating points further from the Waxman geometry the ceiling does bind (see the worked examples in `examples/`). There is currently no experimental data point where the ceiling actually changes the answer. Rather than silently override `m_dot_Dyer` with a value that is theoretically sound but empirically unconfirmed in the regime where it matters, `dyer_mass_flow()` returns both `m_dot_Dyer` (always unchanged) and `m_dot_crit_HF` plus a `choked` boolean, side by side (see `docs/future_work.md`, Priority 1). The interface surfaces this as a warning.

### File location

`src/model/injector_two_phase.py`

---
*Next section: 4.5 `full_system.py` -- orchestrating the full tank to feed line to injector path.*

## 4.5 `full_system.py` -- Full system orchestration and coupled solver

### Purpose

Chains `feed_line.py`, `injector_spi.py`, and `injector_two_phase.py` in the order the sizing problem requires (Section 3.5 synthesis): evaluate the feed line to get conditions at the injector inlet, then automatically decide whether SPI, Dyer, or HEM with a two-phase inlet applies. The self-consistent operating point satisfies both

```math
P_{inlet}^{*} = P_{tank} - \Delta P_{line}(\dot m^{*}), \qquad \dot m^{*} = \dot m_{injector}(P_{inlet}^{*})
```

and is found by damped fixed-point iteration ($`\alpha = 0.5`$, tolerance $`10^{-4}`$ on the relative change of $`\dot m`$). Non-convergence raises a `RuntimeError` with the iteration history rather than returning a wrong number. `m_dot_design` is only the initial guess.

### `design_injector_area()`

The Design-mode logic -- choose SPI or Dyer, and iterate the Dyer area to hit the target flow -- lives in a pure function so it can be unit tested without Streamlit. If the flow stays single-phase through the orifice ($`P_{chamber} \geq P_{sat}(T_{tank})`$) it returns the SPI area with `regime = "SPI"`; otherwise it iterates the Dyer area (all flow models scale linearly with area, so it converges in one or two steps).

### The Dyer-to-HEM discontinuity and solver convergence

Because `hem_mass_flow_two_phase_inlet` and `dyer_mass_flow` are genuinely different models that disagree at the flashing threshold (Section 4.4), the coupled solver can, for a tank pressure sitting close enough to the threshold for a given line geometry, alternate between predicting flashing (moving it toward the HEM branch) and predicting no flashing (moving it toward the Dyer branch) without settling. This was confirmed directly with the real CoolProp backend for the geometry originally used in `examples/example_03_flashing.md`: a narrow tank-pressure band exists there where the damped fixed-point iteration does not converge even at reduced damping, and the solver correctly raises `RuntimeError` rather than returning an arbitrary value from an unconverged state. Just outside that band, on either side, the solver converges cleanly. `examples/example_03_flashing.md` now uses tank pressures chosen to sit clearly outside this band -- see that file and `docs/future_work.md`, Priority 2.

### Validation

Reusing `feed_line.py`'s Case B geometry with a 55 bar tank and $`A = 3.79`$ mm$`^2`$: the feed line loses a small fraction of a bar at the converged flow, no flashing is detected, but `spi_sufficient` correctly returns `False` given the large pressure drop to a 20 bar chamber, and Dyer is selected automatically, predicting substantially less flow than SPI would have at the same geometry -- see `examples/example_01_sizing.md` for the fully worked, CoolProp-confirmed version of this same comparison (386.1 g/s Dyer vs. 515.0 g/s SPI).

### File location

`src/model/full_system.py`

---
*Next section: 4.6 `src/interface/` -- the interactive Streamlit tool.*

## 4.6 `src/interface/` -- Interactive tool

### Purpose

A Streamlit web app (`app.py`) wrapping `full_system.py`: editable tank, feed line (dynamic pipe/fitting segment list), and injector inputs, updating live as inputs change. Reports flashing/SPI-sufficiency status and the real mass flow, and renders interactive diagrams (`plotting.py`) and a one-page PDF (`export.py`).

### Implementation notes

Segment state is kept in `st.session_state`, since Streamlit re-runs the whole script on every interaction; without it, the segment list would reset on every slider move. All UI inputs are in display-friendly units (bar, degC, mm, g/s) and converted to SI at the UI boundary before calling into `full_system.py`, which continues to operate in SI throughout, per the project's units convention.

**Flashing in the feed line.** In Sizing mode the result cards show the HEM two-phase-inlet estimate (labelled as such), and the diagnostics panel -- suggestions to remove the flashing: raise the subcooling margin to at least 5 bar, shorten the line, enlarge the diameter, replace high-K fittings, pre-cool -- is shown alongside, together with a naive-SPI reference value. In Design mode area sizing is not offered while the line flashes; the same panel is shown.

**Henry-Fauske choking diagnostic.** Both modes surface a warning (amber badge and `st.warning`) whenever `dyer_mass_flow()`'s `choked` flag is `True`. Design mode's warning notes the condition is independent of orifice area (both scale linearly with area). `plot_model_comparison()` draws the ceiling as a reference line, red when exceeded, and the PDF export includes the ceiling value and a note when triggered. None of this changes any displayed mass-flow number.

**Design mode, SPI-valid regime.** When the chamber pressure is at or above $`P_{sat}(T_{tank})`$ the result cards show "Recommended (SPI)" and an explanatory note; the model-comparison chart is replaced by the pressure-along-line chart, and the PDF omits the Dyer rows.

**Grain sizing panel.** See Section 4.7.

### Running

From the repository root: `streamlit run src/interface/app.py`. Requires `pip install streamlit plotly matplotlib numpy reportlab CoolProp`.

### File location

`src/interface/app.py`, `src/interface/plotting.py`, `src/interface/export.py`

---
*Next section: 4.7 `grain_sizing.py` -- fuel grain sizing via the Marxman regression rate correlation.*

## 4.7 `grain_sizing.py` -- Fuel grain sizing

### Purpose

Implements Priority 3 of the roadmap: from the oxidiser mass flow already computed by the rest of this tool and a target O/F ratio, size the initial fuel grain geometry (port radius, given a required grain length and number of ports) via the Marxman regression rate correlation. See `docs/03b_grain_sizing.md` for the full theoretical derivation.

### Theory implemented

```math
\dot m_{fuel} = \frac{\dot m_{ox}}{OF}, \qquad \dot r = a\,G_o^n, \qquad G_o = \frac{\dot m_{ox}}{A_{port}}
```

Combining these for a single circular port of radius $`r`$ and length $`L`$ gives $`\dot m_{fuel}(r) = K r^{1-2n}`$ (see the theory doc for $`K`$), solved for $`r_0`$ by a log-spaced scan (1 micron to 1 km) followed by bisection -- the same convention already used elsewhere in this project for transcendental relationships.

**Scope decisions:**
- $`a`$ and $`n`$ are required inputs, not fixed per-fuel defaults (2-3x scatter between independent studies).
- Grain length $`L`$ is a required input, not derived from an unsourced L/D heuristic.
- No specific-impulse output -- requires a chemical equilibrium code (CEA/RPA).
- Only circular ports (single- or multi-port) -- non-circular shapes tracked as Priority 3b.

`FUEL_PROPERTIES` provides density defaults for paraffin wax, HTPB, ABS and PMMA (each cited in `references.md`), alongside a non-authoritative $`(a,n)`$ reference range per fuel.

### Protection against unit errors

Two real test failures (a wrong-unit $`a`$ solving "successfully" to a port radius of several metres) led to two safeguards:

- `a_from_reference_rate(r_dot_ref_mm_s, G_o_ref, n)` converts a regression-rate data point as read off a plot or table (mm/s at a stated $`G_o`$ in kg/(m$`^2`$ s)) into the SI coefficient $`a`$. The interface asks for exactly these three numbers and never for $`a`$ itself.
- `size_grain()` raises a `RuntimeError` when the solved initial radius lies outside `PLAUSIBLE_PORT_RADIUS_RANGE_M` = 5-300 mm, reporting the regression rate the given $`(a,n)`$ imply at the rejected radius. (`solve_initial_port_radius` keeps its very wide search range on purpose, so that it always finds a root when one exists.)

### A physical subtlety worth restating here

For $`n>0.5`$, the sign of the exponent $`1-2n`$ flips negative: a larger target fuel flow requires a smaller port radius, the reverse of naive intuition, because $`G_o \propto 1/r^2`$ falls faster than the burning perimeter $`\propto r`$ grows. At $`n=0.5`$ exactly, fuel flow is independent of port radius entirely. Both directions are explicitly tested in `test_grain_sizing.py`.

### Validation

47 tests in `test_grain_sizing.py`: hand-computed known values; closed-form-vs-bisection cross-check (agreement to within a few parts in $`10^7`$ relative); both directions of the $`n`$ vs. 0.5 radius-flow relationship, verified numerically; the $`n=0.5`$ degenerate case, including its correctly-unreachable-target failure mode; multi-port perimeter scaling (proportional to $`\sqrt N`$ at fixed total area, checked directly via pure geometry); the plausibility guarantee (both bounds, and the diagnostic content of the message); `a_from_reference_rate` (known value, round trip, scaling, and an end-to-end realistic case); and edge cases. There is no experimental dataset for the grain sizing itself.

### Interface integration

`app.py` has a "Grain sizing" expander in both modes (fed with the already-computed $`\dot m_{ox}`$): O/F, fuel dropdown (density auto-filled, editable), the three-number regression-rate data point (with the literature reference range shown as orientation only), grain length, number of ports, optional burn duration. It displays $`\dot m_{fuel}`$, initial port radius, initial $`G_o`$ and $`\dot r`$, and -- when a burn duration is given -- the conservative final-radius / fuel-consumed estimate, captioned as a first-order approximation rather than a transient simulation.

### File location

`src/model/grain_sizing.py`
