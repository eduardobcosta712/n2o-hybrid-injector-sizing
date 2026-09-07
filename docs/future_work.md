# Future Work

Extensions identified during this project's development, ordered by
technical importance per the project roadmap document.

---

## Implemented (removed from future work)

- Two-phase inlet via isenthalpic flash (x_inlet from feed line)  ✅
- HEM with two-phase inlet enthalpy  ✅
- Coupled feed-line / injector solver (iterative, damped fixed-point)  ✅
- Two-phase HEM pressure-drop model in feed line (rho_mix, mu_mix, x(s) per segment)  ✅
- Saturated vapour viscosity mu_v(T) from NIST WebBook (Table A.3, Millat et al. 1991)  ✅
- Saturated liquid Cp, mu_l, entropy (liquid+vapour) from NIST WebBook (Table A.4, Lemmon & Span 2006)  ✅
- HEM isenthalpic critical flow (hem_critical_flow, Waxman Eq.5) -- standalone  ✅
- Definitive validation against Waxman (2013/2014): MAPE = 3.51%, all 4 cases within +/-5%  ✅
- Sensitivity tornado plot  ✅
- Dyer formula weight correction (Solomon 2011)  ✅

---

## Priority 1 — Isentropic choking limit (automatic cap in Dyer)

**Current status (updated September 2026).** `hem_critical_flow()` in `injector_two_phase.py` implements the Waxman (2013) Eq.(5) isenthalpic maximum scan, correctly identifying the physical choking ceiling (41.1 g/s for Waxman conditions). At design-regime conditions, the Dyer predictions are validated within MAPE = 3.51% without any cap.

The isentropic path itself is now **implemented and validated, standalone**:
- `s_liquid_sat(T)`, `s_vapor_sat(T)`, `s_fg(T)` added to `n2o_properties.py`, interpolating Table A.4 (same pattern as `cp_liquid_sat`/`mu_liquid_sat`). Table A.4 is narrower than the correlation's own range (stops at 307.33 K, not 309.52 K) — this is now an explicit, separate range check (`T_MIN_A4`/`T_MAX_A4`, `_check_range_a4`), instead of silently falling through to an unhelpful generic error, which is also what `cp_liquid_sat`/`mu_liquid_sat` now use (previously they incorrectly used the wider, correlation-level range check).
- `vapor_quality_isentropic(s_upstream, T_downstream)` and `hem_critical_flow_isentropic(...)` added to `injector_two_phase.py`, mirroring `vapor_quality_isenthalpic`/`hem_critical_flow` exactly but along $s=$const instead of $h=$const. Kept **side by side** with the isenthalpic versions (not a replacement) — the isenthalpic scan is already validated and cited in `validation/waxman_2013_results.md`, and remains useful as a direct comparison.
- At Waxman conditions (T₁=280 K, P₁=4.36 MPa): isentropic $\dot m_{crit}$ = **41.64 g/s** vs. isenthalpic **41.05 g/s** (+1.43%) — both still below the experimental Dyer-regime range (44.0–48.0 g/s), consistent with the existing interpretation that Dyer's non-equilibrium correction legitimately predicts above either HEM-only ceiling.
- 12 new tests in `test_n2o_properties.py` (`TestEntropyFunctions`) and 14 new tests in `test_injector_two_phase.py` (`TestVaporQualityIsentropic`, `TestHemCriticalFlow` — the isenthalpic scan had no dedicated tests before this — and `TestHemCriticalFlowIsentropic`). All existing tests continue to pass.

**What remains.** Decide and implement **how** `hem_critical_flow_isentropic()` is applied automatically as a cap in `dyer_mass_flow()`: always `min(m_dot_Dyer, m_dot_crit_isentropic)`, or only outside the validated design-ΔP range (to avoid touching the already-validated Waxman result)? Re-run the full Waxman validation afterwards to confirm the MAPE does not change in the design regime (expected: it shouldn't, since the cap should not bind there).

**Domain note.** `hem_critical_flow_isentropic()` raises `ValueError` for $T_{upstream}$ above 307.33 K (Table A.4's limit) — `hem_critical_flow()` (isenthalpic) remains available as a fallback in that narrow band near the critical point.

**Dependency note.** Priority 4 below (CoolProp/REFPROP integration) would give direct access to N₂O entropy from a validated equation of state, instead of linear interpolation of the 27-point NIST Table A.4 used here, and would also remove the 307.33–309.52 K gap above. Since `s_liquid_sat`/`s_vapor_sat` are now small, isolated functions, re-pointing them at CoolProp later (if Priority 4 goes ahead) is a contained change — everything built on top of them (`vapor_quality_isentropic`, `hem_critical_flow_isentropic`, tests) is unaffected by the swap.

## Priority 2 — Full-system experimental validation

**Current state.** The Dyer injector model is validated against Waxman
(2013/2014) with mean error −1.9% at moderate pressure drops. The coupled
solver is new and not yet validated against full-system data.

**What is needed.** Experimental data covering the full path
(tank → line → injector) with known geometry, discharge coefficient,
and measured mass flow. The validation should be re-run after the coupled
solver since the one-pass results used previously may differ slightly.

---

## Priority 3 — OF ratio and fuel grain sizing

**Motivation.** The model currently outputs $\dot{m}_{oxidizer}$ at the converged
operating point. Given a target OF ratio (specified by the user from thermochemical
sizing), the fuel mass flow rate follows directly:

$$\dot{m}_{fuel} = \dot{m}_{oxidizer} / OF$$

With the regression rate correlation for paraffin or HTPB (Marxman):

$$\dot{r} = a\, G_o^n, \qquad G_o = \dot{m}_{oxidizer} / A_{port}$$

the initial grain geometry (port radius $r_0$, length $L$) can be estimated to
deliver the required $\dot{m}_{fuel}$ at the design burn duration.

**Why deferred.** Requires semi-empirical coefficients $a$, $n$ specific to
the fuel (paraffin, HTPB) and oxidiser (N₂O) combination, which are sourced
from static fire data or the literature. The model would expose these as
user inputs — it does not model the combustion or thermal processes that
govern regression rate. This is a post-processing step on the injector sizing
output, not a change to the two-phase flow model.

**Proposed inputs.** OF ratio (from CEA or equivalent), fuel type (dropdown
with literature $a$, $n$ values for common fuels), burn duration, number of
ports. Output: initial port radius, grain length, grain mass, estimated
$I_{sp}$ at design OF.

*Items implemented in this version are listed at the top. New items are added
here as they are identified, in the priority order used throughout this file.*

---

## Reordering note (September 2026)

Priorities 4 and 5 were swapped relative to the original roadmap draft:
**improved N₂O thermophysical properties (CoolProp/REFPROP integration)**
is now Priority 4, ahead of the **Monte Carlo uncertainty analysis**, now
Priority 5.

**Rationale.** CoolProp directly improves the accuracy of every downstream
model in this project — including the entropy data ($s_l$, $s_v$) needed
for the isentropic cap in Priority 1 — whereas running a Monte Carlo
analysis on top of correlations with a known, unquantified ~2–3%
systematic error near the critical point (Perry/McGill vs. REFPROP, see
Priority 4 below) would understate the real output uncertainty: the
quoted confidence interval would reflect input-parameter noise only, not
the model-form error already present in the property correlations.
Fixing the property backbone first is the more useful order of operations.

## Priority 4 — Improved N₂O thermophysical properties

**Motivation.** The Perry/McGill correlations have ~2–3% error near the
critical point vs. REFPROP. CoolProp (open source, no licence) provides
REFPROP-quality properties and would improve accuracy in the critical region.

**Why deferred.** The current property implementation is sufficient for
design-regime accuracy. CoolProp integration would also enable the
isentropic maximum search needed for choking (Priority 1).

## Priority 5 — Quantitative uncertainty analysis (Monte Carlo)

**Motivation.** The current sensitivity tornado shows which input matters
most. A Monte Carlo analysis would quantify the output uncertainty given
realistic input distributions, producing a confidence interval on mass flow,
and later in O/F ratio in the chamber.

**Proposed approach.** Sample input distributions (normal for temperature
and pressure, uniform for Cd and geometry), run the coupled solver for
each sample, and report the 5th–95th percentile of the mass-flow
distribution. Once Priority 4 is done, sample CoolProp-backed properties
directly rather than the Perry/McGill correlations, so the reported
confidence interval reflects only genuine input uncertainty, not
correlation error.

---

## Priority 6 — Improved exporting of the results

**Motivation.** The export.py files only informs the user about results already
presented in the interface. In consequence, the exporting should also be expanded 
to .JSON or .CSV files in order to export the test results to a global motor testing
or CAD modeling.

**Proposed approach.** Using a decoupled, schema-driven export architecture utilizing 
a unified Data Transfer Object (DTO) to seamlessly stream simulation parameters,
transient blowdown datasets, and hole-pattern geometric vectors into standardized 
JSON state files, CSV numerical tables, CAD-compatible coordinate scripts, and automated
PDF engineering reports.

---

## Priority 7 — Tank thermal model

**Motivation.** Tank temperature is currently a direct user input. A
thermal model would predict T_tank from ambient conditions, solar irradiance,
and tank geometry.

**Proposed approach.** Retrieve weather data via Open-Meteo API (free, no
key required) and implement a simplified steady-state tank-ambient heat
balance as a first approximation before a full transient model.

---

## Priority 8 — Transient / blowdown model

**Motivation.** The current model is steady-state. A transient model would
predict how mass flow, tank pressure, and temperature evolve over the burn.

**Why deferred.** Requires a reliable steady-state coupled solver (now
implemented) as the inner loop, plus a tank thermodynamic model (Priority 7).

---
