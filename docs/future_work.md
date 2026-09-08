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
- HEM **isentropic** critical flow (hem_critical_flow_isentropic) -- standalone, kept side-by-side
  with the isenthalpic version for comparison; entropy functions s_liquid_sat/s_vapor_sat/s_fg
  added to n2o_properties.py (September 2026)  ✅
- Henry-Fauske (1971) non-equilibrium critical flow ceiling (henry_fauske_critical_flow) --
  surfaced as a diagnostic (m_dot_crit_HF, choked flag) alongside dyer_mass_flow(), not applied
  as an automatic cap; validated against Waxman (does not perturb MAPE = 3.51%) (September 2026)  ✅
- Definitive validation against Waxman (2013/2014): MAPE = 3.51%, all 4 cases within +/-5%  ✅
- Sensitivity tornado plot  ✅
- Dyer formula weight correction (Solomon 2011)  ✅

---

## Priority 1 — Non-equilibrium critical flow ceiling for Dyer ✅ RESOLVED (September 2026)

**Final state.** `henry_fauske_critical_flow()` was added to
`injector_two_phase.py`, implementing the simplified Henry-Fauske
(1971) non-equilibrium critical mass flux model for saturated/subcooled
liquid discharging through a converging nozzle, transcribed from
Simoneau, Henry, Hendricks & Watterson (1971), NASA TM X-67863, Eqs.
(2)-(5) -- see `references.md`. It reuses the isentropic quality
machinery already built for the (now superseded, but retained)
`hem_critical_flow_isentropic()`.

**Validation.** At Waxman conditions: `m_dot_crit` = **50.67 g/s** --
above the equilibrium HEM ceiling (41.05/41.64 g/s isenthalpic/
isentropic, correct direction: non-equilibrium exceeds equilibrium),
and above all 4 already-validated Dyer predictions (42.25-49.55 g/s),
so it does not cut into any validated result. Re-running the full
Waxman validation through the coupled solver (`evaluate_full_system`)
confirms **MAPE = 3.51%, unchanged**, with `choked = False` at all 4
points -- see `validation/waxman_2013_results.md`.

**Design decision: diagnostic, not automatic override.** At operating
points further from the Waxman geometry -- e.g.
`examples/example_01_sizing.md`'s conditions (20°C, 58→22 bar, 6×1.5mm)
-- the Henry-Fauske ceiling *does* bind, roughly 13-17% below the
uncapped Dyer blend. There is no experimental data point in this
project's validation set where the ceiling actually changes the answer
(all 4 Waxman points sit below it), so applying it as a silent
`min(m_dot_Dyer, m_dot_crit_HF)` override would risk quietly changing
already-published results with no empirical confirmation in the regime
where it matters. Instead, `dyer_mass_flow()` now returns
**both** values side by side -- `m_dot_Dyer` (unchanged) and
`m_dot_crit_HF` plus a `choked` boolean flag -- so calling code can
surface a warning when `choked = True` without silently altering the
headline number. The Streamlit interface (`app.py`) should be updated
to display this flag; not yet done, see below.

**What is genuinely still open.**
1. Wire the `choked`/`m_dot_crit_HF` fields into `app.py`'s result
   cards (a visible warning badge, similar to the existing combustion
   stability indicator) -- currently only available programmatically.
2. If/when experimental data becomes available at conditions where the
   ceiling binds, validate it there and reconsider whether to promote
   it to an automatic cap.
3. `hem_critical_flow_isentropic()` and the entropy functions
   (`s_liquid_sat`, `s_vapor_sat`, `s_fg`) remain in the codebase as a
   correct, useful standalone diagnostic (the equilibrium ceiling), even
   though they turned out not to be the fix for Dyer's extreme-ΔP
   behaviour -- Henry-Fauske was.

**History (kept for context).** This priority went through three
framings before landing here: (1) originally scoped as "build the
isentropic HEM scan and cap Dyer with it" -- built, but (2) shown
numerically to be wrong (would break the validated MAPE, since Dyer
legitimately exceeds the *equilibrium* ceiling by design), which
correctly redirected the search toward a genuine *non-equilibrium*
critical-flow model -- Henry-Fauske -- which (3) was implemented from a
primary source once the equations were sourced (initial web search
attempts returned image-embedded equations, unusable without risking
fabricated coefficients; the user supplied NASA TM X-67863 directly),
validated against Waxman, and found to change already-published example
numbers meaningfully outside the validated regime -- leading to the
side-by-side diagnostic design above rather than an automatic cap.

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
