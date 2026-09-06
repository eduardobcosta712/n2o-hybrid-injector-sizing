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

**Current status.** `hem_critical_flow()` in `injector_two_phase.py` implements the Waxman (2013) Eq.(5) isenthalpic maximum scan, correctly identifying the physical choking ceiling (41.1 g/s for Waxman conditions). It is available as a standalone function. At design-regime conditions, the Dyer predictions are validated within MAPE = 3.51% without any cap.

**What remains.** The isenthalpic scan uses $h = $ const, but the correct choking condition is along $s = $ const (isentropic). Entropy data ($s_l$, $s_v$) is now available in Table A.4. The remaining step is to implement the isentropic quality:

$$x_{is}(P_2) = \frac{s_{up} - s_l(T_{sat}(P_2))}{s_v(T_{sat}(P_2)) - s_l(T_{sat}(P_2))}$$

and find $\max_{P_2} \dot{m}_{HEM}(x_{is})$. This can then be applied as an automatic cap in `dyer_mass_flow()` for the regime where $\Delta P$ exceeds the design range.

**Why not yet.** The isentropic inversion requires interpolating entropy at arbitrary $P_2$, which needs care near the critical point. All required data is in the model — this is the next concrete implementation step.


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
here as they are identified. See `roadmap/Future_Improvements_Roadmap.md`
for the full priority order and rationale.*

---

## Priority 4 — Quantitative uncertainty analysis (Monte Carlo)

**Motivation.** The current sensitivity tornado shows which input matters
most. A Monte Carlo analysis would quantify the output uncertainty given
realistic input distributions, producing a confidence interval on mass flow,
and later in O/F ratio in the chamber.

**Proposed approach.** Sample input distributions (normal for temperature
and pressure, uniform for Cd and geometry), run the coupled solver for
each sample, and report the 5th–95th percentile of the mass-flow distribution.

---

## Priority 5 — Improved N₂O thermophysical properties

**Motivation.** The Perry/McGill correlations have ~2–3% error near the
critical point vs. REFPROP. CoolProp (open source, no licence) provides
REFPROP-quality properties and would improve accuracy in the critical region.

**Why deferred.** The current property implementation is sufficient for
design-regime accuracy. CoolProp integration would also enable the
isentropic maximum search needed for choking (Priority 3).

---

## Priority 6 - Improved exporting of the results

**Motivation.** The export.py files only informs the user about results already
presented in the interface. In consequence, the exporting should also be expanded 
to .JSON or .CSV files in order to export the test results to a global motor testing
or CAD modeling.

**Proposed approach.** Using a decoupled, schema-driven export architecture utilizing 
a unified Data Transfer Object (DTO) to seamlessly stream simulation parameters,
transient blowdown datasets, and hole-pattern geometric vectors into standardized 
JSON state files, CSV numerical tables, CAD-compatible coordinate scripts, and automated
PDF engineering reports.



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


