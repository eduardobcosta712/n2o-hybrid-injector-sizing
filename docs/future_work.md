# Future Work

Extensions identified during this project's development, ordered by
technical importance per the project roadmap document.

---

## Implemented (removed from future work)

- Two-phase inlet via isenthalpic flash (x_inlet from feed line)  
- HEM with two-phase inlet enthalpy  
- Coupled feed-line / injector solver (iterative, damped fixed-point)  
- Two-phase HEM pressure-drop model in feed line (rho_mix, mu_mix, x(s) per segment)  
- Saturated vapour viscosity mu_v(T) from NIST WebBook (Table A.3)  
- Saturated liquid Cp, mu_l, entropy (liquid+vapour) from NIST WebBook (Table A.4, Lemmon & Span 2006)  
- Temperature-dependent liquid viscosity mu_l(T) actually used in the feed line  
- HEM isenthalpic critical flow (hem_critical_flow, Waxman Eq.5) -- standalone  
- HEM **isentropic** critical flow (hem_critical_flow_isentropic) -- standalone, kept side-by-side
  with the isenthalpic version for comparison; entropy functions s_liquid_sat/s_vapor_sat/s_fg
  added to n2o_properties.py (September 2026)  
- Henry-Fauske (1971) non-equilibrium critical flow ceiling (henry_fauske_critical_flow) --
  surfaced as a diagnostic (m_dot_crit_HF, choked flag) alongside dyer_mass_flow(), not applied
  as an automatic cap; validated against Waxman (does not perturb MAPE = 3.51%) (September 2026);
  wired into the interface (warning badge in both modes, comparison-chart line, PDF note)  
- Fuel grain sizing via the Marxman regression rate correlation (grain_sizing.py) -- initial port
  radius, fuel mass flow, multi-port geometry, first-order conservative burnback estimate; a and n
  are required user inputs, entered as a regression-rate data point converted by
  a_from_reference_rate(); results outside a plausible port-radius range are refused; density
  defaults for paraffin/HTPB/ABS/PMMA (September 2026); integrated in the interface 
- Definitive validation against Waxman (2013/2014): MAPE = 3.51%, all 4 cases within +/-5%  
- Sensitivity tornado plot  
- Dyer formula weight correction (Solomon 2011)  
- Design-mode sizing logic extracted into a tested function (full_system.design_injector_area),
  fixing a crash when the flow is single-phase through the orifice   
- Explicit range errors for Table A.3 (mu_vapor_sat) and for T_sat(P) outside the correlation  
- CoolProp implementation and validation, confirming, at least for the waxman results, a more accurate result

---

## Priority 1 — Non-equilibrium critical flow ceiling for Dyer - RESOLVED (September 2026)

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
points further from the Waxman geometry the Henry-Fauske ceiling *does*
bind: 17% below the Dyer prediction in `examples/example_01_sizing.md`
(315 vs 382 g/s), 10% below the target in `example_02_design.md`
(452 vs 500 g/s), and it is also exceeded by the corrected design of
`example_03_flashing.md`. There is no experimental data point in this
project's validation set where the ceiling actually changes the answer,
so applying it as a silent `min(m_dot_Dyer, m_dot_crit_HF)` override
would risk quietly changing published results with no empirical
confirmation in the regime where it matters. Instead, `dyer_mass_flow()`
returns **both** values side by side -- `m_dot_Dyer` (unchanged) and
`m_dot_crit_HF` plus a `choked` boolean flag -- and the interface shows
a warning when `choked = True`.

**What is genuinely still open.**
1. If/when experimental data becomes available at conditions where the
   ceiling binds (roughly 20-50 bar pressure drop), validate it there
   and reconsider whether to promote it to an automatic cap.
2. `hem_critical_flow_isentropic()` and the entropy functions
   (`s_liquid_sat`, `s_vapor_sat`, `s_fg`) remain in the codebase as a
   correct, useful standalone diagnostic (the equilibrium ceiling), even
   though they turned out not to be the fix for Dyer's extreme-ΔP
   behaviour -- Henry-Fauske was.
3. `apply_choking_limit()` (the retracted equilibrium-cap plan) is
   deprecated and unused; delete it once no external script depends on it.

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
(2013/2014) with mean error −0.3% at pressure drops of 8-14 bar, using
a line with negligible losses. Three parts of the model have **no
experimental validation at all**:

1. **Pressure drops of 20-50 bar** (the range of typical motor designs).
   Every worked example lies here, and the Henry-Fauske diagnostic
   flags all of them.
2. **The coupled feed-line / injector solver** with a line that matters
   (the Waxman line is short and wide).
3. **The two-phase inlet path** (flashing in the line): the two-phase
   line model and the HEM two-phase-inlet injector model are implemented
   and unit tested only. In addition, the switch from Dyer to HEM at the
   flashing threshold is discontinuous (HEM at x_inlet -> 0 gives roughly
   half of the Dyer flow at the same conditions -- see
   `examples/example_03_flashing.md`); a smooth transition, or at least
   a quantified uncertainty band around the threshold, is worth
   designing once data exist.

**What is needed.** Experimental data covering the full path
(tank → line → injector) with known geometry, discharge coefficient,
and measured mass flow, ideally including runs at 20-50 bar drop and
runs with and without flashing in the line. A team's own cold-flow
(water or N₂O) data would already help calibrate Cd and check the
line-loss model.

---

## Priority 3 — OF ratio and fuel grain sizing - RESOLVED (September 2026)

**Final state.** `grain_sizing.py` implements the Marxman regression rate correlation ($\dot r = a G_o^n$) for sizing the initial fuel grain geometry from a target O/F ratio and the oxidiser mass flow already computed by the rest of this tool. See `docs/03b_grain_sizing.md` for the full theoretical derivation.

**Scope deliberately corrected relative to the original plan below** (kept for context; do not re-read it as the current spec):

1. **$a$ and $n$ are required user inputs, not fixed per-fuel defaults.** Researching citable coefficients surfaced genuine, large scatter between independent studies of nominally the same fuel/oxidiser pair (paraffin/N₂O regression rates of ≈2, ≈3.5, and 4–5 mm/s each separately reported at comparable oxidiser mass flux — see `references.md`). Unlike density, $a$ and $n$ are test-article-specific empirical fits (injector design, motor scale, chamber pressure range all matter), not universal material constants. Shipping a literature pair as a default would imply false precision. `FUEL_PROPERTIES` in `grain_sizing.py` provides **density only** as a safe default (paraffin, HTPB, ABS, PMMA — all properly sourced, see `references.md`), plus a labelled, non-authoritative reference range for $a, n$ per fuel.
2. **Grain length $L$ is a required input, not a derived output.** The original plan implied deriving $L$ from an L/D heuristic; rather than use an unsourced ratio, $L$ is left as a direct input — most teams already know their available case length as a hard constraint.
3. **No estimated $I_{sp}$ output.** That requires a chemical equilibrium code (CEA/RPA) this project does not implement or wrap. Get $I_{sp}$ at the design O/F from CEA/RPA directly.
4. **Only circular ports (single- or multi-port) are supported.** Non-circular shapes (star, wagon-wheel) are tracked separately below, not delivered here.

**What is delivered.** `size_grain()`: given $\dot m_{ox}$, target OF, $(a, n)$, fuel density, grain length, and number of ports, solves for the initial port radius, reports the initial oxidiser mass flux and regression rate, and — if a burn duration is supplied — a first-order, deliberately **conservative** (over-)estimate of final port radius and fuel mass consumed (using the *initial* regression rate held constant; real regression rate falls as the port opens up, so this over-estimates burnback, appropriate as a safety-margin check, not a transient simulation). Two protections against unit errors: the interface asks for a regression-rate data point (mm/s at a stated $G_o$) converted by `a_from_reference_rate()` instead of asking for $a$, and `size_grain()` refuses results outside a plausible port-radius range (5-300 mm).

**Physical subtlety surfaced during testing, worth flagging prominently:** for $n>0.5$ (e.g. many HTPB fits), a *larger* target fuel flow requires a *smaller* port radius — the reverse of naive intuition — because $G_o \propto 1/r^2$ falls faster than the burning perimeter $\propto r$ grows. At $n=0.5$ exactly, fuel flow is independent of port radius entirely (a well-known, practically valuable property of $n\approx0.5$ fuels like paraffin). Both directions are explicitly tested (not assumed) in `test_grain_sizing.py`.

**Validation.** 47 tests (`test_grain_sizing.py`): known-value hand cross-checks, closed-form-vs-bisection agreement, both directions of the $n$ vs. $0.5$ radius-flow relationship, the $n=0.5$ degenerate case, multi-port perimeter scaling, the plausibility guarantee, `a_from_reference_rate`, and edge cases. There is no experimental dataset for the grain sizing itself.

**What remains open.**
1. **Non-circular port shapes** (star, wagon-wheel) — Priority 3b below.
2. **Full transient burn simulation** (radius, $G_o$, $\dot r$, O/F all evolving over the burn) remains Priority 8, unaffected by this work beyond providing its $t=0$ starting point.

---

*Original plan, superseded by the corrections above — kept only for historical context:*

> **Motivation.** The model currently outputs $\dot{m}_{oxidizer}$ at the converged
> operating point. Given a target OF ratio (specified by the user from thermochemical
> sizing), the fuel mass flow rate follows directly:
>
> $$\dot{m}_{fuel} = \dot{m}_{oxidizer} / OF$$
>
> With the regression rate correlation for paraffin or HTPB (Marxman):
>
> $$\dot{r} = a\, G_o^n, \qquad G_o = \dot{m}_{oxidizer} / A_{port}$$
>
> the initial grain geometry (port radius $r_0$, length $L$) can be estimated to
> deliver the required $\dot{m}_{fuel}$ at the design burn duration.
>
> **Proposed inputs.** OF ratio (from CEA or equivalent), fuel type (dropdown
> with literature $a$, $n$ values for common fuels), burn duration, number of
> ports. Output: initial port radius, grain length, grain mass, estimated
> $I_{sp}$ at design OF.

---

## Reordering note (September 2026)

Priorities 4 and 5 were swapped relative to the original roadmap draft:
**improved N₂O thermophysical properties (CoolProp/REFPROP integration)**
is now Priority 4, ahead of the **Monte Carlo uncertainty analysis**, now
Priority 5. Priority 4 is now implemented (see above), which is what makes
this swap concrete rather than theoretical: Priority 5's sampled property
values are the CoolProp ones as of this section, not the superseded
Perry/McGill correlations.

**Rationale.** CoolProp directly improves the accuracy of every downstream
model in this project — including the entropy data ($s_l$, $s_v$) used by
the choking models — whereas running a Monte Carlo analysis on top of
correlations with a known, unquantified ~2–3% systematic error near the
critical point (Perry/McGill vs. REFPROP, see Priority 4 below) would
understate the real output uncertainty: the quoted confidence interval
would reflect input-parameter noise only, not the model-form error
already present in the property correlations. Fixing the property
backbone first is the more useful order of operations.

## Priority 3b — Non-circular grain port shapes (star, wagon-wheel)

**Motivation.** Split off from Priority 3 above during implementation
(September 2026): `grain_sizing.py` supports only circular ports
(single- or multi-port). Non-circular cross-sections are a real, widely
used technique to increase initial burning perimeter (hence fuel mass
flow) without the geometric complexity of many small circular ports.

**Why deferred separately.** Unlike a circular port, a non-circular
port's burning perimeter does not stay self-similar as it regresses --
sharp concave features (e.g. the inner points of a star) round off at
a different rate than convex ones, and the perimeter-vs-burned-web
relationship for an arbitrary shape has no simple closed form. Rigorous
treatment needs either (a) published closed-form perimeter-vs-web
formulas for specific classical shapes -- solid-rocket grain design
literature has solved this for some standard shapes (star, wagon-wheel,
etc.), analogous in spirit to how this project sourced the Henry-Fauske
equations from a primary source rather than approximating from memory
-- or (b) a numerical burnback simulation (e.g. level-set or polygon
offsetting), which is a substantially larger implementation than
anything else in this priority list.

**Proposed approach.** Source specific classical-shape formulas (star
grain is the most commonly documented case) before writing any code;
implement as an additional port-shape option in `grain_sizing.py`
alongside the existing circular case, not a replacement for it.

---

## Priority 4 — Improved N₂O thermophysical properties -  RESOLVED (September 2026)

**What was done.** `n2o_properties.py`'s thermodynamic functions now call
**CoolProp** (Bell et al., 2014), implementing the Lemmon & Span (2006)
equation of state for N₂O, replacing the closed-form Perry/McGill
correlations and McGill Table A.1.

**Confirmed with the real package (September 2026).** `pip install CoolProp`
succeeded on Eduardo's machine (CoolProp 8.0.0) and was independently
reproduced by Claude in a second sandboxed environment (also 8.0.0,
network access to PyPI). Both runs agree to the displayed precision. The
full test suite (262 tests) passes, `python src/model/n2o_properties.py`'s
self-check matches the literature reference values it targets, and
`validation/waxman_2013_validation.py` has been re-run end to end — see
`validation/waxman_2013_results.md` for the regenerated report.

**One test needed a genuine fix, not just a re-run.**
`tests/test_n2o_properties.py::TestPsat::test_agrees_with_perry_correlation`
failed at `rel_tol=0.04`: the real CoolProp P_sat differs from the Perry
correlation by up to **4.81 %** at T = 230 K (far from the critical
point), decreasing smoothly to about 0.9 % near T_MAX. This is the
opposite trend from what the module's own docstring claimed ("Perry/McGill
carried ~2–3 % error near the critical point") — the real data shows the
Perry correlation is *better* near the critical point and *worse* well
below it. The test's tolerance was widened to `rel_tol=0.06` (with a
comment recording this finding) rather than silently loosened without
explanation; the docstring in `n2o_properties.py` should be corrected to
match this the next time that module is touched (not yet done in this
pass — see "Audit follow-ups").

**Numbers that moved.** Every figure downstream of P_sat/h_fg shifted by
roughly 1–5 %, in the direction the known Perry/McGill error would
predict. Regenerated so far: `validation/waxman_2013_results.md` (full
rewrite), `examples/example_01_sizing.md`, `examples/example_02_design.md`.
`examples/example_03_flashing.md` could **not** be regenerated as-is — see
the Priority 2 finding above. `docs/04_implementation.md`'s per-module
numeric callouts (§4.1's literature self-check numbers, §4.4's isolated
Dyer/HEM validation figures) have not yet been re-verified line by line
in this pass.

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

**Motivation.** `export.py` only reproduces, as a PDF, results already
presented in the interface. The export should also be available as .JSON
or .CSV files, so test results can feed a global motor testing or CAD
modeling workflow.

**Proposed approach.** A decoupled, schema-driven export architecture
using a unified Data Transfer Object (DTO) to stream simulation
parameters, transient blowdown datasets and hole-pattern geometric
vectors into standardized JSON state files, CSV numerical tables,
CAD-compatible coordinate scripts and automated PDF engineering reports.

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

