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
  as an automatic cap; validated against Waxman (does not perturb MAPE = 2.76%) (September 2026);
  wired into the interface (warning badge in both modes, comparison-chart line, PDF note)  
- Fuel grain sizing via the Marxman regression rate correlation (grain_sizing.py) -- initial port
  radius, fuel mass flow, multi-port geometry, first-order conservative burnback estimate; a and n
  are required user inputs, entered as a regression-rate data point converted by
  a_from_reference_rate(); results outside a plausible port-radius range are refused; density
  defaults for paraffin/HTPB/ABS/PMMA (September 2026); integrated in the interface 
- Definitive validation against Waxman (2013/2014): MAPE = 2.76%, all 4 cases within +/-5%  
- Sensitivity tornado plot  
- Dyer formula weight correction (Solomon 2011)  
- Design-mode sizing logic extracted into a tested function (full_system.design_injector_area),
  fixing a crash when the flow is single-phase through the orifice   
- Explicit range errors for Table A.3 (mu_vapor_sat) and for T_sat(P) outside the correlation  
- CoolProp implementation and validation, confirming, at least for the waxman results, a more accurate result
- Extended experimental validation against digitised Waxman (2013) Figs. 11-16 (injector 3, nine
  supercharge levels, dP up to 46 bar) -- see Priority 1 and Priority 2 below for what this changed

---

## Priority 1 — Non-equilibrium critical flow ceiling for Dyer - PARTIALLY RE-OPENED (September 2026)

**Original resolution (kept for context).** `henry_fauske_critical_flow()` was
added to `injector_two_phase.py`, implementing the simplified Henry-Fauske
(1971) non-equilibrium critical mass flux model (Simoneau et al. 1971, Eqs.
2-5). At the four originally-validated Waxman points it sits above every
Dyer prediction (`choked = False` in all four cases), so it was surfaced as
a side-by-side diagnostic (`m_dot_crit_HF`, `choked` flag) rather than an
automatic cap -- see the original reasoning below, which still explains
*why* an automatic cap is risky in general.

**New finding (September 2026, extended validation).** Eduardo digitised
Waxman (2013) Figs. 11-16 directly from the source PDF -- the full
injector-3 mass-flow map across nine supercharge levels, dP up to 46 bar
(104 usable points, `validation/digitized/`, see
`validation/waxman_2013_results.md` Part B). Against this much larger
dataset:

- The Henry-Fauske ceiling is exceeded (`choked = True`) at **31 of 64**
  two-phase points -- previously it had *never* bound at any validated
  point, which is what motivated keeping it diagnostic-only.
- At **30 of those 31** points, capping the Dyer prediction at the ceiling
  (`min(m_dot_Dyer, m_dot_crit_HF)`) moves the prediction *closer* to the
  experimental value, not farther.
- At 28 of the 31, the true experimental value sits *below* the ceiling
  too, meaning even a tighter non-equilibrium correction than
  Henry-Fauske would help further.
- The pattern that controls where this matters is **tank supercharge
  (subcooling margin)**, not injector pressure drop: MAPE is ~2% for
  Dyer alone when supercharge >= 200 psi (~1.38 MPa, ~14 bar), but climbs
  to ~12% (worst single point: 33%) below that, and this is exactly the
  regime where the ceiling starts to bind.

**Why this is not yet promoted to an automatic cap.** The cap does not
close the whole gap (MAPE with capping is 4.9% vs 6.9% uncapped, still
worse than the high-supercharge regime's own 2%), and at low supercharge
even the capped prediction still overshoots most points. This looks like
a genuine partial correction, not the fix. The four originally-validated
points (Part A) are all comfortably at high supercharge and still show
`choked = False`, so nothing published earlier is disturbed by this
finding.

**What is now genuinely still open (revised).**
1. Investigate a supercharge-dependent (rather than binary choked/not)
   correction -- the data suggests the *size* of the needed correction
   scales with how close the tank sits to saturation, not just whether
   the ceiling is crossed.
2. Re-run this same digitisation exercise for injector geometries other
   than injector 3 (only a 4-point check exists for injector 2, the
   square-edge geometry actually used in Part A) before generalising.
3. If/when a supercharge-aware correction is validated, reconsider
   promoting it from diagnostic to (still supercharge-gated, not
   unconditional) automatic cap.
4. `apply_choking_limit()` (the retracted equilibrium-cap plan, distinct
   from Henry-Fauske) remains deprecated and unused; delete it once no
   external script depends on it.

**Original design rationale (kept for context -- still valid in general,
even though the "never binds" premise it was partly based on is now
outdated).** At operating points further from the Waxman geometry the
Henry-Fauske ceiling *does* bind: 17% below the Dyer prediction in
`examples/example_01_sizing.md`, 10% below the target in
`example_02_design.md`. Silently overriding `m_dot_Dyer` with a value
whose net effect on accuracy was, until now, never checked against data
would have been the wrong precedent regardless of how the individual
numbers happened to compare -- the extended validation above is the
"if/when experimental data becomes available" moment this section
originally deferred to, and it says "helps, but not enough to fully
trust yet", not "don't bother".

**History.** This priority went through three framings before the
September 2026 update: (1) originally scoped as "build the isentropic HEM
scan and cap Dyer with it" -- shown numerically to be wrong (Dyer
legitimately exceeds the *equilibrium* ceiling by design); (2) redirected
to the genuine *non-equilibrium* Henry-Fauske model, implemented from a
primary source and validated as a diagnostic against the 4-point Waxman
set (`choked = False` everywhere, so harmless but also untested where it
would matter); (3) the extended digitised dataset above is the first time
the ceiling has actually been checked against data at conditions where it
fires, and shows real (if partial) predictive value.

## Priority 2 — Full-system experimental validation - PARTIALLY RESOLVED (September 2026)

**What changed.** The main gap this priority flagged -- "pressure drops of
20-50 bar are outside the validated band" -- is now partially closed.
The digitised Waxman Fig. 13 dataset (Priority 1, above) covers injector-3
pressure drops from well below 1 bar up to 46 bar, at nine supercharge
levels. Counter-intuitively, the **largest** pressure drops (30-46 bar,
the range every worked example in `examples/` actually uses) give the
**best** Dyer agreement (MAPE 2.5%) once the tank has a reasonable
supercharge -- the 14-30 bar band at low supercharge is where the model
struggles most (see `validation/waxman_2013_results.md`, Part B3/B4). The
controlling variable is tank subcooling margin, not pressure drop, which
was not previously known.

**What is still genuinely unvalidated (unchanged from before).**
1. **The coupled feed-line / injector solver** with a line that matters
   (both the original and the new Waxman data use a short, wide upstream
   line with negligible losses -- the coupling itself is still untested
   against a measured tank-to-chamber flow with real line losses).
2. **The two-phase inlet path** (flashing in the line): the two-phase
   line model and the HEM two-phase-inlet injector model are implemented
   and unit tested only. The switch from Dyer to HEM at the flashing
   threshold remains discontinuous (HEM at x_inlet -> 0 gives roughly
   half of the Dyer flow at the same conditions -- see
   `examples/example_03_flashing.md`); a smooth transition, or at least
   a quantified uncertainty band around the threshold, is worth
   designing once data exist. The digitised Waxman data does not cover
   this regime (all points there are liquid at the injector inlet).
3. **Geometries other than injector 3.** The extended dataset used only
   the rounded-inlet 1.5 mm injector (Waxman's "injector 3"); the
   square-edge injector 2 used in Part A (the original 4-point
   validation) still has only those 4 points at one supercharge level.

**What is needed next.** A team's own cold-flow (water or N2O) data,
specifically spanning a range of tank supercharge at a fixed, large
pressure drop, would directly test the Priority 1 finding above (that
supercharge, not dP, controls Dyer's error) on a different injector and
rig than Waxman's.

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

**Confirmed with the real package (September 2026).** The full test suite (262 tests) 
passes, `python src/model/n2o_properties.py`'s
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
