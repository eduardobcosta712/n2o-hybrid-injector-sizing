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
- Sensitivity tornado plot  ✅
- Dyer formula weight correction (Solomon 2011)  ✅

---

## Priority 1 — Two-phase choking limit

**Motivation.** The Dyer and HEM formulas use
$\dot{m} = C_d A \sqrt{2 \rho \Delta P}$, which grows without bound with
$\Delta P$. Real two-phase critical flow reaches a physical maximum when
the mixture velocity equals the two-phase speed of sound. At large pressure
drops ($\Delta P \gg 50$ bar), the model over-predicts without bound.

**Why not yet.** The HEM critical flow requires finding the maximum of
$\dot{m}_{HEM}(P_2)$ along a constant-entropy path — which requires
entropy data for N₂O not present in the current Perry/McGill dataset.

**Proposed approach.** Integrate CoolProp (open source) for entropy data,
enabling the isentropic maximum search (Waxman 2013, Eq. 5). Alternatively,
implement the Henry-Fauske critical flow model, which does not require
entropy data.

---

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

## Priority 4 — Quantitative uncertainty analysis (Monte Carlo)

**Motivation.** The current sensitivity tornado shows which input matters
most. A Monte Carlo analysis would quantify the output uncertainty given
realistic input distributions, producing a confidence interval on mass flow.

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

## Priority 6 — Tank thermal model

**Motivation.** Tank temperature is currently a direct user input. A
thermal model would predict T_tank from ambient conditions, solar irradiance,
and tank geometry.

**Proposed approach.** Retrieve weather data via Open-Meteo API (free, no
key required) and implement a simplified steady-state tank-ambient heat
balance as a first approximation before a full transient model.

---

## Priority 7 — Transient / blowdown model

**Motivation.** The current model is steady-state. A transient model would
predict how mass flow, tank pressure, and temperature evolve over the burn.

**Why deferred.** Requires a reliable steady-state coupled solver (now
implemented) as the inner loop, plus a tank thermodynamic model (Priority 9).

---


