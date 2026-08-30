# Future Work

Extensions identified during this project's development, ordered by
technical importance per the project roadmap document.

---

## Implemented (removed from future work)

- Two-phase inlet via isenthalpic flash (x_inlet from feed line)  ✅
- HEM with two-phase inlet enthalpy  ✅
- Coupled feed-line / injector solver (iterative, damped fixed-point)  ✅
- Sensitivity tornado plot  ✅
- Dyer formula weight correction (Solomon 2011)  ✅

---

## Priority 2 — Physically consistent two-phase feed-line model

**Motivation.** Once flashing is detected in the feed line, the current model
continues to use liquid-phase properties (density, viscosity) for the
remaining pressure-drop calculation. After the flashing onset, the fluid is
a two-phase mixture with different density and pressure-gradient behaviour.

**Why not yet.** Implementing a correct two-phase pressure-drop model
(HEM homogeneous, or Lockhart-Martinelli type) requires knowing the local
void fraction at each point after the onset — which in turn requires the
coupled solver (now implemented) to track x(s) along the line, not just at
the inlet. This is the natural next step now that the coupled solver exists.

**Proposed approach.** Extend `feed_line.py` to track enthalpy (and
therefore x) segment by segment after the flashing onset, and use the
HEM two-phase friction multiplier (Φ²_lo) for pressure-drop in the two-phase
region.

---

## Priority 3 — Two-phase choking limit

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

## Priority 4 — Full-system experimental validation

**Current state.** The Dyer injector model is validated against Waxman
(2013/2014) with mean error −1.9% at moderate pressure drops. The coupled
solver is new and not yet validated against full-system data.

**What is needed.** Experimental data covering the full path
(tank → line → injector) with known geometry, discharge coefficient,
and measured mass flow. The validation should be re-run after the coupled
solver since the one-pass results used previously may differ slightly.

---

## Priority 7 — Quantitative uncertainty analysis (Monte Carlo)

**Motivation.** The current sensitivity tornado shows which input matters
most. A Monte Carlo analysis would quantify the output uncertainty given
realistic input distributions, producing a confidence interval on mass flow.

**Proposed approach.** Sample input distributions (normal for temperature
and pressure, uniform for Cd and geometry), run the coupled solver for
each sample, and report the 5th–95th percentile of the mass-flow distribution.

---

## Priority 8 — Improved N₂O thermophysical properties

**Motivation.** The Perry/McGill correlations have ~2–3% error near the
critical point vs. REFPROP. CoolProp (open source, no licence) provides
REFPROP-quality properties and would improve accuracy in the critical region.

**Why deferred.** The current property implementation is sufficient for
design-regime accuracy. CoolProp integration would also enable the
isentropic maximum search needed for choking (Priority 3).

---

## Priority 9 — Tank thermal model

**Motivation.** Tank temperature is currently a direct user input. A
thermal model would predict T_tank from ambient conditions, solar irradiance,
and tank geometry.

**Proposed approach.** Retrieve weather data via Open-Meteo API (free, no
key required) and implement a simplified steady-state tank-ambient heat
balance as a first approximation before a full transient model.

---

## Priority 10 — Transient / blowdown model

**Motivation.** The current model is steady-state. A transient model would
predict how mass flow, tank pressure, and temperature evolve over the burn.

**Why deferred.** Requires a reliable steady-state coupled solver (now
implemented) as the inner loop, plus a tank thermodynamic model (Priority 9).

---

*Items implemented in this version are listed at the top. New items are added
here as they are identified. See `roadmap/Future_Improvements_Roadmap.md`
for the full priority order and rationale.*
