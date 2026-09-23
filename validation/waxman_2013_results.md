# Validation Report — Waxman (2013/2014) Dataset (Updated)
## Injector Model: SPI + Dyer/NHNE + Coupled Solver + Henry-Fauske Diagnostic

**Date:** September 2026 (regenerated with the real CoolProp package installed,
by running `validation/waxman_2013_validation.py` on Eduardo's machine and
independently cross-checked by Claude in a second environment with CoolProp
8.0.0 installed — both runs agree to the displayed precision)
**Author:** Eduardo Costa, Instituto Superior Técnico
**Model state:** Full coupled solver, two-phase line model, CoolProp
(Lemmon & Span 2006) properties, Henry-Fauske choking diagnostic

---

> **Property-backend status: CONFIRMED.** Earlier versions of this report
> carried a "Property-backend note" saying the CoolProp integration was
> written and tested against a table-interpolation stand-in only, because
> the real package could not be installed in the development sandbox. That
> is now resolved: `pip install CoolProp` succeeded (CoolProp 8.0.0), the
> full test suite (262 tests) passes, and every number below was produced
> by actually running the code with the real Lemmon & Span (2006) equation
> of state. The figures below **replace** the earlier (Perry/McGill-backed)
> validation; they moved by 1–2 percentage points, in the direction expected
> from the ~2–5 % known error of the old correlations (see Section 8).

## 1. Purpose

Definitive validation of the complete injector model against the Waxman
(2013/2014) dataset — the only open-access experimental dataset in the
correct domain (supercharged N₂O, QF_upstream = 0) found for this project.

**Scope of the evidence, stated up front.** The four operating points
below have injector pressure drops of 8–14 bar. They validate the Dyer
model *in that band, for a short and wide upstream line*. They do **not**
validate the model at the 20–50 bar drops typical of motor designs, nor
the two-phase line model, nor the HEM two-phase-inlet path (Section 9).

---

## 2. References

**Primary dataset:**
Waxman, B. S. (2014). *An Investigation of Injectors for Use with High
Vapour Pressure Propellants with Applications to Hybrid Rockets.*
PhD thesis, Stanford University. See also Waxman, B.S., Zimmerman, J.E.,
Cantwell, B., & Zilliac, G. (2013), AIAA 2013-3636, the source of the
model equations (Eq. 5, Eq. 9) and of the discharge-coefficient data.

**Tabulated operating points:**
Niño, E. V., and Razavi, M. R. (2019). *Design of Two-Phase Injectors
Using Analytical and Numerical Methods with Application to Hybrid Rockets.*
AIAA 2019-4154. (Table 4: four tabulated points; Table 3: model error
summary — Dyer MAPE = 3.91% with Cd = 0.63 from their correlation.)

**Henry-Fauske non-equilibrium critical flow:**
Henry, R.E. & Fauske, H.K. (1971). *The Two-Phase Critical Flow of
One-Component Mixtures in Nozzles, Orifices, and Short Tubes.* ASME J.
Heat Transfer, 93(2), 179-187. Equations transcribed from: Simoneau,
R.J., Henry, R.E., Hendricks, R.C. & Watterson, R. (1971). *Two-Phase
Critical Discharge of High Pressure Liquid Nitrogen.* NASA TM X-67863.

---

## 3. Model State at Validation

| Component | Status |
|---|---|
| Coupled solver | ✅ damped fixed-point, α=0.5, tol=1e-4 |
| Feed-line: single-phase | ✅ Darcy-Weisbach, Swamee-Jain, μ_l(T) from Table A.4 |
| Feed-line: two-phase HEM | ✅ ρ_mix, μ_mix per segment (not exercised here) |
| Injector: SPI | ✅ Bernoulli |
| Injector: Dyer/NHNE | ✅ corrected weights (Solomon 2011) |
| Injector: HEM two-phase inlet | ✅ x_inlet from feed-line flash (not exercised here) |
| N₂O properties | ✅ **CoolProp (Lemmon & Span 2006), real package, confirmed** |
| HEM critical flow (equilibrium) | ✅ isenthalpic + isentropic, standalone |
| **Henry-Fauske critical flow (non-equilibrium)** | ✅ **standalone diagnostic, side-by-side with Dyer, not auto-applied** |

---

## 4. Test Configuration

| Parameter | Value | Source |
|---|---|---|
| D | 1.50 mm | Waxman Table 1 |
| L/D | 12.3 | Waxman Table 1 |
| Inlet | Square edge | Waxman injector no. 2 |
| T₁ | 280 K | Niño & Razavi Table 4 |
| P₁ | 4.36 MPa | Niño & Razavi Table 4 |
| Cd | 0.65 | Conservative literature estimate (Waxman's measured value for this injector: 0.71) |
| Line | ID=25.4mm, L=50mm | Waxman upstream chamber |

**P_sat(T1) with the real CoolProp backend:** 37.068 bar, vs. Niño & Razavi's
37.40 bar — model error **−0.9 %**. (With the old Perry/McGill correlation
this comparison used to read appreciably worse; the CoolProp switch
directly improves the single number that sets the Dyer κ denominator for
every case below.) Model supercharge = 6.53 bar vs. Niño & Razavi's 6.20 bar.

**Note on `waxman_2013_experimental_data.csv`.** That file records a raw
extraction of Waxman's paper (used, for instance, for the measured Cd
values) with reference upstream conditions P₁ = 4.96 MPa and supercharge
1.26 MPa, whereas the operating points above use P₁ = 4.36 MPa and
supercharge 0.62 MPa (Niño & Razavi Table 4). The CSV is **not read by any
code**. This discrepancy is still unresolved (see `docs/future_work.md`,
"Audit follow-ups") — out of scope for this regeneration.

---

## 5. Results (CoolProp-confirmed)

| Case | ΔP [bar] | m_exp [g/s] | m_dot [g/s] | Error | Regime | Choked? |
|---|---|---|---|---|---|---|
| Pre-critical | 8.40 | 44.0 | 42.91 | −2.5% | Dyer | No |
| Critical | 9.80 | 46.5 | 45.34 | −2.5% | Dyer | No |
| Post-critical 1 | 10.90 | 47.5 | 46.97 | −1.1% | Dyer | No |
| Post-critical 2 | 13.70 | 48.0 | 50.37 | +4.9% | Dyer | No |
| **Mean** | — | — | — | **−0.3%** | — | — |

**MAPE = 2.76%** (previous, Perry/McGill-backed figure: 3.51%; Niño & Razavi
reference: 3.91% with Cd=0.63). The injector inlet pressure sits at 43.600 bar
for all four cases (the 50 mm / 25.4 mm ID upstream chamber has negligible
losses, as designed), so the only thing driving the per-case error is the
Dyer model itself at each chamber pressure.

All 4 cases within ±5%. All 4 cases within ±10%. The coupled solver converged
in 10–11 iterations for every case.

The experimental values are read off a graph (Niño & Razavi Fig. 2), with
~±3% read-off uncertainty, which bounds the achievable accuracy of this
comparison.

**What changed vs. the Perry/McGill-backed report.** Every individual case's
predicted mass flow moved by roughly 1–3 g/s (about 1–2 %), and in the same
direction the ~2–5 % documented Perry/McGill P_sat error would predict.
The corrected-weights Dyer model, the coupled solver logic, and the
Henry-Fauske diagnostic are all unchanged — only the property backend
feeding them changed. The net effect here is a *better* MAPE (2.76 % vs.
3.51 %), but that should be read as this specific test case improving, not
as a general claim that CoolProp always reduces error — see Section 8.

---

## 6. Equilibrium HEM Critical Flow Reference

`hem_critical_flow()` (isenthalpic, Waxman Eq. 5) gives **m_crit = 42.48 g/s**
at P₂_crit = 30.5 bar, x_crit = 0.0744. `hem_critical_flow_isentropic()`
gives **43.02 g/s** at P₂_crit = 29.8 bar, x_crit = 0.0774 — a +1.3%
difference, both retained side-by-side in the codebase.

The Dyer predictions (42.91–50.37 g/s) are above **both** equilibrium
ceilings, which is physically correct: the Dyer non-equilibrium correction
(κ weighting toward SPI) accounts for the fact that real injectors do not
reach full thermodynamic equilibrium inside the orifice. The experimental
values (44.0–48.0 g/s) confirm this — but this also means neither
equilibrium ceiling is the right bound to cap Dyer with (see Section 7).

---

## 7. Henry-Fauske Non-Equilibrium Critical Flow

**Result at Waxman conditions (CoolProp-confirmed):** `henry_fauske_critical_flow()`
gives **m_dot_crit = 52.15 g/s** at P2_crit = 31.76 bar, x_crit(equilibrium) =
0.0576, N = 0.412 — above the equilibrium HEM ceiling (correct direction:
non-equilibrium exceeds equilibrium) and above all 4 validated Dyer
predictions (42.91–50.37 g/s), confirmed by the `choked = False` column in
Section 5.

**Design decision — diagnostic, not automatic cap.** Unchanged from the
previous version of this report: at operating points further from the
Waxman geometry the Henry-Fauske ceiling *does* bind, and there is still no
experimental confirmation in this project's validation set at conditions
where the ceiling actually changes the answer. `dyer_mass_flow()` therefore
returns `m_dot_Dyer` (always unchanged) alongside `m_dot_crit_HF` and a
`choked` boolean, rather than silently overriding the headline number.

---

## 8. Error Sources

**(a)** P_sat correlation: with the CoolProp backend this is no longer an
error source at T1 = 280 K (the equation of state now matches Niño &
Razavi's reference value to −0.9 %, see Section 4) — this line item, carried
over from the Perry/McGill era, is now resolved for this specific case.

**(b)** h_fg: now from the same equation of state as P_sat (CoolProp), so
the ~3–5 % Perry/McGill vs. NIST latent-heat gap no longer applies here
either.

**(c)** Cd = 0.65 vs. Waxman measured Cd = 0.71 for this injector. Using
Cd = 0.71 would still shift every case by a similar amount to before; a
team-calibrated Cd from a water cold-flow test is recommended for
production use. Unaffected by the property-backend change.

**(d)** Experimental read-off uncertainty: ~±3%. m_dot values estimated from
Niño & Razavi Fig. 2 (graph). This bounds the achievable accuracy of this
comparison and is now, along with (c), the dominant remaining error source.

**(e)** The Henry-Fauske diagnostic (Section 7) is validated only in the
sense that it does not perturb the 4 known-good points — it has not itself
been checked against experimental data in the regime where it binds.

**Independent numerical cross-check.** These figures were reproduced from
scratch in a second Python environment (CoolProp 8.0.0) by Claude, from the
same source files, and matched Eduardo's own terminal output to the
displayed precision (see the September 2026 audit chat). This is a
software cross-check (same equations, same equation-of-state library, two
independent installs) — it is not a second, independent experimental
dataset, so it does not add new evidence beyond what is in (a)–(d) above.

---

## 9. Conclusions and Limits of the Evidence

The model is validated in its correct domain (supercharged liquid inlet,
short upstream line) **at injector pressure drops of 8–14 bar** with
MAPE = 2.76%, better than both the Niño & Razavi reference (3.91%) and the
earlier Perry/McGill-backed run of this same model (3.51%). The
Henry-Fauske non-equilibrium ceiling is confirmed not to perturb this
validated result and is available as a diagnostic warning for operating
points outside it.

What this report does **not** establish:

1. **Accuracy at 20–50 bar pressure drop.** No experimental data in this
   project covers it. At such conditions the Henry-Fauske ceiling is
   exceeded by the Dyer prediction at some worked examples (see
   `examples/`), i.e. two models disagree.
2. **The two-phase feed-line model and the HEM two-phase-inlet injector
   model.** Not exercised by the Waxman geometry (negligible line losses,
   liquid inlet). Tested only by unit tests. **New finding from this
   CoolProp confirmation run:** with the real (lower) CoolProp P_sat, the
   `examples/example_03_flashing.md` scenario no longer has the same
   subcooling margin it had under Perry/McGill, and the coupled solver
   fails to converge there even at very small damping (α down to 0.05) —
   the fixed-point iteration oscillates between the Dyer (liquid-inlet)
   and HEM (two-phase-inlet) branches rather than settling. This is a
   genuine consequence of the already-documented Dyer→HEM discontinuity at
   the flashing threshold (`docs/future_work.md`, Priority 2, item 3)
   showing up as non-convergence right at a borderline operating point,
   not a bug introduced by the CoolProp switch. Example 3 needs a
   deliberate fix (a different operating point with clearer margin, or a
   convergence strategy for the borderline case) before its numbers can be
   regenerated — left open pending Eduardo's input.
3. **A line with significant losses feeding the injector.** The coupled
   solver's benefit is realised there (see Examples 1–2, regenerated with
   CoolProp, both converge cleanly), but it has not been checked against a
   measured tank-to-chamber flow.

Predictions in these regimes should be read as model estimates with
unquantified error, to be confirmed by a cold-flow or hot-fire measurement.

---

## 10. Files

| File | Location |
|---|---|
| `waxman_2013_results.md` | `validation/` |
| `waxman_2013_validation.py` | `validation/` (runs from any directory) |
| `waxman_2013_experimental_data.csv` | `validation/` (raw extraction; not read by code) |
