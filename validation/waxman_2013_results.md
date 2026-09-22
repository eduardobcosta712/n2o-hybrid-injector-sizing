# Validation Report — Waxman (2013/2014) Dataset (Updated)
## Injector Model: SPI + Dyer/NHNE + Coupled Solver + Henry-Fauske Diagnostic

**Date:** September 2026 (regenerated in the project audit by running `validation/waxman_2013_validation.py`)
**Author:** Eduardo Costa, Instituto Superior Técnico
**Model state:** Full coupled solver, two-phase line model, NIST properties, Henry-Fauske choking diagnostic

---

> **Property-backend note (September 2026).** The MAPE and every figure below were computed with the Perry/McGill + NIST property set used before the CoolProp integration (`docs/future_work.md`, Priority 4). They have not yet been regenerated with the real CoolProp package (unavailable in the development sandbox); expect small shifts once `python validation/waxman_2013_validation.py` is re-run with CoolProp installed.

## 1. Purpose

Definitive validation of the complete injector model against the Waxman
(2013/2014) dataset — the only open-access experimental dataset in the
correct domain (supercharged N₂O, QF_upstream = 0) found for this project.

This report supersedes the earlier validation report which used the
one-pass model without the coupled solver, updated properties, or the
Henry-Fauske choking diagnostic.

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

**Henry-Fauske non-equilibrium critical flow (September 2026 addition):**
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
| N₂O properties | ✅ McGill/Perry (A.1) + NIST A.3 (μ_v) + NIST A.4 (C_pl, μ_l, entropy) |
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

**Note on `waxman_2013_experimental_data.csv`.** That file records a raw
extraction of Waxman's paper (used, for instance, for the measured Cd
values) with reference upstream conditions P₁ = 4.96 MPa and supercharge
1.26 MPa, whereas the operating points above use P₁ = 4.36 MPa and
supercharge 0.62 MPa (Niño & Razavi Table 4). The CSV is **not read by any
code**. Whether the two sets describe different tests of the same paper
or one of them is mis-transcribed has not been resolved (see
`docs/future_work.md`, "Audit follow-ups").

---

## 5. Results

| Case | ΔP [bar] | m_exp [g/s] | m_dot [g/s] | Error | Regime | Choked? |
|---|---|---|---|---|---|---|
| Pre-critical | 8.40 | 44.0 | 42.25 | −4.0% | Dyer | No |
| Critical | 9.80 | 46.5 | 44.60 | −4.1% | Dyer | No |
| Post-critical 1 | 10.90 | 47.5 | 46.20 | −2.7% | Dyer | No |
| Post-critical 2 | 13.70 | 48.0 | 49.55 | +3.2% | Dyer | No |
| **Mean** | — | — | — | **−1.9%** | — | — |

**MAPE = 3.51%** (Niño & Razavi reference: 3.91% with Cd=0.63) — **identical to the
pre-Henry-Fauske result**: the diagnostic ceiling (Section 7) sits above all 4
`m_dot` values, so `choked = False` throughout and nothing changed numerically. It is
also identical before and after the audit's change to a temperature-dependent
liquid viscosity in the line (the Waxman line is short and wide, so friction is negligible).

All 4 cases within ±5%. All 4 cases within ±10%. The coupled solver converged in 10–11 iterations.

The experimental values are read off a graph (Niño & Razavi Fig. 2), with
~±3% read-off uncertainty, which bounds the achievable accuracy of this comparison.

---

## 6. Equilibrium HEM Critical Flow Reference

`hem_critical_flow()` (isenthalpic, Waxman Eq. 5) gives **m_crit = 41.05 g/s**
at P₂_crit = 30.37 bar, x_crit = 0.086. `hem_critical_flow_isentropic()`
(added September 2026, the thermodynamically correct path for an
*equilibrium* choking condition) gives **41.64 g/s** — a +1.43% difference,
both retained side-by-side in the codebase.

The Dyer predictions (42.25–49.55 g/s) are above **both** equilibrium
ceilings, which is physically correct: the Dyer non-equilibrium correction
(κ weighting toward SPI) accounts for the fact that real injectors do not
reach full thermodynamic equilibrium inside the orifice. The experimental
values (44.0–48.0 g/s) confirm this — but this also means neither
equilibrium ceiling is the right bound to cap Dyer with (see Section 7).

---

## 7. Henry-Fauske Non-Equilibrium Critical Flow (added September 2026)

**Motivation.** Section 6 shows Dyer legitimately exceeds the *equilibrium*
HEM ceiling at every validated point — expected, correct behaviour, not an
error to be capped. The physically appropriate ceiling for a *non-equilibrium*
prediction like Dyer is therefore a genuine non-equilibrium critical flow
model: Henry-Fauske (1971).

**Result at Waxman conditions:** `henry_fauske_critical_flow()` gives
**m_dot_crit = 50.67 g/s** at the same geometry — above the equilibrium HEM
ceiling (correct direction: non-equilibrium exceeds equilibrium) and above
all 4 validated Dyer predictions (42.25–49.55 g/s), confirmed by the
`choked = False` column in Section 5.

**Design decision — diagnostic, not automatic cap.** At operating points
further from the Waxman geometry (the worked examples: 17% below the Dyer
prediction in Example 1, 10% below the target in Example 2), the
Henry-Fauske ceiling *does* bind, and there is currently no experimental
confirmation in this project's validation set at conditions where the
ceiling actually changes the answer. `dyer_mass_flow()` therefore returns
`m_dot_Dyer` (always unchanged) alongside `m_dot_crit_HF` and a `choked`
boolean, rather than silently overriding the headline number with a value
that, while theoretically well-founded and correctly sourced, is
empirically unconfirmed outside the Waxman-validated regime.

---

## 8. Error Sources

**(a)** P_sat correlation: +1.5% at 280 K (Perry/McGill vs. NIST). Shifts κ
denominator, small bias in Dyer weights.

**(b)** h_fg from Perry interpolation vs. NIST: ~3–5%. Dominant error source
in Cases III–IV (larger ΔP). Would be resolved by CoolProp integration.

**(c)** Cd = 0.65 vs. Waxman measured Cd = 0.71 for this injector. Using
Cd = 0.71 gives errors of +1.9% to −6.6%. A team-calibrated Cd from a
water cold-flow test is recommended for production use.

**(d)** Experimental read-off uncertainty: ~±3%. m_dot values estimated from
Niño & Razavi Fig. 2 (graph). This bounds the achievable accuracy.

**(e)** The Henry-Fauske diagnostic (Section 7) is validated only in the
sense that it does not perturb the 4 known-good points — it has not itself
been checked against experimental data in the regime where it binds.

---

## 9. Conclusions and Limits of the Evidence

The model is validated in its correct domain (supercharged liquid inlet,
short upstream line) **at injector pressure drops of 8–14 bar** with
MAPE = 3.51%, better than the Niño & Razavi reference (3.91%). The Henry-Fauske
non-equilibrium ceiling is confirmed not to perturb this validated result and is
available as a diagnostic warning for operating points outside it.

What this report does **not** establish:

1. **Accuracy at 20–50 bar pressure drop.** This is the range of typical
   motor designs, and no experimental data in this project covers it. (An
   earlier version of this report asserted accuracy "competitive with the
   state of the art" at these pressures; that claim was not supported by the
   data and has been withdrawn.) At such conditions the Henry-Fauske ceiling
   is exceeded by the Dyer prediction, i.e. two models disagree.
2. **The two-phase feed-line model and the HEM two-phase-inlet injector
   model.** Not exercised by the Waxman geometry (negligible line losses, liquid
   inlet). Tested only by unit tests (`tests/test_feed_line.py::TestTwoPhaseLineModel`,
   `tests/test_full_system.py::TestFlashingInLine`).
3. **A line with significant losses feeding the injector.** The coupled solver's
   benefit is realised there, but it has not been checked against a
   measured tank-to-chamber flow.

Predictions in these regimes should be read as model estimates with unquantified
error, to be confirmed by a cold-flow or hot-fire measurement.

---

## 10. Files

| File | Location |
|---|---|
| `waxman_2013_results.md` | `validation/` |
| `waxman_2013_validation.py` | `validation/` (runs from any directory; also run by the CI workflow) |
| `waxman_2013_experimental_data.csv` | `validation/` (raw extraction; not read by code) |
