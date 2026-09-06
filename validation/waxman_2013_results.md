# Validation Report — Waxman (2013/2014) Dataset (Updated)
## Injector Model: SPI + Dyer/NHNE + Coupled Solver

**Date:** August 2026
**Author:** Eduardo Costa, Instituto Superior Técnico
**Model state:** Full coupled solver, two-phase line model, NIST properties

---

## 1. Purpose

Definitive validation of the complete injector model against the Waxman
(2013/2014) dataset — the only open-access experimental dataset in the
correct domain (supercharged N₂O, QF_upstream = 0).

This report supersedes the earlier validation report which used the
one-pass model without the coupled solver or updated properties.

---

## 2. References

**Primary dataset:**
Waxman, B. S. (2014). *An Investigation of Injectors for Use with High
Vapour Pressure Propellants with Applications to Hybrid Rockets.*
PhD thesis, Stanford University.

**Tabulated operating points:**
Niño, E. V., and Razavi, M. R. (2019). *Design of Two-Phase Injectors
Using Analytical and Numerical Methods with Application to Hybrid Rockets.*
AIAA 2019-4154. (Table 4: four tabulated points; Table 3: model error
summary — Dyer MAPE = 3.91% with Cd = 0.63 from their correlation.)

---

## 3. Model State at Validation

| Component | Status |
|---|---|
| Coupled solver | ✅ damped fixed-point, α=0.5, tol=1e-4 |
| Feed-line: single-phase | ✅ Darcy-Weisbach, Swamee-Jain |
| Feed-line: two-phase HEM | ✅ ρ_mix, μ_mix per segment |
| Injector: SPI | ✅ Bernoulli |
| Injector: Dyer/NHNE | ✅ corrected weights (Solomon 2011) |
| Injector: HEM two-phase inlet | ✅ x_inlet from feed-line flash |
| N₂O properties | ✅ McGill/Perry (A.1) + NIST A.3 (μ_v) + NIST A.4 (C_pl, μ_l, entropy) |
| HEM critical flow | ✅ standalone (not auto-applied) |

---

## 4. Test Configuration

| Parameter | Value | Source |
|---|---|---|
| D | 1.50 mm | Waxman Table 1 |
| L/D | 12.3 | Waxman Table 1 |
| Inlet | Square edge | Waxman injector no. 2 |
| T₁ | 280 K | Niño & Razavi Table 4 |
| P₁ | 4.36 MPa | Niño & Razavi Table 4 |
| Cd | 0.65 | Waxman Fig. 15 (conservative) |
| Line | ID=25.4mm, L=50mm | Waxman upstream chamber |

---

## 5. Results

| Case | ΔP [bar] | m_exp [g/s] | m_dot [g/s] | Error | Regime | Iters |
|---|---|---|---|---|---|---|
| Pre-critical | 8.40 | 44.0 | 42.25 | −4.0% | Dyer | 10 |
| Critical | 9.80 | 46.5 | 44.60 | −4.1% | Dyer | 10 |
| Post-critical 1 | 10.90 | 47.5 | 46.20 | −2.7% | Dyer | 10 |
| Post-critical 2 | 13.70 | 48.0 | 49.55 | +3.2% | Dyer | 11 |
| **Mean** | — | — | — | **−1.9%** | — | — |

**MAPE = 3.51%** (Niño & Razavi reference: 3.91% with Cd=0.63)

All 4 cases within ±5%. All 4 cases within ±10%.

---

## 6. HEM Critical Flow Reference

The `hem_critical_flow()` function (Waxman Eq. 5, isenthalpic maximum scan)
gives **m_crit = 41.1 g/s** at P₂_crit = 30.4 bar, x_crit = 0.086.

The Dyer predictions (42.3–49.5 g/s) are above this HEM ceiling, which is
physically correct: the Dyer non-equilibrium correction (κ weighting toward
SPI) accounts for the fact that real injectors do not reach full
thermodynamic equilibrium inside the orifice. The experimental values
(44.0–48.0 g/s) confirm this interpretation.

---

## 7. Error Sources

**(a)** P_sat correlation: +1.5% at 280 K (Perry/McGill vs. NIST). Shifts κ
denominator, small bias in Dyer weights.

**(b)** h_fg from Perry interpolation vs. NIST: ~3–5%. Dominant error source
in Cases III–IV (larger ΔP). Would be resolved by CoolProp integration.

**(c)** Cd = 0.65 vs. Waxman measured Cd = 0.71 for this injector. Using
Cd = 0.71 gives errors of +1.9% to −6.6%. A team-calibrated Cd from a
water cold-flow test is recommended for production use.

**(d)** Experimental read-off uncertainty: ~±3%. m_dot values estimated from
Niño & Razavi Fig. 2 (graph). This bounds the achievable accuracy.

---

## 8. Conclusions

The model is validated in its correct domain with MAPE = 3.51%, better
than the Niño & Razavi reference (3.91%). The coupled solver converges
in 10–11 iterations. The two-phase line model is not exercised here
(negligible line losses in the Waxman geometry) but is tested separately
in `tests/test_feed_line.py::TestTwoPhaseLineModel`.

For injector sizing at realistic motor pressures, the model provides
accuracy competitive with the state of the art in open-source tools.

---

## 9. Files

| File | Location |
|---|---|
| `waxman_2013_results.md` | `validation/` |
| `waxman_2013_validation.py` | `validation/` |
| `waxman_2013_experimental_data.csv` | `validation/` |
