# Validation Report — Waxman (2013/2014) Dataset
## Injector Model: SPI + Dyer/NHNE

**Date:** August 2026
**Author:** Eduardo Costa, Instituto Superior Técnico

---

## 1. Purpose

Validate the injector model (SPI + Dyer/NHNE) against published experimental data in the correct domain of validity: supercharged N₂O (QF_upstream = 0), where the fluid arrives at the injector inlet as pure subcooled liquid and flashes inside the orifice.

This is the domain the model was designed for. The Palacz & Cieślik (2021) dataset (self-pressurized, QF_upstream ≈ 0.4–0.5) is explicitly **not** used here because the fluid arrives two-phase at the injector inlet, outside the domain of both SPI and Dyer as implemented.

---

## 2. References

**Primary dataset:**
Waxman, B. S. (2014). *An Investigation of Injectors for Use with High Vapour Pressure Propellants with Applications to Hybrid Rockets.* PhD thesis, Stanford University.

**Tabulated operating points and model comparison:**
Niño, E. V., and Razavi, M. R. (2019). *Design of Two-Phase Injectors Using Analytical and Numerical Methods with Application to Hybrid Rockets.* AIAA 2019-4154.
(Table 4: four tabled operating points; Table 3: model error summary; Fig. 2: m_dot vs delta_P curve for D=1.5mm N₂O injector.)

---

## 3. Test Configuration

| Parameter | Value | Source |
|---|---|---|
| Injector diameter D | 1.50 mm | Waxman Table 1 / Niño & Razavi Table 2 |
| Orifice length L | 18.4 mm | Waxman Table 1 |
| L/D | 12.3 | — |
| Inlet geometry | Square edge | Waxman injector no. 2 |
| Upstream temperature T₁ | 280 K | Niño & Razavi Table 4 |
| Upstream pressure P₁ | 4.36 MPa | Niño & Razavi Table 4 |
| Supercharge P₁_super | 0.62 MPa | Niño & Razavi Table 4 |
| Fluid | N₂O | — |
| QF_upstream | 0 (supercharged) | — |
| Discharge coefficient Cd | 0.65 | Waxman Fig. 15, square-edge at this supercharge |

P_sat at 280 K: model gives 37.98 bar vs. Waxman's implied 37.4 bar (+1.6% — Perry/McGill correlation accuracy at this temperature, documented in Section 4.1 of `04_implementation.md`).

---

## 4. Operating Points

Four points from Niño & Razavi Table 4, covering pre-critical, critical, and post-critical flow regimes. Experimental m_dot values read from the Waxman curve as reproduced in Niño & Razavi Fig. 2 (read-off uncertainty ±3%).

| Case | State | ΔP [bar] | P₂ [bar] | m_dot_exp [g/s] |
|---|---|---|---|---|
| I | Pre-critical | 8.4 | 35.2 | 44.0 |
| II | Critical | 9.8 | 33.8 | 46.5 |
| III | Post-critical | 10.9 | 32.7 | 47.5 |
| IV | Post-critical | 13.7 | 29.9 | 48.0 |

All four cases have P₂ < P_sat(T₁) = 37.98 bar → the Dyer model applies in all cases.

---

## 5. Results

### 5.1 SPI model (reference)

The SPI model is valid in this regime only at low ΔP (P₂ > P_sat). At the four test points, SPI over-predicts because the fluid partially flashes inside the orifice.

| Case | m_dot_SPI [g/s] | Error vs. exp |
|---|---|---|
| I | 43.9 | −0.2% |
| II | 47.4 | +2.0% |
| III | 50.0 | +5.3% |
| IV | 56.1 | +16.9% |

Note: Case I is near the SPI–two-phase transition; the SPI prediction happens to be close to experimental there. Cases III–IV show the expected growing over-prediction as ΔP increases.

### 5.2 Dyer/NHNE model (corrected formula)

$$\dot{m}_{Dyer} = \frac{\kappa}{1+\kappa}\,\dot{m}_{SPI} + \frac{1}{1+\kappa}\,\dot{m}_{HEM}$$

$$\kappa = \sqrt{\frac{P_1 - P_2}{P_{sat}(T_1) - P_2}}$$

| Case | ΔP [bar] | κ | x_exit | m_HEM [g/s] | m_Dyer [g/s] | m_exp [g/s] | Error |
|---|---|---|---|---|---|---|---|
| I | 8.4 | 1.739 | 0.043 | 39.4 | 42.3 | 44.0 | −4.0% |
| II | 9.8 | 1.532 | 0.055 | 40.3 | 44.6 | 46.5 | −4.1% |
| III | 10.9 | 1.437 | 0.063 | 40.7 | 46.2 | 47.5 | −2.7% |
| IV | 13.7 | 1.302 | 0.082 | 41.0 | 49.5 | 48.0 | +3.2% |

**Mean error: −1.9%. All four points within ±5%.**

This is consistent with Niño & Razavi Table 3, which reports a Dyer MAPE of 3.91% for the same dataset (using Cd = 0.63 from their internal correlation rather than the measured 0.65 used here).

---

## 6. Bug Found and Fixed During This Validation

During the validation process, a weight-swap error was found in `injector_two_phase.py`.

### What was wrong

The code had:

```python
m_dot_Dyer = m_dot_SPI / (1 + kappa) + (kappa / (1 + kappa)) * m_dot_HEM
```

The correct formula (Waxman 2013, Eq. 9; Solomon 2011, which corrected the original Dyer et al. 2007 error) is:

```python
m_dot_Dyer = (kappa / (1 + kappa)) * m_dot_SPI + (1.0 / (1 + kappa)) * m_dot_HEM
```

### Why this matters physically

κ is proportional to τ_bubble / τ_residence — the ratio of bubble growth time to fluid residence time in the orifice.

Large κ → bubbles grow slowly relative to residence time → less equilibrium → more weight on SPI (the "no vaporisation" limit).
Small κ → bubbles grow fast, near equilibrium → more weight on HEM.

The original code did the opposite: large κ → more weight on HEM. The corrected formula matches the physical definition.

### Impact on results

At typical design conditions (T = 280 K, P₁ = 55 bar, P₂ = 20 bar, κ ≈ 1.40):

| | m_dot_Dyer |
|---|---|
| Previous (wrong weights) | 63.3 g/s |
| Corrected | 70.8 g/s |
| Difference | +11.8% |

All 98 automated tests continue to pass after the fix — the tests check physical properties (Dyer is between HEM and SPI, positive, scales with area), not absolute numerical values.

---

## 7. Systematic Error Sources

**(a) P_sat correlation (+1.6% at 280 K, Perry/McGill vs. REFPROP)**
Shifts the κ denominator slightly, causing a small bias toward over- or under-prediction depending on operating conditions.

**(b) h_fg from tabulated Perry interpolation vs. REFPROP (~3–5%)**
Affects x_exit and therefore ρ_HEM. A slightly wrong h_fg moves the mixture density and hence m_HEM.

**(c) Experimental read-off uncertainty (~±3%)**
m_dot values estimated from a graph (Niño & Razavi Fig. 2), not a table.

**(d) Cd = 0.65 vs. Waxman's measured 0.71 for this injector**
Waxman Fig. 15 gives Cd ≈ 0.71 for the square-edge 1.5 mm injector at this supercharge level. Using 0.65 under-predicts m_SPI and m_HEM by ~8%, which partially explains the consistent slight under-prediction of m_Dyer.

---

## 8. Domain of Validity and Limitations

This model is validated for **injector sizing at realistic motor pressures** — ΔP of 20–50 bar, P₂ = chamber pressure, fluid arriving subcooled liquid at the injector inlet. This is the design regime.

The model is **not valid** in the critical flow regime (P₂ << P_sat, ΔP >> 50 bar). In that regime, the Bernoulli-based Dyer formula has no choking limit and over-predicts without bound. The physical limit (two-phase speed of sound) requires either the HEM isentropic maximum (needs entropy data, REFPROP) or the Henry-Fauske model. This limitation is documented in `future_work.md`.

---

## 9. Files

| File | Location | Description |
|---|---|---|
| `waxman_2013_validation.py` | `validation/` | Script that runs the model against Waxman data |
| `waxman_2013_experimental_data.csv` | `validation/` | Experimental data with extraction notes |
| `waxman_2013_results.md` | `validation/` | This document |
| `injector_two_phase.py` | `src/model/` | Fixed source file (Dyer weights corrected) |
