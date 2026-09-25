# Validation Report — Waxman (2013/2014) Dataset (Updated, extended domain)
## Injector Model: SPI + Dyer/NHNE + Coupled Solver + Henry-Fauske Diagnostic

**Date:** September 2026 (regenerated with the real CoolProp package installed,
by running `validation/waxman_2013_validation.py`, extended with digitised
data from the Stanford AIAA 2013-3636 PDF)

---

## 1. Purpose

Validation of the complete injector model against the Waxman (2013/2014)
dataset. This report now uses **two** independent slices of the same paper:

- **Part A** — the original four Niño & Razavi (2019) tabulated points
  (injector 2 geometry, dP = 8–14 bar). Kept exactly as before.
- **Part B/C/D** — a much larger dataset digitised directly from the
  AIAA 2013-3636 PDF, using WebPlotDigitizer: the full injector-3 mass-flow map
  (Fig. 13, nine supercharge levels, dP up to ~46 bar), the matching
  effective-Cd curves (Fig. 14), the critical-flow-vs-supercharge curve
  (Fig. 16), and the SPI Cd of injectors 1/2/5 (Fig. 15). Digitised CSVs
  live in `validation/digitized/`.

This is the first time the model is checked against **more than four
points** and against pressure drops up to **46 bar** — inside the range
used by every worked example in `examples/`.

---

## 2. References

Same as before (Waxman 2013/2014, Niño & Razavi 2019, Henry & Fauske 1971,
Simoneau et al. 1971) plus the source PDF itself, from which Figs. 11–16
and Table 1 were digitised: Waxman, B.S., Zimmerman, J.E., Cantwell, B., &
Zilliac, G. (2013), AIAA 2013-3636.

---

## 3. Digitised dataset (Part B/C/D)

**Injector 3** (Table 1): straight hole, rounded inlet, D = 1.50 mm,
L = 18.4 mm, A = 1.7671 mm². Nine supercharge levels, legends give
P1 (psig), P1super (psi), and T1 to the nearest kelvin. Because a ±0.5 K
uncertainty on T1 moves P_sat by several psi — large next to a 41 psi
supercharge — **T1 is recovered from P1 and P1super** (`T1 = T_sat(P1 −
P1super)`) rather than taken as the stated integer; Part B5 shows this
barely changes anything (MAPE 6.85 % vs 7.03 %), so the choice is not
load-bearing.

**Calibration.** The discharge coefficient is fit **once**, pooled over
all nine curves, using only the single-phase window of each curve
(30 psi ≤ dP ≤ supercharge — pure SPI there, no two-phase model involved):
**Cd = 0.781** (n = 40 points, individual values 0.755–00.810). Waxman's
own text quotes ≈0.77 for this injector under one specific condition; the
match is consistent. Every two-phase (Dyer) point is then a genuine
prediction — none of them entered the fit.

**Consistency checks (Part E, all pass):**
- Fig. 11 (single test, P1super = 169 psi) matches the corresponding Fig.
  13 curve exactly once its y-axis is rescaled by 0.1 (a digitiser
  calibration factor) — median ratio 1.002 after correction. The script
  raises `RuntimeError` if this check ever fails.
- The Cd implied by Fig. 13 (mass flow ÷ SPI formula) matches the Cd read
  directly off Fig. 14 to within 1 % median, 5th–95th percentile
  0.991–1.020.
- The critical mass flow of Fig. 16 matches the plateau (maximum) of the
  matching Fig. 13 curve to within 1.5 %.

---

## 4. Results — Part A (unchanged)

| Case | ΔP [bar] | m_exp [g/s] | m_dot [g/s] | Error | Regime |
|---|---|---|---|---|---|
| Pre-critical | 8.40 | 44.0 | 42.91 | −2.5% | Dyer |
| Critical | 9.80 | 46.5 | 45.34 | −2.5% | Dyer |
| Post-critical 1 | 10.90 | 47.5 | 46.97 | −1.1% | Dyer |
| Post-critical 2 | 13.70 | 48.0 | 50.37 | +4.9% | Dyer |

**MAPE = 2.76 %**, mean −0.3 %, all 4 within ±5 %. Unchanged from before.
Sensitivity check: re-running these same 4 points with Cd = 0.681 (the
digitised Fig. 15 mean for injector **2**, the actual geometry of these
points, vs the 0.65 assumed by the original report) gives MAPE = 4.41 %
— the Cd choice matters more than anything else at this scale.

---

## 5. Results — Part B: full injector-3 map (Fig. 13), dP up to 46 bar

**B1 — single-phase window (calibration residual, not a prediction):**
MAPE = 1.25 %, all 40/40 points within ±5 %. Confirms SPI + the fitted Cd
reproduce the single-phase branch essentially exactly.

**B2 — two-phase region (Dyer, dP > supercharge), n = 64, genuine
predictions:**

| | mean | MAPE | max\|err\| | within ±5% | within ±10% |
|---|---|---|---|---|---|
| Dyer (as used by the tool) | +6.70% | 6.85% | 33.4% | 41/64 | 48/64 |
| min(Dyer, Henry-Fauske) — diagnostic only | +4.68% | 4.87% | 18.8% | 41/64 | 50/64 |

Dyer **systematically over-predicts** outside the previously-validated
8–14 bar band — always in the same direction, never under. 31 of the 64
points exceed the Henry-Fauske ceiling (`choked = True`); at **30 of
those 31**, capping at the ceiling would have moved the prediction closer
to experiment, and at **28 of the 31** the true experimental value sits
*below* the Henry-Fauske ceiling too (i.e. the ceiling would not even
have had to be perfectly tight to help). This is new evidence — Part 1
of `future_work.md` previously had *zero* points where the ceiling bound
anything; there are now 31, and it looks like a useful (if imperfect)
correction rather than a coincidence.

**B3 — by injector pressure drop:**

| dP band | Dyer MAPE | capped MAPE |
|---|---|---|
| ≤ 14 bar (previously validated) | 5.27% | 5.22% |
| 14–30 bar | 9.18% | 6.14% |
| 30–46 bar | **2.52%** | 1.02% |

Counter-intuitively, the **largest** pressure drops (30–46 bar — the
range every worked example in `examples/` actually uses) give the
**best** agreement, not the worst. The 14–30 bar band is where Dyer
over-predicts most.

**B4 — by supercharge (tank subcooling margin):**

| supercharge | Dyer MAPE | capped MAPE |
|---|---|---|
| 41 psi (2.8 bar) | 17.7% | 12.3% |
| 79 psi (5.5 bar) | 13.5% | 9.8% |
| 115 psi (7.9 bar) | 6.6% | 5.6% |
| 169 psi (11.7 bar) | 5.9% | 4.1% |
| 206–371 psi (14–26 bar) | 1.1–3.5% | 0.6–3.1% |

Aggregated: **supercharge ≥ 200 psi (≥1.38 MPa) → MAPE 1.96 %**;
**supercharge < 200 psi → MAPE 11.74 %**. The controlling variable is not
the pressure drop itself but how subcooled the tank is: Dyer degrades
sharply at low supercharge, regardless of dP.

**B5 — sensitivity check on the two-phase result:**

| variant | mean | MAPE |
|---|---|---|
| baseline (Cd = 0.781 pooled, T1 recovered) | +6.70% | 6.85% |
| Cd = 0.77 (Waxman's own quoted text value) | +5.22% | 6.12% |
| Cd fitted curve-by-curve instead of pooled | +5.96% | 6.30% |
| T1 = stated integer kelvin | +6.93% | 7.03% |

None of these choices explains away the over-prediction at low
supercharge — it is a property of the Dyer model itself in this regime,
not a calibration artefact.

---

## 6. Results — Part C: critical mass flow (Fig. 16)

Using Waxman's own 5 %-deviation criterion (first dP where the effective
Cd drops below 95 % of the SPI value) applied to the model's own Dyer
curve:

| | mean | MAPE | max\|err\| |
|---|---|---|---|
| Model (Dyer, 95% criterion) | −1.0% | **2.6%** | 11.9% |
| Henry-Fauske ceiling | +6.2% | 6.2% | 15.3% |
| HEM (isenthalpic) ceiling | −6.1% | 6.1% | 7.9% |

The model's own critical-flow definition reproduces Fig. 16 well (2.6 %
MAPE, n = 8), confirming the earlier finding that Dyer legitimately sits
above the equilibrium HEM ceiling and below (mostly) the Henry-Fauske
ceiling — both ceilings by construction over/under-shoot the true
critical point, in opposite directions, exactly as the theory in
`docs/03_two_phase_flow.md` predicts.

---

## 7. Interpretation and consequences for the tool

1. **The previously-validated 8–14 bar band remains fine** (Part A,
   unchanged) and is now also confirmed by the larger dataset at
   comparable supercharge (Part B, 169–206 psi curves: MAPE 4–6 %).
2. **The 20–50 bar band used by every worked example is now, for the
   first time, backed by data** — and the news is mixed: at *high*
   supercharge (≥ 1.38 MPa ≈ 14 bar) it is good (MAPE ≈ 2 %); at *low*
   supercharge it is not (MAPE up to 18 %), and Dyer always over-predicts
   there, never under.
3. **The Henry-Fauske ceiling is no longer merely a theoretical,
   never-triggered diagnostic.** In this larger dataset it fires on
   31/64 two-phase points and, when it fires, applying it as a cap
   improves the prediction in 30/64 cases. This is still not strong
   enough evidence to switch it from diagnostic to automatic cap
   (it does not fully close the gap, and 28/31 experimental points sit
   below the ceiling too, meaning a *tighter* correction than Henry-Fauske
   would do even better) — but the negative result recorded before
   ("no data point where it binds") is now out of date and should be
   corrected in `docs/future_work.md`, Priority 1 and 2.
4. **Practical implication for own designs:** trust Dyer most
   when the tank has a solid subcooling margin (≥ 10–14 bar) even at
   large ΔP; be more cautious the closer the tank sits to saturation,
   *regardless* of how large the pressure drop across the injector is.

## 8. What this does **not** establish

- Only injector 3 (rounded inlet, 1.5 mm) was validated at more than one
  supercharge; Part A's injector-2 geometry (square edge) still has only
  4 points.
- Fig. 13/14/16 come from a single facility/rig; no independent
  cross-check dataset exists.
- The two-phase feed-line model and HEM two-phase-inlet path remain
  entirely unvalidated (unchanged from before).
- The digitisation itself carries plot-reading uncertainty (see Part E
  consistency checks, all within ~2 %, which bounds this).

---

## 9. Files

| File | Location |
|---|---|
| `waxman_2013_results.md` | `validation/` (this file) |
| `waxman_2013_validation.py` | `validation/` — Parts A–E, `--no-plot` to skip the figure |
| `waxman_2013_fig13_comparison.png` | `validation/` — 9-panel model-vs-experiment figure |
| `waxman_2013_experimental_data.csv` | `validation/` (unchanged; still not read by code) |
| `digitized/waxman_fig11_mdot_vs_dP_single_test.csv` | `validation/digitized/` |
| `digitized/waxman_fig12_cd_vs_dP_single_test.csv` | `validation/digitized/` |
| `digitized/waxman_fig13_mdot_vs_dP_by_supercharge.csv` | `validation/digitized/` |
| `digitized/waxman_fig14_cd_vs_dP_by_supercharge.csv` | `validation/digitized/` |
| `digitized/waxman_fig15_cd_vs_supercharge_injectors_1_2_5.csv` | `validation/digitized/` |
| `digitized/waxman_fig16_critical_mdot_vs_supercharge.csv` | `validation/digitized/` |
