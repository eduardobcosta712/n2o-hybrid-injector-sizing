# 4. Implementation

This document tracks the computational implementation of the model described in Sections 1–3, module by module, in the order they were built. Each section documents: the theory being translated into code, the implementation itself, its validation, and any issues encountered along the way — including issues that turned out to be informative rather than mere bugs.

## 4.1 `n2o_properties.py` — Saturated N₂O properties

### Purpose

Every other module in this project depends on this one: the subcooling margin (Section 1.4), the flashing criterion (Section 1.6), and the SPI/HEM/Dyer models (Sections 2–3) all require, at minimum, $P_{sat}(T)$ and its inverse $T_{sat}(P)$. This module implements those, plus saturated liquid density.

### Source of the correlations

Rather than deriving an equation of state from scratch (outside the scope of this project, and unnecessary — accurate, experimentally-fitted correlations already exist in the literature), this module uses published correlations for N₂O saturated properties, originally from Perry's Chemical Engineers' Handbook (Green & Perry, 2008), re-transcribed with full coefficients in Jean-Philyppe, J. (2023), *"A computational model for the design of a nitrous oxide-paraffin wax hybrid rocket engine,"* McGill Rocket Team technical report, arXiv:2302.06725, Appendix A.1. Full citation in `references.md`.

The saturation pressure correlation is:

$$P_{sat}(T) = \exp\left[c_1 + \frac{c_2}{T} + c_3 \ln T + c_4 T^{c_5}\right]$$

with T in Kelvin, P in Pa, valid for T ∈ [182.33, 309.52] K (triple point to near the critical point), and coefficients $c_1 = 96.512$, $c_2 = -4045$, $c_3 = -12.277$, $c_4 = 2.886 \times 10^{-5}$, $c_5 = 2$.

The saturated liquid molar volume correlation is:

$$\nu_l(T) = \frac{c_2^{\,1 + (1 - T/c_3)^{c_4}}}{c_1}$$

with coefficients $c_1 = 2.781$, $c_2 = 0.27244$, $c_3 = 309.57$, $c_4 = 0.2882$, converted to density via $\rho_l = M_{N_2O} / \nu_l$, with $M_{N_2O} = 44.013$ kg/kmol.

Both correlations are implemented with an explicit valid-range check (`_check_range`): any call with T outside [182.33, 309.52] K raises an error rather than silently extrapolating, since the fit is not guaranteed valid there.

### Inverting $P_{sat}(T)$: Newton-Raphson with step damping

The design tool needs $T_{sat}(P)$ (saturation temperature given a pressure) as often as $P_{sat}(T)$ — this is the function used to compute the subcooling margin $\Delta T_{sub} = T_{sat}(P) - T$ at any point along the feed system. The correlation above cannot be algebraically inverted for $T$ in closed form ($T$ appears inside a fraction, a logarithm, and raised to a power simultaneously), so $T_{sat}(P)$ is solved numerically via Newton-Raphson, using the analytically-differentiated $dP_{sat}/dT$.

**Issue encountered and resolved.** An initial, undamped Newton-Raphson implementation (fixed default initial guess of 250 K) failed for target pressures corresponding to temperatures far from the guess: for example, inverting $P_{sat}(293.15\ \text{K})$ from a starting guess of 250 K produced a single Newton step of $-67$ K, overshooting directly past the correlation's valid upper bound (309.52 K) in one iteration.

This is a known failure mode of plain Newton-Raphson: $P_{sat}(T)$ becomes increasingly steep (non-linear) approaching the critical point (Section 1.5), so the local linear approximation used by each Newton step can be inaccurate over large distances from the current estimate, producing an overshooting correction.

**Fix.** Each Newton step is clamped to a maximum magnitude (`max_step`, default 20 K), and the resulting temperature estimate is additionally clamped to stay strictly inside $[T_{min}, T_{max}]$ after every iteration. This is a standard *damped* (or *safeguarded*) Newton's method: convergence takes a few more iterations when the initial guess is far from the root, but every intermediate step is guaranteed to remain within the domain where the correlation is valid.

This issue is noted here deliberately, rather than corrected silently, because it is informative: it confirms that care is needed specifically in the region nearest the critical point — which is precisely the region of greatest interest for this project (Section 3.2), since that is where N₂O's flashing behavior is most sensitive.

### Tabulated saturated properties: $\nu_v(T)$, $h_l(T)$, $h_v(T)$, $h_{fg}(T)$, $\mu_v(T)$, $\mu_l(T)$, $c_{pl}(T)$

The HEM/Dyer two-phase injector model (Section 3.4) and the two-phase feed-line model (Section 4.2) require several saturated properties beyond the closed-form correlations. These are stored in `n2o_saturation_table.csv` across four sections:

**Table A.1** (McGill/Perry, arXiv:2302.06725): $\nu_v$, $h_l$, $h_v$ at 27 saturation points (182–309 K). $h_{fg}(T) = h_v(T) - h_l(T)$ directly — not via Clausius-Clapeyron, which amplifies numerical error near the critical point.

**Table A.2** (McGill/Perry): derivatives of Table A.1 (not currently used).

**Table A.3** (NIST WebBook, Millat et al. 1991 viscosity correlation, downloaded August 2026): saturated vapour dynamic viscosity $\mu_v(T)$ in 26 points (182–307 K). Uncertainty ~2% for T > 150 K. Used in `mu_vapor_sat(T)` and `mu_mixture(x, T)` for the two-phase feed-line model.

**Table A.4** (NIST WebBook, Lemmon & Span 2006 EOS + Laesecke & Hafer 1998 liquid viscosity, downloaded August 2026): isobaric heat capacity of saturated liquid $c_{pl}(T)$, saturated liquid dynamic viscosity $\mu_l(T)$, and liquid/vapour entropy $s_l(T)$, $s_v(T)$. Used in `cp_liquid_sat(T)` and `mu_liquid_sat(T)`. Entropy data is available for future use in the isentropic choking limit (see `future_work.md`, Priority 1).

Note: $\mu_l(T)$ replaces the previous constant `MU_LIQUID_N2O = 1.5e-4 Pa·s` in the feed-line model. The constant is retained as a fallback where T is not known.

### Validation

Run via `python n2o_properties.py`, which checks:

1. **$P_{sat}(T)$ against reference values** quoted in `01_n2o_thermodynamics.md`: agreement within $+2.6\%$ at 0 °C, $+0.9\%$ at 20 °C, and $-4.8\%$ at 34 °C. The larger discrepancy near the critical point is consistent with the physical sensitivity discussed in Section 3.2 — small differences between independently-sourced reference values are expected to be amplified there, rather than indicating an implementation error.
2. **$T_{sat}(P)$ as the exact numerical inverse of $P_{sat}(T)$**: recovered temperatures agree with the original inputs to within $10^{-12}$ K across the tested range, confirming the damped Newton-Raphson solver is implemented correctly.
3. **Saturated liquid density at 20 °C**: 784.8 kg/m³ computed, against a typical literature reference value of ≈786 kg/m³ (agreement within 0.15%).

### File location

`src/model/n2o_properties.py`, reading `src/model/n2o_saturation_table.csv`

---
*Next section: 4.2 `feed_line.py` — pressure drop and subcooling margin along the feed line (Darcy-Weisbach friction losses, fitting losses).*

## 4.2 `feed_line.py` — Feed line pressure drop and subcooling margin

**Two-phase HEM model (implemented August 2026).** The feed line model now handles two distinct flow regimes:

**Single-phase region** (`P > P_sat(T_tank)`): Darcy-Weisbach with pure liquid properties — `rho_l(T_tank)` and `mu_l(T_tank)` from `mu_liquid_sat(T)` (Table A.4, NIST). Previously a constant (`MU_LIQUID_N2O = 1.5e-4 Pa·s`); now temperature-dependent.

**Two-phase region** ($P \leq P_{\mathrm{sat}}(T_{\mathrm{tank}})$): once flashing is detected, all subsequent segments use HEM mixture properties updated at each segment's local pressure:

$$x(s) = \frac{h_l(T_{tank}) - h_l(T_{sat}(P(s)))}{h_{fg}(T_{sat}(P(s)))}, \qquad \rho_{mix} = \frac{1}{\dfrac{1-x}{\rho_l} + \dfrac{x}{\rho_v}}, \qquad \mu_{mix} = (1-x)\,\mu_l + x\,\mu_v$$

where $\mu_v$ comes from `mu_vapor_sat(T)` (NIST WebBook / Millat et al. 1991, interpolated from `n2o_saturation_table.csv` Table A.3). The per-segment trace includes `x_quality`, `rho_eff_kg_m3`, and `mu_eff_Pa_s` for inspection and plotting.

**Physical significance.** With typical conditions (N₂O at 20°C, $x \approx 0.05$), $\rho_{mix} \approx 60\text{–}100\,\text{kg/m}^3$ versus $\rho_l \approx 785\,\text{kg/m}^3$ — a factor of 8–13× reduction. Since $\Delta P \propto \dot{m}^2 / (\rho_{mix} A^2)$, the two-phase friction losses in this region are correspondingly 8–13× larger than the single-phase model would predict. This is the dominant effect; the viscosity correction ($\mu_{mix}$ vs $\mu_l$) is secondary.

**`x_inlet` output.** The final `x_inlet` is computed via isenthalpic flash at $P_{final}$, which feeds directly into `hem_mass_flow_two_phase_inlet` in `full_system.py` when `flashing_detected = True`.

### Purpose

Implements Module 1 from the implementation plan (Section 3.5 synthesis): tracks pressure and subcooling margin from the tank exit to the injector inlet, through an arbitrary sequence of straight pipe segments and fittings (valves, elbows), flagging the point (if any) where the fluid crosses the saturation curve.

### Theory implemented

- **Reynolds number**, $Re = \rho v D / \mu$, to classify the flow regime.
- **Darcy friction factor** f: exact laminar solution `f = 64/Re` for Re < 2300; the explicit **Swamee-Jain approximation** to the (implicit) Colebrook equation for Re ≥ 2300,
$$f = \frac{0.25}{\left[\log_{10}\left(\dfrac{\varepsilon/D}{3.7} + \dfrac{5.74}{Re^{0.9}}\right)\right]^2}$$
chosen over solving Colebrook directly because it is explicit (no iteration needed) while remaining within ~1% of the implicit solution for the range of interest. The 2300–4000 transitional regime is conservatively treated with the turbulent formula, documented as a deliberate simplification in the function's docstring.
- **Darcy-Weisbach friction loss**, $\Delta P = f (L/D)(\rho v^2/2)$, and **fitting (minor) losses**, $\Delta P = K(\rho v^2/2)$.
- The line is assumed **adiabatic** (constant temperature, per the Section 1.6/2.1 justification): only pressure is tracked segment by segment; temperature stays fixed at the tank value throughout.

### Validation

Two cases were run with an identical line geometry (2 m of 8 mm ID tubing, one ball valve, one 90° elbow), differing only in initial tank subcooling:

- **Case A — tank exactly at saturation** ($\Delta T_{sub} = 0$ initially): the very first segment already pushes the margin negative (−1.07 K), and the flag `flashing_detected` correctly triggers. This is expected, not a bug: with zero initial margin, *any* pressure drop, however small, crosses the saturation curve — the code is correctly enforcing the definition from Section 1.4 in the least forgiving case.
- **Case B — tank with 5 bar of initial subcooling**: the identical line geometry produces essentially the same pressure drop (≈2.63 bar) but the final margin stays positive (+2.05 K), and no flashing is flagged. Same hardware, different outcome — driven entirely by how much margin the tank started with.

**Known limitation surfaced by this comparison.** The pressure drop is nearly identical between the two cases because `rho_liquid_sat(T_tank)` depends only on temperature, not pressure, in the current implementation — consistent with the incompressibility assumption underlying the SPI model (Section 2.4), but a simplification worth stating explicitly: for feed lines with much larger pressure excursions than in this example, real liquid density would vary somewhat with pressure too, an effect this module does not capture.

### File location

`src/model/feed_line.py`

---
*Next section: 4.3 `injector_spi.py` — the SPI injector model, and the criterion for when it is sufficient on its own.*

## 4.3 `injector_spi.py` — SPI injector model and sufficiency criterion

### Purpose

Implements the SPI model derived in Section 2 — both directions (mass flow from a known orifice area, and required orifice area from a target mass flow) — plus the decision criterion (planned in the pre-implementation discussion, "Module 2" logic) that determines whether SPI alone is valid at a given operating point, or whether two-phase effects must be considered.

### Theory implemented

$$\dot{m}_{SPI} = C_d \, A \sqrt{2 \rho \Delta P}$$

and its algebraic inverse, $A = \dot m_{target} / (C_d \sqrt{2\rho\Delta P})$ — the direction more commonly needed in practice, since the target mass flow is usually fixed by the motor's thermochemical sizing (O/F ratio, chamber pressure) and the orifice area is the unknown being solved for.

**Sufficiency criterion.** Per Section 3.1 (pressure is minimum, and therefore flashing risk is greatest, at the orifice's narrowest cross-section), SPI alone is valid only if even the downstream pressure stays at or above the saturation pressure evaluated at the upstream temperature:

$$\text{SPI sufficient} \iff P_{downstream} \geq P_{sat}(T_{upstream})$$

If this fails, the flow crosses the saturation curve somewhere inside the orifice, and the two-phase model (Section 4.4, Dyer) must be used instead.

### Validation

A deliberately demanding test case was used: N₂O at 20 °C, $P_{upstream} = 50$ bar (already slightly below $P_{sat}(20°C) \approx 51.4$ bar — i.e. entering the orifice already at the edge of saturation), $P_{downstream} = 20$ bar. The `spi_sufficient` criterion correctly returns `False`, and `spi_mass_flow` reports 169.0 g/s — a number now understood to be *not* the physically correct flow rate for this operating point, but rather the reference value SPI would (incorrectly) predict by assuming single-phase liquid throughout. This value is retained as the comparison baseline for validating the Dyer model in Section 4.4, where it is expected to predict a lower, more physically realistic mass flow (per the two-phase choking discussion in Section 3.3).

### File location

`src/model/injector_spi.py`

---
*Next section: 4.4 `injector_two_phase.py` — HEM and Dyer models.*

## 4.4 `injector_two_phase.py` — HEM and Dyer two-phase models

### Purpose

Implements the two-phase injector models from Section 3.4: HEM (full thermodynamic equilibrium) and Dyer (weighted blend of SPI and HEM), used when `spi_sufficient` returns `False`.

### Theory implemented

Vapor quality at the orifice exit is obtained assuming an isenthalpic process and full equilibrium at the exit (exit sits on the saturation curve at T(downstream) = T(sat)(P(downstream))):

$$x = \frac{h_{upstream} - h_l(T_{downstream})}{h_{fg}(T_{downstream})}$$

The HEM mixture density follows from the mass-weighted average of the two phases' specific volumes, and $\dot m_{HEM}$ from the same orifice equation used throughout the project, with $\rho_{HEM}$ in place of the pure-liquid density. Dyer blends $\dot m_{SPI}$ and $\dot m_{HEM}$ via the non-equilibrium parameter $\kappa$:

$$\kappa = \sqrt{\frac{P_{upstream} - P_{downstream}}{P_{sat}(T_{upstream}) - P_{downstream}}}, \qquad \dot m_{Dyer} = \frac{\dot m_{SPI}}{1+\kappa} + \frac{\kappa}{1+\kappa}\dot m_{HEM}$$

**Domain restriction.** $\kappa$ requires $P_{upstream} > P_{sat}(T_{upstream})$ — the fluid must still be liquid at the orifice inlet, consistent with this model addressing vaporization *inside* the orifice (Section 3.1), not an already-two-phase feed line (that case is `feed_line.py`'s `flashing_detected`). Enforced with an explicit `ValueError` rather than an extreme or undefined $\kappa$.

### Validation

An operating point with the inlet modestly subcooled (55 bar at 20 °C, $P_{sat}(20°C) \approx 51.4$ bar) but a large enough pressure drop (to 20 bar) to cross saturation inside the orifice: Dyer predicts 128.8 g/s against SPI's 182.6 g/s (≈30% lower) — consistent with the expected direction of the two-phase correction (Section 3.3). The injector_spi.py example point (50 bar upstream) was found to already violate the domain restriction above and was not reused here.


### HEM critical flow — `hem_critical_flow()`

Available as a standalone function in `injector_two_phase.py`. Implements Waxman (2013) Eq.(5): scans $P_2$ from $P_{sat}(T_{upstream})$ downward and finds the maximum of the isenthalpic HEM mass-flow curve:

$$\dot{m}_{crit} = \max_{P_2 < P_{sat}} \left[ C_d A \sqrt{2\,\rho_{mix}(P_2)\,\Delta P} \right]$$

This maximum is the physical choking limit — the two-phase speed-of-sound condition expressed through the isenthalpic path. At Waxman conditions (`T1 = 280 K`, `P1 = 4.36 MPa`): `m_dot_crit = 41.1 g/s` at `P2_crit = 30.4 bar`.

The function is **not applied automatically** in `dyer_mass_flow()` because the Dyer non-equilibrium correction legitimately predicts above the HEM-only ceiling (confirmed by Waxman experimental data: 44–48 g/s vs. HEM cap of 41.1 g/s).

### Isentropic choking scan — `hem_critical_flow_isentropic()` (added September 2026, Priority 1)

`hem_critical_flow()` above uses the isenthalpic path ($h=$ const), which correctly describes the real thermodynamic *state* of the fluid at the orifice exit (an orifice is adiabatic, so the 1st law gives $h_{up}=h_{down}$ regardless of internal irreversibility), but is only an approximation to the true choking condition. Choking is set by the two-phase speed of sound, $c^2 = (\partial P/\partial\rho)_s$ — a derivative taken at constant **entropy**, because an acoustic disturbance is a small, fast, essentially reversible perturbation on top of the (possibly irreversible) mean flow.

`hem_critical_flow_isentropic()` implements this more rigorous scan, mirroring `hem_critical_flow()` exactly except along $s=$const:

$$x_{is}(P_2) = \frac{s_{up} - s_l(T_{sat}(P_2))}{s_v(T_{sat}(P_2)) - s_l(T_{sat}(P_2))}, \qquad \dot m_{crit,\,is} = \max_{P_2 < P_{sat}} \Big[C_d A \sqrt{2\,\rho_{HEM}(x_{is})\,(P_{up}-P_2)}\Big]$$

The two required entropy functions, `s_liquid_sat(T)` and `s_vapor_sat(T)`, were added to `n2o_properties.py`, interpolating $s_l$, $s_v$ from Table A.4 (NIST WebBook) — the same table already used for `cp_liquid_sat`/`mu_liquid_sat`, previously loaded but with no public accessor. `s_fg(T) = s_v(T) - s_l(T)` is the entropy-domain analogue of `h_fg(T)`.

**Table A.4 range gap.** While implementing this, a pre-existing gap surfaced: Table A.4 only covers 182.33–307.33 K, narrower than the module's main correlation range (up to 309.52 K, the critical point). `cp_liquid_sat`/`mu_liquid_sat` were previously checked against the wider range, so a call between 307.33 K and 309.52 K would silently fall through to `_interp`'s generic "should be unreachable" `RuntimeError`. This is now an explicit, named range check (`T_MIN_A4`, `T_MAX_A4`, `_check_range_a4`), applied consistently to all four Table-A.4 functions plus the two new entropy functions. `hem_critical_flow_isentropic()` checks `T_upstream` against this range up front and raises a `ValueError` pointing to `hem_critical_flow()` (isenthalpic) as a fallback, and to Priority 4 (CoolProp/REFPROP) as the eventual fix.

**Kept side by side, not replaced.** `hem_critical_flow()` (isenthalpic) is unchanged and remains the version cited in `validation/waxman_2013_results.md`. `hem_critical_flow_isentropic()` is additive, for direct comparison and eventual use as the Dyer cap.

**Validation.** At Waxman conditions ($T_1=280$ K, $P_1=4.36$ MPa, $D=1.5$ mm, $C_d=0.65$):

| Path | $\dot m_{crit}$ | $P_{2,crit}$ | $x_{crit}$ |
|---|---|---|---|
| Isenthalpic (`hem_critical_flow`) | 41.05 g/s | 30.37 bar | 0.0856 |
| Isentropic (`hem_critical_flow_isentropic`) | 41.64 g/s | 29.65 bar | 0.0880 |

+1.43% difference — both remain below the experimental Dyer-regime range (44.0–48.0 g/s), consistent with the existing interpretation that Dyer's non-equilibrium correction legitimately predicts above either HEM-only ceiling.

**What remains** (see `future_work.md`, Priority 1): deciding how `hem_critical_flow_isentropic()` should be applied automatically as a cap inside `dyer_mass_flow()`, and re-confirming the Waxman MAPE afterwards.

### File location

`src/model/injector_two_phase.py`

---
*Next section: 4.5 `full_system.py` — orchestrating the full tank -> feed line -> injector path.*

## 4.5 `full_system.py` — Full system orchestration and coupled solver

### Purpose

Chains `feed_line.py`, `injector_spi.py`, and `injector_two_phase.py` in the order the sizing problem requires (Section 3.5 synthesis): evaluate the feed line to get conditions at the injector inlet, then automatically decide via `spi_sufficient` whether SPI alone is valid there or whether Dyer must be used, without manual intervention. Introduces no new physics.

### Validation

Reusing `feed_line.py`'s Case B geometry with a 55 bar tank: the feed line loses 2.63 bar (52.37 bar at the injector inlet, no flashing along the line), but `spi_sufficient` correctly returns `False` at the injector given the large pressure drop to the 20 bar chamber, and Dyer is selected automatically (124.8 g/s, vs. 175.6 g/s SPI would have predicted).

### File location

`src/model/full_system.py`

---
*Next section: 4.6 `src/interface/` — the interactive Streamlit tool.*

## 4.6 `src/interface/` — Interactive tool

### Purpose

A Streamlit web app (`app.py`) wrapping `full_system.py`: editable tank, feed line (dynamic pipe/fitting segment list), and injector inputs, updating live as inputs change. Reports flashing/SPI-sufficiency status and the real mass flow, and renders two diagrams (`plotting.py`): pressure along the feed line against $P_{sat}(T_{tank})$, and a P-T diagram with the saturation curve and the tank/injector-inlet/chamber operating points.

### Implementation notes

Segment state is kept in `st.session_state`, since Streamlit re-runs the whole script on every interaction; without it, the segment list would reset on every slider move. All UI inputs are in display-friendly units (bar, °C, mm, g/s) and converted to SI at the UI boundary before calling into `full_system.py`, which continues to operate in SI throughout, per the project's units convention.

### Running

From the repository root: `streamlit run src/interface/app.py`. Requires `pip install streamlit matplotlib numpy`.

### File location

`src/interface/app.py`, `src/interface/plotting.py`
