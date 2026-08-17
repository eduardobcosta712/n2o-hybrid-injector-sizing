# 4. Implementation

This document tracks the computational implementation of the model described in Sections 1–3, module by module, in the order they were built. Each section documents: the theory being translated into code, the implementation itself, its validation, and any issues encountered along the way — including issues that turned out to be informative rather than mere bugs.

## 4.1 `n2o_properties.py` — Saturated N₂O properties

### Purpose

Every other module in this project depends on this one: the subcooling margin (Section 1.4), the flashing criterion (Section 1.6), and the SPI/HEM/Dyer models (Sections 2–3) all require, at minimum, $P_{sat}(T)$ and its inverse $T_{sat}(P)$. This module implements those, plus saturated liquid density.

### Source of the correlations

Rather than deriving an equation of state from scratch (outside the scope of this project, and unnecessary — accurate, experimentally-fitted correlations already exist in the literature), this module uses published correlations for N₂O saturated properties, originally from Perry's Chemical Engineers' Handbook (Green & Perry, 2008), re-transcribed with full coefficients in Jean-Philyppe, J. (2023), *"A computational model for the design of a nitrous oxide-paraffin wax hybrid rocket engine,"* McGill Rocket Team technical report, arXiv:2302.06725, Appendix A.1. Full citation in `references.md`.

The saturation pressure correlation is:

$$P_{sat}(T) = \exp\left[c_1 + \frac{c_2}{T} + c_3 \ln T + c_4 T^{c_5}\right]$$

with $T$ in Kelvin, $P$ in Pa, valid for $T \in [182.33, 309.52]$ K (triple point to near the critical point), and coefficients $c_1 = 96.512$, $c_2 = -4045$, $c_3 = -12.277$, $c_4 = 2.886 \times 10^{-5}$, $c_5 = 2$.

The saturated liquid molar volume correlation is:

$$\nu_l(T) = \frac{c_2^{\,1 + (1 - T/c_3)^{c_4}}}{c_1}$$

with coefficients $c_1 = 2.781$, $c_2 = 0.27244$, $c_3 = 309.57$, $c_4 = 0.2882$, converted to density via $\rho_l = M_{N_2O} / \nu_l$, with $M_{N_2O} = 44.013$ kg/kmol.

Both correlations are implemented with an explicit valid-range check (`_check_range`): any call with $T$ outside $[182.33, 309.52]$ K raises an error rather than silently extrapolating, since the fit is not guaranteed valid there.

### Inverting $P_{sat}(T)$: Newton-Raphson with step damping

The design tool needs $T_{sat}(P)$ (saturation temperature given a pressure) as often as $P_{sat}(T)$ — this is the function used to compute the subcooling margin $\Delta T_{sub} = T_{sat}(P) - T$ at any point along the feed system. The correlation above cannot be algebraically inverted for $T$ in closed form ($T$ appears inside a fraction, a logarithm, and raised to a power simultaneously), so $T_{sat}(P)$ is solved numerically via Newton-Raphson, using the analytically-differentiated $dP_{sat}/dT$.

**Issue encountered and resolved.** An initial, undamped Newton-Raphson implementation (fixed default initial guess of 250 K) failed for target pressures corresponding to temperatures far from the guess: for example, inverting $P_{sat}(293.15\ \text{K})$ from a starting guess of 250 K produced a single Newton step of $-67$ K, overshooting directly past the correlation's valid upper bound (309.52 K) in one iteration.

This is a known failure mode of plain Newton-Raphson: $P_{sat}(T)$ becomes increasingly steep (non-linear) approaching the critical point (Section 1.5), so the local linear approximation used by each Newton step can be inaccurate over large distances from the current estimate, producing an overshooting correction.

**Fix.** Each Newton step is clamped to a maximum magnitude (`max_step`, default 20 K), and the resulting temperature estimate is additionally clamped to stay strictly inside $[T_{min}, T_{max}]$ after every iteration. This is a standard *damped* (or *safeguarded*) Newton's method: convergence takes a few more iterations when the initial guess is far from the root, but every intermediate step is guaranteed to remain within the domain where the correlation is valid.

This issue is noted here deliberately, rather than corrected silently, because it is informative: it confirms that care is needed specifically in the region nearest the critical point — which is precisely the region of greatest interest for this project (Section 3.2), since that is where N₂O's flashing behavior is most sensitive.

### Tabulated saturated vapor properties: $\nu_v(T)$, $h_l(T)$, $h_v(T)$, $h_{fg}(T)$

The HEM/Dyer two-phase injector model (Section 3.4) additionally requires the saturated vapor molar volume, $\nu_v(T)$, and the latent heat of vaporization, $h_{fg}(T)$. No closed-form correlation for $\nu_v(T)$ was available for this project. An ideal-gas approximation was considered and rejected: it degrades precisely near the critical point, the region of greatest interest for this project (Section 3.2), where intermolecular forces are no longer negligible. A Peng-Robinson equation of state was also considered, but a direct tabulated source was judged preferable — real measured data rather than a general model's prediction, avoiding the systematic error any cubic EOS carries near the critical point.

Table A.1 of the same source already in use for the closed-form correlations (Jean-Philyppe, J., 2023, arXiv:2302.06725, Appendix A) tabulates $\nu_v$, $h_l$, and $h_v$ at 27 saturation points spanning the full valid range ($T \in [182.33, 309.52]$ K). These values are stored in `n2o_saturation_table.csv` and used via linear interpolation between points — the same scheme the source paper itself uses for this quantity. The latent heat is then obtained directly as

$$h_{fg}(T) = h_v(T) - h_l(T)$$

rather than via the Clausius-Clapeyron relation, $h_{fg} = T(\nu_v - \nu_l)\,dP_{sat}/dT$: the direct subtraction avoids amplifying numerical error from $(\nu_v - \nu_l)$, which shrinks toward zero near the critical point.

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

### Purpose

Implements Module 1 from the implementation plan (Section 3.5 synthesis): tracks pressure and subcooling margin from the tank exit to the injector inlet, through an arbitrary sequence of straight pipe segments and fittings (valves, elbows), flagging the point (if any) where the fluid crosses the saturation curve.

### Theory implemented

- **Reynolds number**, $Re = \rho v D / \mu$, to classify the flow regime.
- **Darcy friction factor** $f$: exact laminar solution $f = 64/Re$ for $Re < 2300$; the explicit **Swamee-Jain approximation** to the (implicit) Colebrook equation for $Re \geq 2300$,
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

Vapor quality at the orifice exit is obtained assuming an isenthalpic process and full equilibrium at the exit (exit sits on the saturation curve at $T_{downstream} = T_{sat}(P_{downstream})$):

$$x = \frac{h_{upstream} - h_l(T_{downstream})}{h_{fg}(T_{downstream})}$$

The HEM mixture density follows from the mass-weighted average of the two phases' specific volumes, and $\dot m_{HEM}$ from the same orifice equation used throughout the project, with $\rho_{HEM}$ in place of the pure-liquid density. Dyer blends $\dot m_{SPI}$ and $\dot m_{HEM}$ via the non-equilibrium parameter $\kappa$:

$$\kappa = \sqrt{\frac{P_{upstream} - P_{downstream}}{P_{sat}(T_{upstream}) - P_{downstream}}}, \qquad \dot m_{Dyer} = \frac{\dot m_{SPI}}{1+\kappa} + \frac{\kappa}{1+\kappa}\dot m_{HEM}$$

**Domain restriction.** $\kappa$ requires $P_{upstream} > P_{sat}(T_{upstream})$ — the fluid must still be liquid at the orifice inlet, consistent with this model addressing vaporization *inside* the orifice (Section 3.1), not an already-two-phase feed line (that case is `feed_line.py`'s `flashing_detected`). Enforced with an explicit `ValueError` rather than an extreme or undefined $\kappa$.

### Validation

An operating point with the inlet modestly subcooled (55 bar at 20 °C, $P_{sat}(20°C) \approx 51.4$ bar) but a large enough pressure drop (to 20 bar) to cross saturation inside the orifice: Dyer predicts 128.8 g/s against SPI's 182.6 g/s (≈30% lower) — consistent with the expected direction of the two-phase correction (Section 3.3). The injector_spi.py example point (50 bar upstream) was found to already violate the domain restriction above and was not reused here.

### File location

`src/model/injector_two_phase.py`

---
*Next section: 4.5 `full_system.py` — orchestrating the full tank -> feed line -> injector path.*

## 4.5 `full_system.py` — Full system orchestration

### Purpose

Chains `feed_line.py`, `injector_spi.py`, and `injector_two_phase.py` in the order the sizing problem requires (Section 3.5 synthesis): evaluate the feed line to get conditions at the injector inlet, then automatically decide via `spi_sufficient` whether SPI alone is valid there or whether Dyer must be used, without manual intervention. Introduces no new physics.

### Validation

Reusing `feed_line.py`'s Case B geometry with a 55 bar tank: the feed line loses 2.63 bar (52.37 bar at the injector inlet, no flashing along the line), but `spi_sufficient` correctly returns `False` at the injector given the large pressure drop to the 20 bar chamber, and Dyer is selected automatically (124.8 g/s, vs. 175.6 g/s SPI would have predicted).

### File location

`src/model/full_system.py`
