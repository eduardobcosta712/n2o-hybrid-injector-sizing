# 3. Two-Phase Flow: Where the SPI Model Fails

## 3.1 The vena contracta as the system's critical point

From Section 2.1: velocity is maximum and pressure is minimum at the orifice's narrowest cross-section (*vena contracta*). From Section 1.6: it is the pressure drop that pushes the fluid across the saturation curve. Combining the two: **the point in the system most likely to trigger flashing is inside the injector orifice itself**, at its minimum cross-section — not necessarily the feed line upstream, where velocity (and the associated pressure drop) is comparatively small.

Direct consequence: even if N₂O reaches the injector inlet perfectly liquid and subcooled, it can still begin vaporizing *inside* the orifice, purely from the acceleration occurring there. This is a distinct (though related) flashing mechanism from the one associated with losses along the feed line.

## 3.2 Amplified sensitivity near the critical point

From Section 1.5: T_crit(N₂O) ≈ 36.4 °C. The practical relevance of this proximity to ambient conditions is tied to the behavior of the **latent heat of vaporization**, $h_{fg}$ — the energy required to convert, at the same pressure/temperature, a unit mass of saturated liquid into saturated vapor.

$h_{fg}$ decreases as temperature approaches $T_{crit}$, reaching zero at the critical point itself (where liquid and vapor become indistinguishable). At typical ambient temperatures (15–25 °C), already relatively close to $T_{crit}$, N₂O's $h_{fg}$ is significantly reduced compared to fluids operating well away from their critical point (contrast with water, $T_{crit} \approx 374\ ^\circ\text{C}$). A given pressure drop generates, for N₂O, a disproportionately larger vapor quality x than it would for other, comparatively "well-behaved" fluids — this is the physical reason this effect cannot be treated as negligible for this specific propellant.

## 3.3 Two-phase choking and the drop in mass flow

When a fraction of the liquid vaporizes inside the orifice, the mixture's average density, ρ_mix, drops sharply — the specific volume of saturated N₂O vapor is substantially larger than that of the saturated liquid (potentially, depending on pressure, on the order of tens to hundreds of times larger), so even a small mass-based vapor quality $x$ already occupies most of the available volume in the orifice cross-section.

Two consequences, which reinforce each other:

1. **Direct reduction in $\dot m$** — for the same area A, a mixture with much lower density carries less mass per unit time, even at similar velocity, since ṁ depends directly on ρ.
2. **Two-phase choking** — the liquid-vapor mixture is significantly compressible, unlike the pure liquid; there is a maximum speed at which pressure disturbances can propagate (analogous to a local speed of sound in the mixture). Once that limiting speed is reached in the orifice, further downstream pressure reductions stop increasing the flow rate — the pressure disturbance can no longer propagate upstream past the choking point. Flow rate becomes capped at a maximum value, regardless of further downstream pressure reductions.

The combined effect of these two mechanisms explains why real mass flow, under significant vaporization, can fall far below what the SPI model predicts — a model that always assumes dense, incompressible liquid with no propagation speed limit. Discrepancies of several-fold between predicted and measured flow, documented in the literature (`references.md`), reflect a physical regime change, not a fine calibration error.

## 3.4 Two-phase flow models: HEM and Dyer

### HEM — Homogeneous Equilibrium Model

Assumes that, at every point in the flow, the liquid and vapor phases are in mutual thermal and mechanical equilibrium — same temperature, same pressure, treated as a single effective fluid whose properties (density, among others) depend on the local vapor quality x. It is relatively simple to implement computationally, but assumes instantaneous vaporization as soon as local pressure crosses $P_{sat}(T)$.

This instantaneous-equilibrium assumption is a known limitation in short orifices (like those of a rocket injector): vaporization is not instantaneous — it requires nucleation sites for bubble formation — and the fluid's residence time in a short orifice can be insufficient for full thermodynamic equilibrium to establish.

### Dyer model (NHNE — Non-Homogeneous Non-Equilibrium, in its general formulation)

Directly addresses the HEM limitation identified above. It is formulated as a **weighted combination** between the SPI limit ("no time to vaporise") and the HEM limit ("full thermodynamic equilibrium"):

$$\dot{m}_{Dyer} = \frac{\kappa}{1+\kappa}\,\dot{m}_{SPI} + \frac{1}{1+\kappa}\,\dot{m}_{HEM}$$

where the non-equilibrium parameter $\kappa = \sqrt{(P_{up} - P_{down})\,/\,(P_{sat}(T_{up}) - P_{down})}$ is proportional to $\tau_{bubble}/\tau_{residence}$. Large $\kappa$ → less equilibrium → more weight on SPI. Small $\kappa$ → near-equilibrium → more weight on HEM. The weight on SPI, $\kappa/(1+\kappa)$, increases with $\kappa$ — consistent with this physical interpretation.

> **Implementation note.** The original Dyer et al. (2007) paper contained a sign error in the weighting formula (weights on SPI and HEM were swapped). This was corrected by Solomon (2011) and confirmed by Waxman (2013, Eq. 9). The formula above and the implementation in `injector_two_phase.py` reflect the corrected version.

Because it reasonably captures observed behaviour in short orifices without the added complexity of more general non-homogeneous, non-equilibrium models, the Dyer model is the most widely adopted in the university rocketry literature for this class of problem, and is the reference model for this project. Validation against Waxman (2013/2014) gives MAPE = 3.51% at moderate pressure drops — see `validation/waxman_2013_results.md`.

## 3.5 Synthesis for implementation

From the analysis in Sections 1–3, the sizing problem breaks down into two complementary fronts:

1. **Feed line path** (tank → injector inlet) — tracking pressure drop and subcooling margin $\Delta T_{sub}$ along the line, incorporating friction losses (Darcy-Weisbach formulation) and fitting losses (valves, bends). When flashing occurs in the line, the model uses HEM two-phase properties (mixture density $\rho_{mix}$, mixture viscosity $\mu_{mix}$ via McAdams rule) updated segment-by-segment at the local pressure and isenthalpic vapour quality $x(s)$.

2. **Injector orifice** — computing real mass flow, explicitly accounting for possible partial vaporisation inside the orifice itself, via SPI, HEM, or Dyer, depending on how close pressure is to saturation. When fluid arrives two-phase at the inlet ($x_{inlet} > 0$ from feed-line flashing), the HEM two-phase inlet model is used with upstream enthalpy $h_{up} = h_l(T_{tank}) + x_{inlet} \cdot h_{fg}(T_{tank})$.

The two fronts are coupled via the iterative solver in `full_system.py`, which finds the self-consistent operating point where feed-line losses and injector flow are mutually consistent (see `04_implementation.md`, Section 4.5).


## Summary

Near the critical point, N₂O vaporizes with disproportionate ease compared to other fluids; when that vaporization occurs inside the injector orifice, the mixture's density drop — and possible two-phase choking — reduce real mass flow far more severely than the pure-liquid model can predict. The Dyer model, as a bridge between the SPI and HEM limits, is the approach adopted in this project to capture that behavior.

---
*Previous document: [02_spi_model.md](02_spi_model.md) · Next document: [04_implementation.md](04_implementation.md) — computational implementation of the SPI, HEM, and Dyer models, and of the sizing tool (in progress).*
