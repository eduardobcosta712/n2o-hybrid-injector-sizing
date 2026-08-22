# 2. The SPI (Single Phase Incompressible) Model

## 2.1 Why pressure drops at a restriction

Before deriving the orifice model, it's worth clearing up a point of intuition that is frequently inverted: why does an area restriction (like an injector orifice) cause a pressure **drop**, rather than an increase.

By conservation of mass, in steady state, for a tube with no branches:

$$\dot{m} = \rho \, A \, v = \text{constant along the tube}$$

For a liquid (practically incompressible, ρ ≈ constant), a decrease in area A forces an increase in velocity v, to keep $\dot m$ constant.

From Bernoulli's equation, along a streamline, for frictionless flow with no external work, and neglecting elevation changes:

$$P + \frac{1}{2}\rho v^2 \approx \text{constant}$$

The sum of pressure energy and kinetic energy is what stays constant — not each term on its own. If v increases due to the restriction, the term ½ρv² increases, and P must necessarily drop for the sum to remain constant.

Physically: for the fluid to accelerate as it enters the narrower section, it needs a net force in that direction — provided by the pressure difference between upstream (higher) and inside/immediately downstream of the restriction (lower). Pressure energy converts into kinetic energy. It is at the point of minimum cross-section (the *vena contracta*, for an orifice) that velocity is maximum and pressure is minimum — and therefore the point where the subcooling margin ($\Delta T_{sub}$, Section 1.4) is smallest, and where the risk of crossing the saturation curve is greatest.

## 2.2 Model derivation

Consider the injector orifice: point 1, immediately before the orifice, inside the feed line (wide cross-section, so $v_1 \approx 0$ by mass conservation relative to the orifice's much smaller area); point 2, at the *vena contracta*, immediately after the orifice.

From Bernoulli's equation between the two points:

$$P_1 + \frac{1}{2}\rho v_1^2 = P_2 + \frac{1}{2}\rho v_2^2$$

With $v_1 \approx 0$:

$$\frac{1}{2}\rho v_2^2 = P_1 - P_2 \equiv \Delta P$$

$$v_2 = \sqrt{\frac{2 \Delta P}{\rho}}$$

## 2.3 Mass flow rate and the discharge coefficient

From $\dot{m} = \rho A v$:

$$\dot{m}_{theoretical} = \rho \, A \, v_2 = A\sqrt{2\rho \Delta P}$$

The real measured flow rate is consistently lower than this theoretical value, for two reasons:

1. The jet does not occupy the full geometric area of the orifice as it exits — it contracts (the *vena contracta* effect proper; the fluid doesn't perfectly follow the orifice edge).
2. There is energy loss to friction/turbulence not captured by the ideal, frictionless Bernoulli equation.

A **discharge coefficient**, $C_d$, is introduced (dimensionless, typically 0.6–0.9 for rocket injector orifices, depending on geometry — sharp edges: lower $C_d$; rounded orifices: $C_d$ closer to 1), determined empirically:

$$\boxed{\dot{m}_{SPI} = C_d \, A \sqrt{2 \rho \Delta P}}$$

This is the SPI (*Single Phase Incompressible*) equation — the reference model used by most university teams for injector area sizing. Given a target mass flow (set by the design O/F, from the motor's thermochemical sizing) and a design $\Delta P$ (typically 15–20% of $P_c$, for combustion stability reasons, outside the scope of this document), the equation is inverted to obtain the required orifice area.

## 2.4 Implicit assumptions and where they fail

The derivation above implicitly assumes:

1. **Single phase** — $\rho$ corresponds to a pure liquid, a single well-defined value; no vapor or two-phase mixture present.
2. **Incompressibility** — $\rho$ constant along the path, despite the pressure drop. A good approximation for liquid far from saturation; it degrades even before actual vaporization begins, because liquid density becomes more pressure-sensitive near saturation.
3. **Frictionless flow**, partially corrected by $C_d$ — but the empirical value of $C_d$ is typically calibrated with well-behaved fluids (often water), not necessarily representative of N₂O near saturation.

Assumption (1) is the most critical one in the context of this project. From Section 1.6: along the flow, P decreases. If at any point between the tank and the orifice exit — including inside the orifice itself, where P reaches its minimum (Section 2.1) — pressure crosses $P_{sat}(T)$, the single-phase assumption stops holding. The SPI model keeps computing as if the fluid were dense, incompressible liquid, when in reality a meaningful fraction is already vapor, with much lower density.

This mechanism explains the large-magnitude discrepancies between predicted (SPI) and measured mass flow documented in the literature (see `references.md`) — not a 10–20% calibration error, but a physical regime change the model does not capture.

## Summary

The SPI model, $\dot{m} = C_d A \sqrt{2\rho \Delta P}$, is Bernoulli's equation applied to the injector orifice, empirically corrected by $C_d$. It is widely used because it is simple and works well when the fluid stays pure subcooled liquid throughout the path — the exact assumption violated under flashing conditions, addressed in Section 3.

---
*Previous document: [01_n2o_thermodynamics.md](01_n2o_thermodynamics.md) · Next document: [03_two_phase_flow.md](03_two_phase_flow.md) — the HEM and Dyer models, which relax the single-phase assumption.*
