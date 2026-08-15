# 1. N₂O Basic Thermodynamics

## 1.1 Pure substances and what governs their state

N₂O (nitrous oxide) is, for the purposes of this document, a pure substance — a single chemical compound, not a mixture. This matters because, for a pure substance, the thermodynamic state (liquid, vapor, or a mix of both) is fully determined by **only two independent properties**. Pressure $P$ and temperature $T$ are usually chosen, being the easiest to measure and control in a real system.

There is, however, an important exception to this rule, which is the starting point of this entire document: **when the substance is changing phase** (liquid to vapor, or vice versa), $P$ and $T$ stop being independent of each other.

## 1.2 Vapor pressure / saturation pressure

Consider a sealed, airtight container, partially filled with liquid N₂O, with the remaining volume under vacuum, held at constant temperature. Part of the liquid evaporates: molecules at the surface with enough thermal energy escape into the vapor phase, while vapor molecules collide with the surface and some condense back.

Initially, the evaporation rate exceeds the condensation rate and pressure rises. As vapor pressure increases, the condensation rate rises too, until a **dynamic equilibrium** is reached: evaporation rate equals condensation rate, and pressure stabilizes at a specific value for that temperature.

This equilibrium value is called the **vapor pressure** or **saturation pressure**, $P_{sat}$, and is a function of temperature alone:

$$P_{sat} = f(T) \quad \text{only}$$

Raising the container's temperature increases the thermal energy of the liquid's molecules, makes evaporation easier, and shifts the new equilibrium to a higher pressure. There is thus a unique curve $P_{sat}(T)$ for each substance — the **saturation curve**, or **liquid-vapor coexistence curve**. For N₂O, as an order-of-magnitude reference:

| Temperature | $P_{sat}$ |
|---|---|
| 0 °C | ≈ 31.3 bar |
| 20 °C | ≈ 50.9 bar |
| 36 °C | ≈ 72.5 bar (close to the critical point — see 1.5) |

Note that N₂O's saturation pressure is high even at room temperature — this is why N₂O is *self-pressurizing* in rocket tanks: no external pressurant gas is needed, because N₂O itself generates enough pressure as long as it coexists in liquid and vapor phase inside the tank.

## 1.3 The central practical consequence

> **If the temperature of a liquid in equilibrium with its vapor is known, its pressure is automatically determined. The two quantities cannot be chosen independently under these conditions.**

This has a direct implication for a rocket feed system: if the liquid N₂O in the tank is at, say, 20 °C, the tank pressure **must** be ≈50.9 bar (assuming liquid-vapor coexistence, the normal case for a self-pressurized tank). It is not a possible equilibrium state to have liquid N₂O at 20 °C with the tank at 30 bar.

## 1.4 Subcooled liquid vs. saturated liquid

With the curve $P_{sat}(T)$ defined, the liquid's state falls into three categories:

- **Saturated liquid** — sits exactly on the curve; its pressure is exactly $P_{sat}(T)$ for its temperature. Any small perturbation (a pressure drop or a temperature increase) triggers vaporization.
- **Subcooled liquid** — sits at a pressure **higher** than $P_{sat}(T)$ for its temperature; it has a margin before reaching saturation.
- **Saturated vapor** — the vapor-side equivalent.

The **degree of subcooling** is defined as:

$$\Delta T_{sub} = T_{sat}(P) - T$$

where $T_{sat}(P)$ is the saturation temperature corresponding to the fluid's current pressure (the inverse of $P_{sat}(T)$). If $\Delta T_{sub} > 0$, the liquid is subcooled, with margin; if $\Delta T_{sub} = 0$, it is exactly saturated, with no margin.

This is the central concept of the entire project: the problem of premature vaporization in the feed system essentially comes down to identifying at what point along the path the N₂O, which started out subcooled, lost that margin and crossed the saturation curve.

## 1.5 The critical point

Following the curve $P_{sat}(T)$ to increasing temperatures and pressures, a special point is reached — the **critical point**. For N₂O:

$$T_{crit} \approx 36.4\ ^\circ\text{C}, \qquad P_{crit} \approx 72.5\ \text{bar}$$

Above this temperature, the distinction between liquid and vapor disappears — the substance becomes a single homogeneous supercritical fluid, regardless of applied pressure; the saturation curve terminates at that point.

How close $T_{crit}$ is to ambient temperature (compare with water, $T_{crit} \approx 374\ ^\circ\text{C}$) is the fundamental reason N₂O is significantly more prone to vaporization effects than more conventional liquid propellants: a tank sitting in the sun, or even just on a warm day, can bring N₂O close to its critical point. In that neighborhood, small pressure or temperature changes cause disproportionately large changes in fluid properties (density, latent heat of vaporization — see 1.6 and Section 3).

## 1.6 Flashing

Consider subcooled liquid N₂O ($\Delta T_{sub} > 0$) flowing — for example, along a tube, from the tank toward the injector. In any flow with friction, or through a restriction (a valve, an area reduction, the injector orifice itself), fluid pressure drops along the path (energy conservation — see Section 2 for the formal derivation via Bernoulli's equation).

In a sufficiently fast flow, there is no significant time or mechanism for heat exchange with the surroundings, so fluid temperature stays approximately constant while pressure drops. Since $T_{sat}(P)$ decreases with $P$ (the saturation curve is increasing in $P(T)$), the margin $\Delta T_{sub} = T_{sat}(P) - T$ **shrinks** as $P$ falls. If pressure drops enough, $P = P_{sat}(T)$ is reached — the margin vanishes — and any further pressure drop forces a fraction of the liquid to vaporize instantly, with no external heat input.

This phenomenon is called **flashing** (sudden vaporization). The energy needed to vaporize that fraction (the latent heat) is drawn from the liquid itself, which cools slightly as it gives it up — which explains, for instance, the frost that forms on the outside of N₂O lines during testing (the internal fluid cools enough to condense moisture from the ambient air).

## 1.7 Two-phase flow and vapor quality

Once flashing begins, the flow is no longer single-phase — liquid and vapor coexist simultaneously in the same flow: **two-phase flow**. The vaporized fraction is quantified through the **vapor quality**, $x$:

$$x = \frac{m_{vapor}}{m_{vapor} + m_{liquid}}$$

that is, the mass fraction in the vapor phase. $x = 0$ corresponds to pure saturated liquid; $x = 1$, to pure saturated vapor; intermediate values, to the two-phase mixture.

$x > 0$ at any point in the feed system is the central criterion used in this project as a signal that the pure-liquid model is no longer valid at that point — with the consequences detailed in Section 3 (sharp reduction in mixture density, possible two-phase choking, and a drop in real mass flow relative to what a single-phase model would predict).

## Summary

N₂O only behaves as a predictable pure liquid while its local pressure stays above $P_{sat}(T)$ for its local temperature. How close its critical point is to ambient temperature makes this margin exceptionally easy to lose to flow effects (pressure drop from friction or restrictions) — and once lost, the system enters two-phase operation, whose behavior is covered in the following sections.

---
*Next document: [02_spi_model.md](02_spi_model.md) — the SPI (Single Phase Incompressible) model used for injector sizing, derived from Bernoulli's equation, and the implicit assumptions that make it invalid under flashing conditions.*
