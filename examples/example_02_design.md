# Example 2 — Design Mode: Find Orifice Area for a Target Mass Flow

> **Property-backend status: CONFIRMED (September 2026).** All figures below
> were regenerated after installing the real CoolProp package (8.0.0) and
> confirming the full test suite (262 tests) against it. They replace the
> earlier, Perry/McGill-backed figures (32.7% area increase, 13.638 mm²),
> which are now superseded.

This example shows how to use **Design mode** to find the orifice diameter required to deliver a specific oxidiser mass flow, with the Dyer correction applied so the injector is not under-sized relative to the SPI prediction.

---

## Scenario

A hybrid motor design specifies an oxidiser mass flow of **500 g/s** at a chamber pressure of **20 bar**. The team has:

- N₂O tank at **15 °C**, pressurised to **60 bar**
- Feed line: 1000 mm of 12 mm ID tubing → ball valve → 600 mm of 12 mm ID tubing
- Injector plate: **8 holes** (Cd = 0.65)

The question is: what diameter should each hole be, and how does the SPI-sized injector compare to the Dyer-corrected one?

---

## Inputs

| Parameter | Value |
|---|---|
| Tank temperature | 15 °C |
| Tank pressure | 60 bar |
| P_sat at 15 °C | 45.04 bar |
| Subcooling margin | **14.96 bar** |
| Target mass flow | 500 g/s |
| Chamber pressure | 20 bar |
| Discharge coefficient Cd | 0.65 |
| Number of orifices | 8 |

**Feed line segments:**

| # | Type | Length | ID | K |
|---|---|---|---|---|
| 1 | Pipe | 1000 mm | 12 mm | — |
| 2 | Fitting | — | 12 mm | 0.05 (ball valve, fully open) |
| 3 | Pipe | 600 mm | 12 mm | — |

---

## Results

**Feed line evaluation (at the target flow):**
- No flashing detected — pressure at the injector inlet is **59.77 bar**, well above P_sat(15 °C) = 45.04 bar.
- The 12 mm ID line with a single low-K fitting preserves nearly all the tank pressure.

**Orifice sizing:**

| | SPI | Dyer (recommended) |
|---|---|---|
| Total orifice area | 9.521 mm² | **12.497 mm²** |
| Diameter per hole | 1.231 mm | **1.410 mm** |
| Area increase vs SPI | — | **+31.3%** |

**Dyer model parameters:**
- Non-equilibrium parameter κ = 1.260
- Exit vapour quality x = 0.259 (≈ 26% of the mass exits as vapour)

---

## Non-equilibrium choking diagnostic

Design mode also shows a **choking warning**: at these tank/chamber conditions the Henry-Fauske (1971) non-equilibrium choking ceiling is **456.2 g/s** for the recommended area — below the 500 g/s target itself.

**Important physical point:** this ceiling scales linearly with orifice area, exactly like the Dyer prediction, so the *ratio* between them does not depend on the area chosen. This is not something the Dyer area above can be adjusted to avoid. If the ceiling turns out, with future validation, to be the physically binding one, the 500 g/s target would not be achievable at this ΔP (about 40 bar) regardless of hole sizing — the fix would be tank pressure or chamber pressure, not area. As with Example 1, this ceiling is theoretically sound but not experimentally confirmed at this example's pressure drop (outside the 8–14 bar Waxman-validated band), so the 12.497 mm² / 1.410 mm Dyer sizing above remains the tool's primary recommendation — the warning is a flag to investigate further, not a correction to apply by hand.

---

## Interpretation

The Dyer model requires a total orifice area **31.3% larger** than SPI. The difference comes from two-phase effects inside the orifice: even though the feed line delivers fully-subcooled liquid at the injector inlet, the large pressure drop across the orifice (59.77 → 20 bar) causes partial vaporisation inside the orifice itself, reducing the effective mixture density and therefore the mass flow for a given area.

If the team manufactured the SPI-sized holes (1.231 mm), the Dyer model says the system would deliver only:

```
m_dot_Dyer(A_SPI) ≈ 380.9 g/s
```

instead of the target 500 g/s — a shortfall of about 24%.

**Combustion stability check:**
- Injector pressure drop (59.77 − 20 = 39.77 bar) / P_chamber (20 bar) = **198.8%** — well above the 15% minimum. Stable.

**Recommended manufacture:** 8 holes of **1.410 mm** diameter (or the nearest available drill size — re-run the tool in Sizing mode with that diameter to verify the resulting mass flow). Remember that the Dyer model is validated only at 8–14 bar pressure drops: treat the recommendation as a model estimate and confirm it with a water/cold-flow test before machining the final plate.

---

## How to reproduce in the tool

1. Open **Design mode**.
2. Set sidebar: T = 15 °C, P = 60 bar, Cd = 0.65, P_chamber = 20 bar.
3. Set target: 500 g/s, 8 holes.
4. Add feed line segments as listed above.
5. Read the results — Dyer hole diameter should show **1.410 mm** and area increase **31.3%**. An amber choking-ceiling warning (456.2 g/s) should also appear underneath the area caption.
