# Example 2 — Design Mode: Find Orifice Area for a Target Mass Flow

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
| P_sat at 15 °C | 45.92 bar |
| Subcooling margin | **14.08 bar** |
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

**Feed line evaluation:**
- No flashing detected — pressure at the injector inlet is **59.75 bar**, well above P_sat(15 °C) = 45.92 bar.
- The 12 mm ID line with a single low-K fitting preserves nearly all the tank pressure.

**Orifice sizing:**

| | SPI | Dyer (recommended) |
|---|---|---|
| Total orifice area | 9.525 mm² | **13.711 mm²** |
| Diameter per hole | 1.231 mm | **1.477 mm** |
| Area increase vs SPI | — | **+44.0%** |

**Dyer model parameters:**
- Non-equilibrium parameter κ = 1.238
- Exit vapour quality x = 0.267 (≈ 27% of the mass exits as vapour)

---

## Interpretation

The Dyer model requires a total orifice area **44% larger** than SPI. The difference comes from two-phase effects inside the orifice: even though the feed line delivers fully-subcooled liquid at the injector inlet, the large pressure drop across the orifice (59.75 → 20 bar) causes partial vaporisation inside the orifice itself, reducing the effective mixture density and therefore the mass flow for a given area.

If the team manufactured the SPI-sized holes (1.231 mm), the system would deliver only:

```
m_dot_Dyer(A_SPI) ≈ 500 / 1.44 ≈ 347 g/s
```

instead of the target 500 g/s — a shortfall of approximately 30%.

**Combustion stability check:**
- Injector pressure drop (59.75 − 20 = 39.75 bar) / P_chamber (20 bar) = **199%** — well above the 15% minimum. Stable.

**Recommended manufacture:** 8 holes of **1.477 mm** diameter (or nearest available drill size, e.g. 1.50 mm — re-run the tool in Sizing mode with 1.50 mm to verify the resulting mass flow).

---

## How to reproduce in the tool

1. Open **Design mode**.
2. Set sidebar: T = 15 °C, P = 60 bar, Cd = 0.65, P_chamber = 20 bar.
3. Set target: 500 g/s, 8 holes.
4. Add feed line segments as listed above.
5. Read the results — Dyer hole diameter should show **1.477 mm** and area increase **44.0%**.
