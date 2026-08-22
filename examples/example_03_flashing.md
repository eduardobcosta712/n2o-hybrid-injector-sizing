# Example 3 — Flashing Case: Diagnosis and Design Correction

This example shows how to use the tool to **diagnose a flashing feed line**, understand the cause, and apply corrections until the system operates safely. It illustrates the physical mechanism that caused unexpected results in several documented university hybrid motor tests.

---

## Scenario

A team designs a motor with N₂O at **22 °C**, tank pressure **54 bar**, and the following feed line:

- 2500 mm of 8 mm ID tubing
- Needle valve (fully open, K = 2.0)
- 90° elbow (standard, K = 0.9)
- Injector: 4 holes, 1.8 mm diameter (Cd = 0.65)
- Chamber pressure: 18 bar

The SPI model was used for the initial design. On test, the motor produces less thrust than expected and combustion is unstable.

---

## Initial configuration — what goes wrong

| Parameter | Value |
|---|---|
| Tank temperature | 22 °C |
| Tank pressure | 54 bar |
| P_sat at 22 °C | **53.68 bar** |
| Subcooling margin | **only 0.32 bar** |

With only 0.32 bar of margin, even small friction losses are enough to push the fluid below its saturation pressure.

**Feed line trace:**

| Segment | Type | P after (bar) | ΔT_sub (K) |
|---|---|---|---|
| Tank exit | — | 54.00 | +0.19 |
| After pipe (2.5 m) | Pipe | 51.99 | −1.46 |
| After needle valve | Fitting (K=2.0) | 51.17 | −2.18 |
| After elbow | Fitting (K=0.9) | 50.79 | −2.51 |

**Result:** Flashing detected after the first pipe segment. The needle valve alone generates enough pressure drop to push the fluid 2 °C below its saturation temperature at that pressure. The injector models cannot be evaluated — the fluid arrives two-phase at the orifice inlet.

**SPI reference (not valid here):** The SPI model would predict approximately 370 g/s. The real mass flow is substantially lower — in documented test cases with similar conditions, reductions of 5–10× below the SPI prediction have been reported.

---

## Diagnosis

The tool's diagnostic panel identifies three contributing factors:

1. **Tank pressure critically close to P_sat** — only 0.32 bar of margin. Any friction loss causes flashing.
2. **Long feed line** — 2.5 m of 8 mm ID tubing generates significant friction losses at the design mass flow.
3. **High-K fitting** — the needle valve (K = 2.0) contributes a large local pressure loss. A ball valve (K = 0.05) would reduce this by a factor of 40.

---

## Correction — step by step

### Step 1: Replace the needle valve

Swap the needle valve (K = 2.0) for a ball valve (K = 0.05). This alone reduces fitting losses by ≈ 97%.

Result: still flashing — the pipe friction alone is enough to consume the 0.32 bar margin.

### Step 2: Increase tank pressure

Raise the regulator to **60 bar**, increasing the subcooling margin from 0.32 bar to 6.32 bar.

Result: still flashing — the 2.5 m pipe is too long.

### Step 3: Increase pipe diameter

Replace the 8 mm ID tubing with **12 mm ID** (pressure drop scales as 1/D⁴ — a 50% increase in diameter reduces friction losses by a factor of ≈ 5).

**Combined fix — all three changes applied:**

| Parameter | Original | Fixed |
|---|---|---|
| Tank pressure | 54 bar | **60 bar** |
| Pipe ID | 8 mm | **12 mm** |
| Valve type | Needle (K=2.0) | **Ball (K=0.05)** |
| Flashing detected | **YES** | **NO** |
| Injector inlet pressure | — | 59.73 bar |
| Real mass flow | N/A | **361.2 g/s** |

After the fix, the system delivers 361 g/s with no flashing and a comfortable subcooling margin throughout the line.

---

## Key lessons

1. **Subcooling margin matters more than absolute pressure.** A tank at 54 bar with only 0.32 bar of margin is more vulnerable than a tank at 45 bar with 5 bar of margin.
2. **The needle valve is the worst fitting for N₂O systems.** Its high K value generates large local pressure losses. Use ball valves for isolation and flow control.
3. **Pipe diameter has a disproportionate effect.** Doubling the diameter reduces friction losses by a factor of 16. Undersized tubing is a common and avoidable cause of feed line flashing.
4. **The SPI model cannot predict this.** It has no mechanism to detect or account for flashing — it always produces a number, even when the real flow is a fraction of that number. This tool exists to close that gap.

---

## How to reproduce in the tool

**Original (flashing) case:**
1. Open **Sizing mode**.
2. T = 22 °C, P = 54 bar, Cd = 0.65, P_chamber = 18 bar. Orifices: 4 × 1.8 mm.
3. Segments: pipe 2500 mm / 8 mm → needle valve / 8 mm → elbow 90° / 8 mm.
4. Observe: Flashing YES, diagnostic panel with three suggestions.

**Corrected case:**
1. Change T = 22 °C, P = **60 bar**.
2. Change pipe segments to **12 mm ID**.
3. Change needle valve to **Ball valve (fully open)**.
4. Observe: Flashing NO, Dyer model applied, mass flow **361.2 g/s**.
