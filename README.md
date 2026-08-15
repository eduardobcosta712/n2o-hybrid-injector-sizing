# n2o-hybrid-injector-sizing

A tool and technical report for predicting and mitigating premature vaporization (*flashing*) of N₂O in the feed systems of paraffin/N₂O hybrid rocket motors, along the tank → line → injector orifice path.

## The problem

N₂O is the most common oxidizer in university-built hybrid motors because it is self-pressurizing and comparatively safe to handle. It has, however, a critical point at only ≈36.4 °C — dangerously close to typical ambient launch temperatures. This makes it exceptionally prone to partial vaporization before it reaches the combustion chamber, whenever pressure in the feed system drops — due to line friction, restrictions, or inside the injector orifice itself.

When this happens, the model most commonly used in the literature and by university teams (SPI — *Single Phase Incompressible*, essentially Bernoulli's equation with an empirical discharge coefficient) is no longer valid, because it implicitly assumes the fluid stays liquid from start to finish. Documented results in the literature (see `docs/references.md`) show real mass flow far below the SPI prediction — by a factor of several times in published cases — which shifts the motor's real O/F far from the design point, directly affecting thrust and combustion stability.

This project exists to give a hybrid propulsion team a way to predict, before testing the motor, whether their feed system is at risk of this, and what design margin they actually have.

## What this project delivers

1. **Theoretical foundation** (`docs/`) — N₂O thermodynamics, the SPI model and its assumptions, and the two-phase models (HEM, Dyer) that correct those assumptions, written from first principles to be understandable without prior knowledge of two-phase flow.
2. **Calculation model** (`src/`) — implementation of the SPI, HEM, and Dyer models, applied to the full path (feed line + injector orifice), validated against published data.
3. **Interactive tool** (`src/`) — design inputs (line geometry, tank orientation, N₂O temperature, injector geometry) → subcooling margin along the path and expected real mass flow rate.

## Scope and known limitations

- N₂O tank temperature is, in this version, a **direct user input** — it is not derived from ambient conditions (location, time of day, solar radiation). Modeling the tank's thermal balance from weather data is a natural extension, identified but deliberately left out of this version's scope to ensure a validated deliverable within the available timeframe (see `docs/future_work.md`).
- The model assumes steady-state flow — it does not capture transient effects during motor start-up.
- Discharge coefficients and Dyer model weighting are typical literature values, not calibrated against this team's own tested injector geometry.

## Repository structure

```
n2o-hybrid-injector-sizing/
├── README.md
├── docs/
│   ├── 01_n2o_thermodynamics.md
│   ├── 02_spi_model.md
│   ├── 03_two_phase_flow.md
│   ├── 04_implementation.md        (in progress)
│   ├── future_work.md
│   └── references.md
├── src/
│   ├── model/                      (calculation functions — SPI, HEM, Dyer)
│   └── interface/                  (interactive tool)
└── results/
    └── (plots and validation, as they become available)
```

## Author

Eduardo Costa — Aerospace Engineering student, Técnico Lisboa. Developed independently as a personal project in hybrid rocket propulsion.
