# Future Work

Extensions identified during this project's development, deliberately left out of the current version's scope to guarantee a validated deliverable within the available timeframe.

## Estimating tank temperature from ambient conditions

**Motivation.** N₂O tank temperature is, in the current version, a direct user input. In practice, that temperature results from a thermal balance between the tank and its environment (incident solar radiation, convection with ambient air, insulation, orientation, exposure), which varies with geographic location, time of day, and weather conditions. Since N₂O operates dangerously close to its critical point (≈36.4 °C — see `01_n2o_thermodynamics.md`), this environmental variation could be enough to significantly shift the subcooling margin available at the start of the path, even before any losses along the line or the injector.

**Why it isn't included in this version.** Correctly modeling this thermal balance requires:

- A model of incident solar radiation as a function of location, time, season, and tank orientation;
- A model of convection with ambient air, dependent on wind speed and tank geometry;
- A transient energy balance (dependent on exposure time, not just instantaneous conditions), including the effect of N₂O's own internal phase change on tank temperature.

This is a thermal modeling problem of comparable scope and complexity to the flow model that forms the core of this project, so including it at this stage would risk compromising the validation and robustness of the flow model itself.

**Proposed future implementation.** Retrieval of weather data (air temperature, solar irradiance) via a public API (e.g. Open-Meteo, no API key required), combined with a simplified tank-ambient thermal equilibrium model as a first approximation — explicitly assuming no transient effects and no direct solar radiation in that first iteration — before a full transient model is justified.

## Other identified extensions

- Transient regime at motor start-up (the current model assumes steady-state flow).
- Calibration of $C_d$ and Dyer model parameters against the team's own experimental data, rather than literature reference values.
