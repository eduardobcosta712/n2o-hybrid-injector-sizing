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

## Pressure-dependent liquid density (compressibility correction)

**Motivation.** The current model treats saturated liquid density as a function of temperature only, $\rho_{sat}(T)$, consistent with the incompressible-liquid assumption underlying the SPI model (Section 2.4). In reality, liquid density depends weakly on pressure too, through the isothermal compressibility $\kappa_T = -\frac{1}{\nu}\left(\frac{\partial \nu}{\partial P}\right)_T$. For feed lines with large pressure excursions, a first-order correction

$$\rho(T, P) \approx \rho_{sat}(T)\left[1 + \kappa_T \cdot (P - P_{sat}(T))\right]$$

would let `feed_line.py` recompute density segment by segment using local pressure, rather than a single value fixed at the tank temperature.

**Why it isn't included in this version.** The comparison between Case A and Case B in `04_implementation.md` (Section 4.2) shows the constant-density approximation does not distort the central flashing-detection conclusion for the pressure ranges tested. A correlation for $\kappa_T(T)$ specific to liquid N₂O would also need to be sourced separately from the saturation correlations already in use (see `references.md`). This is judged a second-order refinement relative to the two-phase injector modeling (Sections 3–4), which remains the project's core focus.

## Other identified extensions

- Transient regime at motor start-up (the current model assumes steady-state flow).
- Calibration of Cd and Dyer model parameters against the team's own experimental data, rather than literature reference values.
