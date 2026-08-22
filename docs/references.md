# References

Sources used in the theoretical foundation and computational implementation of this project.

---

## Primary source — N₂O property correlations and two-phase model

**Jean-Philyppe, J. (2023).** *A computational model for the design of a nitrous oxide-paraffin wax hybrid rocket engine.* McGill Rocket Team technical report. arXiv:2302.06725.

Used for:
- Saturation pressure correlation P_sat(T) — Appendix A.1 (coefficients originally from Perry's Chemical Engineers' Handbook, re-transcribed here)
- Saturated liquid molar volume correlation ν_l(T) — Appendix A.1
- Saturated vapour and liquid enthalpies h_v(T), h_l(T), and molar volume ν_v(T) — Table A.1, transcribed manually for use as a look-up table in `n2o_saturation_table.csv`
- Dyer (NHNE) injector model formulation and κ parameter — Section 3
- General methodology for the tank → feed line → injector sizing problem

---

## N₂O thermophysical properties — original source

**Green, D.W. & Perry, R.H. (Eds.) (2008).** *Perry's Chemical Engineers' Handbook*, 8th edition. McGraw-Hill.

The saturation property correlations used in `n2o_properties.py` originate from Perry's Handbook and were re-transcribed via the McGill source above. Perry's is the primary experimental data source underlying the correlations.

---

## Two-phase flow models

**Dyer, R.S. (1976).** *The effect of dissolved gas and submicron particles on cavitation inception in water.* PhD thesis, California Institute of Technology.
*(Original formulation of the NHNE non-equilibrium injector model, later adapted for liquid oxidiser applications.)*

**Waxman, B.S., Zimmerman, J.E., Cantwell, B., & Zilliac, G. (2014).** *Mass flow rate and isolation pressure measurements in nitrous oxide with the Dyer injector model.* AIAA 2014-3834.
*(Validation of the Dyer model against N₂O experimental data; documents the SPI over-prediction effect and the Dyer correction.)*

---

## Feed line pressure drop

**Colebrook, C.F. (1939).** Turbulent flow in pipes with particular reference to the transition region between smooth and rough pipe laws. *Journal of the Institution of Civil Engineers*, 11, 133–156.
*(Original implicit friction factor correlation for turbulent pipe flow.)*

**Swamee, P.K. & Jain, A.K. (1976).** Explicit equations for pipe-flow problems. *Journal of the Hydraulics Division, ASCE*, 102(5), 657–664.
*(Explicit approximation to the Colebrook equation used in `feed_line.py`; within ~1% of the implicit solution for the range of interest.)*

---

## Loss coefficients for pipe fittings

**Idel'chik, I.E. (1994).** *Handbook of Hydraulic Resistance*, 3rd edition. CRC Press.
*(Reference for fitting loss coefficients K used in the fitting dropdown table. Individual values also consistent with:)*

**Crane Technical Paper 410 (2013).** *Flow of Fluids Through Valves, Fittings, and Pipe*. Crane Co.

---

## Typical fitting K values used in this project

| Fitting | K | Source |
|---|---|---|
| Ball valve (fully open) | 0.05 | Crane TP 410 |
| Ball valve (½ open) | 5.5 | Crane TP 410 |
| Needle valve (fully open) | 2.0 | Idel'chik |
| Globe valve (fully open) | 10.0 | Crane TP 410 |
| Check valve | 2.5 | Crane TP 410 |
| Elbow 90° (standard) | 0.9 | Crane TP 410 |
| Elbow 90° (long radius) | 0.4 | Crane TP 410 |
| Elbow 45° | 0.4 | Crane TP 410 |
| Tee (flow-through) | 0.4 | Crane TP 410 |
| Tee (branch) | 1.5 | Crane TP 410 |
| Union / coupling | 0.04 | Idel'chik |

---

## Discharge coefficient reference range

The range Cd = 0.61–0.82 cited in the tool (sharp-edged to well-rounded orifices) is consistent with:

**Lienhard, J.H. & Lienhard, J.H. IV (2020).** *A Heat Transfer Textbook*, 5th edition. Phlogiston Press.
*(Section on orifice flow; also available freely at ahtt.mit.edu.)*

---

## N₂O critical point reference values

T_crit = 36.4 °C (309.52 K), P_crit = 72.45 bar — consistent with:

**NIST WebBook, National Institute of Standards and Technology.** Nitrous oxide (N₂O) thermophysical properties. https://webbook.nist.gov/cgi/cbook.cgi?ID=10024-97-2&Type=SatT&Offset=0

---

## Notes

All numerical coefficients in `n2o_properties.py` are sourced from the McGill/Perry reference above and are documented with their source in the module docstring. No coefficients were derived independently or taken from unverified sources, consistent with the project's policy of not using "magic numbers" without a traceable origin.
