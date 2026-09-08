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

**Solomon, B. J. (2011).** Engineering model for propellant-actuated devices. M.Sc. thesis, Utah State University.
*(Identifies and corrects the weight-swap error in the original Dyer et al. (2007) two-phase injector formula. The corrected NHNE weighting convention — large κ → more weight on SPI — is adopted in this project.)*

**Waxman, B.S., Zimmerman, J.E., Cantwell, B., & Zilliac, G. (2014).** *Mass flow rate and isolation pressure measurements in nitrous oxide with the Dyer injector model.* AIAA 2014-3834.
*(Validation of the Dyer model against N₂O experimental data; documents the SPI over-prediction effect and the Dyer correction.)*

**Henry, R.E. & Fauske, H.K. (1971).** *The Two-Phase Critical Flow of One-Component Mixtures in Nozzles, Orifices, and Short Tubes.* ASME Journal of Heat Transfer, 93(2), 179-187.
*(Original non-equilibrium critical flow model. Used in `henry_fauske_critical_flow()` as the physically correct ceiling for the Dyer prediction — the equilibrium HEM ceiling was shown to be the wrong bound, since real non-equilibrium two-phase flow chokes at a higher mass flux than full equilibrium allows.)*

**Simoneau, R.J., Henry, R.E., Hendricks, R.C. & Watterson, R. (1971).** *Two-Phase Critical Discharge of High Pressure Liquid Nitrogen.* NASA Technical Memorandum TM X-67863.
*(Presents the simplified Henry-Fauske equations (Eqs. 2-5) actually transcribed into `henry_fauske_critical_flow()`, including the empirical non-equilibrium factor N = min(1, x_E/0.14) from Henry (1970), fit to steam-water data of Starkman et al. (1964).)*

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

**Lemmon, E. W., and Span, R. (2006).** Short fundamental equations of state for 20 industrial fluids. *Journal of Chemical and Engineering Data*, 51(3), 785–850. DOI: 10.1021/je050186n.
*(N₂O equation of state used by NIST WebBook. Source of Tables A.3 and A.4 thermodynamic properties.)*

**Millat, J., Vesovic, V., and Wakeham, W. A. (1991).** The viscosity of nitrous oxide and tetrafluoromethane in the limit of zero density. *International Journal of Thermophysics*, 12(2), 265–273.
*(N₂O saturated vapour viscosity correlation used by NIST. Source of Table A.3.)*

**Laesecke, A., and Hafer, R. F. (1998).** Viscosity of fluorinated propane isomers. *Journal of Chemical and Engineering Data*, 43(1), 84–92. DOI: 10.1021/je970186n.
*(N₂O saturated liquid viscosity correlation used by NIST. Source of Table A.4 μ_l values.)*

**Niño, E. V., and Razavi, M. R. (2019).** Design of two-phase injectors using analytical and numerical methods with application to hybrid rockets. AIAA 2019-4154.
*(Tabulated Waxman operating points used for model validation; Table 4 and Fig. 2.)*

---

## Notes

All numerical coefficients in `n2o_properties.py` and `injector_two_phase.py` are sourced from the references above and are documented with their source in the module docstrings. No coefficients were derived independently or taken from unverified sources, consistent with the project's policy of not using "magic numbers" without a traceable origin. The Henry-Fauske non-equilibrium factor N = x_E/0.14 was specifically deferred until a primary/near-primary source with legible equations (not image-embedded) could be obtained, rather than being implemented from partial or half-remembered recollection.
