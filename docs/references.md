# References

Sources used in the theoretical foundation and computational implementation of this project.

Bibliographic entries were re-checked in the September 2026 audit. Entries the audit **could verify** against an independent source are marked ; entries it **could not** verify are collected in "Open bibliographic points" at the end and must not be cited elsewhere without checking.

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

**Dyer, J., Doran, E., Dunn, Z., Lohner, K., et al. (2007).** *Modeling Feed System Flow Physics for Self-Pressurizing Propellants.* AIAA 2007-5702. ✔
*(Original formulation of the non-homogeneous non-equilibrium (NHNE) two-phase injector model used in this project. Its weighting formula contained a swap of the two weights, corrected by Solomon (2011).)*

**Solomon, B. J. (2011).** *Engineering Model to Calculate Mass Flow Rate of a Two-Phase Saturated Fluid through an Injector Orifice.* M.Sc. thesis, Utah State University. ✔
*(Identifies and corrects the weight-swap error in the original Dyer et al. (2007) two-phase injector formula. The corrected NHNE weighting convention — large κ → more weight on SPI — is adopted in this project.)*

**Waxman, B.S., Zimmerman, J.E., Cantwell, B., & Zilliac, G. (2013).** *Mass Flow Rate and Isolation Characteristics of Injectors for Use with Self-Pressurizing Oxidizers in Hybrid Rockets.* AIAA 2013-3636. ✔
*(Source of the equations cited in the code as "Waxman (2013) Eq. (5)" — the HEM critical flow — and "Eq. (9)" — the corrected Dyer weighting — and of the discharge-coefficient measurements; documents the SPI over-prediction effect and the Dyer correction.)*

**Waxman, B.S. (2014).** *An Investigation of Injectors for Use with High Vapour Pressure Propellants with Applications to Hybrid Rockets.* PhD thesis, Stanford University.
*(Original experimental dataset behind the Niño & Razavi operating points used for validation.)*

**Henry, R.E. & Fauske, H.K. (1971).** *The Two-Phase Critical Flow of One-Component Mixtures in Nozzles, Orifices, and Short Tubes.* ASME Journal of Heat Transfer, 93(2), 179-187.
*(Original non-equilibrium critical flow model. Used in `henry_fauske_critical_flow()` as the physically correct ceiling for the Dyer prediction — the equilibrium HEM ceiling was shown to be the wrong bound, since real non-equilibrium two-phase flow chokes at a higher mass flux than full equilibrium allows.)*

**Simoneau, R.J., Henry, R.E., Hendricks, R.C. & Watterson, R. (1971).** *Two-Phase Critical Discharge of High Pressure Liquid Nitrogen.* NASA Technical Memorandum TM X-67863.
*(Presents the simplified Henry-Fauske equations (Eqs. 2-5) actually transcribed into `henry_fauske_critical_flow()`, including the empirical non-equilibrium factor N = min(1, x_E/0.14) from Henry (1970), fit to steam-water data of Starkman et al. (1964).)*

**Niño, E. V., and Razavi, M. R. (2019).** *Design of two-phase injectors using analytical and numerical methods with application to hybrid rockets.* AIAA 2019-4154.
*(Tabulated Waxman operating points used for model validation; Table 4 and Fig. 2. Not independently re-checked in the audit — see "Open bibliographic points".)*

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

## N₂O critical point and equation of state

T_crit = 36.4 °C (309.52 K), P_crit = 72.45 bar.  (consistent with the NIST WebBook and the CoolProp fluid page for nitrous oxide.)

**NIST WebBook, National Institute of Standards and Technology.** Nitrous oxide (N₂O) thermophysical properties. https://webbook.nist.gov/cgi/cbook.cgi?ID=10024-97-2&Type=SatT&Offset=0
*(Source of the numerical data in Tables A.3 and A.4 of `n2o_saturation_table.csv`.)*

**Lemmon, E. W., and Span, R. (2006).** Short fundamental equations of state for 20 industrial fluids. *Journal of Chemical and Engineering Data*, 51(3), 785–850. DOI: 10.1021/je050186n. 
*(N₂O equation of state used by NIST WebBook and CoolProp. Source of the thermodynamic properties — including entropy — in Table A.4.)*

---

## CoolProp and alternatives (added September 2026, Priority 4)

**Bell, I.H., Wronski, J., Quoilin, S., & Lemort, V. (2014).** Pure and pseudo-pure fluid thermophysical property evaluation and the open-source thermophysical property library CoolProp. *Industrial & Engineering Chemistry Research*, 53(6), 2498-2508. DOI: 10.1021/ie4033999. 
*(The open-source library `n2o_properties.py` now uses for all saturation thermodynamics -- P_sat, T_sat, densities, enthalpy, entropy, cp_l -- via its `PropsSI` interface. CoolProp's own N₂O fluid page (coolprop.org/fluid_properties/fluids/NitrousOxide.html) lists only the equation of state (Lemmon & Span 2006, above) and a surface-tension correlation as references so CoolProp does not supply an N₂O viscosity model; this project's own NIST-table-based `mu_liquid_sat`/`mu_vapor_sat` remain the viscosity source.)*

**Possible lead for the unverified NIST viscosity attribution below.** Huber, M.L. (2018). *Models for Viscosity, Thermal Conductivity, and Surface Tension of Selected Pure Fluids as Implemented in REFPROP v10.0.* NIST Interagency/Internal Report NISTIR 8209. DOI: 10.6028/NIST.IR.8209. This report describes *interim* transport-property models NIST built specifically for fluids -- like N₂O -- that lack a published reference-quality viscosity/thermal-conductivity model, which matches the situation this project is trying to source. Not yet confirmed to actually cover N₂O, or to be the correlation behind the WebBook's saturated-viscosity table used in Tables A.3/A.4 of `n2o_saturation_table.csv`.

**If CoolProp cannot be installed** (e.g. a fully offline machine -- `pip install CoolProp` needs network access; it is not otherwise unusual or paid, standard wheels exist for Windows/Linux/macOS and Python 3.8-3.13):
1. **REFPROP** (NIST) — the reference implementation CoolProp itself is validated against, and CoolProp can be configured to call it directly (`AbstractState("REFPROP", "NitrousOxide")`) if a REFPROP licence and installation are available (a Técnico institutional licence, if one exists, would qualify). Likely the single best option if accessible, but it is commercial software requiring a separate purchase/install, not just a `pip install`.
2. **Revert to the closed-form Perry/McGill correlations** this project used before September 2026 (kept in `tests/test_n2o_properties.py` as `_perry_P_sat`/`_perry_rho_l`, and in `n2o_saturation_table.csv` Table A.1) — no dependency, but reintroduces the ~2-3 % P_sat error and 3-5 % h_fg error this integration was meant to remove (see docs/04_implementation.md, "Error sources"), and loses the extended (up to the critical point) entropy range needed by the isentropic and Henry-Fauske choking models.
3. **Implement the Lemmon & Span (2006) Helmholtz-energy equation of state by hand** — the correct long-term fix if neither of the above is acceptable, but a substantial undertaking (the EOS has dozens of fitted terms and needs a numerical solver for phase-equilibrium conditions) that duplicates what CoolProp already provides for free; not recommended unless CoolProp becomes genuinely and permanently unavailable.

## Open bibliographic points

These items were flagged, not resolved, by the September 2026 audit. Each needs a manual check against the original source or the NIST WebBook page before being relied on.

1. **Origin of the NIST viscosity correlations (Tables A.3 and A.4).** Earlier versions of this file, of `n2o_properties.py` and of the CSV header attributed the saturated-vapour viscosity to *Millat, Vesovic & Wakeham (1991), Int. J. Thermophys. 12, 265* ("viscosity of nitrous oxide and tetrafluoromethane in the limit of zero density") and the saturated-liquid viscosity to *Laesecke & Hafer (1998), J. Chem. Eng. Data 43(1), 84* ("Viscosity of fluorinated propane isomers"). The audit could not confirm either attribution, and the second is almost certainly wrong: by its title it concerns fluorinated propanes, not N₂O. The CoolProp documentation for nitrous oxide lists only an equation of state and a surface-tension correlation as references, i.e. it carries no N₂O viscosity model to cross-check against. **Action:** open the WebBook saturation table for N₂O (link above) and note the viscosity reference it lists; until then treat the viscosity numbers as "NIST WebBook data, correlation reference unverified". The numerical data themselves are unaffected.
2. **Niño & Razavi (2019), AIAA 2019-4154** — title and content taken from the project's earlier notes; not re-checked.
3. **Waxman et al. (2014).** Earlier versions listed *"Mass flow rate and isolation pressure measurements in nitrous oxide with the Dyer injector model", AIAA 2014-3834*. The audit could not find such a paper and has removed the entry; the verified Waxman references are the 2013 AIAA paper and the 2014 PhD thesis above.
4. **Author list of Dyer et al. (2007).** The paper's identity (AIAA 2007-5702, title above) was verified; the full author list should be copied from the paper itself.
5. **`waxman_2013_experimental_data.csv`** records upstream conditions P1 = 4.96 MPa, supercharge 1.26 MPa, while the validation uses P1 = 4.36 MPa, supercharge 0.62 MPa (Niño & Razavi Table 4). Both may be correct (different tests in the same paper), but the CSV is not read by any code; see `validation/waxman_2013_results.md`.

---

## Notes

All numerical coefficients in `n2o_properties.py` and `injector_two_phase.py` are sourced from the references above and are documented with their source in the module docstrings. No coefficients were derived independently or taken from unverified sources, consistent with the project's policy of not using "magic numbers" without a traceable origin. The Henry-Fauske non-equilibrium factor N = x_E/0.14 was specifically deferred until a primary/near-primary source with legible equations could be obtained.

---

## Fuel grain sizing (`grain_sizing.py`, Priority 3)

**Marxman regression rate correlation — theoretical basis:**

**Marxman, G.A. (1964).** Combustion in the turbulent boundary layer on a vaporizing surface. *Tenth Symposium (International) on Combustion*, 1337–1349.
*(Original turbulent-boundary-layer derivation of the G_o^n regression-rate scaling; predicts n ≈ 0.8, which real fitted data for most propellant combinations falls below — see below.)*

**On why (a, n) are NOT shipped as fixed per-fuel defaults in this project** — evidence of scatter between independent studies of nominally the same fuel/oxidiser pair, consistent with the general understanding that these are test-article-specific empirical fits, not universal material constants:

- Paraffin/N₂O regression rates of ≈2 mm/s, ≈3.5 mm/s, and 4–5 mm/s have each been separately reported in the literature for broadly comparable oxidiser mass flux ranges, with at least one research group (Libre-Space/ULB paraffin-N₂O studies) explicitly noting they could not fully reconcile their measured rate against other published values for the same propellant combination.
- **Karabeyoglu, M.A., Cantwell, B.J., & Zilliac, G. (2005).** Development of Scalable Space-Time Averaged Regression Rate Expressions for Hybrid Rockets. AIAA 2005-3544. — reports HTPB/N₂O regression rates roughly 3–5× lower than paraffin/N₂O at comparable oxidiser mass flux, illustrating the fuel-to-fuel spread that compounds the study-to-study spread within a single fuel.

**Fuel densities (genuine material properties, used as defaults in `FUEL_PROPERTIES`):**

- **Paraffin wax**: 900.0 kg/m³ (theoretical/bulk), 894.0 kg/m³ (experimental, ~7% void fraction) — both reported in the same characterisation source; 900.0 kg/m³ used as the module default.
- **PMMA**: 1185.2 kg/m³, measured value from **"Characterization of PolyMethylMethAcrylate as a Fuel for Hybrid Rocket Motors"** (AIAA), a specific clear-cast PMMA grain characterisation study — consistent with commercial PMMA sheet stock, typically quoted 1180–1190 kg/m³.
- **HTPB**: 920.0 kg/m³, a standard value in the propulsion literature. (HTPB's *other* thermochemical properties — notably heat of formation — are much less consistently reported across sources; density specifically is not contested in the same way.)
- **ABS**: 1050.0 kg/m³, standard commercial ABS resin bulk density (typically quoted 1.04–1.06 g/cm³); 3D-printed grains may have a somewhat lower effective density depending on infill fraction.

**ABS as a hybrid fuel (comparison to HTPB):**

**Whitmore, S.A. et al.** Analytical and Experimental Comparisons of HTPB and ABS as Hybrid Rocket Fuels. AIAA paper.
*(ABS/N₂O regression rate and combustion performance reported comparable to, slightly below, HTPB/N₂O in direct side-by-side testing at equal grain geometry — cited in `grain_sizing.py`'s FUEL_PROPERTIES reference notes for ABS. Full citation details, e.g. AIAA paper number, still to be added.)*
