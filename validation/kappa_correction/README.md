# validation/kappa_correction/

Exploratory validation of a supercharge-dependent correction of Dyer's kappa
(`dyer_mass_flow_corrected()` in `src/model/injector_two_phase.py`), on the
digitised Waxman (2013) Fig. 13 data (`../digitized/`). Nothing here changes the
model's default behaviour. Written up in `../waxman_2013_results.md`, Sections
10 and 11, and in `docs/future_work.md`, Priority 1.

Run from anywhere, CoolProp installed, in this order:

| Script | Question it answers | Output |
|---|---|---|
| `loocv_kappa_correction.py` | Does the gated correction (beta, supercharge_ref) survive out of sample? Leave-one-curve-out over the 9 supercharge curves; Cd re-calibrated in every fold; Part A (injector 2) never used in a fit; limit check kappa' -> 0 gives Dyer -> HEM. | `loocv_kappa_results.txt` |
| `loocv_functional_forms.py` | Which form is the simplest that works (linear, soft gate, parameter-free)? Is it a disguised Cd shift (Cd fitted jointly in mode B)? | `loocv_functional_forms_results.txt` |
| `scan_threshold_continuity.py` | Does the correction remove the Dyer <-> HEM discontinuity and the non-convergent band of `examples/example_03_flashing.md` in the coupled solver? | `scan_threshold_continuity_results.txt` |

Dependencies:
- `loocv_functional_forms.py` imports `loocv_kappa_correction.py` (same folder);
  both import `../explore_supercharge_correction.py` and
  `../waxman_2013_validation.py`.
- `scan_threshold_continuity.py` needs the optional `kappa_correction` argument of
  `full_system.py` and the geometries of `tests/generate_example_03.py`.
- The unit tests of the option are in `tests/test_kappa_correction_option.py`
  (these run in CI; the scripts here do not, they need CoolProp and take longer).

Status: exploratory. beta, supercharge_ref and S0 are fitted to one injector
geometry at 280-283 K; see `../waxman_2013_results.md`, Section 11.7 for what is
not established.
