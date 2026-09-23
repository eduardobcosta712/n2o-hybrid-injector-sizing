# Future Work

Extensions identified during this project's development, ordered by
technical importance per the project roadmap document.

---

## Implemented (removed from future work)

(unchanged from the previous version of this file — see git history / the
September 2026 audit chat for the full list; only Priority 4's status and
the "Audit follow-ups" section changed in this update.)

---

## Priority 1 — Non-equilibrium critical flow ceiling for Dyer ✅ RESOLVED (September 2026)

(unchanged — see `injector_two_phase.henry_fauske_critical_flow()` and
`validation/waxman_2013_results.md`, Section 7, now regenerated with the
confirmed CoolProp backend: 52.15 g/s at Waxman conditions.)

## Priority 2 — Full-system experimental validation

**Current state.** The Dyer injector model is validated against Waxman
(2013/2014) with mean error −0.3% at pressure drops of 8-14 bar (MAPE
2.76%, CoolProp-confirmed, September 2026), using a line with negligible
losses. Three parts of the model have **no experimental validation at
all** — unchanged from before:

1. **Pressure drops of 20-50 bar** (the range of typical motor designs).
2. **The coupled feed-line / injector solver** with a line that matters.
3. **The two-phase inlet path** (flashing in the line): implemented and
   unit tested only. The switch from Dyer to HEM at the flashing
   threshold is discontinuous.

**New finding (September 2026, CoolProp confirmation run).** Item 3's
discontinuity is not only a "different model, different number" issue —
it can make the coupled fixed-point solver (`full_system.evaluate_full_system`)
**fail to converge entirely** for an operating point close enough to the
flashing threshold. Confirmed directly: with the real CoolProp P_sat
(lower than the old Perry/McGill correlation predicted at 22 °C), the
tank/feed-line geometry used in `examples/example_03_flashing.md`
(54 bar, 22 °C, needle valve, 8 mm line) now sits inside a narrow
pressure band (roughly 54.0–55.5 bar at this geometry) where the damped
fixed-point iteration oscillates between the Dyer (liquid-inlet, ~290 g/s)
and HEM (two-phase-inlet, ~240 g/s) branches indefinitely, even at very
small damping (α down to 0.05, 200 iterations) — it never satisfies the
convergence tolerance and `RuntimeError`s out, correctly (the solver is
not silently returning a wrong answer). Just outside that band it
converges cleanly on either side (53.5 bar → 193.9 g/s, HEM; 56 bar →
344.4 g/s, Dyer). This is a genuine, newly-exposed instance of the
already-documented Dyer→HEM discontinuity, not a bug in the CoolProp
switch itself — but it means `example_03_flashing.md`'s specific
54 bar starting point can no longer be reproduced with the confirmed
property backend and needs a deliberate fix: either move the example's
starting pressure to a point with clearer margin on one side, or give
the coupled solver a dedicated strategy for this borderline regime
(e.g. detect oscillation and report a range instead of failing, or
switch to a slower/bisection-based update near the threshold). Left open
pending Eduardo's decision — see "Audit follow-ups" below.

**What is needed.** Experimental data covering the full path
(tank → line → injector) with known geometry, discharge coefficient,
and measured mass flow, ideally including runs at 20-50 bar drop and
runs with and without flashing in the line.

---

## Priority 3 — OF ratio and fuel grain sizing ✅ RESOLVED (September 2026)

(unchanged)

## Priority 3b — Non-circular grain port shapes (star, wagon-wheel)

(unchanged)

## Priority 4 — Improved N₂O thermophysical properties ✅ RESOLVED AND CONFIRMED (September 2026)

**What was done.** `n2o_properties.py`'s thermodynamic functions now call
**CoolProp** (Bell et al., 2014), implementing the Lemmon & Span (2006)
equation of state for N₂O, replacing the closed-form Perry/McGill
correlations and McGill Table A.1.

**Confirmed with the real package (September 2026).** `pip install CoolProp`
succeeded on Eduardo's machine (CoolProp 8.0.0) and was independently
reproduced by Claude in a second sandboxed environment (also 8.0.0,
network access to PyPI). Both runs agree to the displayed precision. The
full test suite (262 tests) passes, `python src/model/n2o_properties.py`'s
self-check matches the literature reference values it targets, and
`validation/waxman_2013_validation.py` has been re-run end to end — see
`validation/waxman_2013_results.md` for the regenerated report.

**One test needed a genuine fix, not just a re-run.**
`tests/test_n2o_properties.py::TestPsat::test_agrees_with_perry_correlation`
failed at `rel_tol=0.04`: the real CoolProp P_sat differs from the Perry
correlation by up to **4.81 %** at T = 230 K (far from the critical
point), decreasing smoothly to about 0.9 % near T_MAX. This is the
opposite trend from what the module's own docstring claimed ("Perry/McGill
carried ~2–3 % error near the critical point") — the real data shows the
Perry correlation is *better* near the critical point and *worse* well
below it. The test's tolerance was widened to `rel_tol=0.06` (with a
comment recording this finding) rather than silently loosened without
explanation; the docstring in `n2o_properties.py` should be corrected to
match this the next time that module is touched (not yet done in this
pass — see "Audit follow-ups").

**Numbers that moved.** Every figure downstream of P_sat/h_fg shifted by
roughly 1–5 %, in the direction the known Perry/McGill error would
predict. Regenerated so far: `validation/waxman_2013_results.md` (full
rewrite), `examples/example_01_sizing.md`, `examples/example_02_design.md`.
`examples/example_03_flashing.md` could **not** be regenerated as-is — see
the Priority 2 finding above. `docs/04_implementation.md`'s per-module
numeric callouts (§4.1's literature self-check numbers, §4.4's isolated
Dyer/HEM validation figures) have not yet been re-verified line by line
in this pass.

## Priority 5 — Quantitative uncertainty analysis (Monte Carlo)

(unchanged)

## Priority 6 — Improved exporting of the results

(unchanged)

## Priority 7 — Tank thermal model

(unchanged)

## Priority 8 — Transient / blowdown model

(unchanged)

---

## Audit follow-ups (September 2026) — open items needing the author

Carried over from the previous audit, plus new items from the CoolProp
confirmation pass:

1. **NIST viscosity attribution** -- still unresolved, unaffected by this pass.
2. **Table A.1/A.2 data** -- still unresolved, unaffected by this pass.
3. **Waxman CSV vs. validation conditions** -- still unresolved (P1 = 4.96
   MPa in the CSV vs. 4.36 MPa in the validation script); unaffected by
   this pass.
4. **`plotting.py` micro-cleanups** -- still pending, unaffected by this pass.
5. **Dyer → HEM discontinuity** at the flashing threshold — **now shown to
   cause solver non-convergence, not just a value jump** (see Priority 2
   above). Higher priority than before.
6. **`apply_choking_limit()`** -- still pending deletion.
7. **Line pressure reaching zero** -- still pending.
8. **NEW: `n2o_properties.py` module docstring is now factually wrong**
   about where the Perry/McGill correlation's error is largest (claims
   "near the critical point"; the confirmed CoolProp run shows the
   opposite — largest at low T, smallest near critical). Needs a one-line
   correction next time that file is edited.
9. **NEW: README's test-count / module list does not match the repository
   as currently provided.** The README (before this update) claimed 7 test
   modules including `test_audit_regressions.py` (15 tests, 257 total).
   The repository as supplied for this audit has 6 test modules and 262
   tests total, with no `test_audit_regressions.py` file. This update's
   README reflects the 6-module, 262-test reality and flags the
   discrepancy rather than guessing at it — Eduardo should confirm whether
   that 7th file was merged into the others, deleted, or simply not
   included in what was shared for this audit.
10. **NEW: `.github/workflows/tests.yml` does not actually run
    `validation/waxman_2013_validation.py` or grep for a MAPE value.** The
    README and the project's own memory notes both describe a CI step that
    "runs the Waxman validation script and fails if the validated MAPE
    changes" — no such step is present in the workflow file as supplied
    (it only runs `pytest tests/ -v`). Either the workflow file shared for
    this audit is out of date relative to the real repository, or that CI
    step was removed at some point and the documentation was never
    updated. Needs Eduardo to check the actual `.github/workflows/tests.yml`
    on GitHub and reconcile one way or the other — nothing was invented or
    silently added to the workflow file in this pass.
11. **NEW (Priority 2 above): `examples/example_03_flashing.md` needs a new
    starting scenario or a solver-side fix** before it can be regenerated
    with the confirmed CoolProp backend.
