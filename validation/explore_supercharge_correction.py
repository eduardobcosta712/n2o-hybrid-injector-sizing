"""
explore_supercharge_correction.py

*** EXPLORATORY -- does NOT change injector_two_phase.py or full_system.py ***

Tests the hypothesis that a supercharge-dependent correction of Dyer's
parameter kappa,

    kappa' = kappa * (supercharge / supercharge_ref) ** beta

reduces the MAPE observed at low supercharge (docs/future_work.md,
Priority 1, item 1) WITHOUT degrading the already-validated band
(dP = 8-14 bar, Part A).

HOW TO RUN (in your local repo, with CoolProp installed):
  python validation/explore_supercharge_correction.py
  (runs from any directory, like waxman_2013_validation.py). It needs
  dyer_mass_flow_corrected(), which now lives in
  src/model/injector_two_phase.py right after dyer_mass_flow().

It directly reuses the loading and calibration functions of
validation/waxman_2013_validation.py (load_multiseries, curve_state,
spi_cd_samples) -- there is no duplicated CSV-parsing logic.

Methodology (mirrors run_part_b() of waxman_2013_validation.py):
  1. Cd is calibrated ONCE, pooled over the single-phase window of each
     curve (30 psi <= dP <= supercharge) -- kappa does not enter there, so
     the Cd calibration is not affected by the correction.
  2. For each (beta, supercharge_ref) on a grid, the corrected mass flow is
     recomputed at ALL two-phase points (dP > supercharge) and the global
     MAPE and the MAPE per supercharge band are measured.
  3. Hard constraint: the "sacred" band (original Part A, the 4
     Nino & Razavi points, dP 8-14 bar, injector-2 geometry, Cd = 0.65)
     must not degrade significantly -- IMPORTANT: always use
     CD_PART_A = 0.65, NOT the pooled Cd of Fig. 13 (injector 3) -- they
     are different geometries.

Note on small numerical differences against waxman_2013_results.md: this
script evaluates the injector model directly (no upstream line), whereas
validation/waxman_2013_validation.py goes through the coupled solver with a
5 cm line. That gives a baseline Part A MAPE of 2.75 % here vs 2.76 % there.
Likewise CD_PART_A_FIG15 is the ROUNDED value 0.681 (unrounded digitised
mean: 0.6806); with 0.681 the no-correction Part A MAPE is 4.47 %, with
0.6806 it is 4.41 % (Section 4 of the results report).
"""

import math
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "model"))

from n2o_properties import P_sat, T_sat, rho_liquid_sat, nu_vapor_sat, M_N2O
from injector_two_phase import dyer_mass_flow, dyer_mass_flow_corrected

# Reuse the loader and calibration already written and tested in
# waxman_2013_validation.py -- avoids duplicating parsing logic.
sys.path.insert(0, _REPO_ROOT)
from validation.waxman_2013_validation import (
    load_multiseries, curve_state, spi_cd_samples, A_INJ3, MIN_DP_PSI,
    PSI_PA, BAR_PER_PSI, pct_error, summarize,
)

# --- Original Part A (Nino & Razavi 2019, injector-2 geometry) -------------
# DIFFERENT geometry from injector 3 (Fig. 13) -- own Cd, not the pooled
# Cd calibrated below.
T1_A = 280.0
P1_A = 4.36e6
D_A = 0.0015
A_A = math.pi * (D_A / 2.0) ** 2
CD_PART_A = 0.65   # generic value already assumed by the repo for injector 2
CD_PART_A_FIG15 = 0.681   # MEASURED Cd for injector 2 (Waxman Fig. 15,
                          # digitised; rounded from 0.6806). Already known,
                          # see waxman_2013_validation.py print_part_d():
                          # merely swapping 0.65 -> 0.681 raises the
                          # baseline MAPE from 2.76 % to ~4.4 %, with no
                          # kappa correction involved. Used below to
                          # separate how much of the degradation seen
                          # under the "aggressive" correction is really a
                          # kappa effect, and how much is just the assumed
                          # Cd no longer being the right one.
CASES_A = [
    ("Pre-critical",  0.84, 44.0),
    ("Critical",      0.98, 46.5),
    ("Post-critical1", 1.09, 47.5),
    ("Post-critical2", 1.37, 48.0),
]


def _rho_v(T):
    return M_N2O / nu_vapor_sat(T)


def two_phase_points(curves):
    """All two-phase points (dP > supercharge) of injector 3, with the
    upstream state -- same selection as dy_rows in run_part_b()."""
    pts = []
    for c in curves:
        T, P1_pa = curve_state(c)
        rho_l_up = rho_liquid_sat(T)
        for dP, m_exp in zip(c["dP_psi"], c["y"]):
            if dP <= c["super_psi"] or dP < MIN_DP_PSI:
                continue
            P2 = P1_pa - dP * PSI_PA
            if P2 <= 1e4:
                continue
            pts.append({"super_psi": c["super_psi"], "dP_psi": dP,
                        "dP_bar": dP * BAR_PER_PSI,
                        "m_exp": m_exp, "T": T, "P1": P1_pa, "P2": P2,
                        "rho_l_up": rho_l_up})
    return pts


def predict_baseline(pts, cd):
    """Original Dyer (beta = 0), no upstream line -- negligible, as
    assumed in the original script (LINE ~= 5 cm of 25.4 mm ID)."""
    out = []
    for p in pts:
        T_d = T_sat(p["P2"])
        rho_l_d = rho_liquid_sat(T_d)
        rho_v_d = _rho_v(T_d)
        r = dyer_mass_flow(cd, A_INJ3, p["T"], p["P1"], p["P2"],
                           p["rho_l_up"], rho_l_d, rho_v_d)
        out.append({**p, "m_model": r["m_dot_Dyer"],
                    "err": pct_error(r["m_dot_Dyer"], p["m_exp"])})
    return out


def predict_corrected(pts, cd, beta, supercharge_ref_pa):
    out = []
    for p in pts:
        T_d = T_sat(p["P2"])
        rho_l_d = rho_liquid_sat(T_d)
        rho_v_d = _rho_v(T_d)
        r = dyer_mass_flow_corrected(cd, A_INJ3, p["T"], p["P1"], p["P2"],
                                     p["rho_l_up"], rho_l_d, rho_v_d,
                                     beta=beta, supercharge_ref_Pa=supercharge_ref_pa)
        m = r["m_dot_Dyer_corrected"]
        out.append({**p, "m_model": m, "err": pct_error(m, p["m_exp"]),
                    "kappa": r["kappa"], "kappa_corr": r["kappa_corrected"]})
    return out


def fmt(s):
    return f"n={s['n']:3d}  mean={s['mean']:+6.2f}%  MAPE={s['mape']:5.2f}%  max|err|={s['max']:5.1f}%"


def run_part_a_corrected(beta, supercharge_ref_pa, cd=None):
    """
    Checks the sacred band (8-14 bar, Nino & Razavi) with the correction
    applied.

    cd : float or None
        If None (default), uses CD_PART_A (0.65 -- the generic value
        already assumed by the repo). Pass CD_PART_A_FIG15 (0.681, the
        MEASURED value for this geometry) to test whether the degradation
        seen with the "aggressive" correction is a genuine effect of the
        kappa correction, or (partly/entirely) an artefact of the assumed
        Cd already being wrong before any correction comes into play.
    """
    if cd is None:
        cd = CD_PART_A
    rows = []
    rho_l_up = rho_liquid_sat(T1_A)
    for label, dP_MPa, m_exp_gs in CASES_A:
        P2 = P1_A - dP_MPa * 1e6
        T_d = T_sat(P2)
        rho_l_d = rho_liquid_sat(T_d)
        rho_v_d = _rho_v(T_d)
        r = dyer_mass_flow_corrected(cd, A_A, T1_A, P1_A, P2,
                                     rho_l_up, rho_l_d, rho_v_d,
                                     beta=beta, supercharge_ref_Pa=supercharge_ref_pa)
        m = r["m_dot_Dyer_corrected"]
        m_exp_kgs = m_exp_gs * 1e-3
        rows.append({"label": label, "err": pct_error(m, m_exp_kgs),
                     "m_model": m, "m_exp": m_exp_kgs})
    return rows


def main():
    print("=" * 78)
    print("EXPLORATION: GATED supercharge-dependent correction of kappa")
    print("(real CoolProp -- check that 'import CoolProp' works first)")
    print("=" * 78)

    curves, _n_dropped = load_multiseries("waxman_fig13_mdot_vs_dP_by_supercharge.csv")
    samples = [v for c in curves for v in spi_cd_samples(c)]
    cd_pooled = sum(samples) / len(samples)
    print(f"\nPooled Cd (single-phase window, injector 3, n={len(samples)}): {cd_pooled:.4f}")

    pts = two_phase_points(curves)
    print(f"Two-phase points (dP > supercharge): {len(pts)}")

    supercharge_A_psi = (P1_A - P_sat(T1_A)) / PSI_PA
    print(f"\nNOTE: the supercharge of the 4 Part A points is "
          f"{supercharge_A_psi:.1f} psi -- ANY supercharge_ref >= "
          f"{supercharge_A_psi:.0f} psi makes the correction act on "
          f"Part A too (being 'gated' alone does not protect it).")

    # --- Baseline (beta = 0, original Dyer) ---
    base = predict_baseline(pts, cd_pooled)
    print(f"\nBASELINE (original Dyer) -- must match waxman_2013_results.md:")
    print("  Global:", fmt(summarize([r["err"] for r in base])))
    for lo, hi, name in [(0, 100, "supercharge < 100 psi"),
                         (100, 200, "100 <= supercharge < 200 psi"),
                         (200, 1000, "supercharge >= 200 psi")]:
        b = [r["err"] for r in base if lo <= r["super_psi"] < hi]
        print(f"  {name:<32}", fmt(summarize(b)))

    part_a_base = run_part_a_corrected(beta=0.0, supercharge_ref_pa=200 * PSI_PA)
    print("  Part A (4 points, 8-14 bar, should give ~2.75 %):",
          fmt(summarize([r["err"] for r in part_a_base])))

    # --- Grid search (GATED: factor = 1 outside the corrected region) ---
    print(f"\n{'='*78}\nGRID SEARCH (GATED): beta x supercharge_ref\n{'='*78}")
    betas = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]
    # refs_psi includes values BELOW and ABOVE the supercharge of Part A
    # (~95 psi), to make the trade-off visible in the report.
    refs_psi = [50, 60, 70, 80, 100, 150, 200, 250, 300]

    results = []
    for ref_psi in refs_psi:
        ref_pa = ref_psi * PSI_PA
        for beta in betas:
            corr = predict_corrected(pts, cd_pooled, beta, ref_pa)
            s_global = summarize([r["err"] for r in corr])
            part_a = run_part_a_corrected(beta, ref_pa)
            s_a = summarize([r["err"] for r in part_a])
            results.append({"beta": beta, "ref_psi": ref_psi,
                            "mape_global": s_global["mape"], "mean_global": s_global["mean"],
                            "mape_a": s_a["mape"], "mean_a": s_a["mean"],
                            "protects_A": ref_psi <= supercharge_A_psi})

    print(f"\n{'ref[psi]':>8} {'beta':>5} {'global MAPE':>12} {'global mean':>12} "
          f"{'PartA MAPE':>11} {'PartA mean':>11}  {'Protects A?':>11}")
    for r in results:
        flag = "YES (gate)" if r["protects_A"] else ("degrades" if r["mape_a"] > 3.6 else "-")
        print(f"{r['ref_psi']:8.0f} {r['beta']:5.1f} {r['mape_global']:12.2f} "
              f"{r['mean_global']:+12.2f} {r['mape_a']:11.2f} {r['mean_a']:+11.2f}  {flag:>11}")

    # --- Two readings, to make the trade-off explicit ---
    print(f"\n{'='*78}")
    protected = [r for r in results if r["protects_A"]]
    if protected:
        best_protected = min(protected, key=lambda r: r["mape_global"])
        print("CONSERVATIVE OPTION -- Part A protected BY CONSTRUCTION "
              f"(supercharge_ref <= {supercharge_A_psi:.0f} psi):")
        print(f"  beta = {best_protected['beta']}, "
              f"supercharge_ref = {best_protected['ref_psi']} psi")
        print(f"  Global MAPE: {best_protected['mape_global']:.2f}% "
              f"(baseline: {summarize([r['err'] for r in base])['mape']:.2f}%)")
        print(f"  Part A MAPE: {best_protected['mape_a']:.2f}% "
              f"(identical to the baseline, by construction)")

    print()
    unrestricted = [r for r in results if r["mape_a"] <= 5.0]
    if unrestricted:
        best_unrestricted = min(unrestricted, key=lambda r: r["mape_global"])
        print("AGGRESSIVE OPTION -- only requires Part A MAPE <= 5% (Part A "
              "may degrade a little):")
        print(f"  beta = {best_unrestricted['beta']}, "
              f"supercharge_ref = {best_unrestricted['ref_psi']} psi")
        print(f"  Global MAPE: {best_unrestricted['mape_global']:.2f}% "
              f"(baseline: {summarize([r['err'] for r in base])['mape']:.2f}%)")
        print(f"  Part A MAPE: {best_unrestricted['mape_a']:.2f}% "
              f"(baseline: {summarize([r['err'] for r in part_a_base])['mape']:.2f}%)")

        best_corr = predict_corrected(pts, cd_pooled, best_unrestricted["beta"],
                                       best_unrestricted["ref_psi"] * PSI_PA)
        print(f"\n  Breakdown by supercharge (aggressive option):")
        for lo, hi, name in [(0, 100, "< 100 psi"), (100, 200, "100-200 psi"),
                             (200, 1000, ">= 200 psi")]:
            b = [r["err"] for r in best_corr if lo <= r["super_psi"] < hi]
            print(f"    {name:<15}", fmt(summarize(b)))

        # -------------------------------------------------------------
        # TEST: how much of the Part A degradation is a kappa effect vs.
        # an effect of the assumed Cd (0.65) already not being the
        # correct one to begin with? CD_PART_A_FIG15 (0.681) is the value
        # MEASURED directly for this injector (Waxman Fig. 15,
        # digitised) -- it was already known, before this correction,
        # that merely swapping the Cd raised the baseline MAPE from
        # 2.76 % to ~4.4 % (see waxman_2013_validation.py
        # print_part_d() and the sensitivity done there).
        # -------------------------------------------------------------
        print(f"\n{'='*78}")
        print("SEPARATING THE EFFECT OF Cd FROM THE EFFECT OF KAPPA on Part A")
        print(f"{'='*78}")

        beta_best = best_unrestricted["beta"]
        ref_best_pa = best_unrestricted["ref_psi"] * PSI_PA

        combos = [
            ("Cd=0.65 (assumed), beta=0 (no correction)", CD_PART_A, 0.0, 200 * PSI_PA),
            ("Cd=0.65 (assumed), beta=best (corrected)", CD_PART_A, beta_best, ref_best_pa),
            ("Cd=0.681 (measured, Fig.15), beta=0 (no correction)", CD_PART_A_FIG15, 0.0, 200 * PSI_PA),
            ("Cd=0.681 (measured, Fig.15), beta=best (corrected)", CD_PART_A_FIG15, beta_best, ref_best_pa),
        ]
        for label, cd_test, beta_test, ref_test in combos:
            rows = run_part_a_corrected(beta_test, ref_test, cd=cd_test)
            s = summarize([r["err"] for r in rows])
            print(f"  {label:<52}", fmt(s))

        print(f"\n  Reading: compare line 1 with line 3 (effect of swapping only the")
        print(f"  Cd, WITHOUT the kappa correction) against line 2 vs line 4 (the same")
        print(f"  effect, BUT with the correction already applied). If the 3-1")
        print(f"  difference is similar to the 4-2 difference, the Cd explains most")
        print(f"  of the MAPE 'increase' we attributed to the kappa correction -- not")
        print(f"  the correction itself.")


if __name__ == "__main__":
    main()
