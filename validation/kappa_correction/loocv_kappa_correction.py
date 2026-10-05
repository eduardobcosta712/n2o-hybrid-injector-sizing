"""
loocv_kappa_correction.py

Leave-one-curve-out cross-validation (LOOCV) of the supercharge-gated kappa
correction (dyer_mass_flow_corrected), on the 9 supercharge curves of Waxman
(2013) Fig. 13, injector 3.

*** EXPLORATORY -- does NOT change any model file. ***

QUESTION
    The 1.86 % global MAPE of the gated correction is IN-SAMPLE: beta and
    supercharge_ref were fitted on the same 64 two-phase points that the
    error is measured on. Does the correction survive out-of-sample?

PROTOCOL (for each of the 9 folds, k = one supercharge curve held out)
    1. Cd is re-calibrated on the single-phase windows of the 8 OTHER curves
       (no information from the held-out curve, not even through Cd).
    2. (beta, supercharge_ref) are fitted on the two-phase points of the 8
       other curves, by minimising their MAPE on a fixed grid, with that Cd.
    3. The held-out curve's two-phase points are predicted with the fitted
       (beta, ref) and the fold Cd. Nothing of them entered steps 1-2.
    Two fitting variants are run:
       FREE        any (beta, ref) on the grid
       GATED<=A    supercharge_ref restricted to <= Part A's own supercharge
                   (~95 psi), i.e. the variant that leaves the original four
                   Part A points untouched by construction.
    Baseline = original Dyer (beta = 0) with the same fold Cd.

    Part A (4 points, injector 2) never enters any fit; it is evaluated
    afterwards with every fold's parameters, as an independent check.

WHAT TO LOOK AT
    * pooled out-of-sample MAPE vs baseline (6.85 %) and vs in-sample (~1.9 %)
    * the 41 psi fold: the only curve close to saturation, hence the test of
      extrapolation towards the flashing threshold
    * the (beta, ref) chosen in each fold: if they jump around, the
      two-parameter form is under-determined by the data
    * the limit check at the end (does the corrected Dyer tend to HEM as the
      supercharge -> 0? That is what would make Dyer -> HEM continuous at
      the flashing threshold)

Usage (from anywhere, CoolProp installed):
    python validation/kappa_correction/loocv_kappa_correction.py
The output is also written to validation/kappa_correction/loocv_kappa_results.txt.
"""

import os
import sys

import numpy as np

# this file lives in <repo>/validation/kappa_correction/
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Search grid (fixed BEFORE looking at any out-of-sample result).
BETAS = np.round(np.arange(0.0, 3.0001, 0.1), 10)
REFS_PSI = np.arange(20.0, 400.0001, 5.0)


# ---------------------------------------------------------------------------
# Pure logic (no CoolProp needed -- operates on precomputed arrays)
# ---------------------------------------------------------------------------
def corrected_unit_flow(kappa, m_spi, m_hem, S_psi, beta, ref_psi):
    """
    Corrected Dyer flow for Cd = 1 (all flows scale linearly with Cd, and kappa
    does not depend on it). Same formula as dyer_mass_flow_corrected().
    """
    if beta > 0:
        factor = np.where(S_psi < ref_psi, (S_psi / ref_psi) ** beta, 1.0)
    else:
        factor = np.ones_like(S_psi)
    kc = kappa * factor
    return (kc / (1.0 + kc)) * m_spi + (1.0 / (1.0 + kc)) * m_hem


def build_grid(pts):
    """U[i_beta, i_ref, n_points]: unit-Cd corrected flows on the grid."""
    U = np.empty((len(BETAS), len(REFS_PSI), len(pts["m_exp"])))
    for i, b in enumerate(BETAS):
        for j, r in enumerate(REFS_PSI):
            U[i, j, :] = corrected_unit_flow(pts["kappa"], pts["m_spi"],
                                             pts["m_hem"], pts["S_psi"], b, r)
    return U


def _mape(Cd, U, m_exp, mask):
    err = np.abs(Cd * U[..., mask] - m_exp[mask]) / m_exp[mask]
    return 100.0 * err.mean(axis=-1)


def fit_params(U, pts, Cd, train_mask, ref_max_psi=None):
    """(beta, ref_psi, train MAPE %) minimising the training MAPE on the grid."""
    m = _mape(Cd, U, pts["m_exp"], train_mask)
    if ref_max_psi is not None:
        m = np.where(REFS_PSI[None, :] <= ref_max_psi, m, np.inf)
    i, j = np.unravel_index(np.argmin(m), m.shape)
    return float(BETAS[i]), float(REFS_PSI[j]), float(m[i, j])


def signed_errors(Cd, pts, mask, beta, ref_psi):
    m = Cd * corrected_unit_flow(pts["kappa"][mask], pts["m_spi"][mask],
                                 pts["m_hem"][mask], pts["S_psi"][mask],
                                 beta, ref_psi)
    return 100.0 * (m - pts["m_exp"][mask]) / pts["m_exp"][mask]


def run_loocv(pts, cd_samples, part_a_limit_psi):
    """
    pts        : dict of numpy arrays (kappa, m_spi, m_hem, S_psi, m_exp, sup)
    cd_samples : {super_psi: [implied Cd values of the single-phase window]}
    Returns (folds, cd_all, pooled) -- see main() for how they are printed.
    """
    U = build_grid(pts)
    sups = sorted(cd_samples)
    folds = []
    pooled = {"base": [], "free": [], "gated": [], "sup": []}

    for k in sups:
        test = pts["sup"] == k
        if not test.any():
            continue
        train = ~test
        cal = [v for s in sups if s != k for v in cd_samples[s]]
        cd_f = sum(cal) / len(cal)

        b_f, r_f, mt_f = fit_params(U, pts, cd_f, train)
        b_g, r_g, mt_g = fit_params(U, pts, cd_f, train, ref_max_psi=part_a_limit_psi)

        e_base = signed_errors(cd_f, pts, test, 0.0, 200.0)
        e_free = signed_errors(cd_f, pts, test, b_f, r_f)
        e_gate = signed_errors(cd_f, pts, test, b_g, r_g)
        folds.append({
            "sup": k, "n": int(test.sum()), "cd": cd_f,
            "free": (b_f, r_f), "gated": (b_g, r_g),
            "train_free": mt_f, "train_gated": mt_g,
            "base": e_base, "e_free": e_free, "e_gate": e_gate})
        pooled["base"].extend(e_base)
        pooled["free"].extend(e_free)
        pooled["gated"].extend(e_gate)
        pooled["sup"].extend([k] * int(test.sum()))

    all_cd = [v for s in sups for v in cd_samples[s]]
    cd_all = sum(all_cd) / len(all_cd)
    in_sample = {
        "free": fit_params(U, pts, cd_all, np.ones(len(pts["m_exp"]), bool)),
        "gated": fit_params(U, pts, cd_all, np.ones(len(pts["m_exp"]), bool),
                            ref_max_psi=part_a_limit_psi),
    }
    return folds, cd_all, pooled, in_sample


def stats(errs):
    e = np.asarray(errs, dtype=float)
    if e.size == 0:
        return "n=  0"
    return (f"n={e.size:3d}  mean={e.mean():+6.2f}%  MAPE={np.abs(e).mean():5.2f}%  "
            f"max|err|={np.abs(e).max():5.1f}%")


# ---------------------------------------------------------------------------
# Project-dependent part (needs CoolProp)
# ---------------------------------------------------------------------------
def _load_project():
    sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "model"))
    sys.path.insert(0, _REPO_ROOT)
    sys.path.insert(0, os.path.join(_REPO_ROOT, "validation"))
    global n2o, inj, wv, ex
    import n2o_properties as n2o
    import injector_two_phase as inj
    import validation.waxman_2013_validation as wv
    import explore_supercharge_correction as ex
    return n2o, inj, wv, ex


def build_points(curves):
    """Two-phase points (same selection as the exploration script), with the
    Cd-independent pieces of the corrected Dyer flow computed at Cd = 1."""
    pts_list = ex.two_phase_points(curves)
    cols = {k: [] for k in ("kappa", "m_spi", "m_hem", "S_psi", "m_exp", "sup")}
    for p in pts_list:
        T_d = n2o.T_sat(p["P2"])
        r = inj.dyer_mass_flow_corrected(
            1.0, wv.A_INJ3, p["T"], p["P1"], p["P2"], p["rho_l_up"],
            n2o.rho_liquid_sat(T_d), n2o.M_N2O / n2o.nu_vapor_sat(T_d),
            beta=0.0, supercharge_ref_Pa=1.0)
        cols["kappa"].append(r["kappa"])
        cols["m_spi"].append(r["m_dot_SPI"])
        cols["m_hem"].append(r["m_dot_HEM"])
        cols["S_psi"].append((p["P1"] - n2o.P_sat(p["T"])) / wv.PSI_PA)
        cols["m_exp"].append(p["m_exp"])
        cols["sup"].append(p["super_psi"])
        p["_T_d"] = T_d
    return {k: np.array(v, dtype=float) for k, v in cols.items()}, pts_list


def self_check(pts, pts_list):
    """The vectorised formula must reproduce dyer_mass_flow_corrected()."""
    for idx in (0, len(pts_list) // 2, len(pts_list) - 1):
        p = pts_list[idx]
        for beta, ref in ((1.0, 200.0), (2.5, 80.0), (0.0, 200.0)):
            r = inj.dyer_mass_flow_corrected(
                1.0, wv.A_INJ3, p["T"], p["P1"], p["P2"], p["rho_l_up"],
                n2o.rho_liquid_sat(p["_T_d"]),
                n2o.M_N2O / n2o.nu_vapor_sat(p["_T_d"]),
                beta=beta, supercharge_ref_Pa=ref * wv.PSI_PA)
            mine = corrected_unit_flow(pts["kappa"][idx:idx + 1], pts["m_spi"][idx:idx + 1],
                                       pts["m_hem"][idx:idx + 1], pts["S_psi"][idx:idx + 1],
                                       beta, ref)[0]
            if abs(mine - r["m_dot_Dyer_corrected"]) > 1e-9 * abs(mine):
                raise RuntimeError(
                    f"self-check failed at point {idx}, beta={beta}, ref={ref}: "
                    f"{mine} vs {r['m_dot_Dyer_corrected']}")


def part_a_stats(beta, ref_psi, cd):
    rows = ex.run_part_a_corrected(beta, ref_psi * wv.PSI_PA, cd=cd)
    return wv.summarize([r["err"] for r in rows])


def limit_check(beta, ref_psi, cd, out):
    """Does the corrected Dyer tend to HEM as supercharge -> 0?  (T = 280 K,
    P_down = 30 bar, as in the Waxman conditions.)"""
    T, P2 = 280.0, 30e5
    T_d = n2o.T_sat(P2)
    rl_d, rv_d = n2o.rho_liquid_sat(T_d), n2o.M_N2O / n2o.nu_vapor_sat(T_d)
    rl_u = n2o.rho_liquid_sat(T)
    out(f"  beta = {beta}, ref = {ref_psi:.0f} psi, Cd = {cd:.3f}")
    out(f"  {'S [psi]':>8} {'Dyer':>9} {'corrected':>10} {'HEM':>9} "
        f"{'corr/HEM':>9} {'corr/Dyer':>10}")
    for S in (0.05, 0.5, 2.0, 5.0, 10.0, 20.0, 41.0, 95.0, 200.0):
        Pu = n2o.P_sat(T) + S * wv.PSI_PA
        r = inj.dyer_mass_flow_corrected(cd, wv.A_INJ3, T, Pu, P2, rl_u, rl_d, rv_d,
                                         beta=beta, supercharge_ref_Pa=ref_psi * wv.PSI_PA)
        d0 = r["m_dot_SPI"] * r["kappa"] / (1 + r["kappa"]) + r["m_dot_HEM"] / (1 + r["kappa"])
        out(f"  {S:8.2f} {d0*1000:9.2f} {r['m_dot_Dyer_corrected']*1000:10.2f} "
            f"{r['m_dot_HEM']*1000:9.2f} {r['m_dot_Dyer_corrected']/r['m_dot_HEM']:9.3f} "
            f"{r['m_dot_Dyer_corrected']/d0:10.3f}")


def main():
    _load_project()
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    bar = "=" * 78
    curves, _ = wv.load_multiseries("waxman_fig13_mdot_vs_dP_by_supercharge.csv")
    cd_samples = {c["super_psi"]: wv.spi_cd_samples(c) for c in curves}
    pts, pts_list = build_points(curves)
    self_check(pts, pts_list)

    sup_A_psi = (ex.P1_A - n2o.P_sat(ex.T1_A)) / wv.PSI_PA
    folds, cd_all, pooled, ins = run_loocv(pts, cd_samples, sup_A_psi)

    out(bar)
    out("LOOCV of the gated kappa correction -- Waxman Fig. 13, injector 3")
    out(bar)
    out(f"Two-phase points: {len(pts['m_exp'])} on {len(folds)} curves | "
        f"pooled Cd (all curves) = {cd_all:.4f} | grid: {len(BETAS)} betas x "
        f"{len(REFS_PSI)} refs")
    out(f"Part A supercharge = {sup_A_psi:.1f} psi -> GATED<=A restricts ref to <= that")
    out("Vectorised formula verified against dyer_mass_flow_corrected(): OK")
    out()

    out("IN-SAMPLE fit on all 64 points (for comparison; this is what 1.86 % was):")
    for name in ("free", "gated"):
        b, r, m = ins[name]
        out(f"  {name:6s}: beta={b:.1f} ref={r:.0f} psi -> in-sample MAPE {m:.2f} %")
    out()

    out("PER FOLD (held-out curve; parameters fitted on the other 8 curves)")
    out(f"{'held-out':>9} {'n':>3} {'Cd_fold':>8} | {'baseline':>9} | "
        f"{'FREE b/ref':>11} {'train':>6} {'OOS MAPE':>9} | "
        f"{'GATED b/ref':>11} {'train':>6} {'OOS MAPE':>9}")
    for f in folds:
        out(f"{f['sup']:7.0f}psi {f['n']:3d} {f['cd']:8.4f} | "
            f"{np.abs(f['base']).mean():8.2f}% | "
            f"{f['free'][0]:4.1f}/{f['free'][1]:5.0f} {f['train_free']:5.2f}% "
            f"{np.abs(f['e_free']).mean():8.2f}% | "
            f"{f['gated'][0]:4.1f}/{f['gated'][1]:5.0f} {f['train_gated']:5.2f}% "
            f"{np.abs(f['e_gate']).mean():8.2f}%")
    out()

    out("POOLED OUT-OF-SAMPLE (every point predicted by a model that never saw its curve)")
    out(f"  baseline (Dyer)     : {stats(pooled['base'])}")
    out(f"  FREE                : {stats(pooled['free'])}")
    out(f"  GATED<=A            : {stats(pooled['gated'])}")
    out()

    sup = np.array(pooled["sup"])
    out("OUT-OF-SAMPLE by supercharge band:")
    for lo, hi, name in ((0, 100, "< 100 psi"), (100, 200, "100-200 psi"), (200, 1000, ">= 200 psi")):
        m = (sup >= lo) & (sup < hi)
        for key, lab in (("base", "baseline"), ("free", "FREE    "), ("gated", "GATED<=A")):
            out(f"  {name:<12} {lab}: {stats(np.array(pooled[key])[m])}")
    out()

    f41 = next((f for f in folds if abs(f["sup"] - 41) < 1), None)
    if f41:
        out("THE 41 psi FOLD (closest to saturation = extrapolation towards the threshold):")
        out(f"  baseline : {stats(f41['base'])}")
        out(f"  FREE     : {stats(f41['e_free'])}   (beta={f41['free'][0]}, ref={f41['free'][1]:.0f})")
        out(f"  GATED<=A : {stats(f41['e_gate'])}   (beta={f41['gated'][0]}, ref={f41['gated'][1]:.0f})")
        out()

    bs = [f["free"][0] for f in folds]; rs = [f["free"][1] for f in folds]
    out("STABILITY of the FREE fit across folds:")
    out(f"  beta: min {min(bs):.1f}  max {max(bs):.1f} | ref: min {min(rs):.0f}  max {max(rs):.0f} psi")
    out(f"  (in-sample optimum: beta={ins['free'][0]:.1f}, ref={ins['free'][1]:.0f})")
    out()

    out("PART A (4 points, injector 2; never used in any fit)")
    out(f"  {'parameters':<34} {'Cd=0.65':>30} | {'Cd=0.681 (measured)':>30}")
    sets = [("baseline (beta=0)", 0.0, 200.0)]
    sets.append((f"in-sample FREE ({ins['free'][0]:.1f}/{ins['free'][1]:.0f})", *ins["free"][:2]))
    sets.append((f"in-sample GATED<=A ({ins['gated'][0]:.1f}/{ins['gated'][1]:.0f})", *ins["gated"][:2]))
    for name, b, r in sets:
        a1 = part_a_stats(b, r, ex.CD_PART_A)
        a2 = part_a_stats(b, r, ex.CD_PART_A_FIG15)
        out(f"  {name:<34} mean={a1['mean']:+6.2f}% MAPE={a1['mape']:5.2f}% | "
            f"mean={a2['mean']:+6.2f}% MAPE={a2['mape']:5.2f}%")
    mp = [part_a_stats(*f["free"], ex.CD_PART_A)["mape"] for f in folds]
    mp2 = [part_a_stats(*f["free"], ex.CD_PART_A_FIG15)["mape"] for f in folds]
    out(f"  FREE params of the 9 folds: Part A MAPE range {min(mp):.2f}-{max(mp):.2f} % "
        f"(Cd 0.65) | {min(mp2):.2f}-{max(mp2):.2f} % (Cd 0.681)")
    out()

    out("LIMIT CHECK: corrected Dyer vs HEM as supercharge -> 0 (T=280 K, P_down=30 bar)")
    out("  (corr/HEM -> 1 means Dyer joins HEM at the flashing threshold; the HEM")
    out("   two-phase-inlet path at x_inlet -> 0 is exactly this HEM flow.)")
    limit_check(*ins["free"][:2], cd_all, out)
    if ins["gated"][:2] != ins["free"][:2]:
        out()
        limit_check(*ins["gated"][:2], cd_all, out)

    path = os.path.join(_REPO_ROOT, "validation", "kappa_correction", "loocv_kappa_results.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n(output also written to {path})")


if __name__ == "__main__":
    main()
