"""
loocv_functional_forms.py

Leave-one-curve-out comparison of candidate functional forms for a
supercharge-dependent correction of Dyer's kappa.  Follow-up of
loocv_kappa_correction.py (results: validation/kappa_correction/loocv_kappa_results.txt).

*** EXPLORATORY -- does NOT change any model file. ***

WHY
    The first LOOCV showed the correction generalises (pooled out-of-sample
    MAPE 6.88 % -> 1.36 %), but the optimum sat at the edge of the grid
    (supercharge_ref = 400 psi, data only reach 371 psi): the data do not
    identify a gate, they identify an (almost) linear rescaling
    kappa' ~ kappa * S / S0.  Here we ask which form is the simplest one
    that does the job, and whether the Cd interacts with it.

CANDIDATES  (S = supercharge, kappa0 = original Dyer kappa)
    F0  baseline Dyer                                   0 parameters
    F1  gated power   kappa0 * (S/ref)^b  for S < ref   2 (wide grid, ref up to 3000 psi)
    F2  linear        kappa0 * S / S0                   1
    F4  kappa' = sqrt(kappa0^2 - 1)  (= sqrt(S/D), D = P_sat - P_down)   0
    F5  kappa' = kappa0^2 - 1        (= S/D)                            0
    F6  soft gate     kappa0 * S / (S + S1)             1  (bounded: -> kappa0 for large S)
    F7  soft gate     kappa0 * (S/(S+S1))^b             2
    Identity used for F4/F5:  kappa0^2 = (P_up - P_down)/(P_sat - P_down)
                                       = 1 + S/D,  so kappa0^2 - 1 = S/D.
    Every form with factor -> 0 as S -> 0 makes Dyer -> HEM at the flashing
    threshold (column "f(S->0)" is the numerical check).

TWO WAYS TO TREAT THE DISCHARGE COEFFICIENT
    mode A  Cd from the single-phase window of the training curves (as before)
    mode B  Cd fitted JOINTLY with the form parameters on the two-phase
            training points. If mode B returns Cd ~ 0.78 and the same
            ranking, the correction is not a disguised Cd shift; if the
            baseline F0 in mode B picks a much lower Cd and closes most of
            the gap, part of the "correction" was Cd.

PROTOCOL: identical to loocv_kappa_correction.py (9 folds, everything fitted
on the 8 other curves, prediction of the held-out curve). Part A (injector 2,
4 points) is evaluated afterwards with the parameters fitted on ALL
injector-3 data (it never enters a fit).

Usage:  python validation/kappa_correction/loocv_functional_forms.py
Output also written to validation/kappa_correction/loocv_functional_forms_results.txt
"""

import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import loocv_kappa_correction as lk   # noqa: E402  (reuses loaders / builders)

CDS = np.round(np.arange(0.60, 0.9001, 0.005), 6)     # mode-B Cd grid
S_PROBE_PSI = (41.0, 100.0, 200.0, 371.0, 600.0, 1000.0)


# ---------------------------------------------------------------------------
# Candidate forms. factors(S_psi, kappa, params) -> (n_params, n_points):
# the multiplier applied to kappa0.
# ---------------------------------------------------------------------------
def make_forms():
    forms = []

    forms.append(dict(
        name="F0 baseline Dyer", npar=0, S_only=True,
        params=np.zeros((1, 0)),
        factors=lambda S, k, P: np.ones((P.shape[0], len(S))),
        describe=lambda p: "-"))

    b = np.round(np.arange(0.1, 3.0001, 0.1), 10)
    r = np.geomspace(20.0, 3000.0, 45)
    P1 = np.array([(bb, rr) for bb in b for rr in r])
    forms.append(dict(
        name="F1 gated power (b, ref)", npar=2, S_only=True, params=P1,
        factors=lambda S, k, P: np.where(S[None, :] < P[:, 1:2],
                                         (S[None, :] / P[:, 1:2]) ** P[:, 0:1], 1.0),
        describe=lambda p: f"b={p[0]:.1f} ref={p[1]:.0f}"))

    S0 = np.geomspace(50.0, 5000.0, 60)[:, None]
    forms.append(dict(
        name="F2 linear k*S/S0", npar=1, S_only=True, params=S0,
        factors=lambda S, k, P: S[None, :] / P[:, 0:1],
        describe=lambda p: f"S0={p[0]:.0f}"))

    forms.append(dict(
        name="F4 sqrt(k0^2-1)", npar=0, S_only=False, params=np.zeros((1, 0)),
        factors=lambda S, k, P: np.tile(np.sqrt(np.maximum(k ** 2 - 1.0, 0.0)) / k,
                                        (P.shape[0], 1)),
        describe=lambda p: "-"))

    forms.append(dict(
        name="F5 k0^2-1", npar=0, S_only=False, params=np.zeros((1, 0)),
        factors=lambda S, k, P: np.tile((k ** 2 - 1.0) / k, (P.shape[0], 1)),
        describe=lambda p: "-"))

    S1 = np.geomspace(10.0, 5000.0, 60)[:, None]
    forms.append(dict(
        name="F6 soft gate k*S/(S+S1)", npar=1, S_only=True, params=S1,
        factors=lambda S, k, P: S[None, :] / (S[None, :] + P[:, 0:1]),
        describe=lambda p: f"S1={p[0]:.0f}"))

    b7 = np.round(np.arange(0.5, 3.0001, 0.5), 10)
    s7 = np.geomspace(10.0, 5000.0, 40)
    P7 = np.array([(bb, ss) for bb in b7 for ss in s7])
    forms.append(dict(
        name="F7 soft gate (b, S1)", npar=2, S_only=True, params=P7,
        factors=lambda S, k, P: (S[None, :] / (S[None, :] + P[:, 1:2])) ** P[:, 0:1],
        describe=lambda p: f"b={p[0]:.1f} S1={p[1]:.0f}"))
    return forms


def unit_flows(form, S, k, ms, mh):
    """Cd = 1 flows, shape (n_params, n_points)."""
    kp = k[None, :] * form["factors"](S, k, form["params"])
    return (kp / (1.0 + kp)) * ms[None, :] + mh[None, :] / (1.0 + kp)


def _mape_rows(Cd, U, m_exp, mask):
    return 100.0 * (np.abs(Cd * U[:, mask] - m_exp[mask]) / m_exp[mask]).mean(axis=1)


def fit(U, m_exp, mask, cd=None):
    """Mode A (cd given) or mode B (cd=None -> joint fit over CDS).
    Returns (param index, Cd, training MAPE %)."""
    if cd is not None:
        m = _mape_rows(cd, U, m_exp, mask)
        i = int(np.argmin(m))
        return i, cd, float(m[i])
    best = (None, None, np.inf)
    for c in CDS:
        m = _mape_rows(c, U, m_exp, mask)
        i = int(np.argmin(m))
        if m[i] < best[2]:
            best = (i, float(c), float(m[i]))
    return best


def run_forms(pts, partA, cd_samples):
    """Pure logic. Returns {form name: {mode: result dict}} and bookkeeping."""
    forms = make_forms()
    sups = sorted(cd_samples)
    n = len(pts["m_exp"])
    allmask = np.ones(n, bool)
    all_cd = [v for s in sups for v in cd_samples[s]]
    cd_all = sum(all_cd) / len(all_cd)
    res = {}

    for form in forms:
        U = unit_flows(form, pts["S_psi"], pts["kappa"], pts["m_spi"], pts["m_hem"])
        UA = unit_flows(form, partA["S_psi"], partA["kappa"], partA["m_spi"], partA["m_hem"])
        res[form["name"]] = {"form": form, "f_zero": None}

        for mode in ("A", "B"):
            pooled, sup_list, chosen, cd_chosen, f41 = [], [], [], [], None
            for k in sups:
                test = pts["sup"] == k
                if not test.any():
                    continue
                train = ~test
                if mode == "A":
                    cal = [v for s in sups if s != k for v in cd_samples[s]]
                    i, cd, _ = fit(U, pts["m_exp"], train, cd=sum(cal) / len(cal))
                else:
                    i, cd, _ = fit(U, pts["m_exp"], train)
                e = 100.0 * (cd * U[i, test] - pts["m_exp"][test]) / pts["m_exp"][test]
                pooled.extend(e)
                sup_list.extend([k] * int(test.sum()))
                chosen.append(i)
                cd_chosen.append(cd)
                if abs(k - 41.0) < 1.0:
                    f41 = e
            i_all, cd_in, m_in = fit(U, pts["m_exp"], allmask,
                                     cd=cd_all if mode == "A" else None)
            if mode == "A":   # factor kappa'/kappa0 at S -> 0 with the fitted in-sample parameters
                res[form["name"]]["f_zero"] = float(form["factors"](
                    np.array([0.05]), np.array([1.0001]),
                    form["params"][i_all:i_all + 1])[0, 0])
            partA_cds = (0.65, 0.681) if mode == "A" else (0.65, 0.681, cd_in)
            pa = {}
            for c in partA_cds:
                ea = 100.0 * (c * UA[i_all] - partA["m_exp"]) / partA["m_exp"]
                pa[c] = float(np.abs(ea).mean())
            res[form["name"]][mode] = {
                "pooled": np.array(pooled), "sup": np.array(sup_list),
                "chosen": chosen, "cd_chosen": cd_chosen, "f41": f41,
                "in_idx": i_all, "in_cd": cd_in, "in_mape": m_in, "partA": pa}
    return res, cd_all


def s_probe_table(res, out):
    out("Behaviour of the in-sample (mode A) factor kappa'/kappa0 outside the data (data end at 371 psi):")
    out(f"  {'form':<26} {'fitted':<18}" + "".join(f"{s:>8.0f}" for s in S_PROBE_PSI) + "   [psi]")
    for name, r in res.items():
        form = r["form"]
        if form["npar"] == 0 and not form["S_only"]:
            out(f"  {name:<26} {'-':<18}  (depends on kappa0, not on S alone)")
            continue
        p = form["params"][r["A"]["in_idx"]:r["A"]["in_idx"] + 1]
        f = form["factors"](np.array(S_PROBE_PSI), np.ones(len(S_PROBE_PSI)), p)[0]
        out(f"  {name:<26} {form['describe'](form['params'][r['A']['in_idx']]):<18}"
            + "".join(f"{v:8.2f}" for v in f))


def report(res, cd_all, out):
    bar = "=" * 100
    out(bar)
    out("LOOCV of candidate functional forms -- Waxman Fig. 13 (64 two-phase points, 9 curves)")
    out(bar)
    out(f"pooled single-phase Cd (all curves) = {cd_all:.4f}")
    out("mode A: Cd from single-phase windows of training curves | mode B: Cd fitted jointly with the form")
    out()
    for mode in ("A", "B"):
        out(f"MODE {mode}")
        out(f"  {'form':<26} {'npar':>4} {'f(S->0)':>8} | {'pooled OOS':<46} | {'<100psi':>8} {'41psi':>7} | "
            f"{'PartA MAPE':>22} | in-sample")
        for name, r in res.items():
            m = r[mode]
            e, s = m["pooled"], m["sup"]
            low = e[s < 100]
            f41 = m["f41"]
            pa = "  ".join(f"{c:.3f}:{v:5.2f}%" for c, v in m["partA"].items())
            out(f"  {name:<26} {r['form']['npar']:>4} {r['f_zero']:8.3f} | {lk.stats(e):<46} | "
                f"{np.abs(low).mean():7.2f}% {np.abs(f41).mean():6.2f}% | {pa:>22} | {m['in_mape']:5.2f}%")
        out()
    out("Cd chosen by the joint fit (mode B), range over the 9 folds, and in-sample:")
    for name, r in res.items():
        c = r["B"]["cd_chosen"]
        out(f"  {name:<26} folds {min(c):.3f}-{max(c):.3f} | in-sample {r['B']['in_cd']:.3f}")
    out()
    out("Parameter stability across folds (mode A):")
    for name, r in res.items():
        form = r["form"]
        if form["npar"] == 0:
            continue
        P = form["params"][r["A"]["chosen"]]
        rng = " ".join(f"p{j}: {P[:, j].min():.3g}-{P[:, j].max():.3g}" for j in range(P.shape[1]))
        out(f"  {name:<26} {rng}   (in-sample: {form['describe'](form['params'][r['A']['in_idx']])})")
    out()
    s_probe_table(res, out)


def main():
    lk._load_project()
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    curves, _ = lk.wv.load_multiseries("waxman_fig13_mdot_vs_dP_by_supercharge.csv")
    cd_samples = {c["super_psi"]: lk.wv.spi_cd_samples(c) for c in curves}
    pts, pts_list = lk.build_points(curves)
    lk.self_check(pts, pts_list)

    # Part A points (injector 2), unit-Cd pieces
    ex, n2o, inj = lk.ex, lk.n2o, lk.inj
    S_A = (ex.P1_A - n2o.P_sat(ex.T1_A)) / lk.wv.PSI_PA
    cols = {k: [] for k in ("kappa", "m_spi", "m_hem", "S_psi", "m_exp")}
    for _label, dP_MPa, m_gs in ex.CASES_A:
        P2 = ex.P1_A - dP_MPa * 1e6
        T_d = n2o.T_sat(P2)
        r = inj.dyer_mass_flow_corrected(
            1.0, ex.A_A, ex.T1_A, ex.P1_A, P2, n2o.rho_liquid_sat(ex.T1_A),
            n2o.rho_liquid_sat(T_d), n2o.M_N2O / n2o.nu_vapor_sat(T_d),
            beta=0.0, supercharge_ref_Pa=1.0)
        cols["kappa"].append(r["kappa"]); cols["m_spi"].append(r["m_dot_SPI"])
        cols["m_hem"].append(r["m_dot_HEM"]); cols["S_psi"].append(S_A)
        cols["m_exp"].append(m_gs * 1e-3)
    partA = {k: np.array(v) for k, v in cols.items()}

    # consistency: F1 factors must reproduce the reference formula of loocv_kappa_correction
    f1 = next(f for f in make_forms() if f["name"].startswith("F1"))
    U1 = unit_flows(f1, pts["S_psi"], pts["kappa"], pts["m_spi"], pts["m_hem"])
    for idx in (0, len(f1["params"]) // 2, len(f1["params"]) - 1):
        b, r_ = f1["params"][idx]
        ref = lk.corrected_unit_flow(pts["kappa"], pts["m_spi"], pts["m_hem"], pts["S_psi"], b, r_)
        if not np.allclose(U1[idx], ref, rtol=1e-12):
            raise RuntimeError("F1 form disagrees with the reference formula")

    res, cd_all = run_forms(pts, partA, cd_samples)
    report(res, cd_all, out)

    path = os.path.join(lk._REPO_ROOT, "validation", "kappa_correction", "loocv_functional_forms_results.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n(output also written to {path})")


if __name__ == "__main__":
    main()
