"""
waxman_2013_validation.py

Validation of the coupled injector model against the Waxman experimental
programme (Stanford / NASA Ames):

    Waxman, Zimmerman, Cantwell & Zilliac (2013), AIAA 2013-3636,
    "Mass Flow Rate and Isolation Characteristics of Injectors for Use with
    Self-Pressurizing Oxidizers in Hybrid Rockets".

Two datasets are used.

PART A -- the original four operating points (Nino & Razavi 2019, Table 4)
    Injector pressure drops of 8-14 bar at a single upstream state. Kept
    exactly as before so that the earlier reported figures stay reproducible.

PART B -- the full injector-3 mass-flow map, Waxman (2013) Fig. 13
    Nine supercharge levels (41-371 psi), pressure drops up to ~46 bar.
    Digitised by the project author from the paper's figures
    (validation/digitized/). This extends the validated domain from
    8-14 bar to 0-46 bar, and to a wide range of subcooling margins.

    Methodology (fixed BEFORE looking at the two-phase results):
      * The discharge coefficient is calibrated ONLY on the single-phase
        window of each curve (30 psi <= dP <= supercharge, where the model
        is pure SPI), pooled over all nine curves -- Cd of an injector is
        not expected to depend on supercharge (Waxman Fig. 15).
      * The Dyer model is then evaluated, with that one Cd, on every point
        that lies beyond the saturation threshold (dP > supercharge).
        Those points are genuine predictions: none of them entered the
        calibration.
      * The Henry-Fauske non-equilibrium ceiling is reported next to the
        Dyer prediction, exactly as in the tool (diagnostic, NOT applied).
        A "capped" error column, min(Dyer, HF), is computed only to
        quantify how the ceiling would have performed against data.

PART C -- critical mass flow, Waxman Fig. 16 (same 5 % criterion as the paper)
PART D -- SPI discharge coefficient of injectors 1, 2, 5, Waxman Fig. 15
PART E -- internal consistency checks of the digitised data (fail loudly)

Digitised-data notes (see validation/waxman_2013_results.md, Section 4):
  * Fig. 11 y-values in the raw file are 10x too large (axis calibration in
    the digitiser); FIG11_Y_SCALE corrects it and PART E verifies the
    correction against Fig. 13.
  * Fig. 13 contains one stray calibration point (negative dP); it is dropped.
  * Points with dP < 30 psi are ignored, as in the paper (Cd_eff scatters
    strongly as dP and mass flow tend to zero).

Usage (from any directory):
    python validation/waxman_2013_validation.py            # tables + figure
    python validation/waxman_2013_validation.py --no-plot  # tables only
"""

import math
import os
import re
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "model"))

from full_system import evaluate_full_system
from injector_spi import spi_mass_flow
from injector_two_phase import (hem_critical_flow, hem_critical_flow_isentropic,
                                 henry_fauske_critical_flow, dyer_mass_flow)
from n2o_properties import P_sat, T_sat, rho_liquid_sat, nu_vapor_sat, M_N2O

DIGITIZED_DIR = os.path.join(_REPO_ROOT, "validation", "digitized")
FIGURE_PATH = os.path.join(_REPO_ROOT, "validation", "waxman_2013_fig13_comparison.png")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
PSI_PA = 6894.757            # Pa per psi
P_ATM_PSI = 14.696           # psig -> psia
BAR_PER_PSI = PSI_PA / 1e5

# Waxman injector 3: straight hole, rounded inlet, D = 1.50 mm, L = 18.4 mm (Table 1)
D_INJ3 = 1.50e-3
A_INJ3 = math.pi * (D_INJ3 / 2.0) ** 2
CD_INJ3_PAPER = 0.77         # Waxman p.14 (Fig. 11/12 conditions), used as a sensitivity only

# Upstream volume of the rig (25.4 mm ID, 50 mm): losses are negligible, as in Part A.
LINE = [{"type": "pipe", "L": 0.05, "D": 0.0254}]

MIN_DP_PSI = 30.0            # ignore the noisy low-dP region (paper, p.14)
CRIT_RATIO = 0.95            # Waxman "critical dP": Cd_eff first below 95 % of SPI
FIG11_Y_SCALE = 0.1          # raw Fig. 11 y-values are 10x too large (see docstring)

# ---------------------------------------------------------------------------
# PART A -- original four points (Nino & Razavi 2019, Table 4). UNCHANGED.
# ---------------------------------------------------------------------------
T1 = 280.0
P1 = 4.36e6
Cd = 0.65
D = 0.0015
A = math.pi * (D / 2.0) ** 2

SEGMENTS = [{"type": "pipe", "L": 0.05, "D": 0.0254}]

CASES = [
    ("Pre-critical",  0.84, 44.0),
    ("Critical",      0.98, 46.5),
    ("Post-critical1", 1.09, 47.5),
    ("Post-critical2", 1.37, 48.0),
]


def pct_error(pred, exp):
    return 100.0 * (pred - exp) / exp


def run_part_a(cd=Cd, verbose=True):
    """Original 4-point validation. Returns (list of errors [%], list of results)."""
    if verbose:
        Psat = P_sat(T1)
        print(f"P_sat(T1) model = {Psat/1e5:.3f} bar  (Nino&Razavi: {(P1-0.62e6)/1e5:.2f} bar,  "
              f"model err = {100*(Psat-(P1-0.62e6))/(P1-0.62e6):+.1f}%)")
        print(f"Supercharge: model = {(P1-Psat)/1e5:.2f} bar  |  Nino&Razavi = 6.20 bar")
        print()
    results = []
    for label, dP_MPa, m_exp in CASES:
        P_chamber = P1 - dP_MPa * 1e6
        r = evaluate_full_system(m_exp * 0.9e-3, T1, P1, SEGMENTS, cd, A, P_chamber)
        m = r["m_dot_real"]
        si = r["solver_info"]
        err = pct_error(m, m_exp * 1e-3)
        ir = r.get("injector_result") or {}
        choked_flag = "CHOKED" if ir.get("choked") else "-"
        results.append((label, dP_MPa, m_exp, m, err, r))
        if verbose:
            print(f"  {label:<16} dP={dP_MPa*10:6.2f} bar  m_exp={m_exp:5.1f} g/s  "
                  f"m_dot={m*1000:7.3f} g/s  err={err:+.2f}%  "
                  f"P_in={r['P_injector_inlet']/1e5:.3f} bar  regime={r['regime']}  "
                  f"iters={si['iterations']}  {choked_flag}")
    errs = [r[4] for r in results]
    if verbose:
        print()
        print(f"Mean error: {sum(errs)/len(errs):+.2f}%")
        print(f"MAPE:       {sum(abs(e) for e in errs)/len(errs):.2f}%")
        print(f"Within +/-5%: {sum(1 for e in errs if abs(e)<=5)}/{len(errs)}")
        print(f"Within +/-10%: {sum(1 for e in errs if abs(e)<=10)}/{len(errs)}")
        print()
        crit = hem_critical_flow(cd, A, T1, P1)
        crit_s = hem_critical_flow_isentropic(cd, A, T1, P1)
        print(f"HEM isenthalpic: m_dot_crit={crit['m_dot_crit']*1000:.2f} g/s @ "
              f"P2={crit['P2_crit']/1e5:.2f} bar x={crit['x_crit']:.4f}")
        print(f"HEM isentropic:  m_dot_crit={crit_s['m_dot_crit']*1000:.2f} g/s @ "
              f"P2={crit_s['P2_crit']/1e5:.2f} bar x={crit_s['x_crit']:.4f}")
        hf = henry_fauske_critical_flow(cd, A, T1, P1)
        print(f"Henry-Fauske:    m_dot_crit={hf['m_dot_crit']*1000:.2f} g/s @ "
              f"P2={hf['P2_crit']/1e5:.2f} bar x_E={hf['x_crit']:.4f} N={hf['N']:.3f}")
    return errs, results


# ---------------------------------------------------------------------------
# Loaders for the digitised figures
# ---------------------------------------------------------------------------
_TITLE_RE = re.compile(
    r"P1\s*super\s*=\s*(\d+)\s*psi\s*,\s*T1\s*=\s*(\d+)\s*K\s*,\s*P1\s*=\s*(\d+)\s*psig")


def _lines(filename):
    path = os.path.join(DIGITIZED_DIR, filename)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Digitised data file not found: {path}")
    with open(path, encoding="utf-8") as f:
        return f.read().splitlines()


def load_multiseries(filename):
    """
    Load Fig. 13 / Fig. 14 style files: nine (X, Y) series, one per supercharge
    level, with the legend text in one header line. The legend contains commas
    itself, so it is parsed with a regex on the whole line, NOT column-wise.

    Returns
    -------
    curves : list of dict
        {"super_psi", "T1_stated_K", "P1_psig", "dP_psi": [...], "y": [...]}
        sorted by decreasing supercharge order as in the file; points sorted
        by dP. Points with dP <= 0 (stray calibration points) are dropped.
    n_dropped : int
    """
    lines = _lines(filename)
    title_line = next((l for l in lines if "P1super" in l or "P1 super" in l), None)
    if title_line is None:
        raise RuntimeError(f"{filename}: legend line with 'P1super' not found")
    titles = _TITLE_RE.findall(title_line)
    if not titles:
        raise RuntimeError(f"{filename}: could not parse the series legend")
    start = next(i for i, l in enumerate(lines) if l.startswith("X,Y"))
    curves = [{"super_psi": float(s), "T1_stated_K": float(t), "P1_psig": float(p),
               "dP_psi": [], "y": []} for s, t, p in titles]
    n_dropped = 0
    for line in lines[start + 1:]:
        if not line.strip():
            continue
        f = line.split(",")
        for k, c in enumerate(curves):
            if 2 * k + 1 < len(f) and f[2 * k].strip() and f[2 * k + 1].strip():
                x, y = float(f[2 * k]), float(f[2 * k + 1])
                if x <= 0:
                    n_dropped += 1
                    continue
                c["dP_psi"].append(x)
                c["y"].append(y)
    for c in curves:
        order = sorted(range(len(c["dP_psi"])), key=lambda i: c["dP_psi"][i])
        c["dP_psi"] = [c["dP_psi"][i] for i in order]
        c["y"] = [c["y"][i] for i in order]
    return curves, n_dropped


def load_xy(filename):
    """Load a plain 'x, y' file with '#' comment lines."""
    pts = []
    for l in _lines(filename):
        if not l.strip() or l.lstrip().startswith("#"):
            continue
        x, y = l.split(",")[:2]
        pts.append((float(x), float(y)))
    return sorted(pts)


def load_fig15():
    """Fig. 15: Cd vs supercharge for injectors 1 (0.79 mm), 2 (1.50 mm), 5 (1.93 mm)."""
    lines = _lines("waxman_fig15_cd_vs_supercharge_injectors_1_2_5.csv")
    start = next(i for i, l in enumerate(lines) if l.startswith("X,Y"))
    series = {1: [], 2: [], 5: []}
    for l in lines[start + 1:]:
        if not l.strip():
            continue
        f = l.split(",")
        for k, inj in enumerate((1, 2, 5)):
            if 2 * k + 1 < len(f) and f[2 * k].strip() and f[2 * k + 1].strip():
                series[inj].append((float(f[2 * k]), float(f[2 * k + 1])))
    return series


# ---------------------------------------------------------------------------
# Physics helpers for Part B / C
# ---------------------------------------------------------------------------
def curve_state(c, T_mode="implied"):
    """
    Upstream (T1, P1 [Pa]) of a curve.

    "implied": T1 = T_sat(P1_abs - P1super). The paper defines the supercharge
        as P1 - P_sat(T1); the legend gives T1 only to the nearest kelvin
        (+/-0.5 K is +/-~5 psi of P_sat, large next to a 41 psi supercharge),
        while P1 (psig) and P1super (psi) are given to ~1 psi. Using them to
        recover T1 keeps the supercharge -- the quantity that controls the
        flashing threshold -- exactly as reported.
    "stated": T1 = the integer kelvin from the legend.
    """
    P1_pa = (c["P1_psig"] + P_ATM_PSI) * PSI_PA
    if T_mode == "implied":
        T = T_sat(P1_pa - c["super_psi"] * PSI_PA)
    elif T_mode == "stated":
        T = c["T1_stated_K"]
    else:
        raise ValueError(f"unknown T_mode {T_mode!r}")
    return T, P1_pa


def spi_cd_samples(c, T_mode="implied"):
    """Implied Cd of the single-phase window points (30 psi <= dP <= supercharge)."""
    T, _ = curve_state(c, T_mode)
    rho = rho_liquid_sat(T)
    out = []
    for dP, m in zip(c["dP_psi"], c["y"]):
        if MIN_DP_PSI <= dP <= c["super_psi"]:
            out.append(m / (A_INJ3 * math.sqrt(2.0 * rho * dP * PSI_PA)))
    return out


def model_flow_direct(T, P_up, P_down, cd):
    """Injector model without the (negligible) upstream line: SPI or Dyer."""
    rho_up = rho_liquid_sat(T)
    if P_down >= P_sat(T):
        return spi_mass_flow(cd, A_INJ3, rho_up, P_up - P_down), None
    Td = T_sat(P_down)
    d = dyer_mass_flow(cd, A_INJ3, T, P_up, P_down, rho_up,
                       rho_liquid_sat(Td), M_N2O / nu_vapor_sat(Td))
    return d["m_dot_Dyer"], d


def predict_curve(c, cd_for_curve, T_mode="implied"):
    """Coupled-solver prediction at every usable point of one curve."""
    T, P1_pa = curve_state(c, T_mode)
    m_hf = henry_fauske_critical_flow(cd_for_curve, A_INJ3, T, P1_pa)["m_dot_crit"]
    rows = []
    for dP, m_exp in zip(c["dP_psi"], c["y"]):
        if dP < MIN_DP_PSI:
            continue
        P2 = P1_pa - dP * PSI_PA
        r = evaluate_full_system(m_exp * 0.9, T, P1_pa, LINE, cd_for_curve, A_INJ3, P2)
        m = r["m_dot_real"]
        m_cap = min(m, m_hf) if r["regime"] == "Dyer" else m
        rows.append({
            "super_psi": c["super_psi"], "dP_psi": dP, "dP_bar": dP * BAR_PER_PSI,
            "m_exp": m_exp, "m_model": m, "regime": r["regime"],
            "err": pct_error(m, m_exp), "m_hf": m_hf,
            "m_capped": m_cap, "err_capped": pct_error(m_cap, m_exp),
            "choked": bool((r.get("injector_result") or {}).get("choked", False)),
        })
    return rows


def summarize(errs):
    n = len(errs)
    if n == 0:
        return {"n": 0, "mean": float("nan"), "mape": float("nan"), "max": float("nan"),
                "w5": 0, "w10": 0}
    return {"n": n, "mean": sum(errs) / n, "mape": sum(abs(e) for e in errs) / n,
            "max": max(abs(e) for e in errs),
            "w5": sum(1 for e in errs if abs(e) <= 5), "w10": sum(1 for e in errs if abs(e) <= 10)}


def _fmt(s):
    return (f"n={s['n']:3d}  mean={s['mean']:+6.2f}%  MAPE={s['mape']:5.2f}%  "
            f"max|err|={s['max']:5.1f}%  within5={s['w5']}/{s['n']}  within10={s['w10']}/{s['n']}")


def run_part_b(cd_mode="pooled", T_mode="implied", verbose=True):
    """
    Full Fig. 13 validation. cd_mode: "pooled" (default), "per_curve", or "paper".
    Returns dict with rows, Cd values and per-curve state.
    """
    curves, n_dropped = load_multiseries("waxman_fig13_mdot_vs_dP_by_supercharge.csv")
    samples = {c["super_psi"]: spi_cd_samples(c, T_mode) for c in curves}
    pooled_all = [v for s in samples.values() for v in s]
    if not pooled_all:
        raise RuntimeError("no single-phase calibration points found")
    cd_pooled = sum(pooled_all) / len(pooled_all)

    rows, cd_used = [], {}
    for c in curves:
        if cd_mode == "pooled":
            cd_c = cd_pooled
        elif cd_mode == "paper":
            cd_c = CD_INJ3_PAPER
        elif cd_mode == "per_curve":
            s = samples[c["super_psi"]]
            cd_c = sum(s) / len(s) if s else cd_pooled
        else:
            raise ValueError(f"unknown cd_mode {cd_mode!r}")
        cd_used[c["super_psi"]] = cd_c
        rows.extend(predict_curve(c, cd_c, T_mode))
    return {"rows": rows, "curves": curves, "cd_pooled": cd_pooled, "cd_used": cd_used,
            "n_cal": len(pooled_all), "cd_cal_range": (min(pooled_all), max(pooled_all)),
            "n_dropped_stray": n_dropped, "T_mode": T_mode, "cd_mode": cd_mode,
            "n_ignored_lowdp": sum(1 for c in curves for d in c["dP_psi"] if d < MIN_DP_PSI)}


def print_part_b():
    res = run_part_b()
    rows, curves = res["rows"], res["curves"]
    print(f"Injector 3: D = {D_INJ3*1e3:.2f} mm, A = {A_INJ3*1e6:.4f} mm2, rounded inlet")
    print(f"Digitised points: {sum(len(c['dP_psi']) for c in curves)} in Fig. 13 "
          f"(+{res['n_dropped_stray']} stray point dropped); "
          f"{res['n_ignored_lowdp']} with dP < {MIN_DP_PSI:.0f} psi ignored; "
          f"{len(rows)} used.")
    print(f"Calibration: Cd pooled over the single-phase window "
          f"({MIN_DP_PSI:.0f} psi <= dP <= supercharge), n = {res['n_cal']} points: "
          f"Cd = {res['cd_pooled']:.4f} (range of individual values "
          f"{res['cd_cal_range'][0]:.3f}-{res['cd_cal_range'][1]:.3f}); "
          f"Waxman reports ~{CD_INJ3_PAPER}.")
    print()
    print("Upstream states (T1 recovered from P1 and P1super; stated integer T1 in brackets):")
    print(f"  {'super[psi]':>10} {'super[bar]':>10} {'P1[bar]':>8} {'T1[K]':>8} {'(stated)':>9} "
          f"{'Psat[bar]':>9} {'HF ceiling[g/s]':>16} {'n(Dyer pts)':>11}")
    for c in sorted(curves, key=lambda c: c["super_psi"]):
        T, P1_pa = curve_state(c)
        hf = rows and next(r["m_hf"] for r in rows if r["super_psi"] == c["super_psi"])
        n_dy = sum(1 for r in rows if r["super_psi"] == c["super_psi"] and r["regime"] == "Dyer")
        print(f"  {c['super_psi']:10.0f} {c['super_psi']*BAR_PER_PSI:10.2f} {P1_pa/1e5:8.2f} "
              f"{T:8.2f} {c['T1_stated_K']:9.0f} {P_sat(T)/1e5:9.2f} {hf*1000:16.1f} {n_dy:11d}")
    print()

    spi_rows = [r for r in rows if r["regime"] == "SPI"]
    dy_rows = [r for r in rows if r["regime"] == "Dyer"]
    print("B1. Single-phase region (model regime SPI) -- CALIBRATION RESIDUAL, not a prediction:")
    print("   ", _fmt(summarize([r["err"] for r in spi_rows])))
    print("B2. Two-phase region (model regime Dyer, dP > supercharge) -- PREDICTION:")
    print("   Dyer (as in the tool):         ", _fmt(summarize([r["err"] for r in dy_rows])))
    print("   min(Dyer, Henry-Fauske) [diag]:", _fmt(summarize([r["err_capped"] for r in dy_rows])))
    print(f"   Dyer prediction above HF ceiling ('choked'): "
          f"{sum(1 for r in dy_rows if r['choked'])}/{len(dy_rows)} points")
    print()

    print("B3. Two-phase region by injector pressure drop:")
    bands = [("dP <= 14 bar  (previously validated band)", 0.0, 14.0),
             ("14 < dP <= 30 bar", 14.0, 30.0),
             ("30 < dP <= 46 bar", 30.0, 100.0)]
    for name, lo, hi in bands:
        b = [r for r in dy_rows if lo < r["dP_bar"] <= hi]
        print(f"   {name:<42} Dyer:   {_fmt(summarize([r['err'] for r in b]))}")
        print(f"   {'':<42} capped: {_fmt(summarize([r['err_capped'] for r in b]))}")
    print()

    print("B4. Two-phase region by supercharge (subcooling margin at the tank):")
    print(f"   {'super[psi]':>10} {'super[bar]':>10}   Dyer                                    "
          f"min(Dyer,HF)")
    for c in sorted(curves, key=lambda c: c["super_psi"]):
        b = [r for r in dy_rows if r["super_psi"] == c["super_psi"]]
        sd, sc = summarize([r["err"] for r in b]), summarize([r["err_capped"] for r in b])
        print(f"   {c['super_psi']:10.0f} {c['super_psi']*BAR_PER_PSI:10.2f}   "
              f"n={sd['n']:2d} mean={sd['mean']:+6.1f}% MAPE={sd['mape']:5.1f}% max={sd['max']:5.1f}%   "
              f"mean={sc['mean']:+6.1f}% MAPE={sc['mape']:5.1f}%")
    hi = [r for r in dy_rows if r["super_psi"] >= 200]
    lo = [r for r in dy_rows if r["super_psi"] < 200]
    print(f"   supercharge >= 200 psi (>= 1.38 MPa): Dyer {_fmt(summarize([r['err'] for r in hi]))}")
    print(f"   supercharge <  200 psi              : Dyer {_fmt(summarize([r['err'] for r in lo]))}")
    print()

    print("B5. Sensitivity of the two-phase-region result (Dyer, as in the tool):")
    print(f"   {'variant':<52} {'mean':>8} {'MAPE':>7} {'max|err|':>9}")
    variants = [("baseline: pooled Cd, T1 recovered from P1 and P1super", "pooled", "implied"),
                (f"Cd = {CD_INJ3_PAPER} (Waxman text) instead of pooled fit", "paper", "implied"),
                ("Cd fitted curve by curve (fallback: pooled)", "per_curve", "implied"),
                ("T1 = integer kelvin from the legend", "pooled", "stated")]
    for name, cm, tm in variants:
        rr = res if (cm, tm) == ("pooled", "implied") else run_part_b(cm, tm)
        s = summarize([r["err"] for r in rr["rows"] if r["regime"] == "Dyer"])
        print(f"   {name:<52} {s['mean']:+7.2f}% {s['mape']:6.2f}% {s['max']:8.1f}%   (n={s['n']})")
    print()

    print("B6. Point-by-point table (two-phase region):")
    print(f"   {'super':>5} {'dP[bar]':>8} {'m_exp[g/s]':>11} {'m_Dyer':>8} {'err%':>7} "
          f"{'HF ceil':>8} {'min(D,HF)':>10} {'err%':>7}")
    for r in sorted(dy_rows, key=lambda r: (r["super_psi"], r["dP_psi"])):
        print(f"   {r['super_psi']:5.0f} {r['dP_bar']:8.2f} {r['m_exp']*1000:11.1f} "
              f"{r['m_model']*1000:8.1f} {r['err']:+7.1f} {r['m_hf']*1000:8.1f} "
              f"{r['m_capped']*1000:10.1f} {r['err_capped']:+7.1f}")
    return res


# ---------------------------------------------------------------------------
# PART C -- critical mass flow (Fig. 16)
# ---------------------------------------------------------------------------
def model_critical_point(T, P1_pa, cd, super_psi):
    """
    Same criterion as Waxman: the first dP for which the effective discharge
    coefficient drops below 95 % of the SPI value, i.e. m_Dyer < 0.95 m_SPI.
    Scans upward from the saturation threshold in 0.25 psi steps.
    Returns (dP_crit [psi], m_crit [kg/s]) or (None, None).
    """
    dP = super_psi * PSI_PA + 0.25 * PSI_PA
    while P1_pa - dP > 1.5e5:
        m, d = model_flow_direct(T, P1_pa, P1_pa - dP, cd)
        if d is not None and d["m_dot_Dyer"] < CRIT_RATIO * d["m_dot_SPI"]:
            return dP / PSI_PA, m
        dP += 0.25 * PSI_PA
    return None, None


def run_part_c(res_b):
    curves = res_b["curves"]
    fig16 = load_xy("waxman_fig16_critical_mdot_vs_supercharge.csv")
    cd = res_b["cd_pooled"]
    out = []
    for sup16, m16 in fig16:
        c = min(curves, key=lambda c: abs(c["super_psi"] - sup16))
        T, P1_pa = curve_state(c)
        m_exp13 = max(c["y"])
        dPc, m_mod = model_critical_point(T, P1_pa, cd, c["super_psi"])
        hf = henry_fauske_critical_flow(cd, A_INJ3, T, P1_pa)["m_dot_crit"]
        hem = hem_critical_flow(cd, A_INJ3, T, P1_pa)["m_dot_crit"]
        out.append({"super16": sup16, "super_curve": c["super_psi"], "m16": m16,
                    "m_exp13_max": m_exp13, "dP_crit_model_psi": dPc, "m_model95": m_mod,
                    "m_hf": hf, "m_hem": hem})
    return out


def print_part_c(res_b):
    out = run_part_c(res_b)
    print("Experimental critical flow: Waxman Fig. 16 (injector 3, N2O). Model: Cd pooled, Dyer,")
    print("first dP where m_Dyer < 0.95 m_SPI (the paper's criterion). HF and HEM are the standalone")
    print("ceilings at the same upstream state (for reference).")
    print()
    print(f"  {'super[psi]':>10} {'m_exp(F16)':>11} {'max(F13)':>9} {'m_Dyer95':>9} {'err vs F16':>11} "
          f"{'dP95[bar]':>10} {'HF':>7} {'err':>7} {'HEM':>7} {'err':>7}")
    errs95, errshf, errshem = [], [], []
    for o in out:
        e95 = pct_error(o["m_model95"], o["m16"]) if o["m_model95"] else float("nan")
        ehf = pct_error(o["m_hf"], o["m16"]); ehem = pct_error(o["m_hem"], o["m16"])
        errs95.append(e95); errshf.append(ehf); errshem.append(ehem)
        print(f"  {o['super16']:10.1f} {o['m16']*1000:11.1f} {o['m_exp13_max']*1000:9.1f} "
              f"{(o['m_model95'] or float('nan'))*1000:9.1f} {e95:+10.1f}% "
              f"{(o['dP_crit_model_psi'] or float('nan'))*BAR_PER_PSI:10.1f} "
              f"{o['m_hf']*1000:7.1f} {ehf:+6.1f}% {o['m_hem']*1000:7.1f} {ehem:+6.1f}%")
    print()
    print("  Dyer@95%  :", _fmt(summarize(errs95)))
    print("  Henry-Fauske ceiling:", _fmt(summarize(errshf)))
    print("  HEM (isenthalpic) ceiling:", _fmt(summarize(errshem)))
    return out


# ---------------------------------------------------------------------------
# PART D -- Fig. 15
# ---------------------------------------------------------------------------
def print_part_d():
    s = load_fig15()
    names = {1: "0.79 mm, square edge", 2: "1.50 mm, square edge", 5: "1.93 mm, square edge"}
    means = {}
    print("  Injector  geometry                 n   mean Cd   range")
    for inj in (1, 2, 5):
        v = [y for _, y in s[inj]]
        means[inj] = sum(v) / len(v)
        print(f"  {inj:^8}  {names[inj]:<22} {len(v):3d}   {means[inj]:.3f}    {min(v):.3f}-{max(v):.3f}")
    return means


# ---------------------------------------------------------------------------
# PART E -- consistency checks of the digitised data
# ---------------------------------------------------------------------------
def _interp(x, pts):
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        if x0 <= x <= x1:
            return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    return None


def run_part_e(res_b, tol_fig11=0.10):
    curves = res_b["curves"]
    c169 = next(c for c in curves if c["super_psi"] == 169.0)
    ref = list(zip(c169["dP_psi"], c169["y"]))
    raw11 = load_xy("waxman_fig11_mdot_vs_dP_single_test.csv")
    ratios_raw, ratios_fix = [], []
    for x, y in raw11:
        yi = _interp(x, ref)
        if yi:
            ratios_raw.append(y / yi); ratios_fix.append(y * FIG11_Y_SCALE / yi)
    med_raw = sorted(ratios_raw)[len(ratios_raw) // 2]
    med_fix = sorted(ratios_fix)[len(ratios_fix) // 2]
    print(f"E1. Fig. 11 (single test, P1super = 169 psi) vs the 169 psi curve of Fig. 13:")
    print(f"    raw file / Fig. 13:                 median ratio = {med_raw:6.2f}  (n={len(ratios_raw)})")
    print(f"    raw x {FIG11_Y_SCALE} (FIG11_Y_SCALE) / Fig. 13:   median ratio = {med_fix:6.3f}")
    if abs(med_fix - 1.0) > tol_fig11:
        raise RuntimeError(
            f"Fig. 11 scale check failed (median ratio {med_fix:.3f} after applying "
            f"FIG11_Y_SCALE={FIG11_Y_SCALE}). The digitised file has changed or the "
            "scale factor is wrong -- do not use Fig. 11 until this is resolved.")

    # Fig. 14 (Cd_eff digitised) vs Cd implied by Fig. 13
    c14, _ = load_multiseries("waxman_fig14_cd_vs_dP_by_supercharge.csv")
    diffs = []
    for a, b in zip(curves, c14):
        if a["super_psi"] != b["super_psi"]:
            raise RuntimeError("Fig. 13 / Fig. 14 curve order mismatch")
        T, _ = curve_state(a)
        rho = rho_liquid_sat(T)
        pts13 = list(zip(a["dP_psi"], a["y"]))
        for dP, cd14 in zip(b["dP_psi"], b["y"]):
            if dP < MIN_DP_PSI:
                continue
            m = _interp(dP, pts13)
            if m:
                diffs.append(m / (A_INJ3 * math.sqrt(2 * rho * dP * PSI_PA)) / cd14)
    diffs.sort()
    print(f"E2. Cd_eff implied by Fig. 13 (model saturated-liquid density) / Cd_eff read from Fig. 14:")
    print(f"    n={len(diffs)}  median={diffs[len(diffs)//2]:.3f}  "
          f"5th-95th percentile={diffs[int(0.05*len(diffs))]:.3f}-{diffs[int(0.95*len(diffs))-1]:.3f}")
    # Fig. 16 vs max of Fig. 13
    f16 = load_xy("waxman_fig16_critical_mdot_vs_supercharge.csv")
    d = []
    for sup16, m16 in f16:
        c = min(curves, key=lambda c: abs(c["super_psi"] - sup16))
        d.append(m16 / max(c["y"]))
    print(f"E3. Fig. 16 critical flow / maximum of the matching Fig. 13 curve: "
          f"min={min(d):.3f} max={max(d):.3f}")
    return {"fig11_med_raw": med_raw, "fig11_med_fixed": med_fix,
            "e2_median": diffs[len(diffs) // 2],
            "e2_p05": diffs[int(0.05 * len(diffs))], "e2_p95": diffs[int(0.95 * len(diffs)) - 1],
            "e3_min": min(d), "e3_max": max(d)}


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_figure(res_b, path=FIGURE_PATH):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cd = res_b["cd_pooled"]
    curves = sorted(res_b["curves"], key=lambda c: c["super_psi"])
    fig, axes = plt.subplots(3, 3, figsize=(12.5, 9.5), sharey=False)
    for ax, c in zip(axes.ravel(), curves):
        T, P1_pa = curve_state(c)
        hf = henry_fauske_critical_flow(cd, A_INJ3, T, P1_pa)["m_dot_crit"] * 1000
        dP_max = max(c["dP_psi"]) * 1.05
        grid = [MIN_DP_PSI + i * (dP_max - MIN_DP_PSI) / 90 for i in range(91)]
        m_line, m_spi = [], []
        rho = rho_liquid_sat(T)
        for dP in grid:
            m, _ = model_flow_direct(T, P1_pa, P1_pa - dP * PSI_PA, cd)
            m_line.append(m * 1000)
            m_spi.append(spi_mass_flow(cd, A_INJ3, rho, dP * PSI_PA) * 1000)
        gbar = [g * BAR_PER_PSI for g in grid]
        ax.plot(gbar, m_spi, ":", color="#94A3B8", lw=1.3, label="SPI")
        ax.plot(gbar, m_line, "-", color="#2563EB", lw=2, label="SPI / Dyer (model)")
        ax.axhline(hf, ls="-.", color="#EA580C", lw=1.2, label="Henry-Fauske ceiling")
        ax.axvline(c["super_psi"] * BAR_PER_PSI, ls="--", color="#DC2626", lw=1,
                   label="P2 = P_sat(T1)")
        cal = [(d * BAR_PER_PSI, m * 1000) for d, m in zip(c["dP_psi"], c["y"])
               if MIN_DP_PSI <= d <= c["super_psi"]]
        val = [(d * BAR_PER_PSI, m * 1000) for d, m in zip(c["dP_psi"], c["y"])
               if d > c["super_psi"]]
        if cal:
            ax.plot(*zip(*cal), "o", ms=5, mfc="white", mec="#334155", label="experiment (calibration)")
        if val:
            ax.plot(*zip(*val), "o", ms=5, color="#0F172A", label="experiment (prediction)")
        ax.set_title(f"P1super = {c['super_psi']:.0f} psi ({c['super_psi']*BAR_PER_PSI:.1f} bar), "
                     f"T1 = {T:.1f} K", fontsize=9)
        ax.set_xlabel("Injector pressure drop (bar)", fontsize=8)
        ax.set_ylabel("Mass flow (g/s)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=6, fontsize=8, frameon=False)
    fig.suptitle(f"Waxman (2013) injector 3, N2O: model vs. digitised Fig. 13   (Cd = {cd:.3f})",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
def main():
    no_plot = "--no-plot" in sys.argv
    bar = "=" * 78
    print(bar); print("PART A -- Nino & Razavi (2019) Table 4 points, dP = 8-14 bar (unchanged)"); print(bar)
    run_part_a()
    print()
    means15 = None
    print(bar); print("PART D -- SPI Cd of injectors 1, 2, 5 (Waxman Fig. 15, digitised)"); print(bar)
    means15 = print_part_d()
    print()
    print(f"PART A sensitivity: same four points with Cd = {means15[2]:.3f} "
          "(mean digitised Fig. 15, injector 2 = the same 1.50 mm square-edge geometry)")
    errs, _ = run_part_a(cd=means15[2], verbose=False)
    s = summarize(errs)
    print("   ", _fmt(s))
    print()
    print(bar); print("PART B -- full injector-3 map, Waxman Fig. 13 (dP up to ~46 bar)"); print(bar)
    res_b = print_part_b()
    print()
    print(bar); print("PART C -- critical mass flow, Waxman Fig. 16"); print(bar)
    print_part_c(res_b)
    print()
    print(bar); print("PART E -- consistency checks of the digitised data"); print(bar)
    run_part_e(res_b)
    if not no_plot:
        try:
            print()
            print(f"Figure written to {make_figure(res_b)}")
        except ImportError:
            print("\nmatplotlib not installed -- figure skipped (pip install matplotlib).")


if __name__ == "__main__":
    main()
