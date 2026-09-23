import math
import sys
import os

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "model"))

from full_system import evaluate_full_system
from injector_two_phase import (hem_critical_flow, hem_critical_flow_isentropic,
                                 henry_fauske_critical_flow)
from n2o_properties import P_sat, rho_liquid_sat

T1   = 280.0
P1   = 4.36e6
Cd   = 0.65
D    = 0.0015
A    = math.pi * (D / 2.0) ** 2

SEGMENTS = [{"type": "pipe", "L": 0.05, "D": 0.0254}]

CASES = [
    ("Pre-critical",  0.84, 44.0),
    ("Critical",      0.98, 46.5),
    ("Post-critical1",1.09, 47.5),
    ("Post-critical2",1.37, 48.0),
]


def pct_error(pred, exp):
    return 100.0 * (pred - exp) / exp


def run_validation():
    Psat = P_sat(T1)
    print(f"P_sat(T1) model = {Psat/1e5:.3f} bar  (Nino&Razavi: {(P1-0.62e6)/1e5:.2f} bar,  "
          f"model err = {100*(Psat-(P1-0.62e6))/(P1-0.62e6):+.1f}%)")
    print(f"Supercharge: model = {(P1-Psat)/1e5:.2f} bar  |  Nino&Razavi = 6.20 bar")
    print()

    results = []
    for label, dP_MPa, m_exp in CASES:
        P_chamber = P1 - dP_MPa * 1e6
        r = evaluate_full_system(m_exp * 0.9e-3, T1, P1, SEGMENTS, Cd, A, P_chamber)
        m = r["m_dot_real"]
        si = r["solver_info"]
        err = pct_error(m, m_exp * 1e-3)
        ir = r.get("injector_result") or {}
        choked_flag = "CHOKED" if ir.get("choked") else "-"
        results.append((label, dP_MPa, m_exp, m, err, r))
        print(f"  {label:<16} dP={dP_MPa*10:6.2f} bar  m_exp={m_exp:5.1f} g/s  "
              f"m_dot={m*1000:7.3f} g/s  err={err:+.2f}%  "
              f"P_in={r['P_injector_inlet']/1e5:.3f} bar  regime={r['regime']}  "
              f"iters={si['iterations']}  {choked_flag}")

    errs = [r[4] for r in results]
    print()
    print(f"Mean error: {sum(errs)/len(errs):+.2f}%")
    print(f"MAPE:       {sum(abs(e) for e in errs)/len(errs):.2f}%")
    print(f"Within +/-5%: {sum(1 for e in errs if abs(e)<=5)}/{len(errs)}")
    print(f"Within +/-10%: {sum(1 for e in errs if abs(e)<=10)}/{len(errs)}")

    print()
    crit = hem_critical_flow(Cd, A, T1, P1)
    crit_s = hem_critical_flow_isentropic(Cd, A, T1, P1)
    print(f"HEM isenthalpic: m_dot_crit={crit['m_dot_crit']*1000:.2f} g/s @ "
          f"P2={crit['P2_crit']/1e5:.2f} bar x={crit['x_crit']:.4f}")
    print(f"HEM isentropic:  m_dot_crit={crit_s['m_dot_crit']*1000:.2f} g/s @ "
          f"P2={crit_s['P2_crit']/1e5:.2f} bar x={crit_s['x_crit']:.4f}")

    hf = henry_fauske_critical_flow(Cd, A, T1, P1)
    print(f"Henry-Fauske:    m_dot_crit={hf['m_dot_crit']*1000:.2f} g/s @ "
          f"P2={hf['P2_crit']/1e5:.2f} bar x_E={hf['x_crit']:.4f} N={hf['N']:.3f}")


if __name__ == "__main__":
    run_validation()
