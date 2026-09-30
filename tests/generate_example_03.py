"""
generate_example_03.py

Regenerates every number quoted in examples/example_03_flashing.md with the
real CoolProp backend, so the example can be closed without any hand-copied
or legacy-backend figures.

Run from anywhere:
    python examples/generate_example_03.py

It prints, for each step of the correction table:
  - feed-line trace (P after each segment, vapour quality x, effective density)
  - flashing yes/no, inlet vapour quality
  - regime, converged mass flow, solver iterations
  - for the Dyer regime: kappa, exit quality, SPI/HEM/Dyer flows,
    Henry-Fauske ceiling and the `choked` flag
  - for the HEM two-phase-inlet regime: exit quality
  - the naive SPI reference flow (full tank-to-chamber drop)

If the coupled solver does not converge at a step (the Dyer<->HEM
discontinuity, docs/future_work.md Priority 2), the failure is reported
instead of aborting, together with a scan that locates the offending
tank-pressure band, so the choice of operating points can be justified.

Paste the output back into the conversation (or straight into the .md).
"""

import math
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "model"))

from full_system import evaluate_full_system
from n2o_properties import P_sat, rho_liquid_sat

T_TANK = 22.0 + 273.15
P_CHAMBER = 18e5
CD = 0.65
A_TOTAL = 4 * math.pi * (0.9e-3) ** 2          # 4 holes of 1.8 mm
M_GUESS = 0.35                                  # initial guess only, kg/s


def seg_list(pipe_D_mm, valve_K):
    D = pipe_D_mm / 1000.0
    return [
        {"type": "pipe", "L": 2.5, "D": D},
        {"type": "fitting", "D": D, "K": valve_K},   # needle valve or ball valve
        {"type": "fitting", "D": D, "K": 0.9},       # 90 deg elbow
    ]


STEPS = [
    ("Step 0: original (53.5 bar, needle valve K=2.0, 8 mm)", 53.5e5, seg_list(8.0, 2.0)),
    ("Step 1: needle valve -> ball valve K=0.05 (53.5 bar, 8 mm)", 53.5e5, seg_list(8.0, 0.05)),
    ("Step 2: raise tank pressure to 56 bar (ball valve, 8 mm)", 56.0e5, seg_list(8.0, 0.05)),
    ("Step 3: pipe ID 8 -> 12 mm (56 bar, ball valve)", 56.0e5, seg_list(12.0, 0.05)),
]


def run_step(label, P_tank, segs):
    print("=" * 74)
    print(label)
    print("=" * 74)
    print(f"P_sat(22 degC) = {P_sat(T_TANK)/1e5:.3f} bar | tank {P_tank/1e5:.2f} bar | "
          f"margin at tank exit {(P_tank - P_sat(T_TANK))/1e5:+.3f} bar")
    rho_l = rho_liquid_sat(T_TANK)
    m_spi_naive = CD * A_TOTAL * math.sqrt(2 * rho_l * (P_tank - P_CHAMBER))
    print(f"Naive SPI (full tank-to-chamber drop): {m_spi_naive*1000:.1f} g/s")
    try:
        r = evaluate_full_system(M_GUESS, T_TANK, P_tank, segs, CD, A_TOTAL, P_CHAMBER)
    except RuntimeError as exc:
        print("SOLVER DID NOT CONVERGE at this operating point:")
        print("  " + str(exc).replace("\n", "\n  "))
        return None

    fl = r["feed_line_result"]
    si = r["solver_info"]
    print(f"Converged in {si['iterations']} iterations (rel. err {si['final_rel_err']:.2e})")
    print("Feed-line trace (at the converged flow):")
    for s in fl["trace"]:
        print(f"  seg {s['segment_index']} {s['segment_type']:7s}: "
              f"P_after = {s['pressure_after_Pa']/1e5:7.3f} bar | "
              f"dP = {s['pressure_drop_Pa']/1e5:6.3f} bar | "
              f"x(start) = {s['x_quality']:.4f} | "
              f"rho_eff = {s['rho_eff_kg_m3']:6.1f} kg/m3")
    print(f"Flashing detected: {fl['flashing_detected']} | "
          f"x_inlet = {fl['x_inlet']:.4f} | "
          f"P_injector_inlet = {r['P_injector_inlet']/1e5:.3f} bar")
    print(f"Regime: {r['regime']} | REAL MASS FLOW = {r['m_dot_real']*1000:.1f} g/s")
    ir = r["injector_result"]
    if r["regime"] == "Dyer":
        print(f"  kappa = {ir['kappa']:.3f} | exit quality x = {ir['x_exit']:.3f}")
        print(f"  SPI = {ir['m_dot_SPI']*1000:.1f} | HEM = {ir['m_dot_HEM']*1000:.1f} | "
              f"Dyer = {ir['m_dot_Dyer']*1000:.1f} g/s")
        if ir["m_dot_crit_HF"] is not None:
            print(f"  Henry-Fauske ceiling = {ir['m_dot_crit_HF']*1000:.1f} g/s | "
                  f"choked = {ir['choked']}")
        else:
            print(f"  Henry-Fauske ceiling unavailable: {ir['HF_unavailable_reason']}")
        dP_inj = (r["P_injector_inlet"] - P_CHAMBER) / 1e5
        print(f"  injector dP = {dP_inj:.2f} bar ({100*dP_inj/(P_CHAMBER/1e5):.0f} % of P_chamber)")
    elif r["regime"] == "HEM_two_phase_inlet":
        print(f"  exit quality x = {ir['x_exit']:.3f} | x_inlet = {ir['x_inlet']:.4f}")
    return r


def scan_convergence_band():
    """Locate tank pressures where the solver fails for the ORIGINAL geometry."""
    print("=" * 74)
    print("Convergence scan (original geometry: needle valve, 8 mm)")
    print("=" * 74)
    segs = seg_list(8.0, 2.0)
    P = 52.0
    while P <= 58.0 + 1e-9:
        try:
            r = evaluate_full_system(M_GUESS, T_TANK, P * 1e5, segs, CD, A_TOTAL, P_CHAMBER)
            print(f"  {P:5.2f} bar: OK  -> {r['m_dot_real']*1000:7.1f} g/s  ({r['regime']})")
        except RuntimeError:
            print(f"  {P:5.2f} bar: NO CONVERGENCE")
        P = round(P + 0.25, 2)


if __name__ == "__main__":
    for lab, P, segs in STEPS:
        run_step(lab, P, segs)
    scan_convergence_band()
    print()
    print("Compare with the values currently in examples/example_03_flashing.md:")
    print("  Step 0 = 193.9 g/s (HEM two-phase inlet), Step 2 = 344.4 g/s (Dyer)")
    print("  (headline figures already CoolProp-confirmed; everything else was legacy.)")
