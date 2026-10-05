"""
scan_threshold_continuity.py

Does the supercharge-gated kappa correction remove the Dyer <-> HEM discontinuity
at the flashing threshold, and with it the band of tank pressures where the
coupled solver does not converge (examples/example_03_flashing.md)?

*** EXPLORATORY -- needs the optional `kappa_correction` argument of
full_system.py (apply_kappa_option.py).  Changes no model file. ***

For the two line geometries of Example 3 that matter
    (A) original: needle valve K = 2.0, 8 mm line   (band at 53.75-55.00 bar)
    (B) ball valve K = 0.05, 8 mm line              (Step 1: fails at 53.5 bar)
the tank pressure is scanned from 52 to 58 bar (0.25 bar steps), with
    - the baseline model (no correction),
    - the correction at the LOOCV optimum (beta = 1.0, ref = 400 psi),
    - the supercharge-protected variant (beta = 2.5, ref = 90 psi),
and (optionally) a custom variant:  --beta X --ref-psi Y.

Per variant it reports: number of non-converging pressures, the largest
relative jump of the converged flow between adjacent pressures (this is the
discontinuity when the two neighbours are on different branches), and the
flows at the ends of the scan.  Finally the four Steps of Example 3 are
recomputed with each variant, because if the correction is adopted these are
the numbers that change.

PREDICTION to check (stated BEFORE running): the corrected variants converge
everywhere and the jump disappears; the corrected Dyer flows in the scan are
much lower than the baseline ones (supercharge 3-5 bar ~ 45-70 psi is deep in
the regime where the correction pulls Dyer towards HEM).

Usage:  python validation/kappa_correction/scan_threshold_continuity.py [--beta 1.0 --ref-psi 400]
Output also written to validation/kappa_correction/scan_threshold_continuity_results.txt
"""

import argparse
import os
import sys

# this file lives in <repo>/validation/kappa_correction/
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "model"))
sys.path.insert(0, os.path.join(_REPO_ROOT, "tests"))

from full_system import evaluate_full_system                       # noqa: E402
from generate_example_03 import (seg_list, T_TANK, P_CHAMBER, CD,   # noqa: E402
                                 A_TOTAL, M_GUESS, STEPS)

PSI_PA = 6894.757


def variants(args):
    v = [("baseline", None),
         ("LOOCV optimum (b=1.0, ref=400 psi)",
          {"beta": 1.0, "supercharge_ref_Pa": 400.0 * PSI_PA}),
         ("Part-A-protected (b=2.5, ref=90 psi)",
          {"beta": 2.5, "supercharge_ref_Pa": 90.0 * PSI_PA})]
    if args.beta is not None and args.ref_psi is not None:
        v.append((f"custom (b={args.beta}, ref={args.ref_psi:g} psi)",
                  {"beta": args.beta, "supercharge_ref_Pa": args.ref_psi * PSI_PA}))
    return v


def run(P_tank, segs, kc):
    try:
        r = evaluate_full_system(M_GUESS, T_TANK, P_tank, segs, CD, A_TOTAL,
                                 P_CHAMBER, kappa_correction=kc)
    except RuntimeError:
        return None
    return r


def scan(label, segs, vlist, out):
    pressures = [52.0 + 0.25 * i for i in range(int(round((58.0 - 52.0) / 0.25)) + 1)]
    out("=" * 96)
    out(f"Tank-pressure scan -- {label}")
    out("=" * 96)
    table = {name: [run(P * 1e5, segs, kc) for P in pressures] for name, kc in vlist}

    head = f"{'P [bar]':>8} " + "".join(f"| {name[:30]:<30} " for name, _ in vlist)
    out(head)
    for i, P in enumerate(pressures):
        row = f"{P:8.2f} "
        for name, _ in vlist:
            r = table[name][i]
            if r is None:
                cell = "NO CONVERGENCE"
            else:
                tag = {"Dyer": "Dyer", "HEM_two_phase_inlet": "HEM2ph", "SPI": "SPI"}[r["regime"]]
                cell = f"{r['m_dot_real']*1000:7.1f} g/s {tag:<7}"
            row += f"| {cell:<30} "
        out(row)
    out()
    out("Summary")
    for name, _ in vlist:
        rs = table[name]
        nc = [P for P, r in zip(pressures, rs) if r is None]
        jumps = []
        for i in range(len(pressures) - 1):
            a, b = rs[i], rs[i + 1]
            if a is not None and b is not None:
                jumps.append((abs(b["m_dot_real"] - a["m_dot_real"]) / a["m_dot_real"],
                              pressures[i], a["regime"], b["regime"]))
        mj = max(jumps) if jumps else None
        first = next((r for r in rs if r is not None), None)
        last = next((r for r in reversed(rs) if r is not None), None)
        out(f"  {name}")
        out(f"    non-converging: {len(nc)}/{len(pressures)}"
            + (f"  at {nc[0]:.2f}..{nc[-1]:.2f} bar" if nc else ""))
        if mj:
            out(f"    largest adjacent jump: {100*mj[0]:.1f} % between {mj[1]:.2f} and "
                f"{mj[1]+0.25:.2f} bar ({mj[2]} -> {mj[3]})")
        if first and last:
            out(f"    flow at the ends of the converged range: {first['m_dot_real']*1000:.1f} -> "
                f"{last['m_dot_real']*1000:.1f} g/s")
    out()


def steps(vlist, out):
    out("=" * 96)
    out("Example 3 steps recomputed with each variant")
    out("=" * 96)
    for label, P, segs in STEPS:
        out(label)
        for name, kc in vlist:
            r = run(P, segs, kc)
            if r is None:
                out(f"  {name:<40} NO CONVERGENCE")
                continue
            ir = r["injector_result"] or {}
            extra = ""
            if r["regime"] == "Dyer":
                hf = ir.get("m_dot_crit_HF")
                extra = (f" kappa={ir['kappa']:.3f}"
                         + (f"->{ir['kappa_corrected']:.3f}" if "kappa_corrected" in ir else "")
                         + (f" | HF ceiling {hf*1000:.1f}, choked={ir['choked']}" if hf else ""))
            out(f"  {name:<40} {r['m_dot_real']*1000:7.1f} g/s  {r['regime']:<20}"
                f" flashing={r['feed_line_result']['flashing_detected']}{extra}")
        out()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta", type=float)
    ap.add_argument("--ref-psi", type=float)
    args = ap.parse_args()
    vlist = variants(args)
    lines = []

    def out(s=""):
        print(s)
        lines.append(s)

    scan("(A) original: needle valve K=2.0, 8 mm line", seg_list(8.0, 2.0), vlist, out)
    scan("(B) ball valve K=0.05, 8 mm line", seg_list(8.0, 0.05), vlist, out)
    steps(vlist, out)

    path = os.path.join(_REPO_ROOT, "validation", "kappa_correction", "scan_threshold_continuity_results.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"(output also written to {path})")


if __name__ == "__main__":
    main()
