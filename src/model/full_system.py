"""
full_system.py

Orchestrates the full tank -> feed line -> injector path, solving for the
self-consistent operating point where feed-line losses and injector flow
are mutually consistent.

PREVIOUS BEHAVIOUR (one-pass):
    A fixed design mass flow was used to compute feed-line losses, giving
    an injector inlet pressure. The injector model then predicted a real
    mass flow from that pressure. If m_dot_real != m_dot_design, the
    feed-line losses used were evaluated at the wrong flow rate.

CURRENT BEHAVIOUR (coupled solver):
    The system is solved iteratively. At each iteration, the feed-line
    losses are evaluated at the current mass-flow estimate, and the
    injector model predicts a new mass flow from the resulting inlet
    pressure. This is repeated until the mass flow converges.

    The operating point (m_dot*, P_inlet*) satisfies both:
        P_inlet* = P_tank - dP_line(m_dot*)      [feed line]
        m_dot*   = m_dot_injector(P_inlet*)       [injector]

    This intersection is found via damped fixed-point iteration:
        m_dot_{k+1} = m_dot_k + alpha * (m_dot_injector(P_inlet_k) - m_dot_k)

    where alpha in (0, 1] is a damping factor that prevents oscillation.

CONVERGENCE:
    Convergence criterion: relative change in m_dot < tolerance.
    Default tolerance: 1e-4 (0.01%). Default max iterations: 50.
    Default damping: alpha = 0.5.

    Non-convergence is reported as a RuntimeError with the iteration
    history attached, rather than silently returning an incorrect result.

BACKWARD COMPATIBILITY:
    evaluate_full_system() retains the same signature. The m_dot_design
    parameter is now used as the initial guess for the iteration, not as
    the fixed evaluation flow rate. Results include a new
    "solver_info" key with convergence details.

Units: SI throughout (Pa, K, kg/m^3, m^2, kg/s).
"""

from n2o_properties import rho_liquid_sat, nu_vapor_sat, T_sat
from feed_line import evaluate_feed_line
from injector_spi import spi_mass_flow, spi_sufficient
from injector_two_phase import dyer_mass_flow, hem_mass_flow_two_phase_inlet

M_N2O = 44.013  # kg/kmol, molar mass of N2O


# ---------------------------------------------------------------------------
# Internal: single-pass injector evaluation at a given m_dot and inlet P
# ---------------------------------------------------------------------------

def _evaluate_injector(Cd, A_injector, T_tank, P_injector_inlet,
                        P_chamber, feed_line_result):
    """
    Given the feed-line result at a specific mass flow, evaluate the
    injector model and return (m_dot_injector, regime, injector_result).

    This is the inner step of the coupled solver loop -- it does not
    re-evaluate the feed line.

    Parameters
    ----------
    Cd : float
        Injector discharge coefficient.
    A_injector : float
        Total orifice area, m^2.
    T_tank : float
        Tank temperature, K.
    P_injector_inlet : float
        Pressure at the injector inlet (= feed line P_final), Pa.
    P_chamber : float
        Downstream chamber pressure, Pa.
    feed_line_result : dict
        Return value of evaluate_feed_line() at the current m_dot.

    Returns
    -------
    tuple: (m_dot, regime_str, injector_result_dict or None)
    """
    if feed_line_result["flashing_detected"]:
        # Two-phase inlet: use HEM with isenthalpic x_inlet from feed line.
        # Dyer is not applied here because P_upstream ~ P_sat implies
        # kappa -> inf, which would collapse Dyer to SPI -- physically wrong
        # for a two-phase inlet (see injector_two_phase.py docstring).
        x_inlet    = feed_line_result["x_inlet"]
        T_down     = T_sat(P_chamber)
        rho_l_down = rho_liquid_sat(T_down)
        rho_v_down = M_N2O / nu_vapor_sat(T_down)

        ir = hem_mass_flow_two_phase_inlet(
            Cd, A_injector, T_tank, x_inlet,
            P_injector_inlet, P_chamber,
            rho_l_down, rho_v_down)
        return ir["m_dot_HEM_2phase"], "HEM_two_phase_inlet", ir

    sufficient = spi_sufficient(P_injector_inlet, T_tank, P_chamber)

    if sufficient:
        # Single-phase throughout the orifice: pure Bernoulli.
        rho_l    = rho_liquid_sat(T_tank)
        delta_P  = P_injector_inlet - P_chamber
        m_dot    = spi_mass_flow(Cd, A_injector, rho_l, delta_P)
        return m_dot, "SPI", None

    # Two-phase inside orifice, liquid at inlet: Dyer/NHNE.
    rho_l_up   = rho_liquid_sat(T_tank)
    T_down     = T_sat(P_chamber)
    rho_l_down = rho_liquid_sat(T_down)
    rho_v_down = M_N2O / nu_vapor_sat(T_down)

    ir = dyer_mass_flow(
        Cd, A_injector, T_tank, P_injector_inlet, P_chamber,
        rho_l_up, rho_l_down, rho_v_down)
    return ir["m_dot_Dyer"], "Dyer", ir


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def evaluate_full_system(m_dot_design, T_tank, P_tank, segments,
                          Cd, A_injector, P_chamber,
                          roughness=1.5e-6,
                          tol=1e-4, max_iter=50, alpha=0.5):
    """
    Solve for the self-consistent operating point of the full
    tank -> feed line -> injector system using damped fixed-point iteration.

    Parameters
    ----------
    m_dot_design : float
        Initial mass-flow guess for the iteration, kg/s. Typically the
        target design mass flow. Does not affect the converged result,
        only convergence speed.
    T_tank : float
        N2O tank temperature, K (assumed constant along the line).
    P_tank : float
        N2O tank pressure, Pa.
    segments : list of dict
        Feed line geometry (see feed_line.py for format).
    Cd : float
        Injector discharge coefficient (dimensionless).
    A_injector : float
        Total injector orifice area, m^2.
    P_chamber : float
        Chamber pressure (injector downstream pressure), Pa.
    roughness : float
        Feed line absolute pipe wall roughness, m.
    tol : float
        Convergence tolerance on relative change in m_dot (dimensionless).
        Default 1e-4 (0.01%).
    max_iter : int
        Maximum number of iterations before raising RuntimeError.
    alpha : float
        Damping factor in (0, 1]. Lower values converge more slowly but
        more robustly. Default 0.5.

    Returns
    -------
    dict with keys:
        "feed_line_result"  : evaluate_feed_line() result at converged m_dot.
        "P_injector_inlet"  : converged injector inlet pressure, Pa.
        "spi_sufficient"    : bool (False if two-phase regime was used).
        "m_dot_real"        : converged mass flow, kg/s.
        "injector_result"   : injector model detail dict, or None (SPI case).
        "regime"            : "SPI", "Dyer", or "HEM_two_phase_inlet".
        "solver_info"       : dict with convergence details:
            "converged"     : bool
            "iterations"    : number of iterations taken
            "history"       : list of m_dot values per iteration, kg/s
            "final_rel_err" : relative change at convergence

    Raises
    ------
    RuntimeError
        If the solver does not converge within max_iter iterations.
        The exception message includes the iteration history.
    ValueError
        If alpha is not in (0, 1] or tol <= 0.
    """
    if not (0 < alpha <= 1):
        raise ValueError(
            f"alpha = {alpha} must be in (0, 1]. "
            "Lower values are more stable; 1.0 is undamped (fastest but may oscillate)."
        )
    if tol <= 0:
        raise ValueError(f"tol = {tol} must be positive.")

    # --- Iteration ---
    # Start from the design guess.
    m_dot = m_dot_design
    history = [m_dot]

    for iteration in range(max_iter):

        # Step 1: evaluate feed-line losses at current m_dot estimate.
        # This gives P_inlet as a function of m_dot.
        fl = evaluate_feed_line(m_dot, T_tank, P_tank, segments, roughness)
        P_inlet = fl["P_final"]

        # Step 2: evaluate injector model at this P_inlet.
        # This gives the injector's predicted m_dot as a function of P_inlet.
        m_dot_new, regime, injector_result = _evaluate_injector(
            Cd, A_injector, T_tank, P_inlet, P_chamber, fl)

        # Step 3: damped update.
        # Without damping (alpha=1): m_dot = m_dot_new directly.
        # With damping: move only a fraction alpha toward the new estimate.
        # This prevents oscillation when the two curves have similar slopes.
        m_dot_next = m_dot + alpha * (m_dot_new - m_dot)

        # Step 4: convergence check.
        # Use relative change to be scale-independent (works for both small
        # and large mass flows without changing the tolerance).
        rel_err = abs(m_dot_next - m_dot) / max(abs(m_dot), 1e-12)
        history.append(m_dot_next)

        m_dot = m_dot_next

        if rel_err < tol:
            # Converged. Run one final evaluation at the converged m_dot
            # so the returned feed_line_result and injector_result are
            # consistent with the converged flow -- not the second-to-last.
            fl_final = evaluate_feed_line(m_dot, T_tank, P_tank,
                                           segments, roughness)
            P_inlet_final = fl_final["P_final"]
            _, regime_final, ir_final = _evaluate_injector(
                Cd, A_injector, T_tank, P_inlet_final, P_chamber, fl_final)

            return {
                "feed_line_result": fl_final,
                "P_injector_inlet": P_inlet_final,
                "spi_sufficient":   regime_final == "SPI",
                "m_dot_real":       m_dot,
                "injector_result":  ir_final,
                "regime":           regime_final,
                "solver_info": {
                    "converged":     True,
                    "iterations":    iteration + 1,
                    "history":       history,
                    "final_rel_err": rel_err,
                },
            }

    # Non-convergence: fail loudly with diagnostics.
    raise RuntimeError(
        f"Coupled solver did not converge in {max_iter} iterations.\n"
        f"  Tolerance requested: {tol:.2e}\n"
        f"  Last relative error: {rel_err:.4e}\n"
        f"  Last m_dot estimate: {m_dot*1000:.3f} g/s\n"
        f"  Iteration history (g/s): "
        f"{[f'{m*1000:.2f}' for m in history]}\n"
        "Possible causes: tank pressure too low for any stable flow, "
        "chamber pressure above tank pressure, or alpha too large "
        "(try alpha=0.3)."
    )


if __name__ == "__main__":
    # -----------------------------------------------------------------------
    # Validation: compare one-pass result vs. coupled solver result for a
    # case where there IS a meaningful discrepancy (long line, high flow).
    # For the Waxman geometry (very short upstream chamber), the discrepancy
    # is negligible -- so we use a longer line to make the effect visible.
    # -----------------------------------------------------------------------
    import math

    T_tank    = 293.15   # K, 20 degC
    P_tank    = 58e5     # Pa, 58 bar
    P_chamber = 22e5     # Pa, 22 bar
    Cd        = 0.65
    # 6 holes of 1.5 mm -- Example 1 geometry from examples/
    A_injector = 6 * math.pi * (0.00075) ** 2

    # Longer line to make the coupling effect visible
    segments = [
        {"type": "pipe",    "L": 2.0, "D": 0.008},
        {"type": "fitting", "D": 0.008, "K": 0.05},
        {"type": "pipe",    "L": 1.5, "D": 0.008},
        {"type": "fitting", "D": 0.008, "K": 0.90},
    ]

    # One-pass result (old behaviour): use target m_dot to size the line
    m_dot_target = 0.368   # kg/s -- the Dyer prediction from Example 1
    from feed_line import evaluate_feed_line
    from injector_spi import spi_mass_flow, spi_sufficient
    from injector_two_phase import dyer_mass_flow
    from n2o_properties import rho_liquid_sat, nu_vapor_sat, T_sat as T_sat_f

    fl_onepass   = evaluate_feed_line(m_dot_target, T_tank, P_tank, segments)
    P_in_onepass = fl_onepass["P_final"]
    rho_l_up     = rho_liquid_sat(T_tank)
    T_down       = T_sat_f(P_chamber)
    rho_l_down   = rho_liquid_sat(T_down)
    rho_v_down   = M_N2O / nu_vapor_sat(T_down)
    ir_onepass   = dyer_mass_flow(Cd, A_injector, T_tank, P_in_onepass,
                                   P_chamber, rho_l_up, rho_l_down, rho_v_down)
    m_onepass    = ir_onepass["m_dot_Dyer"]

    # Coupled solver result (new behaviour)
    result = evaluate_full_system(
        m_dot_target, T_tank, P_tank, segments, Cd, A_injector, P_chamber)

    print("Coupled solver validation")
    print("=" * 60)
    print(f"Tank:          {P_tank/1e5:.1f} bar,  {T_tank-273.15:.1f} degC")
    print(f"Chamber:       {P_chamber/1e5:.1f} bar")
    print()
    print(f"ONE-PASS (old):")
    print(f"  Feed-line evaluated at m_dot = {m_dot_target*1000:.1f} g/s (design guess)")
    print(f"  P_inlet = {P_in_onepass/1e5:.3f} bar")
    print(f"  m_dot_real = {m_onepass*1000:.1f} g/s")
    print(f"  Inconsistency: line was evaluated at {m_dot_target*1000:.1f} g/s, "
          f"but injector gives {m_onepass*1000:.1f} g/s")
    print()
    print(f"COUPLED SOLVER (new):")
    si = result["solver_info"]
    print(f"  Converged in {si['iterations']} iterations "
          f"(rel_err = {si['final_rel_err']:.2e})")
    print(f"  P_inlet = {result['P_injector_inlet']/1e5:.3f} bar")
    print(f"  m_dot_real = {result['m_dot_real']*1000:.1f} g/s  (self-consistent)")
    print(f"  Regime: {result['regime']}")
    print()
    delta = abs(m_onepass - result["m_dot_real"]) / result["m_dot_real"] * 100
    print(f"  Difference one-pass vs coupled: {delta:.2f}%")
    print()
    print("  Iteration history (g/s):")
    for i, m in enumerate(si["history"]):
        print(f"    iter {i:2d}: {m*1000:.3f} g/s")
