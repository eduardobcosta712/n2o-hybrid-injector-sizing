"""
feed_line.py

Pressure drop and subcooling margin evolution along the N2O feed line,
from the tank to the injector inlet.

Implements:
    - Reynolds number
    - Darcy friction factor (laminar exact; turbulent via Swamee-Jain
      explicit approximation to the Colebrook equation)
    - Darcy-Weisbach major (friction) losses
    - Minor (fitting) losses via loss coefficients K
    - A single function that steps through a feed line made of segments
      and fittings, tracking pressure and subcooling margin along the way

Assumption (declared per README.md "Scope and known limitations"):
the feed line is treated as adiabatic — no heat exchange with the
environment, so fluid temperature is assumed constant along the line.
Only pressure changes; see Section 1.6 of the theory docs for why this
is a reasonable approximation for fast flows.

Two-phase pressure drop (Priority 2 of roadmap):
After the onset of flashing, the fluid is a liquid-vapour mixture. The
model uses the Homogeneous Equilibrium Model (HEM) for the two-phase
region: both phases share the same velocity and temperature (consistent
with the HEM used in injector_two_phase.py). The mixture properties are:
    rho_mix = 1 / ((1-x)/rho_l + x/rho_v)   [HEM density]
    mu_mix  = (1-x)*mu_l + x*mu_v             [McAdams mixing rule]
    x(s)    = (h_l(T_tank) - h_l(T_sat(P(s)))) / h_fg(T_sat(P(s)))
The vapour quality x is updated at the start of each segment in the
two-phase region, using the local pressure at that segment. This is
more accurate than a single end-of-line flash because it accounts for
the pressure evolution along the two-phase zone.

Units: SI throughout (Pa, kg/m^3, m, m/s, Pa.s, kg/s).
"""

import math
from n2o_properties import (P_sat, T_sat, rho_liquid_sat, degree_of_subcooling,
                             h_liquid_sat, h_fg, nu_vapor_sat, mu_vapor_sat,
                             mu_mixture, MU_LIQUID_N2O)

M_N2O = 44.013  # kg/kmol, molar mass of N2O (for rho_v = M_N2O / nu_vapor_sat)

# Approximate dynamic viscosity of saturated liquid N2O near room
# temperature (Pa.s). Treated as a constant for this version of the model;
# in reality it varies (mildly) with temperature -- see docs/future_work.md
# for possible refinements.
MU_LIQUID_N2O = 1.5e-4  # Pa.s, order-of-magnitude reference value


def reynolds_number(rho, v, D, mu=MU_LIQUID_N2O):
    """
    Reynolds number, Re = rho * v * D / mu.

    Parameters
    ----------
    rho : float
        Fluid density, kg/m^3.
    v : float
        Mean flow velocity, m/s.
    D : float
        Pipe internal diameter, m.
    mu : float
        Dynamic viscosity, Pa.s.

    Returns
    -------
    float
        Dimensionless Reynolds number.
    """
    return rho * v * D / mu


def darcy_friction_factor(Re, roughness, D):
    """
    Darcy friction factor f, per Darcy-Weisbach.

    - Laminar (Re < 2300): exact solution f = 64 / Re.
    - Turbulent (Re > 4000): Swamee-Jain explicit approximation to the
      implicit Colebrook equation (typically <1% error vs. Colebrook,
      and does not require iteration).
    - Transitional (2300 <= Re <= 4000): flow is physically unstable and
      not well characterized by either formula. This model conservatively
      uses the turbulent (Swamee-Jain) value in this range, which is
      documented here as a deliberate simplification.

    Parameters
    ----------
    Re : float
        Reynolds number.
    roughness : float
        Absolute pipe wall roughness, m (e.g. ~1.5e-6 m for smooth
        stainless steel tubing).
    D : float
        Pipe internal diameter, m.

    Returns
    -------
    float
        Dimensionless Darcy friction factor.
    """
    if Re < 2300:
        return 64.0 / Re

    # Turbulent or transitional: Swamee-Jain approximation.
    relative_roughness = roughness / D
    denominator = math.log10(relative_roughness / 3.7 + 5.74 / Re ** 0.9)
    return 0.25 / denominator ** 2


def friction_pressure_drop(rho, v, f, L, D):
    """
    Darcy-Weisbach major (friction) pressure loss along a straight
    pipe segment of length L.

    delta_P = f * (L / D) * (rho * v^2 / 2)

    Returns
    -------
    float
        Pressure drop in Pa (always >= 0; represents a loss).
    """
    return f * (L / D) * (rho * v ** 2 / 2.0)


def fitting_pressure_drop(rho, v, K):
    """
    Minor (fitting) pressure loss, e.g. for a valve, elbow, or union.

    delta_P = K * (rho * v^2 / 2)

    Parameters
    ----------
    K : float
        Fitting loss coefficient (dimensionless, from reference tables --
        see references.md for typical values for valves/elbows).

    Returns
    -------
    float
        Pressure drop in Pa.
    """
    return K * (rho * v ** 2 / 2.0)


def velocity_from_mass_flow(m_dot, rho, D):
    """
    Mean flow velocity from mass flow rate, via m_dot = rho * A * v.

    Parameters
    ----------
    m_dot : float
        Mass flow rate, kg/s.
    rho : float
        Fluid density, kg/m^3.
    D : float
        Pipe internal diameter, m.

    Returns
    -------
    float
        Mean velocity, m/s.
    """
    A = math.pi * (D / 2.0) ** 2
    return m_dot / (rho * A)


def evaluate_feed_line(m_dot, T_tank, P_tank, segments, roughness=1.5e-6,
                        mu=MU_LIQUID_N2O):
    """
    Step through a feed line made of straight segments and fittings,
    tracking pressure and subcooling margin from the tank to the
    injector inlet.

    The fluid is assumed adiabatic (constant T = T_tank) along the line
    (see module docstring). Density is re-evaluated as the saturated
    liquid density at T_tank -- this is treated as constant along the
    line too, since the model assumes the liquid stays subcooled liquid
    until proven otherwise at each step (see return value below for how
    that is flagged).

    Parameters
    ----------
    m_dot : float
        N2O mass flow rate through the line, kg/s.
    T_tank : float
        N2O temperature at the tank (assumed constant along the line), K.
    P_tank : float
        N2O pressure at the tank exit, Pa.
    segments : list of dict
        Each dict describes one element of the line, in order from tank
        to injector inlet. Two kinds of entries:
          {"type": "pipe", "L": <length, m>, "D": <diameter, m>}
          {"type": "fitting", "D": <diameter, m>, "K": <loss coeff.>}
    roughness : float
        Absolute pipe wall roughness, m. Applied to all "pipe" segments.
    mu : float
        Dynamic viscosity of the liquid, Pa.s.

    Returns
    -------
    dict with keys:
        "P_final": pressure at the injector inlet, Pa
        "delta_T_sub_final": subcooling margin at the injector inlet, K
        "flashing_detected": True if pressure dropped below P_sat(T_tank)
                              at any point along the line
        "trace": list of dicts, one per segment, with the running P and
                 delta_T_sub after that segment -- useful for plotting
                 and for the interactive tool.
    """
    # Liquid-phase properties (constant along the single-phase region)
    rho_l   = rho_liquid_sat(T_tank)
    h_up    = h_liquid_sat(T_tank)   # upstream enthalpy, kJ/kmol (conserved)
    Psat_T  = P_sat(T_tank)          # saturation pressure at tank temperature

    P                = P_tank
    trace            = []
    flashing_detected = False
    x_current        = 0.0    # vapour quality; 0 until flashing onset

    for i, seg in enumerate(segments):
        D = seg["D"]

        # -----------------------------------------------------------------
        # Determine mixture properties at the START of this segment.
        # Before flashing onset: pure liquid (x=0, rho=rho_l, mu=mu_l).
        # After flashing onset: two-phase HEM (x>0, rho=rho_mix, mu=mu_mix).
        # x is updated segment-by-segment using the local pressure so the
        # density and viscosity reflect the actual mixture state here,
        # not just at the final inlet point.
        # -----------------------------------------------------------------
        if P < Psat_T and P > 0:
            # Two-phase region: compute x isenthalpically at local P
            T_sat_local = T_sat(P)
            hfg_local   = h_fg(T_sat_local)
            hl_local    = h_liquid_sat(T_sat_local)
            if hfg_local > 0:
                x_current = max(0.0, min(1.0, (h_up - hl_local) / hfg_local))
            # HEM mixture density: 1 / ((1-x)/rho_l + x/rho_v)
            rho_v_local = M_N2O / nu_vapor_sat(T_sat_local)
            nu_mix      = (1.0 - x_current) / rho_l + x_current / rho_v_local
            rho_eff     = 1.0 / nu_mix
            # McAdams mixture viscosity: (1-x)*mu_l + x*mu_v
            mu_eff      = mu_mixture(x_current, T=T_sat_local, mu_l=mu)
        else:
            # Single-phase liquid: use pure liquid properties
            x_current = 0.0
            rho_eff   = rho_l
            mu_eff    = mu

        # -----------------------------------------------------------------
        # Pressure drop calculation (Darcy-Weisbach with mixture properties)
        # The formula is the same for single-phase and two-phase; only the
        # effective density and viscosity differ.
        # -----------------------------------------------------------------
        v = velocity_from_mass_flow(m_dot, rho_eff, D)

        if seg["type"] == "pipe":
            Re  = reynolds_number(rho_eff, v, D, mu_eff)
            f   = darcy_friction_factor(Re, roughness, D)
            dP  = friction_pressure_drop(rho_eff, v, f, seg["L"], D)
        elif seg["type"] == "fitting":
            dP  = fitting_pressure_drop(rho_eff, v, seg["K"])
        else:
            raise ValueError(f"Unknown segment type: {seg['type']!r}")

        P = P - dP

        # Flashing check (first time P crosses P_sat)
        if P < Psat_T:
            flashing_detected = True

        delta_T_sub = degree_of_subcooling(T_tank, P) if P > 0 else float("nan")

        trace.append({
            "segment_index":    i,
            "segment_type":     seg["type"],
            "velocity_m_s":     v,
            "pressure_drop_Pa": dP,
            "pressure_after_Pa": P,
            "delta_T_sub_K":    delta_T_sub,
            "x_quality":        x_current,    # vapour quality at segment start
            "rho_eff_kg_m3":    rho_eff,      # effective density used for dP
            "mu_eff_Pa_s":      mu_eff,       # effective viscosity used for dP
        })

    # Final x_inlet: quality at the injector inlet (end of last segment)
    # Re-compute at final P for accuracy rather than using x_current from
    # the last segment START -- the pressure may have dropped further.
    if flashing_detected and P > 0 and P < Psat_T:
        T_sat_final = T_sat(P)
        hfg_final   = h_fg(T_sat_final)
        hl_final    = h_liquid_sat(T_sat_final)
        x_inlet     = max(0.0, min(1.0, (h_up - hl_final) / hfg_final
                                   if hfg_final > 0 else 0.0))
    else:
        x_inlet = 0.0

    return {
        "P_final":           P,
        "delta_T_sub_final": degree_of_subcooling(T_tank, P) if P > 0 else float("nan"),
        "flashing_detected": flashing_detected,
        "x_inlet":           x_inlet,
        "trace":             trace,
    }


if __name__ == "__main__":
    # --- Example / sanity check: a simple feed line ---
    # 2 m of 8 mm ID tubing, one ball valve (K ~ 0.05), one 90-degree
    # elbow (K ~ 0.3), feeding N2O at 20 degC.
    #
    # Two cases are shown, deliberately, to illustrate the flashing
    # criterion behaving correctly at both extremes:
    #   Case A: tank exactly at saturation (zero initial margin) --
    #           any pressure drop at all should trigger flashing.
    #   Case B: tank with a modest 5 bar of subcooling margin, more
    #           representative of a real self-pressurized tank.
    T_tank = 293.15  # K (20 degC)
    m_dot = 0.5  # kg/s, example oxidizer mass flow rate

    segments = [
        {"type": "pipe", "L": 1.0, "D": 0.008},
        {"type": "fitting", "D": 0.008, "K": 0.05},   # ball valve
        {"type": "pipe", "L": 1.0, "D": 0.008},
        {"type": "fitting", "D": 0.008, "K": 0.3},    # 90-degree elbow
    ]

    for label, P_tank in [
        ("Case A: tank exactly at saturation (0 K initial margin)", P_sat(T_tank)),
        ("Case B: tank with 5 bar of subcooling margin", P_sat(T_tank) + 5e5),
    ]:
        print("=" * 60)
        print(label)
        print("=" * 60)
        result = evaluate_feed_line(m_dot, T_tank, P_tank, segments)

        print(f"Tank temperature:      {T_tank:.2f} K ({T_tank - 273.15:.1f} degC)")
        print(f"Tank pressure:         {P_tank/1e5:.2f} bar")
        print(f"Mass flow rate:        {m_dot:.3f} kg/s")
        print("-" * 60)
        for seg in result["trace"]:
            print(f"  segment {seg['segment_index']} ({seg['segment_type']:7s}): "
                  f"v={seg['velocity_m_s']:.2f} m/s, "
                  f"dP={seg['pressure_drop_Pa']/1e5:.4f} bar, "
                  f"P_after={seg['pressure_after_Pa']/1e5:.3f} bar, "
                  f"delta_T_sub={seg['delta_T_sub_K']:.2f} K")
        print("-" * 60)
        print(f"Final pressure at injector inlet: {result['P_final']/1e5:.3f} bar")
        print(f"Final subcooling margin:          {result['delta_T_sub_final']:.2f} K")
        print(f"Flashing detected along the line: {result['flashing_detected']}")
        print()
