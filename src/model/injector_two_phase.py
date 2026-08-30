"""
injector_two_phase.py

Two-phase injector models: HEM (Homogeneous Equilibrium Model) and Dyer
(NHNE). See docs/03_two_phase_flow.md for the full theoretical derivation
and physical motivation.

HEM assumes the liquid-vapor mixture inside the orifice is in full
thermodynamic equilibrium at every point, and treats it as a single
effective fluid with a mixture density that depends on the local vapor
quality x. Dyer corrects HEM's "instantaneous equilibrium" assumption
(unrealistic in a short orifice, where the fluid's residence time may be
too short for true equilibrium to establish) by blending it with the pure
SPI prediction (injector_spi.py), weighted by how close the upstream
pressure already is to saturation.

    m_dot_HEM  = Cd * A * sqrt(2 * rho_HEM * delta_P)
    m_dot_Dyer = (kappa / (1 + kappa)) * m_dot_SPI + (1 / (1 + kappa)) * m_dot_HEM

Units: SI throughout (Pa, K, kg/m^3, m^2, kg/s), except vapor quality x
and the Dyer weighting parameter kappa, which are dimensionless.
"""

import math
from n2o_properties import P_sat, T_sat, h_liquid_sat, h_fg
from injector_spi import spi_mass_flow


def vapor_quality_isenthalpic(h_upstream, T_downstream):
    """
    Vapor quality x at the orifice exit, assuming an isenthalpic
    (constant-enthalpy) process across the orifice -- no heat exchange
    with the surroundings, per the same fast-flow justification used for
    the feed line (Section 1.6) -- and full thermodynamic equilibrium at
    the exit (the HEM assumption): the exiting mixture sits exactly on
    the saturation curve at T_downstream = T_sat(P_downstream).

    Derived from h_upstream = h_l(T_downstream) + x * h_fg(T_downstream):

        x = (h_upstream - h_l(T_downstream)) / h_fg(T_downstream)

    Physically: h_upstream is the energy entering the orifice; h_l(T_downstream)
    is what a pure saturated liquid at the exit pressure would carry; the
    difference is the excess energy that must go into vaporizing a
    fraction x of the mass, at a "price" of h_fg(T_downstream) per unit mass.

    Parameters
    ----------
    h_upstream : float
        Upstream fluid molar enthalpy entering the orifice, kJ/kmol.
        For liquid entering subcooled or saturated, this is h_liquid_sat
        evaluated at the upstream temperature (see dyer_mass_flow below).
    T_downstream : float
        Saturation temperature corresponding to the downstream pressure,
        T_sat(P_downstream), K.

    Returns
    -------
    float
        Vapor quality x (dimensionless). Clamped to [0, 1]: values
        slightly outside this range can occur from the linear-interpolation
        and isenthalpic idealizations right at the boundary of validity,
        and are physically meaningless outside [0, 1] (Section 1.7).
    """
    x = (h_upstream - h_liquid_sat(T_downstream)) / h_fg(T_downstream)
    return max(0.0, min(1.0, x))


def hem_mixture_density(x, rho_l, rho_v):
    """
    HEM mixture density, from the mass-weighted average of the two
    phases' specific volumes (Section 4.1 of the theory docs):

        nu_mix = (1 - x) / rho_l + x / rho_v
        rho_HEM = 1 / nu_mix

    Parameters
    ----------
    x : float
        Vapor quality (dimensionless, in [0, 1]).
    rho_l : float
        Saturated liquid density at the mixture's local temperature,
        kg/m^3.
    rho_v : float
        Saturated vapor density at the mixture's local temperature,
        kg/m^3.

    Returns
    -------
    float
        Mixture density, kg/m^3.
    """
    nu_mix = (1.0 - x) / rho_l + x / rho_v
    return 1.0 / nu_mix


def hem_mass_flow(Cd, A, T_upstream, P_upstream, P_downstream,
                   rho_l_upstream, rho_l_downstream, rho_v_downstream):
    """
    HEM mass flow rate through the orifice.

    Combines vapor_quality_isenthalpic and hem_mixture_density with the
    same Bernoulli-derived orifice equation used throughout this project
    (Section 2.2), but with the mixture density rho_HEM in place of the
    pure-liquid density:

        m_dot_HEM = Cd * A * sqrt(2 * rho_HEM * delta_P)

    Parameters
    ----------
    Cd : float
        Discharge coefficient (dimensionless).
    A : float
        Total orifice area, m^2.
    T_upstream : float
        Upstream temperature, K.
    P_upstream : float
        Upstream pressure, Pa.
    P_downstream : float
        Downstream pressure, Pa.
    rho_l_upstream : float
        Saturated liquid density at T_upstream, kg/m^3 -- used only to
        evaluate the upstream enthalpy via h_liquid_sat(T_upstream)
        (liquid entering the orifice is assumed saturated or subcooled
        liquid, per Section 1.4; its enthalpy is well approximated by
        the saturated-liquid value at its own temperature).
    rho_l_downstream : float
        Saturated liquid density at T_downstream = T_sat(P_downstream),
        kg/m^3.
    rho_v_downstream : float
        Saturated vapor density at T_downstream, kg/m^3 (from
        nu_vapor_sat via rho = M_N2O / nu_vapor_sat(T_downstream)).

    Returns
    -------
    dict
        "m_dot_HEM": HEM-predicted mass flow rate, kg/s.
        "x_exit": vapor quality at the orifice exit (dimensionless).
        "rho_HEM": mixture density at the orifice exit, kg/m^3.
        "T_downstream": saturation temperature at P_downstream, K.
    """
    delta_P = P_upstream - P_downstream
    if delta_P < 0:
        raise ValueError(
            f"delta_P = {delta_P:.1f} Pa is negative -- this would imply "
            "reverse flow through the orifice. Check that P_upstream and "
            "P_downstream were passed in the correct order."
        )

    T_downstream = T_sat(P_downstream)
    h_upstream = h_liquid_sat(T_upstream)
    x_exit = vapor_quality_isenthalpic(h_upstream, T_downstream)
    rho_HEM = hem_mixture_density(x_exit, rho_l_downstream, rho_v_downstream)

    m_dot_HEM = Cd * A * math.sqrt(2.0 * rho_HEM * delta_P)

    return {
        "m_dot_HEM": m_dot_HEM,
        "x_exit": x_exit,
        "rho_HEM": rho_HEM,
        "T_downstream": T_downstream,
    }


def hem_mass_flow_two_phase_inlet(Cd, A, T_tank, x_inlet,
                                   P_upstream, P_downstream,
                                   rho_l_downstream, rho_v_downstream):
    """
    HEM mass flow rate when the fluid arrives at the injector inlet already
    partially vaporised (x_inlet > 0), as happens when flashing occurs in
    the feed line before the injector.

    This extends hem_mass_flow() by replacing the pure-liquid upstream
    enthalpy with the enthalpy of a two-phase mixture at the feed line exit:

        h_upstream = h_l(T_tank) + x_inlet * h_fg(T_tank)

    Everything downstream of the inlet (x_exit, rho_HEM, m_dot_HEM) is
    computed identically to hem_mass_flow().

    Note on the Dyer model in this regime: when x_inlet > 0 the fluid is
    at saturation at the orifice inlet, so P_upstream ~ P_sat(T_tank) and
    the Dyer kappa denominator (P_sat - P_downstream) -> 0, giving
    kappa -> infinity and w_SPI -> 1 (Dyer collapses to SPI). That is
    physically wrong: a saturated two-phase inlet means there is NO
    non-equilibrium margin -- HEM is the appropriate model, not SPI.
    This function therefore returns HEM only, and full_system.py uses it
    directly without applying the Dyer blend.

    Parameters
    ----------
    Cd : float
        Discharge coefficient (dimensionless).
    A : float
        Total orifice area, m^2.
    T_tank : float
        Tank temperature (= temperature at feed line inlet), K.
        Used to evaluate h_l(T_tank) and h_fg(T_tank).
    x_inlet : float
        Vapour quality at the injector inlet (from feed_line.py), [0, 1].
    P_upstream : float
        Pressure at the injector inlet (= feed line exit pressure), Pa.
    P_downstream : float
        Downstream pressure (chamber pressure), Pa.
    rho_l_downstream : float
        Saturated liquid density at T_sat(P_downstream), kg/m^3.
    rho_v_downstream : float
        Saturated vapour density at T_sat(P_downstream), kg/m^3.

    Returns
    -------
    dict
        "m_dot_HEM_2phase": HEM mass flow rate with two-phase inlet, kg/s.
        "x_exit": vapour quality at the orifice exit, dimensionless.
        "rho_HEM": effective mixture density at exit, kg/m^3.
        "x_inlet": vapour quality at the inlet (passed through for reporting).
    """
    delta_P = P_upstream - P_downstream
    if delta_P < 0:
        raise ValueError(
            f"delta_P = {delta_P:.1f} Pa is negative -- P_upstream must "
            "exceed P_downstream."
        )

    # Two-phase upstream enthalpy: liquid enthalpy + vapour fraction contribution
    h_up = h_liquid_sat(T_tank) + x_inlet * h_fg(T_tank)

    # Downstream state: T_sat(P_downstream)
    T_downstream = T_sat(P_downstream)

    # Exit quality: isenthalpic expansion from h_up to downstream saturation
    x_exit = vapor_quality_isenthalpic(h_up, T_downstream)

    # HEM mixture density at exit
    rho_HEM = hem_mixture_density(x_exit, rho_l_downstream, rho_v_downstream)

    m_dot = Cd * A * math.sqrt(2.0 * rho_HEM * delta_P)

    return {
        "m_dot_HEM_2phase": m_dot,
        "x_exit": x_exit,
        "rho_HEM": rho_HEM,
        "x_inlet": x_inlet,
    }


def dyer_non_equilibrium_parameter(P_upstream, T_upstream, P_downstream):
    """
    Dyer's non-equilibrium weighting parameter, kappa (Section 3.4):

    kappa = sqrt[(P_upstream - P_downstream) / (P_sat(T_upstream) - P_downstream)]

The numerator is the total pressure drop across the orifice; the
denominator is the pressure drop required to reach saturation. 

A large kappa means the inlet was comfortably subcooled (large margin 
above P_sat), so the fluid spends little time in a flashing state inside 
the orifice. Thus, the flow behaves closer to the single-phase SPI limit 
("no time to vaporize"). 

A small kappa (close to 1) means the inlet was already close to saturation 
(little subcooling margin), so vaporization occurs rapidly and the flow 
behaves closer to the full-equilibrium HEM limit.

    Parameters
    ----------
    P_upstream : float
        Pressure just before the orifice, Pa.
    T_upstream : float
        Temperature just before the orifice, K.
    P_downstream : float
        Pressure just after the orifice, Pa.

    Returns
    -------
    float
        Dyer's kappa parameter (dimensionless, >= 0).
    """
    P_sat_upstream = P_sat(T_upstream)
    if P_upstream <= P_sat_upstream:
        raise ValueError(
            f"P_upstream = {P_upstream/1e5:.2f} bar is at or below "
            f"P_sat(T_upstream) = {P_sat_upstream/1e5:.2f} bar -- the fluid "
            "is already two-phase (or exactly saturated) BEFORE reaching "
            "the orifice inlet. The Dyer model as implemented here assumes "
            "liquid (saturated or subcooled) at the orifice inlet, with "
            "vaporization occurring inside the orifice (Section 3.1); a "
            "two-phase feed line is a separate problem (see feed_line.py's "
            "flashing_detected flag)."
        )
    return math.sqrt((P_upstream - P_downstream) / (P_sat_upstream - P_downstream))


def dyer_mass_flow(Cd, A, T_upstream, P_upstream, P_downstream,
                    rho_l_upstream, rho_l_downstream, rho_v_downstream):
    """
    Dyer (NHNE) mass flow rate: the weighted combination of the SPI and
    HEM limits (Section 3.4), using dyer_non_equilibrium_parameter as the
    weight:

        m_dot_Dyer = (kappa / (1 + kappa)) * m_dot_SPI + (1 / (1 + kappa)) * m_dot_HEM

    Physical interpretation of the weights (Waxman 2013, p.6):
        kappa ~ tau_bubble / tau_residence.
        Large kappa: bubbles grow slowly relative to fluid residence time
            -> less equilibrium -> more weight on SPI (the "no vaporisation"
            limit). w_SPI = kappa/(1+kappa) increases with kappa. Correct.
        Small kappa: full bubble growth, near-equilibrium -> weight on HEM.

    Reference: Waxman (2013) Eq. (9); Solomon (2011); corrects the sign
    error in the original Dyer et al. (2007) formulation.

    This is the reference model adopted for injector sizing in this
    project (Section 3.4/3.5).

    Parameters
    ----------
    Cd, A, T_upstream, P_upstream, P_downstream, rho_l_upstream,
    rho_l_downstream, rho_v_downstream : see hem_mass_flow docstring;
        identical roles here.

    Returns
    -------
    dict
        "m_dot_Dyer": Dyer-predicted mass flow rate, kg/s.
        "m_dot_SPI": SPI-only prediction at the same operating point, kg/s
                     (the "no time to vaporize" limit, for comparison).
        "m_dot_HEM": HEM-only prediction, kg/s (the full-equilibrium
                     limit, for comparison).
        "kappa": Dyer's non-equilibrium weighting parameter.
        "x_exit": vapor quality at the orifice exit, from the HEM
                  sub-calculation (dimensionless).
    """
    kappa = dyer_non_equilibrium_parameter(P_upstream, T_upstream, P_downstream)

    delta_P = P_upstream - P_downstream
    m_dot_SPI = spi_mass_flow(Cd, A, rho_l_upstream, delta_P)

    hem_result = hem_mass_flow(Cd, A, T_upstream, P_upstream, P_downstream,
                                rho_l_upstream, rho_l_downstream,
                                rho_v_downstream)
    m_dot_HEM = hem_result["m_dot_HEM"]

    # Correct NHNE formula (Waxman 2013 Eq.9 / Solomon 2011):
    # large kappa -> more weight on SPI (less equilibrium); small kappa -> HEM.
    m_dot_Dyer = (kappa / (1.0 + kappa)) * m_dot_SPI + (1.0 / (1.0 + kappa)) * m_dot_HEM

    return {
        "m_dot_Dyer": m_dot_Dyer,
        "m_dot_SPI": m_dot_SPI,
        "m_dot_HEM": m_dot_HEM,
        "kappa": kappa,
        "x_exit": hem_result["x_exit"],
    }


if __name__ == "__main__":
    # --- Validation case ---
    #
    #   The injector_spi.py example (50 bar upstream / 20 bar downstream at
    #   20 degC) is NOT usable here: at 20 degC, P_sat is ~51.4 bar, so that
    #   example's 50 bar upstream is already below saturation -- valid for
    #   demonstrating spi_sufficient() == False, but outside this module's
    #   domain (Dyer assumes liquid AT the orifice inlet; see
    #   dyer_non_equilibrium_parameter's docstring). Confirmed directly: an
    #   initial attempt to reuse that exact example raised the expected
    #   ValueError here, which is itself useful confirmation that the
    #   domain check is doing its job.
    #
    #   Instead: N2O at 20 degC upstream, with the tank/run-line pressure
    #   raised to 55 bar (a modest ~3.6 bar of subcooling margin above
    #   P_sat(20 degC) ~= 51.4 bar -- comparable to feed_line.py's own
    #   Case B), and 20 bar downstream. This keeps the inlet liquid (valid
    #   for Dyer) while still crossing saturation inside the orifice
    #   (so SPI is still expected to over-predict, per Section 2.4/3.3).

    from n2o_properties import rho_liquid_sat, nu_vapor_sat

    M_N2O = 44.013  # kg/kmol, molar mass of N2O (n2o_properties.py convention)

    T_upstream = 293.15   # K, 20 degC
    P_upstream = 55e5     # Pa, 55 bar (subcooled: P_sat(20 degC) ~= 51.4 bar)
    P_downstream = 20e5   # Pa, 20 bar chamber pressure
    Cd = 0.65
    A = 3.79e-6            # m^2, same example orifice area as injector_spi.py

    rho_l_upstream = rho_liquid_sat(T_upstream)

    T_downstream = T_sat(P_downstream)
    rho_l_downstream = rho_liquid_sat(T_downstream)
    rho_v_downstream = M_N2O / nu_vapor_sat(T_downstream)

    result = dyer_mass_flow(Cd, A, T_upstream, P_upstream, P_downstream,
                             rho_l_upstream, rho_l_downstream, rho_v_downstream)

    print("Two-phase injector model (HEM + Dyer) -- example evaluation")
    print("-" * 60)
    print(f"Upstream temperature:     {T_upstream:.2f} K ({T_upstream-273.15:.1f} degC)")
    print(f"Upstream pressure:        {P_upstream/1e5:.1f} bar")
    print(f"Downstream pressure:      {P_downstream/1e5:.1f} bar")
    print(f"P_sat at T_upstream:      {P_sat(T_upstream)/1e5:.2f} bar")
    print(f"T_downstream (=T_sat(P_downstream)): {T_downstream:.2f} K "
          f"({T_downstream-273.15:.1f} degC)")
    print(f"rho_l upstream:           {rho_l_upstream:.1f} kg/m^3")
    print(f"rho_l downstream:         {rho_l_downstream:.1f} kg/m^3")
    print(f"rho_v downstream:         {rho_v_downstream:.2f} kg/m^3")
    print("-" * 60)
    print(f"Vapor quality at exit, x: {result['x_exit']:.4f}")
    print(f"Dyer kappa:               {result['kappa']:.3f}")
    print("-" * 60)
    print(f"SPI-predicted mass flow:  {result['m_dot_SPI']*1000:.1f} g/s")
    print(f"HEM-predicted mass flow:  {result['m_dot_HEM']*1000:.1f} g/s")
    print(f"Dyer-predicted mass flow: {result['m_dot_Dyer']*1000:.1f} g/s")
    print("-" * 60)
    reduction_pct = 100.0 * (1.0 - result['m_dot_Dyer'] / result['m_dot_SPI'])
    print(f"Dyer vs SPI reduction:    {reduction_pct:.1f}%")
    print("  -> As expected (Section 3.3), the two-phase-aware models predict")
    print("     a lower mass flow than pure SPI, which over-predicts by")
    print("     assuming single-phase liquid throughout the orifice.")
