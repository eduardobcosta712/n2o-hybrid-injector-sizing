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
    m_dot_Dyer = m_dot_SPI / (1 + kappa) + (kappa / (1 + kappa)) * m_dot_HEM

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


def dyer_non_equilibrium_parameter(P_upstream, T_upstream, P_downstream):
    """
    Dyer's non-equilibrium weighting parameter, kappa (Section 3.4):

        kappa = sqrt[(P_upstream - P_downstream) / (P_sat(T_upstream) - P_downstream)]

    The numerator is the total pressure drop across the orifice; the
    denominator is how much subcooling margin (in pressure terms) the
    fluid had at the orifice inlet before the flow even begins. A large
    kappa means the inlet was already close to saturation (little margin
    to lose), so the flow behaves closer to the full-equilibrium HEM
    limit; a small kappa means the inlet was comfortably subcooled, so
    the flow behaves closer to the "no time to vaporize" SPI limit.

    This function requires P_upstream > P_sat(T_upstream) -- i.e. the
    fluid must still be liquid (saturated or subcooled) AT the orifice
    inlet, per Section 1.4. If P_upstream <= P_sat(T_upstream), the fluid
    has already crossed the saturation curve before reaching the orifice
    at all: this is a modeling error (the two-phase feed line problem,
    not the two-phase orifice problem this module addresses), so it is
    flagged loudly rather than producing a meaningless or infinite kappa.

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

        m_dot_Dyer = m_dot_SPI / (1 + kappa) + kappa / (1 + kappa) * m_dot_HEM

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

    m_dot_Dyer = m_dot_SPI / (1.0 + kappa) + (kappa / (1.0 + kappa)) * m_dot_HEM

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
