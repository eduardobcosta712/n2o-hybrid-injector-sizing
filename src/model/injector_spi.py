"""
injector_spi.py

SPI (Single Phase Incompressible) injector model — Bernoulli's equation
applied to the injector orifice, empirically corrected by a discharge
coefficient. See docs/02_spi_model.md for the full derivation.

    m_dot_SPI = Cd * A * sqrt(2 * rho * delta_P)

This module also provides the criterion used to decide whether SPI alone
is sufficient at a given injector operating point, or whether two-phase
effects (injector_two_phase.py) must be considered -- per the decision
logic laid out in the pre-implementation plan (Section "Module 2" of the
justification) and Section 3.1 of the theory docs.

Units: SI throughout (Pa, kg/m^3, m^2, kg/s).
"""

import math
from n2o_properties import P_sat


def spi_mass_flow(Cd, A, rho, delta_P):
    """
    SPI mass flow rate through an orifice.

    m_dot = Cd * A * sqrt(2 * rho * delta_P)

    Parameters
    ----------
    Cd : float
        Discharge coefficient (dimensionless, typically 0.6-0.9).
    A : float
        Total orifice area, m^2 (if multiple holes, pass N_holes * A_hole).
    rho : float
        Upstream liquid density, kg/m^3.
    delta_P : float
        Pressure drop across the orifice, Pa (P_upstream - P_downstream).
        Must be >= 0; a negative value means the flow direction assumed
        is wrong (fluid would flow backward), which is a modeling error,
        not a physical result -- so this is checked explicitly below.

    Returns
    -------
    float
        Mass flow rate, kg/s.
    """
    if delta_P < 0:
        raise ValueError(
            f"delta_P = {delta_P:.1f} Pa is negative -- this would imply "
            "reverse flow through the orifice. Check that P_upstream and "
            "P_downstream were passed in the correct order."
        )
    return Cd * A * math.sqrt(2.0 * rho * delta_P)


def orifice_area_from_target_flow(m_dot_target, Cd, rho, delta_P):
    """
    Inverse of spi_mass_flow: required orifice area to deliver a target
    mass flow rate, given Cd, upstream density, and design delta_P.

    Obtained by algebraically inverting m_dot = Cd * A * sqrt(2 rho dP):
        A = m_dot / (Cd * sqrt(2 * rho * delta_P))

    This is the sizing direction most commonly used in practice: the
    target m_dot comes from the motor's thermochemical sizing (O/F,
    chamber pressure), and this function returns the orifice area needed
    to deliver it.

    Returns
    -------
    float
        Required total orifice area, m^2.
    """
    return m_dot_target / (Cd * math.sqrt(2.0 * rho * delta_P))


def spi_sufficient(P_upstream, T_upstream, P_downstream):
    """
    Decision criterion (Section 3.1 / implementation plan "Module 2"):
    is the pure SPI model sufficient on its own at this operating point,
    or does the pressure drop across the orifice cross the saturation
    curve, requiring a two-phase model (HEM/Dyer, injector_two_phase.py)?

    SPI alone is sufficient when the fluid stays liquid the entire way
    through the orifice -- i.e. even the downstream pressure (where P is
    lowest, at the vena contracta, per Section 2.1) stays at or above the
    saturation pressure at the upstream temperature.

    Parameters
    ----------
    P_upstream : float
        Pressure just before the orifice, Pa.
    T_upstream : float
        Temperature just before the orifice, K (assumed to not change
        meaningfully across the short orifice itself, per Section 1.6).
    P_downstream : float
        Pressure just after the orifice (e.g. chamber pressure), Pa.

    Returns
    -------
    bool
        True if SPI alone is a valid model here; False if two-phase
        effects must be considered.
    """
    return P_downstream >= P_sat(T_upstream)


if __name__ == "__main__":
    # --- Validation case, referenced in docs/references.md: a published
    #     hybrid motor injector operating point, comparing the SPI
    #     prediction against the documented real (flashing-affected) flow.
    #
    # Reference case parameters (representative of the documented
    # discrepancy discussed in the theory docs -- see references.md for
    # the exact source and figures):
    #   - N2O at ~20 degC upstream of the injector
    #   - Upstream (tank/run-line) pressure and downstream (chamber)
    #     pressure such that the pressure drop crosses saturation
    #   - SPI over-predicts flow because it assumes single-phase liquid
    #     throughout

    from n2o_properties import rho_liquid_sat

    T_upstream = 293.15       # K, 20 degC
    P_upstream = 50e5         # Pa, 50 bar upstream
    P_downstream = 20e5       # Pa, 20 bar chamber pressure (example)
    Cd = 0.65                 # typical injector discharge coefficient
    A = 3.79e-6               # m^2, example total orifice area (from a
                               # published design point, see references.md)

    rho = rho_liquid_sat(T_upstream)
    delta_P = P_upstream - P_downstream

    m_dot = spi_mass_flow(Cd, A, rho, delta_P)
    sufficient = spi_sufficient(P_upstream, T_upstream, P_downstream)

    print("SPI injector model -- example evaluation")
    print("-" * 60)
    print(f"Upstream temperature:   {T_upstream:.2f} K ({T_upstream-273.15:.1f} degC)")
    print(f"Upstream density:       {rho:.1f} kg/m^3")
    print(f"Upstream pressure:      {P_upstream/1e5:.1f} bar")
    print(f"Downstream pressure:    {P_downstream/1e5:.1f} bar")
    print(f"Pressure drop:          {delta_P/1e5:.1f} bar")
    print(f"P_sat at T_upstream:    {P_sat(T_upstream)/1e5:.1f} bar")
    print("-" * 60)
    print(f"SPI-predicted mass flow: {m_dot*1000:.1f} g/s")
    print(f"SPI alone sufficient?    {sufficient}")
    if not sufficient:
        print("  -> Downstream pressure is below P_sat(T_upstream): the flow")
        print("     crosses the saturation curve inside the orifice.")
        print("     SPI will OVER-predict the real mass flow rate here --")
        print("     see injector_two_phase.py for the Dyer model correction.")
