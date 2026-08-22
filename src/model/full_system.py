"""
full_system.py

Orchestrates the full tank -> feed line -> injector path in a single call,
per the Section 3.5 synthesis: run the feed line to get conditions at the
injector inlet, then automatically decide (via spi_sufficient) whether SPI
alone is valid or whether the two-phase Dyer model must be used, and
return the resulting mass flow along with the full trace for inspection.

This module introduces no new physics -- it only sequences the existing
modules (feed_line.py, injector_spi.py, injector_two_phase.py) in the
order the sizing problem requires.

Units: SI throughout (Pa, K, kg/m^3, m^2, kg/s).
"""

from n2o_properties import rho_liquid_sat, nu_vapor_sat, T_sat
from feed_line import evaluate_feed_line
from injector_spi import spi_mass_flow, spi_sufficient
from injector_two_phase import dyer_mass_flow

M_N2O = 44.013  # kg/kmol, molar mass of N2O (n2o_properties.py convention)


def evaluate_full_system(m_dot_design, T_tank, P_tank, segments,
                          Cd, A_injector, P_chamber,
                          roughness=1.5e-6):
    """
    Run the full tank -> feed line -> injector path for a given design
    mass flow and injector area, and report the real (model-corrected)
    mass flow the system would actually deliver.

    Parameters
    ----------
    m_dot_design : float
        Mass flow rate used to evaluate the feed line's pressure drop
        (Section 4.2), kg/s. Typically the target design mass flow.
    T_tank : float
        N2O tank temperature, K (assumed constant along the line).
    P_tank : float
        N2O tank pressure, Pa.
    segments : list of dict
        Feed line geometry, in the format evaluate_feed_line expects
        (see feed_line.py).
    Cd : float
        Injector discharge coefficient (dimensionless).
    A_injector : float
        Total injector orifice area, m^2.
    P_chamber : float
        Chamber pressure (injector downstream pressure), Pa.
    roughness : float
        Feed line absolute pipe wall roughness, m.

    Returns
    -------
    dict
        "feed_line_result": full return value of evaluate_feed_line.
        "P_injector_inlet": pressure at the injector inlet, Pa (== feed
            line's P_final).
        "spi_sufficient": bool, whether SPI alone was valid at the
            injector inlet conditions.
        "m_dot_real": the model-selected real mass flow, kg/s -- equal to
            the SPI prediction if spi_sufficient is True, otherwise the
            Dyer prediction.
        "injector_result": None if spi_sufficient is True; otherwise the
            full dict returned by dyer_mass_flow (m_dot_Dyer, m_dot_SPI,
            m_dot_HEM, kappa, x_exit) for full inspection.
    """
    feed_line_result = evaluate_feed_line(m_dot_design, T_tank, P_tank,
                                           segments, roughness)
    P_injector_inlet = feed_line_result["P_final"]

    # If the feed line itself produced flashing, the fluid arriving at
    # the injector inlet is already two-phase -- outside the domain of
    # both SPI and Dyer (both assume liquid at the orifice inlet, per
    # Section 3.1/injector_two_phase.py's domain restriction). Return
    # early with a clear flag; the injector model is not run.
    if feed_line_result["flashing_detected"]:
        return {
            "feed_line_result": feed_line_result,
            "P_injector_inlet": P_injector_inlet,
            "spi_sufficient": None,
            "m_dot_real": None,
            "injector_result": None,
        }

    sufficient = spi_sufficient(P_injector_inlet, T_tank, P_chamber)

    if sufficient:
        rho_l = rho_liquid_sat(T_tank)
        delta_P = P_injector_inlet - P_chamber
        m_dot_real = spi_mass_flow(Cd, A_injector, rho_l, delta_P)
        injector_result = None
    else:
        rho_l_upstream = rho_liquid_sat(T_tank)
        T_downstream = T_sat(P_chamber)
        rho_l_downstream = rho_liquid_sat(T_downstream)
        rho_v_downstream = M_N2O / nu_vapor_sat(T_downstream)

        injector_result = dyer_mass_flow(
            Cd, A_injector, T_tank, P_injector_inlet, P_chamber,
            rho_l_upstream, rho_l_downstream, rho_v_downstream)
        m_dot_real = injector_result["m_dot_Dyer"]

    return {
        "feed_line_result": feed_line_result,
        "P_injector_inlet": P_injector_inlet,
        "spi_sufficient": sufficient,
        "m_dot_real": m_dot_real,
        "injector_result": injector_result,
    }


if __name__ == "__main__":
    # --- Validation: reuse feed_line.py's own Case B geometry (2 m of
    #     8 mm ID tubing, one ball valve, one 90-degree elbow), now
    #     chained automatically into the injector decision, instead of
    #     evaluating the feed line and the injector as separate,
    #     manually-connected steps.
    T_tank = 293.15   # K, 20 degC
    P_tank = 55e5      # Pa, 55 bar (same subcooling margin used to
                       # validate injector_two_phase.py)
    m_dot_design = 0.5  # kg/s

    segments = [
        {"type": "pipe", "L": 1.0, "D": 0.008},
        {"type": "fitting", "D": 0.008, "K": 0.05},
        {"type": "pipe", "L": 1.0, "D": 0.008},
        {"type": "fitting", "D": 0.008, "K": 0.3},
    ]

    Cd = 0.65
    A_injector = 3.79e-6  # m^2
    P_chamber = 20e5       # Pa, 20 bar

    result = evaluate_full_system(m_dot_design, T_tank, P_tank, segments,
                                   Cd, A_injector, P_chamber)

    print("Full system evaluation (tank -> feed line -> injector)")
    print("-" * 60)
    print(f"Tank:                     {P_tank/1e5:.1f} bar, {T_tank-273.15:.1f} degC")
    print(f"Injector inlet pressure:  {result['P_injector_inlet']/1e5:.2f} bar")
    print(f"Feed line flashing?       {result['feed_line_result']['flashing_detected']}")
    print(f"SPI alone sufficient?     {result['spi_sufficient']}")
    print("-" * 60)
    print(f"Real mass flow (model-selected): {result['m_dot_real']*1000:.1f} g/s")
    if result["injector_result"] is not None:
        ir = result["injector_result"]
        print(f"  (Dyer used -- SPI would have predicted {ir['m_dot_SPI']*1000:.1f} g/s, "
              f"HEM {ir['m_dot_HEM']*1000:.1f} g/s, kappa={ir['kappa']:.2f})")
