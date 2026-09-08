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

Two vapour-quality paths are used in this module, for two different
physical questions (see docs/future_work.md, Priority 1):

    - ISENTHALPIC (h = const): describes the REAL thermodynamic state of
      the fluid at a given point (e.g. the orifice exit). An orifice is
      adiabatic but not reversible, so the 1st law gives h_up = h_down
      regardless of internal losses -- this is the correct path for
      vapor_quality_isenthalpic, hem_mass_flow, dyer_mass_flow, and
      hem_critical_flow (the mass-flow curve itself is still evaluated
      along real, physically-realized states).
    - ISENTROPIC (s = const): describes the speed at which a pressure
      disturbance propagates, c^2 = (dP/drho)_s -- a thermodynamic
      derivative taken at constant entropy, because an acoustic wave is a
      small, fast, essentially reversible perturbation on top of the
      (possibly irreversible) mean flow. This is the physically correct
      path for locating the two-phase CHOKING condition, implemented in
      vapor_quality_isentropic and hem_critical_flow_isentropic.

hem_critical_flow() (isenthalpic) is kept unchanged and side-by-side with
hem_critical_flow_isentropic(): it is already validated and cited in
validation/waxman_2013_results.md, and remains useful as a direct
comparison against the more rigorous isentropic scan.

Units: SI throughout (Pa, K, kg/m^3, m^2, kg/s), except vapor quality x
and the Dyer weighting parameter kappa, which are dimensionless.
"""

import math
from n2o_properties import (P_sat, T_sat, dP_sat_dT, h_liquid_sat, h_fg,
                             rho_liquid_sat, nu_vapor_sat,
                             s_liquid_sat, s_vapor_sat, s_fg,
                             T_MIN, T_MIN_A4, T_MAX_A4)
from injector_spi import spi_mass_flow

M_N2O = 44.013  # kg/kmol, molar mass of N2O


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


def vapor_quality_isentropic(s_upstream, T_downstream):
    """
    Vapor quality x at a given downstream state, assuming an ISENTROPIC
    (constant-entropy) process -- the path that correctly locates the
    two-phase CHOKING condition (see module docstring), as opposed to the
    isenthalpic path (vapor_quality_isenthalpic above) that correctly
    describes the real thermodynamic state at the orifice exit.

    Derived exactly analogously to vapor_quality_isenthalpic, replacing
    enthalpy with entropy: s_upstream = s_l(T_downstream) + x * s_fg(T_downstream)

        x = (s_upstream - s_l(T_downstream)) / s_fg(T_downstream)

    Added September 2026 for the isentropic choking limit (Priority 1,
    docs/future_work.md).

    Parameters
    ----------
    s_upstream : float
        Upstream fluid molar entropy, kJ/(kmol.K). For liquid entering
        subcooled or saturated, this is s_liquid_sat(T_upstream); for a
        two-phase inlet (feed-line flashing), this is the quality-weighted
        mixture entropy -- see hem_critical_flow_isentropic below.
    T_downstream : float
        Saturation temperature corresponding to the downstream pressure
        being scanned, T_sat(P_downstream), K. Must lie within Table
        A.4's range (see n2o_properties.T_MIN_A4 / T_MAX_A4) -- narrower
        than the enthalpy table's range, since entropy data comes from a
        separate, shorter NIST table.

    Returns
    -------
    float
        Vapor quality x (dimensionless), clamped to [0, 1] for the same
        reason as vapor_quality_isenthalpic.
    """
    from n2o_properties import s_liquid_sat, s_fg
    x = (s_upstream - s_liquid_sat(T_downstream)) / s_fg(T_downstream)
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
        "x_exit":           x_exit,
        "rho_HEM":          rho_HEM,
        "x_inlet":          x_inlet,
    }




def hem_critical_flow(Cd, A, T_upstream, P_upstream, x_inlet=0.0,
                       n_steps=200):
    """
    HEM isenthalpic maximum mass flow rate (two-phase critical flow), kg/s.

    Finds the choking limit by locating the maximum of the HEM mass-flow
    rate as a function of downstream pressure P2, scanning P2 from P_sat
    (onset of two-phase flow) down to a minimum pressure. The maximum is
    the physical critical (choked) mass flow -- the ceiling that the
    Bernoulli-based Dyer formula cannot exceed.

    Reference:
        Waxman (2013), Eq. (5):
            m_dot_HEM_crit = Cd * A * max_{P2 < P_sat} sqrt(2 * rho_HEM * dP)
        where rho_HEM and dP are evaluated along the isenthalpic path.

    This approach does not require entropy data or the Henry-Fauske
    approximations. It uses only the enthalpy-based properties already
    available in n2o_properties.py (h_l, h_fg, rho_l, rho_v), making it
    fully consistent with the rest of the model.

    Physical interpretation: the HEM mass flow first increases with
    delta_P (more driving pressure) but then decreases because rho_HEM
    falls rapidly as more vapour forms. The maximum occurs at the pressure
    where these two effects balance -- this is the two-phase speed-of-sound
    condition (the choked state), expressed through the isenthalpic path.

    Approximation note (see hem_critical_flow_isentropic for the more
    rigorous version, added September 2026, Priority 1): the isenthalpic
    path correctly describes the real thermodynamic STATE of the fluid at
    each P2, but the true choking condition -- the propagation speed of a
    pressure disturbance -- is an isentropic derivative, c^2 = (dP/drho)_s.
    This function is kept as-is, unchanged, because it is already
    validated (validation/waxman_2013_results.md) and remains a useful,
    entropy-data-free reference point.

    Parameters
    ----------
    Cd : float
        Discharge coefficient.
    A : float
        Total orifice area, m^2.
    T_upstream : float
        Upstream temperature (= tank temperature), K.
    P_upstream : float
        Upstream pressure, Pa.
    x_inlet : float, optional
        Vapour quality at the orifice inlet (from feed-line flash). Default 0.
    n_steps : int, optional
        Number of P2 values to scan between P_sat and P_min. Default 200.
        Higher values give more accurate peak location.

    Returns
    -------
    dict
        "m_dot_crit"  : critical (maximum) HEM mass flow rate, kg/s
        "P2_crit"     : downstream pressure at the critical condition, Pa
        "x_crit"      : vapour quality at the critical condition
        "rho_crit"    : HEM mixture density at the critical condition, kg/m^3
    """
    import math
    from n2o_properties import (P_sat, T_sat as T_sat_f, rho_liquid_sat,
                                 nu_vapor_sat, h_liquid_sat, h_fg)

    M_N2O = 44.013

    # Upstream enthalpy (conserved along isenthalpic path)
    h_up = h_liquid_sat(T_upstream) + x_inlet * h_fg(T_upstream)

    # Saturation pressure at upstream temperature -- onset of two-phase flow
    P_sat_up = P_sat(T_upstream)

    # Scan range: from just below P_sat down to 5% of P_sat
    P_min   = max(0.05 * P_sat_up, 1e5)   # never below 1 bar
    P_start = min(P_sat_up * 0.999, P_upstream - 1e3)

    best_m = 0.0
    best_P2 = P_start
    best_x  = 0.0
    best_rho = rho_liquid_sat(T_upstream)

    dP_step = (P_start - P_min) / max(n_steps, 1)
    P2 = P_start

    while P2 >= P_min:
        T2      = T_sat_f(P2)
        rho_l2  = rho_liquid_sat(T2)
        rho_v2  = M_N2O / nu_vapor_sat(T2)
        hl2     = h_liquid_sat(T2)
        hfg2    = h_fg(T2)

        if hfg2 <= 0:
            P2 -= dP_step; continue

        x = (h_up - hl2) / hfg2
        x = max(0.0, min(1.0, x))

        nu_mix  = (1.0 - x) / rho_l2 + x / rho_v2
        rho_mix = 1.0 / nu_mix
        dP      = P_upstream - P2

        m = Cd * A * math.sqrt(2.0 * rho_mix * dP) if dP > 0 else 0.0

        if m > best_m:
            best_m   = m
            best_P2  = P2
            best_x   = x
            best_rho = rho_mix

        P2 -= dP_step

    return {
        "m_dot_crit": best_m,
        "P2_crit":    best_P2,
        "x_crit":     best_x,
        "rho_crit":   best_rho,
    }


def hem_critical_flow_isentropic(Cd, A, T_upstream, P_upstream, x_inlet=0.0,
                                  n_steps=200):
    """
    HEM isentropic maximum mass flow rate (two-phase critical flow), kg/s.

    Isentropic-path analogue of hem_critical_flow(): finds the maximum of
    the HEM mass-flow curve along a path of constant molar entropy rather
    than constant molar enthalpy.

    Physical motivation (Section 3.3 / docs/future_work.md, Priority 1):
    the two-phase speed of sound -- and therefore the true choking
    condition -- is a thermodynamic derivative taken at constant entropy,
        c^2 = (dP/drho)_s,
    because an acoustic disturbance is, by definition, a small, fast,
    essentially reversible perturbation superimposed on the (possibly
    irreversible) mean flow. hem_critical_flow() uses the isenthalpic
    path instead, which correctly describes the real thermodynamic STATE
    of the fluid at the orifice exit (energy balance across an adiabatic,
    lossy orifice) but is only an approximation to the propagation-speed
    condition that actually defines choking. This function implements the
    more rigorous isentropic scan; hem_critical_flow() is kept unchanged,
    side-by-side, for direct comparison -- it is already validated and
    cited in validation/waxman_2013_results.md.

    Domain restriction: entropy data (Table A.4, NIST WebBook) only
    covers T in [n2o_properties.T_MIN_A4, n2o_properties.T_MAX_A4], i.e.
    up to 307.33 K -- short of the 309.52 K critical point used
    elsewhere in this project. Since the scan starts at T_sat(P_upstream)
    ~= T_upstream and moves to lower T (lower P2), T_upstream itself is
    the binding constraint: this function raises ValueError immediately
    if T_upstream exceeds Table A.4's range, rather than let the scan
    fail partway through with a less legible error. hem_critical_flow()
    (isenthalpic) remains available as a fallback this close to the
    critical point; CoolProp integration (Priority 4) would remove this
    limitation entirely.

    Parameters
    ----------
    Cd : float
        Discharge coefficient.
    A : float
        Total orifice area, m^2.
    T_upstream : float
        Upstream temperature (= tank temperature), K. Must be
        <= n2o_properties.T_MAX_A4 (307.33 K) -- see domain restriction.
    P_upstream : float
        Upstream pressure, Pa.
    x_inlet : float, optional
        Vapour quality at the orifice inlet (from feed-line flash).
        Default 0.
    n_steps : int, optional
        Number of P2 values to scan between P_sat and P_min. Default 200.

    Returns
    -------
    dict
        "m_dot_crit"  : critical (maximum) HEM mass flow rate, kg/s
        "P2_crit"     : downstream pressure at the critical condition, Pa
        "x_crit"      : vapour quality at the critical condition
        "rho_crit"    : HEM mixture density at the critical condition, kg/m^3

    Raises
    ------
    ValueError
        If T_upstream is outside Table A.4's valid range.
    """
    import math
    from n2o_properties import (P_sat, T_sat as T_sat_f, rho_liquid_sat,
                                 nu_vapor_sat, s_liquid_sat, s_vapor_sat,
                                 s_fg, T_MIN_A4, T_MAX_A4)

    M_N2O = 44.013

    if not (T_MIN_A4 <= T_upstream <= T_MAX_A4):
        raise ValueError(
            f"T_upstream = {T_upstream:.2f} K is outside Table A.4's valid "
            f"range [{T_MIN_A4}, {T_MAX_A4}] K (NIST WebBook, Lemmon & Span "
            "2006), the only source of entropy data in this project. The "
            "isentropic choking scan cannot be evaluated this close to the "
            "critical point with the current data. hem_critical_flow() "
            "(isenthalpic) remains available as a fallback in this regime "
            "-- see also docs/future_work.md, Priority 4 (CoolProp/REFPROP "
            "integration), which would remove this limitation."
        )

    # Upstream entropy (conserved along isentropic path)
    s_up = s_liquid_sat(T_upstream) + x_inlet * s_fg(T_upstream)

    # Saturation pressure at upstream temperature -- onset of two-phase flow
    P_sat_up = P_sat(T_upstream)

    # Scan range: from just below P_sat down to 5% of P_sat (mirrors
    # hem_critical_flow exactly, for a like-for-like comparison)
    P_min   = max(0.05 * P_sat_up, 1e5)   # never below 1 bar
    P_start = min(P_sat_up * 0.999, P_upstream - 1e3)

    best_m = 0.0
    best_P2 = P_start
    best_x  = 0.0
    best_rho = rho_liquid_sat(T_upstream)

    dP_step = (P_start - P_min) / max(n_steps, 1)
    P2 = P_start

    while P2 >= P_min:
        T2 = T_sat_f(P2)
        # T2 <= T_upstream <= T_MAX_A4 always holds since P2 <= P_start <
        # P_sat(T_upstream), so this guard should never trigger -- kept
        # anyway per the project's "fail loudly, never silently" convention,
        # mirroring the hfg2 <= 0 guard in hem_critical_flow.
        if T2 > T_MAX_A4:
            P2 -= dP_step; continue

        rho_l2 = rho_liquid_sat(T2)
        rho_v2 = M_N2O / nu_vapor_sat(T2)
        sfg2   = s_fg(T2)

        if sfg2 <= 0:
            P2 -= dP_step; continue

        x = (s_up - s_liquid_sat(T2)) / sfg2
        x = max(0.0, min(1.0, x))

        nu_mix  = (1.0 - x) / rho_l2 + x / rho_v2
        rho_mix = 1.0 / nu_mix
        dP      = P_upstream - P2

        m = Cd * A * math.sqrt(2.0 * rho_mix * dP) if dP > 0 else 0.0

        if m > best_m:
            best_m   = m
            best_P2  = P2
            best_x   = x
            best_rho = rho_mix

        P2 -= dP_step

    return {
        "m_dot_crit": best_m,
        "P2_crit":    best_P2,
        "x_crit":     best_x,
        "rho_crit":   best_rho,
    }


def henry_fauske_critical_flow(Cd, A, T_upstream, P_upstream, x_inlet=0.0,
                                tol_bar=1e-4, max_iter=100):
    """
    Henry-Fauske (1971) non-equilibrium critical mass flux for saturated
    or subcooled liquid discharging through a converging nozzle/orifice.

    Added September 2026 (docs/future_work.md, Priority 1) after
    `hem_critical_flow()`/`hem_critical_flow_isentropic()` were shown
    NOT to be the right ceiling for a non-equilibrium (Dyer) prediction:
    all 4 validated Waxman operating points already sit 1.03x-1.21x
    ABOVE the equilibrium HEM ceiling, which is correct, validated
    behaviour (real non-equilibrium two-phase flow chokes at a HIGHER
    mass flux than the full-equilibrium limit). Henry-Fauske models that
    non-equilibrium choking directly, instead of assuming full
    equilibrium.

    Source
    ------
    Henry, R.E. & Fauske, H.K. (1971), "The Two-Phase Critical Flow of
    One-Component Mixtures in Nozzles, Orifices, and Short Tubes", ASME
    J. Heat Transfer, 93(2), 179-187.

    Equations transcribed here from the simplified form given in:
    Simoneau, R.J., Henry, R.E., Hendricks, R.C. & Watterson, R. (1971),
    "Two-Phase Critical Discharge of High Pressure Liquid Nitrogen",
    NASA TM X-67863, Eqs. (2)-(5).

    Theory
    ------
    For saturated/subcooled liquid at the nozzle inlet (stagnation
    conditions, subscript 0), the simplified Henry-Fauske model makes
    five assumptions valid for P0/Pc > 0.05 (comfortably true for this
    project's operating range):
        1. Negligible vapour forms before the throat (x_t ~ 0), so the
           inlet-to-throat momentum balance is single-phase Bernoulli:
               eta = P_t/P0 = 1 - v_l0 * G_c^2 / (2 * P0)      (Eq. 2)
        2. The liquid phase is incompressible.
        3. Vapour forming AT the throat is in equilibrium at the local
           (throat) pressure.
        4. Liquid and vapour velocities are equal at the throat (k~1).
        5. The rate of mass transfer (vaporisation) at the throat is a
           fraction N of the full-equilibrium rate:
               N = min(1, x_E / 0.14)
           where x_E is the EQUILIBRIUM quality at the throat -- the
           same isentropic quality already computed by
           vapor_quality_isentropic() for the HEM-equilibrium ceiling.
           N = 0.14 is Henry's (1970) empirical fit to the Starkman et
           al. steam-water low-quality expansion data; capped at 1 once
           the flow is close enough to equilibrium that no further
           non-equilibrium enhancement applies.

    Combining these gives the critical mass flux directly:

        G_c^2 = [ N*(v_gE - v_l0) / (s_gE - s_lE) * ds_lE/dP ]^-1      (Eq. 5)

    evaluated AT the (unknown) throat state -- v_gE, s_gE, s_lE are
    saturated vapour/liquid specific volume and entropy at T_sat(P_t),
    and ds_lE/dP is the derivative of saturated liquid specific entropy
    along the saturation curve, obtained via the chain rule
    ds_lE/dP = (ds_l/dT) / (dP_sat/dT), with ds_l/dT from a small
    central finite difference on s_liquid_sat(T) (no closed-form
    derivative is available, since s_l comes from table interpolation,
    not a fitted correlation like P_sat) and dP_sat/dT from the existing
    analytical dP_sat_dT(T).

    Eqs. (2) and (5) are coupled -- G_c depends on properties evaluated
    AT the throat pressure P_t, and P_t (via Eq. 2) depends on G_c -- so
    they are solved simultaneously here via bisection on P_t (found to
    be far more robust than naive fixed-point iteration, which diverges
    for this system: G_c(P_t) from Eq. 5 blows up as P_t -> P_sat(T0)
    from below, since N -> 0 there).

    Units note: Eq. (5) requires SPECIFIC (per unit mass) volume and
    entropy throughout for the units to work out to a mass flux. This
    project's n2o_properties.py functions are molar (v via
    nu_vapor_sat/M_N2O and rho_liquid_sat already give specific volume
    directly; s_liquid_sat/s_vapor_sat are converted from
    kJ/(kmol.K) to J/(kg.K) internally here via *1000/M_N2O).

    Parameters
    ----------
    Cd : float
        Discharge coefficient (dimensionless). Applied multiplicatively
        to the ideal nozzle mass flux, consistent with how Cd is used
        throughout this project (SPI, HEM, Dyer) -- the original paper's
        G_c itself has no Cd (derived for a characterised venturi nozzle).
    A : float
        Total orifice area, m^2.
    T_upstream : float
        Stagnation (upstream/tank) temperature, K. Must lie within
        Table A.4's range (n2o_properties.T_MIN_A4 to T_MAX_A4) --
        entropy is required.
    P_upstream : float
        Stagnation (upstream/tank) pressure, Pa. Must exceed
        P_sat(T_upstream) -- the model assumes liquid (saturated or
        subcooled) at the nozzle inlet, mirroring
        dyer_non_equilibrium_parameter's domain restriction.
    x_inlet : float, optional
        Vapour quality at the nozzle inlet, for a two-phase (already
        flashing) inlet. Default 0 (pure liquid stagnation).
    tol_bar : float, optional
        Bisection convergence tolerance on the throat pressure, bar.
        Default 1e-4 bar -- tight, converges in ~40-60 iterations
        regardless.
    max_iter : int, optional
        Maximum bisection iterations before returning the best estimate
        found so far without erroring (bisection with a valid bracket
        always converges; this is a safety cap, not expected to bind).

    Returns
    -------
    dict
        "m_dot_crit" : critical (choked) mass flow rate, kg/s.
        "G_c"        : critical mass flux, kg/(m^2.s).
        "P2_crit"    : throat pressure at the critical condition, Pa.
        "x_crit"     : equilibrium quality at the critical throat state
                       (x_E, NOT the actual non-equilibrium quality N*x_E).
        "N"          : Henry's non-equilibrium factor at the critical
                       condition, in [0, 1].

    Raises
    ------
    ValueError
        If T_upstream is outside Table A.4's range, or P_upstream does
        not exceed P_sat(T_upstream).
    RuntimeError
        If a valid bracket for the bisection cannot be found (should not
        happen for physically sensible inputs within this project's
        design range; would indicate a genuinely pathological operating
        point).

    Validation
    -----------
    At Waxman conditions (T=280 K, P=4.36 MPa, D=1.5 mm, Cd=0.65):
    m_dot_crit = 50.67 g/s -- above hem_critical_flow_isentropic's
    41.64 g/s (correct: non-equilibrium exceeds equilibrium), and above
    all 4 validated Dyer predictions (42.25-49.55 g/s), so it does not
    cut into any already-validated result. See
    validation/waxman_2013_results.md for the full comparison.

    IMPORTANT CAVEAT (see docs/future_work.md, Priority 1). At OTHER
    operating points further from the Waxman geometry -- e.g. the
    conditions in examples/example_01_sizing.md (20 degC, 58->22 bar,
    6x1.5mm holes) -- this ceiling DOES bind, sitting roughly 13-17%
    below the uncapped Dyer blend. There is currently NO experimental
    data point in this project's validation set where the Henry-Fauske
    ceiling actually changes the answer (all 4 Waxman points sit below
    it) -- so while the model is theoretically sound and correctly
    implemented from a primary source, it is NOT empirically confirmed
    in the regime where it matters. For this reason it is surfaced as a
    side-by-side diagnostic value with a `choked` flag in
    dyer_mass_flow(), NOT applied as an automatic override of
    `m_dot_Dyer` -- see that function's docstring.
    """
    if not (T_MIN_A4 <= T_upstream <= T_MAX_A4):
        raise ValueError(
            f"T_upstream = {T_upstream:.2f} K is outside Table A.4's valid "
            f"range [{T_MIN_A4}, {T_MAX_A4}] K -- entropy data is required "
            "for the Henry-Fauske non-equilibrium closure. See "
            "docs/future_work.md, Priority 4 (CoolProp/REFPROP)."
        )

    P_sat_up = P_sat(T_upstream)
    if P_upstream <= P_sat_up:
        raise ValueError(
            f"P_upstream = {P_upstream/1e5:.2f} bar is at or below "
            f"P_sat(T_upstream) = {P_sat_up/1e5:.2f} bar -- the fluid must "
            "be liquid (saturated or subcooled) at the nozzle inlet for "
            "this model, mirroring dyer_non_equilibrium_parameter's "
            "domain restriction."
        )

    v_l0 = 1.0 / rho_liquid_sat(T_upstream)
    s_up_molar = s_liquid_sat(T_upstream) + x_inlet * s_fg(T_upstream)  # kJ/(kmol.K)

    h_T = 0.05  # K, finite-difference step for ds_l/dT

    def _G_c_from_eq5(P_t):
        """Eq. (5): the mass-transfer-rate closure evaluated at throat
        pressure P_t. Returns (G_c, x_E, N) or (None, x_E, N) if the
        bracket term is non-positive (no real G_c at this P_t)."""
        T_t = T_sat(P_t)
        v_gE = nu_vapor_sat(T_t) / M_N2O                      # m^3/kg
        s_lE = s_liquid_sat(T_t) * 1000.0 / M_N2O              # J/(kg.K)
        s_gE = s_vapor_sat(T_t) * 1000.0 / M_N2O               # J/(kg.K)

        x_E = vapor_quality_isentropic(s_up_molar, T_t)
        # N = x_E/0.14 formally divides by zero as x_E -> 0 (no vapour
        # yet, no non-equilibrium enhancement needed); use a small floor
        # instead of exactly zero so the bracket term stays finite and
        # bisection can still evaluate the residual there.
        N = min(1.0, x_E / 0.14) if x_E > 1e-9 else 1e-9

        T_plus = min(T_t + h_T, T_MAX_A4)
        T_minus = max(T_t - h_T, T_MIN_A4)
        s_l_plus = s_liquid_sat(T_plus) * 1000.0 / M_N2O
        s_l_minus = s_liquid_sat(T_minus) * 1000.0 / M_N2O
        ds_l_dT = (s_l_plus - s_l_minus) / (T_plus - T_minus)
        ds_lE_dP = ds_l_dT / dP_sat_dT(T_t)

        denom = N * (v_gE - v_l0) / (s_gE - s_lE) * ds_lE_dP
        if denom <= 0:
            return None, x_E, N
        return math.sqrt(1.0 / denom), x_E, N

    def _residual(P_t):
        """Residual between the throat pressure implied by momentum
        (Eq. 2, using G_c from Eq. 5) and the P_t plugged in. Root of
        this function is the self-consistent (P_t, G_c) solution."""
        G_c, x_E, N = _G_c_from_eq5(P_t)
        if G_c is None:
            return None, None, None, None
        P_t_from_momentum = P_upstream - v_l0 * G_c ** 2 / 2.0
        return P_t_from_momentum - P_t, G_c, x_E, N

    # Bracket the root. Just below P_sat(T_upstream), N -> 0 so G_c(Eq.5)
    # blows up, making the residual strongly negative (momentum can't
    # support that flux without P_t going deeply negative). At a low
    # pressure floor, N has saturated at 1 and G_c(Eq.5) has dropped
    # enough that the residual turns positive -- found empirically to
    # bracket reliably across this project's operating range.
    P_hi = P_sat_up * 0.999
    P_lo = max(P_sat(T_MIN + 1.0) * 1.01, 1e5)

    r_hi, *_ = _residual(P_hi)
    r_lo, *_ = _residual(P_lo)
    if r_hi is None or r_lo is None or r_hi * r_lo > 0:
        raise RuntimeError(
            f"henry_fauske_critical_flow: could not bracket a root for "
            f"T_upstream={T_upstream:.2f} K, P_upstream={P_upstream/1e5:.2f} bar "
            f"(r_lo={r_lo} at {P_lo/1e5:.2f} bar, r_hi={r_hi} at {P_hi/1e5:.2f} bar). "
            "This should not happen for physically sensible inputs within "
            "this project's design range -- please report this operating "
            "point."
        )

    P_t = 0.5 * (P_hi + P_lo)
    G_c = x_E = N = None
    for _ in range(max_iter):
        P_t = 0.5 * (P_hi + P_lo)
        r_mid, G_c, x_E, N = _residual(P_t)
        if r_mid is None:
            P_lo = P_t
            continue
        if r_mid > 0:
            P_lo = P_t
        else:
            P_hi = P_t
        if (P_hi - P_lo) < tol_bar * 1e5:
            break

    return {
        "m_dot_crit": Cd * A * G_c,
        "G_c":        G_c,
        "P2_crit":    P_t,
        "x_crit":     x_E,
        "N":          N,
    }


def apply_choking_limit(m_dot_model, Cd, A, T_upstream, P_upstream,
                         x_inlet=0.0):
    """
    Apply the HEM isenthalpic choking limit to a model-predicted mass flow.

    Returns the physically realizable mass flow:
        m_dot_real = min(m_dot_model, m_dot_crit)

    The choking limit is computed by hem_critical_flow() -- the maximum
    of the HEM isenthalpic mass flow curve (Waxman 2013, Eq. 5). This is
    the physical upper bound set by the two-phase speed of sound.

    Parameters
    ----------
    m_dot_model : float
        Mass flow predicted by Dyer or HEM, kg/s.
    Cd : float
        Discharge coefficient.
    A : float
        Total orifice area, m^2.
    T_upstream : float
        Tank temperature, K.
    P_upstream : float
        Upstream pressure, Pa.
    x_inlet : float, optional
        Vapour quality at the orifice inlet, [0, 1].

    Returns
    -------
    dict
        "m_dot_real"  : physically realizable mass flow, kg/s
        "choked"      : True if choking limit was applied
        "m_dot_model" : original model prediction, kg/s
        "m_dot_crit"  : HEM critical flow limit, kg/s
        "crit_result" : full hem_critical_flow() output
    """
    crit = hem_critical_flow(Cd, A, T_upstream, P_upstream, x_inlet)
    m_crit    = crit["m_dot_crit"]
    choked    = m_dot_model > m_crit
    m_dot_out = m_crit if choked else m_dot_model
    return {
        "m_dot_real":  m_dot_out,
        "choked":      choked,
        "m_dot_model": m_dot_model,
        "m_dot_crit":  m_crit,
        "crit_result": crit,
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
        "m_dot_Dyer": Dyer-predicted mass flow rate, kg/s. UNCHANGED by
                     the Henry-Fauske ceiling below -- see "choked".
        "m_dot_SPI": SPI-only prediction at the same operating point, kg/s
                     (the "no time to vaporize" limit, for comparison).
        "m_dot_HEM": HEM-only prediction, kg/s (the full-equilibrium
                     limit, for comparison).
        "kappa": Dyer's non-equilibrium weighting parameter.
        "x_exit": vapor quality at the orifice exit, from the HEM
                  sub-calculation (dimensionless).
        "m_dot_crit_HF": Henry-Fauske non-equilibrium critical mass flow
                  ceiling, kg/s, or None if it could not be computed
                  (see "HF_unavailable_reason"). Surfaced as a DIAGNOSTIC
                  value, side by side with "m_dot_Dyer" -- NOT applied
                  automatically, per the design decision documented
                  below and in docs/future_work.md, Priority 1.
        "choked": True if "m_dot_Dyer" exceeds "m_dot_crit_HF" (i.e. the
                  Dyer blend predicts more flow than the non-equilibrium
                  choking ceiling allows). False if it does not exceed
                  it, or if the ceiling could not be computed.
        "HF_unavailable_reason": None if "m_dot_crit_HF" was computed
                  successfully; otherwise a short string explaining why
                  (e.g. T_upstream outside Table A.4's range).

    Note on why "m_dot_Dyer" is NOT automatically capped at
    "m_dot_crit_HF" (decided September 2026, docs/future_work.md,
    Priority 1). henry_fauske_critical_flow() is correctly implemented
    from a primary source (Henry & Fauske 1971, via Simoneau et al.
    1971) and is internally consistent (always sits above the
    equilibrium HEM ceiling; does not perturb any of the 4 validated
    Waxman operating points, since all 4 already sit below it). BUT at
    other operating points -- e.g. examples/example_01_sizing.md's
    conditions -- it DOES bind, cutting the predicted flow by roughly
    13-17%, and there is currently no experimental data point in this
    project's validation set where the ceiling actually changes the
    answer (all 4 Waxman points sit below it). Silently overriding
    "m_dot_Dyer" with a value that is theoretically well-founded but
    empirically unconfirmed in the regime where it matters would risk
    quietly changing already-published results without evidence.
    Instead, both values are returned, so calling code (full_system.py,
    the Streamlit interface) can choose to display a warning when
    "choked" is True, without silently changing the headline number.
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

    # Henry-Fauske non-equilibrium choking ceiling -- diagnostic only,
    # not applied automatically (see docstring above). Computed
    # best-effort: if it cannot be evaluated (e.g. T_upstream above
    # Table A.4's 307.33 K limit), Dyer's own result is still returned
    # unaffected, with the reason recorded rather than silently ignored.
    try:
        hf = henry_fauske_critical_flow(Cd, A, T_upstream, P_upstream)
        m_dot_crit_HF = hf["m_dot_crit"]
        choked = m_dot_Dyer > m_dot_crit_HF
        hf_unavailable_reason = None
    except (ValueError, RuntimeError) as e:
        m_dot_crit_HF = None
        choked = False
        hf_unavailable_reason = str(e)

    return {
        "m_dot_Dyer": m_dot_Dyer,
        "m_dot_SPI":  m_dot_SPI,
        "m_dot_HEM":  m_dot_HEM,
        "kappa":      kappa,
        "x_exit":     hem_result["x_exit"],
        "m_dot_crit_HF":         m_dot_crit_HF,
        "choked":                choked,
        "HF_unavailable_reason": hf_unavailable_reason,
    }


if __name__ == "__main__":
    # --- Validation case ---
    from n2o_properties import rho_liquid_sat, nu_vapor_sat

    M_N2O = 44.013  # kg/kmol

    T_upstream = 293.15   # K, 20 degC
    P_upstream = 55e5     # Pa, 55 bar (subcooled: P_sat(20 degC) ~= 51.4 bar)
    P_downstream = 20e5   # Pa, 20 bar chamber pressure
    Cd = 0.65
    A = 3.79e-6            # m^2

    rho_l_upstream = rho_liquid_sat(T_upstream)
    T_downstream = T_sat(P_downstream)
    rho_l_downstream = rho_liquid_sat(T_downstream)
    rho_v_downstream = M_N2O / nu_vapor_sat(T_downstream)

    result = dyer_mass_flow(Cd, A, T_upstream, P_upstream, P_downstream,
                             rho_l_upstream, rho_l_downstream, rho_v_downstream)

    print("Two-phase injector model (HEM + Dyer) -- example evaluation")
    print("-" * 60)
    print(f"Dyer-predicted mass flow: {result['m_dot_Dyer']*1000:.1f} g/s")
    print(f"SPI-predicted mass flow:  {result['m_dot_SPI']*1000:.1f} g/s")
    print(f"HEM-predicted mass flow:  {result['m_dot_HEM']*1000:.1f} g/s")

    # --- Waxman validation conditions, for the isenthalpic vs isentropic
    #     critical-flow comparison (added September 2026, Priority 1) ---
    print()
    print("=" * 60)
    print("Isenthalpic vs isentropic critical flow -- Waxman conditions")
    print("(T1 = 280 K, P1 = 4.36 MPa, D = 1.5 mm, Cd = 0.65)")
    print("=" * 60)

    T1 = 280.0
    P1 = 4.36e6
    Cd_w = 0.65
    D_w = 0.0015
    A_w = math.pi * (D_w / 2.0) ** 2

    crit_h = hem_critical_flow(Cd_w, A_w, T1, P1)
    crit_s = hem_critical_flow_isentropic(Cd_w, A_w, T1, P1)

    print(f"  Isenthalpic : m_dot_crit = {crit_h['m_dot_crit']*1000:6.2f} g/s  "
          f"at P2_crit = {crit_h['P2_crit']/1e5:5.2f} bar, "
          f"x_crit = {crit_h['x_crit']:.4f}")
    print(f"  Isentropic  : m_dot_crit = {crit_s['m_dot_crit']*1000:6.2f} g/s  "
          f"at P2_crit = {crit_s['P2_crit']/1e5:5.2f} bar, "
          f"x_crit = {crit_s['x_crit']:.4f}")
    diff_pct = 100.0 * (crit_s['m_dot_crit'] - crit_h['m_dot_crit']) / crit_h['m_dot_crit']
    print(f"  Difference  : {diff_pct:+.2f}%")
    print()
    print("  Reference (Waxman experimental, moderate dP): 44.0-48.0 g/s")
    print("  Both critical-flow ceilings should sit below the experimental")
    print("  Dyer-regime values -- Dyer's non-equilibrium correction")
    print("  legitimately predicts above either HEM-only ceiling.")

    # --- Henry-Fauske diagnostic (added September 2026, Priority 1) ---
    # Two contrasting cases: the Waxman validated point (should NOT be
    # flagged as choked -- the whole point of the earlier validation),
    # and a case further from Waxman (examples/example_01_sizing.md's
    # conditions) where the ceiling DOES bind, demonstrating the
    # side-by-side diagnostic rather than a silent override.
    print()
    print("=" * 60)
    print("Henry-Fauske diagnostic -- side by side with Dyer, not applied")
    print("=" * 60)

    cases = [
        ("Waxman (validated)", 280.0, 4.36e6, 30.4e5, math.pi * (0.0015/2)**2),
        ("Example-1-like (NOT choking-validated)", 293.15, 58e5, 22e5,
         6 * math.pi * (0.75e-3)**2),
    ]
    for label, T_u, P_u, P_d, A_c in cases:
        rho_l_u = rho_liquid_sat(T_u)
        T_d = T_sat(P_d)
        rho_l_d = rho_liquid_sat(T_d)
        rho_v_d = M_N2O / nu_vapor_sat(T_d)
        r = dyer_mass_flow(Cd, A_c, T_u, P_u, P_d, rho_l_u, rho_l_d, rho_v_d)
        print(f"\n  {label}:")
        print(f"    m_dot_Dyer      = {r['m_dot_Dyer']*1000:7.2f} g/s")
        if r["m_dot_crit_HF"] is not None:
            print(f"    m_dot_crit_HF   = {r['m_dot_crit_HF']*1000:7.2f} g/s")
        else:
            print(f"    m_dot_crit_HF   = unavailable ({r['HF_unavailable_reason']})")
        print(f"    choked          = {r['choked']}")
