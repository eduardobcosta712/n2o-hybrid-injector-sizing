"""
grain_sizing.py

OF ratio -> fuel mass flow -> initial fuel grain sizing, via the Marxman
regression rate correlation. See docs/03b_grain_sizing.md for the full
theoretical derivation and docs/references.md for sources.

SCOPE (deliberately limited -- read before extending):

  - This module sizes the INITIAL (t=0) operating point of a circular,
    single- or multi-port grain that delivers a target O/F ratio at a
    given (already-computed) oxidiser mass flow. It does NOT simulate
    the full transient burn (port radius, O/F, and thrust all drift as
    the grain regresses) -- that is future_work.md Priority 8.
  - Only CIRCULAR ports are supported. Non-circular grain cross-sections
    (star, wagon-wheel) are NOT implemented: unlike a circular port,
    their burning perimeter does not stay self-similar as they regress,
    and a rigorous treatment needs either published closed-form
    geometry for specific classical shapes or a numerical burnback
    simulation -- neither is in scope here. See future_work.md for this
    as a tracked, separate extension.
  - The regression-rate coefficients (a, n) are NOT shipped as fixed
    per-fuel defaults. Unlike density (a genuine material property),
    a and n vary by a factor of 2-3x between published studies of
    nominally the same fuel/oxidiser pair (different injector designs,
    scales, and test conditions) -- see docs/references.md for the
    specific numbers this claim is based on. Treating a literature (a, n)
    pair as a reliable design default would be actively misleading.
    a and n are therefore REQUIRED inputs to every function in this
    module; FUEL_PROPERTIES below provides density (a real material
    property, safe to default) and a *reference range* for (a, n) from
    the literature, clearly labelled as non-authoritative.
  - No specific-impulse estimate is provided: that requires a chemical
    equilibrium code (CEA, RPA, or similar), which this project does not
    implement or wrap. Obtain I_sp from CEA/RPA separately.

Units: SI throughout (Pa, K, kg, m, m^2, kg/s, m/s) except where the
Marxman correlation's own units are noted explicitly (mm/s is the
literature's conventional unit for regression rate; this module works
internally in m/s and converts at the boundary).
"""

import math


# ---------------------------------------------------------------------------
# Fuel properties -- density is a genuine material property (safe default).
# a, n reference ranges are for ORIENTATION ONLY -- see module docstring.
# All figures below are cited in docs/references.md.
# ---------------------------------------------------------------------------

FUEL_PROPERTIES = {
    "Paraffin wax": {
        "rho_kg_m3": 900.0,
        # Green & Perry-style "theoretical" (bulk, non-porous) density,
        # 900.0 kg/m^3; a same-source experimental (7% porosity) value of
        # 894.0 kg/m^3 is also reported -- see references.md. 900.0 used
        # as the default; override if a specific formulation is known.
        "a_n_reference_range": {
            "a_note": "Reported G_o-normalised regression rates for "
                      "paraffin/N2O span roughly 2-5 mm/s at G_o on the "
                      "order of 100-400 kg/(m^2.s) across independent "
                      "studies -- NOT a single (a, n) pair. Do not use "
                      "without checking against your own or closely "
                      "matched test data.",
            "n_typical": "~0.5 (near-universal empirical fit for paraffin; "
                         "classical Marxman theory predicts 0.8, but "
                         "liquid-layer entrainment lowers this in practice)",
        },
    },
    "HTPB": {
        "rho_kg_m3": 920.0,
        # Standard propulsion-literature value; HTPB's OTHER thermochemical
        # properties (heat of formation, etc.) are much less consistently
        # reported, but density is not contested in the same way.
        "a_n_reference_range": {
            "a_note": "Reported regression rates for HTPB/N2O are roughly "
                      "3-5x lower than paraffin at comparable G_o (e.g. "
                      "~0.6 mm/s at G_o = 350 kg/(m^2.s) in one study, vs. "
                      "~2 mm/s for solid paraffin at the same G_o) -- "
                      "again NOT a single (a, n) pair, see references.md.",
            "n_typical": "~0.5-0.8 depending on source",
        },
    },
    "ABS": {
        "rho_kg_m3": 1050.0,
        # Standard commercial ABS resin bulk density (~1.04-1.06 g/cm^3);
        # 3D-printed grains may have somewhat lower effective density
        # depending on infill -- override if measured.
        "a_n_reference_range": {
            "a_note": "ABS/N2O regression rate reported comparable to but "
                      "slightly below HTPB/N2O in direct side-by-side "
                      "testing at equal grain geometry -- treat literature "
                      "HTPB values as a rough starting point only, still "
                      "requiring your own confirmation.",
            "n_typical": "~0.5 (reported similar to HTPB)",
        },
    },
    "PMMA": {
        "rho_kg_m3": 1185.2,
        # Measured value from a specific characterisation study (clear
        # cast PMMA grain) -- see references.md. Commercial PMMA sheet is
        # typically quoted 1180-1190 kg/m^3, consistent with this.
        "a_n_reference_range": {
            "a_note": "PMMA is a common laboratory reference fuel "
                      "precisely because of low regression-rate scatter "
                      "between labs relative to paraffin -- but published "
                      "(a, n) pairs still differ meaningfully by source "
                      "and oxidiser; still requires your own confirmation.",
            "n_typical": "~0.5-0.6 depending on source",
        },
    },
}


def regression_rate(a, n, G_o):
    """
    Marxman regression rate correlation.

        r_dot = a * G_o^n

    Parameters
    ----------
    a : float
        Empirical regression-rate coefficient. UNITS DEPEND ON n and on
        the units G_o was fitted in in the source literature -- this is
        a well-known trap (Marxman coefficients are frequently quoted in
        mixed unit systems, e.g. G_o in kg/(m^2.s) with r_dot in mm/s).
        This function is unit-agnostic: pass a, G_o, and interpret the
        returned r_dot consistently with whatever unit system a was
        fitted in. This project's own usage (grain_sizing.py) always
        passes SI G_o [kg/(m^2.s)] and expects SI r_dot [m/s] -- convert
        a literature coefficient to that convention before calling this,
        and record the conversion, since silently mixing unit systems
        here is exactly the kind of error this project's "no magic
        numbers without a traceable source" convention exists to catch.
    n : float
        Regression rate exponent, dimensionless. Typically 0.4-0.8
        depending on fuel (see FUEL_PROPERTIES reference notes).
    G_o : float
        Oxidiser mass flux through the port, kg/(m^2.s).

    Returns
    -------
    float
        Regression rate, in the units implied by `a` (see above -- this
        project's own callers use m/s).

    Raises
    ------
    ValueError
        If G_o <= 0 (undefined for a fractional power with G_o <= 0 in
        general, and physically meaningless -- no flow, no regression).
    """
    if G_o <= 0:
        raise ValueError(
            f"G_o = {G_o} must be positive -- the regression rate "
            "correlation is undefined (and physically meaningless) for "
            "zero or negative oxidiser mass flux."
        )
    return a * G_o ** n


def oxidizer_mass_flux(m_dot_ox_per_port, r_port):
    """
    Oxidiser mass flux through a single circular port.

        G_o = m_dot_ox_per_port / (pi * r_port^2)

    Parameters
    ----------
    m_dot_ox_per_port : float
        Oxidiser mass flow rate through this one port, kg/s (total
        oxidiser mass flow divided by the number of ports, if more than
        one -- see fuel_mass_flow_rate below for the multi-port case).
    r_port : float
        Port radius, m.

    Returns
    -------
    float
        Oxidiser mass flux, kg/(m^2.s).
    """
    if r_port <= 0:
        raise ValueError(f"r_port = {r_port} must be positive.")
    A_port = math.pi * r_port ** 2
    return m_dot_ox_per_port / A_port


def fuel_mass_flow_rate_per_port(a, n, rho_fuel, L, r_port, m_dot_ox_per_port):
    """
    Fuel mass flow rate from a single circular port, combining the
    regression rate with the burning surface geometry:

        m_dot_fuel = rho_fuel * (2 * pi * r_port) * L * r_dot(G_o)

    i.e. (burning perimeter) x (grain length) x (regression rate) x
    (fuel density) -- the volumetric burn rate of a thin annular shell
    at the port surface, converted to mass via density.

    Parameters
    ----------
    a, n : float
        Marxman correlation coefficients (SI: G_o in kg/(m^2.s), regression
        rate in m/s -- see regression_rate()'s docstring on unit care).
    rho_fuel : float
        Solid fuel density, kg/m^3.
    L : float
        Grain length (this port's burning length), m.
    r_port : float
        Current port radius, m.
    m_dot_ox_per_port : float
        Oxidiser mass flow rate through this one port, kg/s.

    Returns
    -------
    float
        Fuel mass flow rate from this port, kg/s.
    """
    G_o = oxidizer_mass_flux(m_dot_ox_per_port, r_port)
    r_dot = regression_rate(a, n, G_o)
    perimeter = 2.0 * math.pi * r_port
    return rho_fuel * perimeter * L * r_dot


def solve_initial_port_radius(m_dot_fuel_target_per_port, a, n, rho_fuel, L,
                               m_dot_ox_per_port,
                               r_min=1e-4, r_max=1.0, tol=1e-9, max_iter=200):
    """
    Solve for the initial port radius r_0 that delivers a target fuel
    mass flow rate (per port), given the Marxman correlation and grain
    geometry -- the core sizing step of this module.

    Physical structure of the equation being solved (see
    docs/03b_grain_sizing.md for the full derivation): substituting
    G_o = m_dot_ox_per_port / (pi r^2) into Marxman and combining terms,

        m_dot_fuel_per_port(r) = K * r^(1 - 2n),
        K = 2 * pi * rho_fuel * L * a * pi^(-n) * m_dot_ox_per_port^n

    For n != 0.5 this inverts in closed form,
        r_0 = (m_dot_fuel_target_per_port / K) ** (1 / (1 - 2n)),
    but this function solves it by bisection regardless (matching this
    project's existing convention for transcendental equations -- see
    n2o_properties.T_sat and injector_two_phase.henry_fauske_critical_flow),
    since it needs no special-casing at n = 0.5 (where the closed form's
    exponent divides by zero) and provides a natural cross-check against
    the closed-form result in testing.

    IMPORTANT, EASY-TO-MISS SUBTLETY: the sign of the exponent (1 - 2n)
    flips at n = 0.5. For n < 0.5 (exponent positive), a larger target
    fuel flow requires a LARGER radius -- the intuitive direction. For
    n > 0.5 (exponent negative, e.g. many HTPB/N2O fits), a larger
    target flow requires a SMALLER radius: increasing r reduces G_o
    faster (as 1/r^2) than it increases the burning perimeter (as r),
    so the net effect on regression rate -- and hence mass flow -- is a
    decrease. This is verified directly in
    test_grain_sizing.py::TestSolveInitialPortRadius (both directions,
    not just assumed) precisely because it is easy to get backwards.

    Special case n = 0.5: m_dot_fuel_per_port becomes INDEPENDENT of
    r_port (the exponent 1 - 2n = 0) -- the target flow rate either
    matches K exactly (any r_0 works; this function returns the
    midpoint of the search bracket) or is unreachable at ANY port radius
    for the given a, L (would require changing L or a instead). This
    project's bisection will simply fail to bracket a root in the
    unreachable case, raising RuntimeError with a clear explanation
    rather than returning a meaningless value.

    Parameters
    ----------
    m_dot_fuel_target_per_port : float
        Required fuel mass flow rate from this one port, kg/s (total
        target fuel mass flow divided by the number of ports).
    a, n : float
        Marxman correlation coefficients (SI convention -- see
        regression_rate()'s docstring).
    rho_fuel : float
        Solid fuel density, kg/m^3.
    L : float
        Grain length (this port's burning length), m.
    m_dot_ox_per_port : float
        Oxidiser mass flow rate through this one port, kg/s.
    r_min, r_max : float, optional
        Search bracket for the port radius, m. Defaults span 0.1 mm to
        1 m, which comfortably covers realistic hybrid motor scales;
        widen if sizing an unusually large or small motor.
    tol : float, optional
        Bisection convergence tolerance on r_port, m.
    max_iter : int, optional
        Maximum bisection iterations.

    Returns
    -------
    float
        Initial port radius r_0, m.

    Raises
    ------
    ValueError
        If m_dot_fuel_target_per_port <= 0, or a, rho_fuel, L,
        m_dot_ox_per_port are not positive.
    RuntimeError
        If no root is bracketed in [r_min, r_max] -- the target fuel
        flow rate is not achievable anywhere in this radius range with
        the given a, n, rho_fuel, L (e.g. n very close to 0.5 and the
        target not equal to the (nearly) r-independent achievable rate;
        or L too short/long for any reasonable radius to reach the
        target). The message reports the achievable range found at the
        bracket endpoints so the next parameter to adjust is obvious.
    """
    if m_dot_fuel_target_per_port <= 0:
        raise ValueError(
            f"m_dot_fuel_target_per_port = {m_dot_fuel_target_per_port} "
            "must be positive."
        )
    for name, val in (("a", a), ("rho_fuel", rho_fuel), ("L", L),
                      ("m_dot_ox_per_port", m_dot_ox_per_port)):
        if val <= 0:
            raise ValueError(f"{name} = {val} must be positive.")

    def residual(r):
        return (fuel_mass_flow_rate_per_port(a, n, rho_fuel, L, r,
                                              m_dot_ox_per_port)
                - m_dot_fuel_target_per_port)

    r_lo, r_hi = r_min, r_max
    f_lo, f_hi = residual(r_lo), residual(r_hi)

    if f_lo * f_hi > 0:
        raise RuntimeError(
            f"solve_initial_port_radius: no root bracketed in "
            f"[{r_min:.4g}, {r_max:.4g}] m. Achievable fuel mass flow at "
            f"the bracket endpoints: {r_min:.4g} m -> "
            f"{fuel_mass_flow_rate_per_port(a, n, rho_fuel, L, r_min, m_dot_ox_per_port)*1000:.3f} g/s, "
            f"{r_max:.4g} m -> "
            f"{fuel_mass_flow_rate_per_port(a, n, rho_fuel, L, r_max, m_dot_ox_per_port)*1000:.3f} g/s "
            f"(target: {m_dot_fuel_target_per_port*1000:.3f} g/s). "
            "If n is close to 0.5, fuel mass flow is nearly independent "
            "of radius -- adjust L or a instead of radius. Otherwise, "
            "widen r_min/r_max."
        )

    for _ in range(max_iter):
        r_mid = 0.5 * (r_lo + r_hi)
        f_mid = residual(r_mid)
        if f_lo * f_mid <= 0:
            r_hi, f_hi = r_mid, f_mid
        else:
            r_lo, f_lo = r_mid, f_mid
        if (r_hi - r_lo) < tol:
            break

    return 0.5 * (r_lo + r_hi)


def size_grain(m_dot_ox, OF, a, n, rho_fuel, L, N_ports=1, burn_time=None):
    """
    Full initial-point grain sizing: from the already-computed oxidiser
    mass flow and a target O/F ratio, find the required fuel mass flow,
    the initial port radius (assuming N_ports identical circular ports,
    each carrying an equal share of the oxidiser flow), and -- if a burn
    duration is given -- a first-order (deliberately conservative)
    estimate of how far the port grows over the burn.

    Parameters
    ----------
    m_dot_ox : float
        Total oxidiser mass flow rate, kg/s (from the existing
        full_system.py solver).
    OF : float
        Target oxidiser-to-fuel mass ratio (from thermochemical sizing,
        e.g. CEA/RPA, external to this project).
    a, n : float
        Marxman correlation coefficients (SI convention -- REQUIRED, see
        module docstring on why no default is provided).
    rho_fuel : float
        Solid fuel density, kg/m^3 (see FUEL_PROPERTIES for reference
        values by fuel type).
    L : float
        Grain length, m -- a direct INPUT in this implementation (see
        module docstring: deliberately not derived from an unsourced
        L/D heuristic).
    N_ports : int, optional
        Number of identical circular ports, sharing the oxidiser flow
        equally. Default 1 (single port). More ports increase total
        burning perimeter for the same total port area (perimeter
        scales as sqrt(N) at fixed total area), increasing fuel mass
        flow for the same regression rate -- a real, commonly used
        design lever, distinct from (and much simpler than) non-circular
        port shapes (see module docstring on why those are out of scope).
    burn_time : float or None, optional
        Target burn duration, s. If provided, a FIRST-ORDER estimate of
        the port radius at the end of the burn is included in the
        result, using the INITIAL regression rate held constant
        (r_final ~= r_0 + r_dot_0 * burn_time). This deliberately
        OVER-estimates how far the port grows: G_o falls as the port
        grows (same m_dot_ox over a larger area), and since n > 0 the
        true regression rate falls over the burn too -- so the real
        final radius is smaller than this estimate. Treat it as a
        conservative (safety-margin) sizing check for the grain's outer
        radius, NOT a transient burn simulation (that is
        future_work.md Priority 8). If None, no burnback estimate is
        computed.

    Returns
    -------
    dict
        "m_dot_fuel"      : total fuel mass flow rate, kg/s
        "m_dot_fuel_per_port" : fuel mass flow rate per port, kg/s
        "m_dot_ox_per_port"   : oxidiser mass flow rate per port, kg/s
        "r_0"             : initial port radius (each port), m
        "G_o_0"           : initial oxidiser mass flux, kg/(m^2.s)
        "r_dot_0"         : initial regression rate, m/s
        "N_ports"         : number of ports (passed through)
        "L"               : grain length (passed through), m
        "r_final_estimate": first-order conservative estimate of final
                             port radius after burn_time, m, or None if
                             burn_time was not provided
        "fuel_mass_consumed_estimate": total fuel mass burned over
                             burn_time (all ports), kg, using the same
                             first-order estimate, or None
    """
    if m_dot_ox <= 0:
        raise ValueError(f"m_dot_ox = {m_dot_ox} must be positive.")
    if OF <= 0:
        raise ValueError(f"OF = {OF} must be positive.")
    if N_ports < 1:
        raise ValueError(f"N_ports = {N_ports} must be >= 1.")

    m_dot_fuel = m_dot_ox / OF
    m_dot_fuel_per_port = m_dot_fuel / N_ports
    m_dot_ox_per_port = m_dot_ox / N_ports

    r_0 = solve_initial_port_radius(
        m_dot_fuel_per_port, a, n, rho_fuel, L, m_dot_ox_per_port)

    G_o_0 = oxidizer_mass_flux(m_dot_ox_per_port, r_0)
    r_dot_0 = regression_rate(a, n, G_o_0)

    r_final_estimate = None
    fuel_mass_consumed_estimate = None
    if burn_time is not None:
        if burn_time <= 0:
            raise ValueError(f"burn_time = {burn_time} must be positive.")
        r_final_estimate = r_0 + r_dot_0 * burn_time
        fuel_mass_consumed_estimate = (
            rho_fuel * math.pi * (r_final_estimate ** 2 - r_0 ** 2)
            * L * N_ports
        )

    return {
        "m_dot_fuel":            m_dot_fuel,
        "m_dot_fuel_per_port":   m_dot_fuel_per_port,
        "m_dot_ox_per_port":     m_dot_ox_per_port,
        "r_0":                   r_0,
        "G_o_0":                 G_o_0,
        "r_dot_0":                r_dot_0,
        "N_ports":               N_ports,
        "L":                     L,
        "r_final_estimate":      r_final_estimate,
        "fuel_mass_consumed_estimate": fuel_mass_consumed_estimate,
    }


if __name__ == "__main__":
    # --- Validation / sanity checks ---
    print("grain_sizing.py -- validation")
    print("=" * 60)

    # Example: paraffin/N2O, illustrative (a, n) chosen to land near a
    # realistic port radius (~2 cm) and G_o (~O(100-400) kg/(m^2.s)) for
    # a 500 g/s-class motor, at n = 0.6 -- deliberately NOT n = 0.5 (the
    # degenerate case is checked separately below). NOT a design value
    # -- see module docstring.
    a_illustrative = 8.24e-5
    n_illustrative = 0.6
    rho_paraffin = FUEL_PROPERTIES["Paraffin wax"]["rho_kg_m3"]

    m_dot_ox = 0.5       # kg/s, 500 g/s, matching several examples/ cases
    OF = 6.0              # typical paraffin/N2O design O/F is O(5-8)
    L = 0.25               # m, 250 mm grain -- illustrative
    N_ports = 1

    result = size_grain(m_dot_ox, OF, a_illustrative, n_illustrative,
                        rho_paraffin, L, N_ports=N_ports, burn_time=8.0)

    print(f"Inputs: m_dot_ox = {m_dot_ox*1000:.0f} g/s, OF = {OF}, "
          f"a = {a_illustrative:.4e} (SI), n = {n_illustrative}, "
          f"rho_fuel = {rho_paraffin:.0f} kg/m^3, L = {L*1000:.0f} mm, "
          f"N_ports = {N_ports}")
    print("-" * 60)
    print(f"m_dot_fuel        = {result['m_dot_fuel']*1000:.2f} g/s")
    print(f"r_0 (per port)    = {result['r_0']*1000:.2f} mm")
    print(f"G_o_0             = {result['G_o_0']:.1f} kg/(m^2.s)")
    print(f"r_dot_0           = {result['r_dot_0']*1000:.3f} mm/s")
    if result["r_final_estimate"] is not None:
        print(f"r_final (8s, 1st-order conservative estimate) = "
              f"{result['r_final_estimate']*1000:.2f} mm")
        print(f"Fuel mass consumed (estimate) = "
              f"{result['fuel_mass_consumed_estimate']*1000:.1f} g")

    # Sanity check: closed-form inversion should match the bisection
    # result for n != 0.5 (use n=0.6 here specifically to exercise the
    # non-degenerate closed form)
    print()
    print("Cross-check vs. closed-form solution (n != 0.5):")
    n_cf = 0.6
    a_cf = a_illustrative
    m_dot_fuel_pp = (m_dot_ox / OF) / N_ports
    m_dot_ox_pp = m_dot_ox / N_ports
    K = 2 * math.pi * rho_paraffin * L * a_cf * math.pi ** (-n_cf) * m_dot_ox_pp ** n_cf
    r0_closed_form = (m_dot_fuel_pp / K) ** (1.0 / (1.0 - 2 * n_cf))
    r0_bisection = solve_initial_port_radius(
        m_dot_fuel_pp, a_cf, n_cf, rho_paraffin, L, m_dot_ox_pp)
    print(f"  Closed form: r_0 = {r0_closed_form*1000:.4f} mm")
    print(f"  Bisection:   r_0 = {r0_bisection*1000:.4f} mm")
    print(f"  Difference:  {abs(r0_closed_form-r0_bisection)/r0_closed_form*100:.4e} %")

    # Sanity check: n = 0.5 special case -- fuel flow independent of radius
    print()
    print("Special case n = 0.5 (fuel flow independent of port radius):")
    n_half = 0.5
    m1 = fuel_mass_flow_rate_per_port(a_illustrative, n_half, rho_paraffin,
                                       L, 0.01, m_dot_ox_pp)
    m2 = fuel_mass_flow_rate_per_port(a_illustrative, n_half, rho_paraffin,
                                       L, 0.03, m_dot_ox_pp)
    print(f"  m_dot_fuel at r=10mm: {m1*1000:.4f} g/s")
    print(f"  m_dot_fuel at r=30mm: {m2*1000:.4f} g/s")
    print(f"  Ratio: {m2/m1:.6f} (should be ~1.0, independent of r)")
