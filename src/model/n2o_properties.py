"""
n2o_properties.py

Saturated thermophysical properties of N2O as a function of temperature.

This module has two sources of data, kept clearly separated below:

1. CLOSED-FORM CORRELATIONS (P_sat, T_sat, rho_liquid_sat) -- originally
   from Perry's Chemical Engineers' Handbook (Green & Perry, 2008),
   re-transcribed in:
       Jean-Philyppe, J. (2023). "A computational model for the design of
       a nitrous oxide-paraffin wax hybrid rocket engine." McGill Rocket
       Team. arXiv:2302.06725, Appendix A.1.

2. TABULATED SATURATED VAPOR/LIQUID PROPERTIES (nu_vapor_sat, h_liquid_sat,
   h_vapor_sat, h_fg) -- needed for the HEM/Dyer two-phase injector model
   (docs/03_two_phase_flow.md), where no simple closed-form correlation
   for the saturated vapor molar volume was available. Read from
   n2o_saturation_table.csv, transcribed by hand from the same source
   (arXiv:2302.06725, Appendix A, Table A.1), via linear interpolation --
   mirroring the interpolation scheme the source paper itself uses.
   See that CSV's header for the full transcription/unit notes.

Valid temperature range: 182.33 K - 309.52 K (triple point to near critical point),
for BOTH the correlations and the tabulated data (the table's own T range).
All functions raise ValueError outside this range rather than silently
extrapolating, since neither the fit nor the table is guaranteed valid there.

IMPORTANT EXCEPTION -- Table A.4 (cp_l, mu_l, s_l, s_v; NIST WebBook) does
NOT cover the full range above: it stops at 307.33 K, short of the
309.52 K critical-point limit used everywhere else in this module. This
is a real, separate constraint (see T_MIN_A4/T_MAX_A4 and
_check_range_a4 below), not an oversight -- it reflects the actual extent
of the downloaded NIST data. Functions built on Table A.4
(cp_liquid_sat, mu_liquid_sat, s_liquid_sat, s_vapor_sat, s_fg) enforce
this narrower range explicitly.

Units: T in Kelvin, P in Pa, molar volumes in m^3/kmol, molar enthalpies
in kJ/kmol, molar entropies in kJ/(kmol.K) [numerically identical to
J/(mol.K), the unit the source table uses -- no factor-of-1000 conversion
needed, unlike h_l/h_v, since kJ/mol and J/mol.K/1000 mol/kmol cancel
exactly] throughout this module. See the CSV header for the raw-file unit
trap on enthalpy.
"""

import csv
import math
import os

# Valid temperature range for all correlations in this module (Kelvin)
T_MIN = 182.33
T_MAX = 309.52

# Valid temperature range for functions built on Table A.4 only (NIST
# WebBook, Lemmon & Span 2006): cp_liquid_sat, mu_liquid_sat,
# s_liquid_sat, s_vapor_sat, s_fg. Narrower than [T_MIN, T_MAX] above --
# Table A.4's last row is 307.33 K, not 309.52 K. Added September 2026
# while implementing the isentropic choking limit (Priority 1), which is
# the first place this distinction actually matters: the isentropic scan
# evaluates entropy at temperatures up to T_upstream, and a design point
# with T_upstream between 307.33 K and 309.52 K would previously have
# passed the (too permissive) _check_range(T) check and then hit an
# opaque "should be unreachable" RuntimeError inside _interp.
T_MIN_A4 = 182.33
T_MAX_A4 = 307.33

# Critical point (reference values, for sanity checks / user-facing warnings)
T_CRIT = 309.52    # K  (~36.4 degC)
P_CRIT = 7.245e6   # Pa (~72.45 bar) -- consistent with the correlation at T_CRIT

# Path to the saturation table CSV, relative to this file's own location
# (not the current working directory), so the module works regardless of
# where the caller's script is run from.
_TABLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "n2o_saturation_table.csv")

# Module-level cache for the parsed table, populated on first use by
# _load_saturation_table(). Avoids re-reading and re-parsing the CSV file
# from disk on every single function call.
_table_cache = None


def _check_range(T):
    if not (T_MIN <= T <= T_MAX):
        raise ValueError(
            f"Temperature {T:.2f} K is outside the correlation's valid range "
            f"[{T_MIN}, {T_MAX}] K. Results outside this range are not reliable."
        )


def _check_range_a4(T):
    """
    Range check for functions based on Table A.4 (NIST WebBook), which is
    narrower than the module's main [T_MIN, T_MAX] range -- see the
    T_MIN_A4/T_MAX_A4 module docstring note above. Raising a specific,
    named error here (rather than falling through to the correlation's
    wider _check_range and then hitting _interp's generic "should be
    unreachable" RuntimeError) is the "fail loudly" convention used
    throughout this project: the caller should immediately understand
    which table ran out of data, not just that something, somewhere,
    was inconsistent.
    """
    if not (T_MIN_A4 <= T <= T_MAX_A4):
        raise ValueError(
            f"Temperature {T:.2f} K is outside Table A.4's valid range "
            f"[{T_MIN_A4}, {T_MAX_A4}] K (NIST WebBook, Lemmon & Span 2006). "
            f"This table is narrower than the main correlation's range "
            f"([{T_MIN}, {T_MAX}] K): it stops short of the critical point. "
            "Affects cp_liquid_sat, mu_liquid_sat, s_liquid_sat, "
            "s_vapor_sat, and s_fg."
        )


def _load_saturation_table():
    """
    Parse n2o_saturation_table.csv and return Table A.1 and Table A.3
    as parallel lists, ready for interpolation. Result is cached in
    _table_cache after the first call.

    The CSV has three sections:
      [TABLE_A1] — thermodynamic properties (McGill/Perry)
      [TABLE_A2] — derivatives of Table A.1 (not parsed)
      [TABLE_A3] — saturated vapour dynamic viscosity mu_v (NIST WebBook,
                   Lemmon & Span 2006 EOS + Millat et al. 1991 viscosity)

    Returns
    -------
    dict
        "T_K"   : list of float, K  (Table A.1 temperatures)
        "nu_v"  : list of float, m^3/kmol
        "h_l"   : list of float, kJ/kmol
        "h_v"   : list of float, kJ/kmol
        "T_K_muv" : list of float, K  (Table A.3 temperatures)
        "mu_v"  : list of float, Pa·s (converted from μPa·s)
    """
    global _table_cache
    if _table_cache is not None:
        return _table_cache

    T_list, nu_v_list, h_l_list, h_v_list = [], [], [], []
    T_muv_list, mu_v_list = [], []

    with open(_TABLE_PATH, "r") as f:
        lines = f.readlines()

    # --- Parse TABLE_A1 ---
    in_a1 = False
    header_skipped = False
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s == "[TABLE_A1]":
            in_a1 = True; header_skipped = False; continue
        if s in ("[TABLE_A2]", "[TABLE_A3]"):
            in_a1 = False; continue
        if in_a1:
            if not header_skipped:
                header_skipped = True; continue
            fields = next(csv.reader([s]))
            T_list.append(float(fields[0]))
            nu_v_list.append(float(fields[1]))
            # h_l and h_v are in kJ/mol in source; convert to kJ/kmol (×1000)
            h_l_list.append(float(fields[3]) * 1000.0)
            h_v_list.append(float(fields[5]) * 1000.0)

    # --- Parse TABLE_A3 (mu_v from NIST) ---
    in_a3 = False
    header_skipped_a3 = False
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s == "[TABLE_A3]":
            in_a3 = True; header_skipped_a3 = False; continue
        if s == "[TABLE_A4]":
            in_a3 = False; continue
        if in_a3:
            if not header_skipped_a3:
                header_skipped_a3 = True; continue
            fields = next(csv.reader([s]))
            T_muv_list.append(float(fields[0]))
            mu_v_list.append(float(fields[1]) * 1e-6)

    # --- Parse TABLE_A4 (cp_l, mu_l, s_l, s_v from NIST) ---
    T_a4, cp_l_list, mu_l_list, s_l_list, s_v_list = [], [], [], [], []
    in_a4 = False
    header_a4 = False
    for line in lines:
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        if s == "[TABLE_A4]":
            in_a4 = True; header_a4 = False; continue
        if in_a4:
            if not header_a4:
                header_a4 = True; continue
            fields = next(csv.reader([s]))
            T_a4.append(float(fields[0]))
            cp_l_list.append(float(fields[1]))   # J/(kg.K)
            mu_l_list.append(float(fields[2]))   # Pa.s
            s_l_list.append(float(fields[3]))    # J/(mol.K)
            s_v_list.append(float(fields[4]))    # J/(mol.K)

    _table_cache = {
        "T_K":     T_list,
        "nu_v":    nu_v_list,
        "h_l":     h_l_list,
        "h_v":     h_v_list,
        "T_K_muv": T_muv_list,
        "mu_v":    mu_v_list,
        "T_K_a4":  T_a4,
        "cp_l":    cp_l_list,   # J/(kg.K)
        "mu_l":    mu_l_list,   # Pa.s
        "s_l":     s_l_list,    # J/(mol.K)
        "s_v":     s_v_list,    # J/(mol.K)
    }
    return _table_cache


def _interp(T, x_values, y_values):
    """
    Linear interpolation of y(T), given the table's x_values (temperatures,
    strictly ascending) and y_values (the property at each x_values[i]).

    This is a thin wrapper around the interpolation formula itself:
        y(T) = y_i + (y_(i+1) - y_i) / (x_(i+1) - x_i) * (T - x_i)
    for the bracketing pair x_i <= T <= x_(i+1). Implemented by hand
    (rather than e.g. numpy.interp) to keep this module dependency-free,
    consistent with the rest of the project so far.

    Parameters
    ----------
    T : float
        Temperature at which to evaluate y, in Kelvin. Must already be
        range-checked by the caller (via _check_range) before calling this.
    x_values : list of float
        Table's temperature column, strictly ascending.
    y_values : list of float
        Table's property column, same length and order as x_values.

    Returns
    -------
    float
        Linearly interpolated y(T).
    """
    # Find the bracketing interval [x_values[i], x_values[i+1]] containing T.
    # A linear scan is used (not binary search) since the table only has
    # ~27 rows -- performance is a non-issue at this size, and a plain
    # scan is easier to read and verify by hand than a bisection routine.
    for i in range(len(x_values) - 1):
        x_i, x_ip1 = x_values[i], x_values[i + 1]
        if x_i <= T <= x_ip1:
            y_i, y_ip1 = y_values[i], y_values[i + 1]
            fraction = (T - x_i) / (x_ip1 - x_i)
            return y_i + fraction * (y_ip1 - y_i)

    # Should be unreachable if _check_range already confirmed T is inside
    # [T_MIN, T_MAX], since T_MIN/T_MAX are themselves the table's first
    # and last rows -- but raise loudly rather than silently returning a
    # wrong value if this invariant is ever broken.
    raise RuntimeError(
        f"T = {T:.2f} K fell outside the saturation table's rows despite "
        "passing the module's range check -- this indicates a mismatch "
        "between T_MIN/T_MAX and the table's actual first/last row."
    )


def P_sat(T):
    """
    Saturation (vapor) pressure of N2O at temperature T.

    P(T) = exp(c1 + c2/T + c3*ln(T) + c4*T^c5)

    Parameters
    ----------
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Saturation pressure in Pa.
    """
    _check_range(T)
    c1, c2, c3, c4, c5 = 96.512, -4045.0, -12.277, 2.886e-5, 2.0
    return math.exp(c1 + c2 / T + c3 * math.log(T) + c4 * T ** c5)


def dP_sat_dT(T):
    """
    Derivative dP_sat/dT, analytically differentiated from P_sat(T).

    Useful later for e.g. relating small temperature changes to pressure
    margin changes without finite-differencing P_sat numerically.

    Returns
    -------
    float
        dP/dT in Pa/K.
    """
    _check_range(T)
    c2, c3, c4, c5 = -4045.0, -12.277, 2.886e-5, 2.0
    return P_sat(T) * (-c2 / T ** 2 + c3 / T + c4 * c5 * T ** (c5 - 1))


def T_sat(P, T_guess=250.0, tol=1e-6, max_iter=100, max_step=20.0):
    """
    Saturation temperature corresponding to a given pressure P — the
    inverse of P_sat(T). Solved numerically (Newton's method) since the
    correlation is not analytically invertible in closed form.

    This is the function used to compute the subcooling margin
    delta_T_sub = T_sat(P) - T at any point along the feed system,
    per Section 1.4 of the theory documents.

    Robustness note: P_sat(T) becomes strongly non-linear (steep) near
    the critical point, so a plain Newton step can occasionally overshoot
    past the correlation's valid range when the initial guess is far from
    the true root. To guard against this, each step is clamped to at most
    `max_step` Kelvin, and clamped again to stay within [T_MIN, T_MAX].
    This trades a few extra iterations for guaranteed robustness.

    Parameters
    ----------
    P : float
        Pressure in Pa.
    T_guess : float
        Initial guess for the iteration, in Kelvin.
    tol : float
        Convergence tolerance on pressure, in Pa.
    max_iter : int
        Maximum Newton iterations before raising an error.
    max_step : float
        Maximum allowed change in T per iteration, in Kelvin.

    Returns
    -------
    float
        Saturation temperature in Kelvin.
    """
    T = T_guess
    for _ in range(max_iter):
        _check_range(T)
        residual = P_sat(T) - P
        if abs(residual) < tol:
            return T
        derivative = dP_sat_dT(T)
        step = residual / derivative
        # Clamp the step so a single iteration can't overshoot far outside
        # the valid range, even where the curve is steep (near T_crit).
        step = max(-max_step, min(max_step, step))
        T_new = T - step
        # Clamp the resulting T to stay strictly inside the valid range.
        T_new = max(T_MIN + 1e-6, min(T_MAX - 1e-6, T_new))
        T = T_new
    raise RuntimeError(
        f"T_sat(P={P:.1f} Pa) did not converge after {max_iter} iterations."
    )


def rho_liquid_sat(T):
    """
    Saturated liquid density of N2O at temperature T.

    The correlation gives molar volume, nu_l(T):
        nu_l(T) = c2 ** (1 + (1 - T/c3)^c4) / c1        [m^3/kmol]
    Density is then rho = M / nu_l, with M the molar mass of N2O.

    Parameters
    ----------
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Saturated liquid density in kg/m^3.
    """
    _check_range(T)
    c1, c2, c3, c4 = 2.781, 0.27244, 309.57, 0.2882
    nu_l = c2 ** (1 + (1 - T / c3) ** c4) / c1  # m^3/kmol
    M_N2O = 44.013  # kg/kmol
    return M_N2O / nu_l


def nu_vapor_sat(T):
    """
    Saturated vapor molar volume of N2O at temperature T.

    Unlike rho_liquid_sat, this is NOT a closed-form correlation: no
    reliable closed-form fit for N2O's saturated vapor molar volume was
    found in the sources available for this project (an ideal-gas
    approximation was considered and rejected -- see docs/04_implementation.md
    -- because it degrades precisely near the critical point, which is
    this project's region of greatest interest). Instead, this value is
    obtained by linear interpolation of Table A.1 (see module docstring
    and n2o_saturation_table.csv), the same scheme the source paper
    itself uses.

    Parameters
    ----------
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Saturated vapor molar volume in m^3/kmol.
    """
    _check_range(T)
    table = _load_saturation_table()
    return _interp(T, table["T_K"], table["nu_v"])


def h_liquid_sat(T):
    """
    Saturated liquid molar enthalpy of N2O at temperature T, from Table A.1
    via linear interpolation. See module docstring for the source and the
    kJ/mol -> kJ/kmol conversion applied when the table is loaded.

    Parameters
    ----------
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Saturated liquid molar enthalpy, h_l, in kJ/kmol.
    """
    _check_range(T)
    table = _load_saturation_table()
    return _interp(T, table["T_K"], table["h_l"])


def h_vapor_sat(T):
    """
    Saturated vapor molar enthalpy of N2O at temperature T, from Table A.1
    via linear interpolation. See module docstring for the source and the
    kJ/mol -> kJ/kmol conversion applied when the table is loaded.

    Parameters
    ----------
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Saturated vapor molar enthalpy, h_v, in kJ/kmol.
    """
    _check_range(T)
    table = _load_saturation_table()
    return _interp(T, table["T_K"], table["h_v"])


def h_fg(T):
    """
    Latent heat of vaporization of N2O at temperature T,
        h_fg(T) = h_v(T) - h_l(T)
    i.e. the enthalpy difference between saturated vapor and saturated
    liquid at the same temperature -- the energy required to vaporize a
    unit (molar) amount of liquid at that temperature (Section 1.6/3.2 of
    the theory docs). Computed directly from the two tabulated enthalpies
    (not via Clausius-Clapeyron), which avoids amplifying numerical noise
    from (nu_v - nu_l) near the critical point, where that difference
    shrinks toward zero -- see docs/04_implementation.md for the
    Clausius-Clapeyron alternative kept as a cross-validation check.

    Parameters
    ----------
    T : float
        Temperature in Kelvin.

    Returns
    -------
    float
        Latent heat of vaporization in kJ/kmol.
    """
    return h_vapor_sat(T) - h_liquid_sat(T)


# ---------------------------------------------------------------------------
# Dynamic viscosity of saturated N2O (liquid and vapour)
# ---------------------------------------------------------------------------

# MU_LIQUID_N2O: fallback constant for backwards compatibility.
# Now that Table A.4 provides temperature-dependent mu_l(T), use
# mu_liquid_sat(T) wherever T is known. This constant is kept for
# cases where T is not available (e.g. module-level defaults).
MU_LIQUID_N2O = 1.5e-4   # Pa.s, mid-range estimate (conservative)



def cp_liquid_sat(T):
    """
    Isobaric heat capacity of saturated liquid N2O, J/(kg.K).

    Interpolated from Table A.4 (NIST WebBook, Lemmon & Span 2006).
    Used in the Henry-Fauske critical flow model (injector_two_phase.py).

    Physical note: cp_l increases sharply near the critical point
    (309.52 K), diverging to infinity at T_crit. Values above 300 K
    should be used with caution in engineering models.

    Range note: Table A.4 stops at 307.33 K (T_MAX_A4), not the module's
    overall T_MAX = 309.52 K -- see _check_range_a4.

    Parameters
    ----------
    T : float
        Temperature, K.

    Returns
    -------
    float
        Liquid heat capacity, J/(kg.K).
    """
    _check_range_a4(T)
    tbl = _load_saturation_table()
    return _interp(T, tbl["T_K_a4"], tbl["cp_l"])


def mu_liquid_sat(T):
    """
    Dynamic viscosity of saturated liquid N2O, Pa.s.

    Interpolated from Table A.4 (NIST WebBook, Laesecke & Hafer 1998).
    Replaces the constant MU_LIQUID_N2O where T is known.

    At design conditions (220-300 K): 60-180 uPa.s.
    Decreases strongly with temperature (liquid viscosity typical behaviour).

    Range note: Table A.4 stops at 307.33 K (T_MAX_A4), not the module's
    overall T_MAX = 309.52 K -- see _check_range_a4.

    Parameters
    ----------
    T : float
        Temperature, K.

    Returns
    -------
    float
        Dynamic viscosity of saturated liquid, Pa.s.
    """
    _check_range_a4(T)
    tbl = _load_saturation_table()
    return _interp(T, tbl["T_K_a4"], tbl["mu_l"])


def s_liquid_sat(T):
    """
    Saturated liquid molar entropy of N2O at temperature T, from Table A.4
    (NIST WebBook, Lemmon & Span 2006) via linear interpolation.

    Added September 2026 for the isentropic choking limit (Priority 1,
    future_work.md): the true two-phase speed-of-sound condition that
    defines choking is a constant-ENTROPY (isentropic) derivative,
        c^2 = (dP/drho)_s,
    because an acoustic disturbance is a small, fast, essentially
    reversible perturbation on top of the (possibly irreversible) mean
    flow -- unlike the real thermodynamic state at an orifice exit, which
    is correctly described by the isenthalpic path (see
    vapor_quality_isenthalpic in injector_two_phase.py) because the
    orifice itself is adiabatic but not reversible.

    Units: kJ/(kmol.K), numerically identical to the source table's
    J/(mol.K) -- no conversion factor needed (1 J/mol.K = 1 kJ/kmol.K),
    unlike h_l/h_v, which the source stores in kJ/mol and this module
    converts to kJ/kmol.

    Range note: Table A.4 stops at 307.33 K (T_MAX_A4), narrower than the
    module's overall T_MAX = 309.52 K -- see _check_range_a4.

    Parameters
    ----------
    T : float
        Temperature, K.

    Returns
    -------
    float
        Saturated liquid molar entropy, s_l, in kJ/(kmol.K).
    """
    _check_range_a4(T)
    tbl = _load_saturation_table()
    return _interp(T, tbl["T_K_a4"], tbl["s_l"])


def s_vapor_sat(T):
    """
    Saturated vapour molar entropy of N2O at temperature T, from Table A.4
    (NIST WebBook, Lemmon & Span 2006) via linear interpolation.

    See s_liquid_sat docstring for the physical motivation, unit note,
    and range note -- identical here.

    Parameters
    ----------
    T : float
        Temperature, K.

    Returns
    -------
    float
        Saturated vapour molar entropy, s_v, in kJ/(kmol.K).
    """
    _check_range_a4(T)
    tbl = _load_saturation_table()
    return _interp(T, tbl["T_K_a4"], tbl["s_v"])


def s_fg(T):
    """
    Entropy of vaporization of N2O at temperature T,
        s_fg(T) = s_v(T) - s_l(T)
    the molar entropy increase on fully vaporizing saturated liquid at
    temperature T -- the entropy-domain analogue of h_fg(T), used as the
    denominator of the isentropic vapour-quality formula (Priority 1):

        x_is(P2) = (s_up - s_l(T_sat(P2))) / s_fg(T_sat(P2))

    Parameters
    ----------
    T : float
        Temperature, K.

    Returns
    -------
    float
        Entropy of vaporization, kJ/(kmol.K). Positive throughout the
        valid range (vapour has higher entropy than liquid at the same
        T, away from the critical point where the two converge -- same
        trend as h_fg, see docs/01_n2o_thermodynamics.md Section 1.5).
    """
    return s_vapor_sat(T) - s_liquid_sat(T)


def mu_mixture(x, T=None, mu_l=None):
    """
    Dynamic viscosity of a two-phase liquid-vapour N2O mixture, Pa.s.

    Uses the McAdams mixing rule (linear in mass quality):
        mu_mix = (1 - x) * mu_l + x * mu_v

    Parameters
    ----------
    x : float
        Vapour quality (mass fraction of vapour), [0, 1].
    T : float or None
        Temperature, K. If provided, mu_l = mu_liquid_sat(T) and
        mu_v = mu_vapor_sat(T). If None, uses constant MU_LIQUID_N2O
        and a mid-range mu_v estimate of 13e-6 Pa.s (~250 K).
    mu_l : float or None
        Override for liquid viscosity, Pa.s. If None, computed from T
        via mu_liquid_sat(T) or defaulted to MU_LIQUID_N2O.

    Returns
    -------
    float
        Mixture dynamic viscosity, Pa.s.
    """
    if T is not None:
        mu_l_val = mu_liquid_sat(T) if mu_l is None else mu_l
        mu_v_val = mu_vapor_sat(T)
    else:
        mu_l_val = MU_LIQUID_N2O if mu_l is None else mu_l
        mu_v_val = 13e-6
    return (1.0 - x) * mu_l_val + x * mu_v_val

def mu_vapor_sat(T):
    """
    Dynamic viscosity of saturated N2O vapour at temperature T, Pa.s.

    Interpolated linearly from Table A.3 of n2o_saturation_table.csv,
    which contains NIST WebBook data (Lemmon & Span 2006 equation of state;
    viscosity from Millat, Vesovic & Wakeham 1991).

    NIST uncertainty: ~2% at T > 150 K (dilute gas limit); higher near
    the critical point (T > 295 K).

    Valid range: 182.33 K (triple point) to 307.33 K (near critical).

    Physical note: vapour viscosity INCREASES with T (unlike liquids).
    At design conditions (220-300 K): 11-18 uPa.s. For comparison,
    MU_LIQUID_N2O ~ 60-130 uPa.s -- roughly 10x larger. In the two-phase
    feed-line model, mu_mix is dominated by the liquid fraction except
    at very high vapour quality.

    Parameters
    ----------
    T : float
        Temperature, K.

    Returns
    -------
    float
        Dynamic viscosity of saturated vapour, Pa.s.
    """
    _check_range(T)
    tbl = _load_saturation_table()
    return _interp(T, tbl["T_K_muv"], tbl["mu_v"])


def degree_of_subcooling(T, P):
    """
    Degree of subcooling, delta_T_sub = T_sat(P) - T, per Section 1.4.

    Positive: liquid is subcooled (has margin before saturation).
    Zero: liquid is exactly saturated.
    Negative: not physically meaningful for a pure liquid state — signals
    that, at this pressure, the fluid at temperature T would already be
    at or past saturation (i.e. two-phase or superheated vapor).

    Parameters
    ----------
    T : float
        Fluid temperature in Kelvin.
    P : float
        Fluid pressure in Pa.

    Returns
    -------
    float
        Degree of subcooling in Kelvin.
    """
    return T_sat(P) - T


if __name__ == "__main__":
    # --- Self-validation against the reference values quoted in
    #     docs/01_n2o_thermodynamics.md, Section 1.2 ---
    print("Validating P_sat(T) against reference values from docs/01_n2o_thermodynamics.md")
    print("-" * 75)

    references = [
        (273.15, 31.3e5, "0 degC"),
        (293.15, 50.9e5, "20 degC"),
        (307.0, 72.5e5, "34 degC (near critical, safely inside valid range)"),
    ]

    for T_K, P_ref, label in references:
        P_calc = P_sat(T_K)
        error_pct = 100.0 * (P_calc - P_ref) / P_ref
        print(f"  T = {label:28s}: P_sat = {P_calc/1e5:6.2f} bar   "
              f"(reference: {P_ref/1e5:.1f} bar, error: {error_pct:+.2f}%)")

    print("-" * 75)
    print(f"Critical point (correlation limit): T = {T_CRIT} K = {T_CRIT - 273.15:.1f} degC, "
          f"P_sat = {P_sat(T_MAX)/1e5:.2f} bar")

    # Sanity check: T_sat should invert P_sat correctly
    print("-" * 75)
    print("Validating T_sat(P) as the inverse of P_sat(T):")
    for T_K, _, label in references:
        P_check = P_sat(T_K)
        T_recovered = T_sat(P_check, T_guess=250.0)
        print(f"  T = {label:28s}: recovered T = {T_recovered:.4f} K "
              f"(original: {T_K:.4f} K, diff: {T_recovered - T_K:+.2e} K)")

    # Sanity check: liquid density at 20 degC (typical reference: ~786 kg/m^3)
    print("-" * 75)
    rho_20C = rho_liquid_sat(293.15)
    print(f"Saturated liquid density at 20 degC: {rho_20C:.1f} kg/m^3 "
          f"(typical reference value: ~786 kg/m^3)")

    # --- Self-validation of the tabulated properties (nu_vapor_sat,
    #     h_liquid_sat, h_vapor_sat, h_fg), added for the HEM/Dyer model ---
    print("-" * 75)
    print("Validating tabulated properties (interpolated from Table A.1):")
    print("-" * 75)

    # Check 1: evaluating exactly AT a tabulated row should reproduce that
    # row's value with no interpolation error (T=290K is an exact row in
    # n2o_saturation_table.csv).
    T_exact = 290.0
    print(f"Exact-row check at T = {T_exact} K (table row, no interpolation):")
    print(f"  nu_vapor_sat = {nu_vapor_sat(T_exact):.5f} m^3/kmol "
          f"(table: 0.30912)")
    print(f"  h_liquid_sat = {h_liquid_sat(T_exact):.1f} kJ/kmol "
          f"(table: {9.0600*1000:.1f})")
    print(f"  h_vapor_sat  = {h_vapor_sat(T_exact):.1f} kJ/kmol "
          f"(table: {17.071*1000:.1f})")

    # Check 2: evaluating BETWEEN two tabulated rows (T=292.5K, halfway
    # between the 290K and 295K rows) should give a value halfway between
    # the two table rows, confirming the interpolation logic itself.
    T_mid = 292.5
    nu_v_mid = nu_vapor_sat(T_mid)
    nu_v_expected_mid = (0.30912 + 0.26142) / 2.0
    print(f"\nMidpoint interpolation check at T = {T_mid} K (halfway between "
          f"290K and 295K rows):")
    print(f"  nu_vapor_sat = {nu_v_mid:.5f} m^3/kmol "
          f"(expected, simple average: {nu_v_expected_mid:.5f})")

    # Check 3: h_fg should be positive throughout the valid range (vapor
    # always has higher enthalpy than liquid at the same T, away from the
    # critical point), and should shrink toward zero approaching T_crit,
    # per Section 1.5/3.2 of the theory docs (latent heat vanishes at the
    # critical point, where liquid and vapor become indistinguishable).
    print(f"\nLatent heat h_fg(T) trend approaching the critical point:")
    for T_check in [220.0, 260.0, 290.0, 305.0, 309.5]:
        print(f"  h_fg({T_check:.1f} K) = {h_fg(T_check)/1000:.2f} kJ/mol")

    # --- Self-validation of the new entropy functions (s_liquid_sat,
    #     s_vapor_sat, s_fg), added September 2026 for Priority 1 ---
    print("-" * 75)
    print("Validating entropy functions (interpolated from Table A.4):")
    print("-" * 75)

    # Table A.4 rows are at 182.33 + 5*n K (triple point + 5 K steps),
    # NOT at round numbers -- 252.33 K is an exact row, 250.0 K is not.
    T_exact_a4 = 252.33  # exact row in Table A.4
    print(f"Exact-row check at T = {T_exact_a4} K:")
    print(f"  s_liquid_sat = {s_liquid_sat(T_exact_a4):.4f} kJ/(kmol.K) "
          f"(table s_l: 24.4610)")
    print(f"  s_vapor_sat  = {s_vapor_sat(T_exact_a4):.4f} kJ/(kmol.K) "
          f"(table s_v: 72.9670)")
    print(f"  s_fg         = {s_fg(T_exact_a4):.4f} kJ/(kmol.K)")

    print(f"\ns_fg(T) trend approaching the (Table A.4) high-T limit:")
    for T_check in [200.0, 240.0, 270.0, 300.0, 307.0]:
        print(f"  s_fg({T_check:.1f} K) = {s_fg(T_check):.3f} kJ/(kmol.K)")

    print(f"\nOut-of-range check (Table A.4 stops at {T_MAX_A4} K, "
          f"below the correlation's own {T_MAX} K):")
    try:
        s_liquid_sat(308.0)
        print("  UNEXPECTED: no error raised")
    except ValueError as e:
        print(f"  Correctly raised ValueError: {e}")
