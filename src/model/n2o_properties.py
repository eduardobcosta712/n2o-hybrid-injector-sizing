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

Units: T in Kelvin, P in Pa, molar volumes in m^3/kmol, molar enthalpies
in kJ/kmol throughout this module (see CSV header for the raw-file unit
trap: the source table stores enthalpies per MOL, not per kmol).
"""

import csv
import math
import os

# Valid temperature range for all correlations in this module (Kelvin)
T_MIN = 182.33
T_MAX = 309.52

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


def _load_saturation_table():
    """
    Parse n2o_saturation_table.csv and return Table A.1 as parallel lists,
    ready for interpolation. Result is cached in _table_cache after the
    first call, so the file is only read from disk once per program run.

    The CSV file has a non-trivial structure: comment lines starting with
    '#', then a '[TABLE_A1]' section header, then the Table A.1 header row
    and data rows, then a '[TABLE_A2]' section (not parsed here -- Table
    A.2 is not used by the current model, see the CSV header for why).

    Returns
    -------
    dict
        Keys "T_K", "nu_v", "h_l", "h_v", each mapping to a list of floats,
        in the same row order as the source table (ascending T).
    """
    global _table_cache
    if _table_cache is not None:
        return _table_cache

    T_list, nu_v_list, h_l_list, h_v_list = [], [], [], []

    with open(_TABLE_PATH, "r") as f:
        lines = f.readlines()

    # Find where Table A.1's data starts: the line right after the
    # "[TABLE_A1]" marker is the header row (column names); every line
    # after that, until "[TABLE_A2]", is a data row.
    in_table_a1 = False
    header_skipped = False
    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line == "[TABLE_A1]":
            in_table_a1 = True
            header_skipped = False
            continue
        if line == "[TABLE_A2]":
            break  # Table A.1 block is over; A.2 is not needed here.
        if in_table_a1:
            if not header_skipped:
                header_skipped = True  # this line is the column-name row
                continue
            fields = next(csv.reader([line]))
            T_list.append(float(fields[0]))
            nu_v_list.append(float(fields[1]))
            # h_l and h_v are columns 3 and 5 (0-indexed: 3 and 5), in
            # kJ/mol in the source file -- convert to kJ/kmol (x1000) to
            # stay consistent with this module's SI-with-kmol convention.
            h_l_list.append(float(fields[3]) * 1000.0)
            h_v_list.append(float(fields[5]) * 1000.0)

    _table_cache = {
        "T_K": T_list,
        "nu_v": nu_v_list,
        "h_l": h_l_list,
        "h_v": h_v_list,
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
