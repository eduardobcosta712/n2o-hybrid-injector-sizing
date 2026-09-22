"""
n2o_properties.py

Saturated thermophysical properties of N2O as a function of temperature.

DATA SOURCES (changed in September 2026, Priority 4 of the roadmap)
-------------------------------------------------------------------
1. THERMODYNAMIC PROPERTIES -- CoolProp (Bell et al., 2014), which
   implements the Lemmon & Span (2006) fundamental equation of state for
   nitrous oxide (fluid "NitrousOxide"). Used for: saturation pressure
   and its inverse, saturated liquid/vapour density, enthalpy, entropy,
   and saturated-liquid heat capacity. These replace the previous Perry /
   McGill correlations and tables (Table A.1), which carried ~2-3 % error
   near the critical point and a 3-5 % error in the latent heat.

2. VISCOSITY -- NIST WebBook tables (Tables A.3 and A.4 of
   n2o_saturation_table.csv), linearly interpolated. CoolProp does NOT
   provide a viscosity model for N2O (its documentation for the fluid
   lists only the equation of state and a surface-tension correlation),
   so these tables remain the viscosity source. They stop at 307.33 K,
   short of the critical point; see T_MAX_A3 / T_MAX_A4 below. The
   literature attribution of the underlying NIST viscosity correlations
   has not been independently verified (docs/references.md, "Open
   bibliographic points"); the numerical values are taken directly from
   the WebBook tables.

The older Perry/McGill data (closed-form P_sat and rho_l correlations,
Table A.1 enthalpies and vapour volumes) and the NIST cp/entropy columns
of Table A.4 are no longer used by the model. They stay in the CSV, and
the test suite (tests/test_n2o_properties.py) uses them as INDEPENDENT
cross-checks of the CoolProp values.

Requirement: `pip install CoolProp`. Importing this module without CoolProp
raises an ImportError with that instruction.

Valid temperature range for the thermodynamic functions:
[T_MIN, T_MAX] = [Ttriple + 0.01 K, Tcrit - 0.02 K] (about 182.34 to
309.50 K). The small margins keep clear of the two end points, where
saturation calls of an equation of state are ill-conditioned. Viscosity
functions are limited to Tables A.3/A.4 (182.33-307.33 K). All functions
raise ValueError outside their range rather than silently extrapolating.

Reference states: CoolProp's zero of enthalpy and entropy is not the one
used in the NIST tables of the CSV. The model only ever uses enthalpy and
entropy DIFFERENCES within one source (x = (h_up - h_l)/h_fg, ...), so the
reference cancels; the tests compare only reference-independent quantities.

Units: T in Kelvin, P in Pa, molar volumes in m^3/kmol, molar enthalpies
in kJ/kmol, molar entropies in kJ/(kmol.K) [numerically identical to
J/(mol.K)], viscosity in Pa.s, cp in J/(kg.K). The molar API of the
previous versions is kept unchanged, so the injector and feed-line modules
needed no change of formulas; conversions from CoolProp's mass-based SI
output use M_N2O.
"""

import csv
import functools
import os

try:
    from CoolProp.CoolProp import PropsSI
except ImportError as exc:  # pragma: no cover - environment problem
    raise ImportError(
        "n2o_properties.py needs the CoolProp package for the N2O equation "
        "of state (Lemmon & Span 2006). Install it with:  pip install CoolProp"
    ) from exc

_FLUID = "NitrousOxide"

# Molar mass of N2O, kg/kmol -- taken from CoolProp so that every molar <->
# mass conversion in the project uses the same value as the equation of
# state. Single definition for the whole project: feed_line.py,
# injector_two_phase.py and full_system.py import it from here.
M_N2O = PropsSI("M", _FLUID) * 1000.0

# Critical point (CoolProp / Lemmon & Span 2006): 309.52 K, 7.245 MPa.
T_CRIT = PropsSI("Tcrit", _FLUID)   # K  (~36.4 degC)
P_CRIT = PropsSI("pcrit", _FLUID)   # Pa (~72.45 bar)
T_TRIPLE = PropsSI("Ttriple", _FLUID)

# Valid temperature range of the thermodynamic functions (see module docstring)
T_MIN = T_TRIPLE + 0.01
T_MAX = T_CRIT - 0.02

# Valid range of the viscosity tables. Table A.4 (liquid viscosity) and
# Table A.3 (vapour viscosity) were both downloaded from the NIST WebBook at
# 5 K spacing starting at the triple point and end at 307.33 K. Kept as
# separate constants because they are separate tables that could be
# extended independently.
T_MIN_A3 = 182.33
T_MAX_A3 = 307.33
T_MIN_A4 = 182.33
T_MAX_A4 = 307.33

# Path to the CSV with the NIST viscosity tables (and the legacy tables used
# as test cross-checks), relative to this file's own location.
_TABLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "n2o_saturation_table.csv")

# Module-level cache for the parsed table, populated on first use.
_table_cache = None


def _check_range(T):
    if not (T_MIN <= T <= T_MAX):
        raise ValueError(
            f"Temperature {T:.2f} K is outside the valid range of the "
            f"saturation properties [{T_MIN:.2f}, {T_MAX:.2f}] K "
            "(triple point to critical point)."
        )


def _check_range_a3(T):
    """
    Range check for functions based on Table A.3 (NIST WebBook saturated
    vapour viscosity), which is narrower than the thermodynamic range: it
    stops at 307.33 K (T_MAX_A3). Raising a specific, named error (instead
    of an opaque failure inside the interpolation) is the "fail loudly"
    convention used throughout this project.
    """
    if not (T_MIN_A3 <= T <= T_MAX_A3):
        raise ValueError(
            f"Temperature {T:.2f} K is outside Table A.3's valid range "
            f"[{T_MIN_A3}, {T_MAX_A3}] K (NIST WebBook, saturated vapour "
            "viscosity), which is narrower than the thermodynamic range "
            f"([{T_MIN:.2f}, {T_MAX:.2f}] K): it stops short of the critical "
            "point. Affects mu_vapor_sat and mu_mixture."
        )


def _check_range_a4(T):
    """Same as _check_range_a3, for Table A.4 (saturated liquid viscosity)."""
    if not (T_MIN_A4 <= T <= T_MAX_A4):
        raise ValueError(
            f"Temperature {T:.2f} K is outside Table A.4's valid range "
            f"[{T_MIN_A4}, {T_MAX_A4}] K (NIST WebBook, saturated liquid "
            "viscosity), which is narrower than the thermodynamic range "
            f"([{T_MIN:.2f}, {T_MAX:.2f}] K): it stops short of the critical "
            "point. Affects mu_liquid_sat."
        )


# ---------------------------------------------------------------------------
# CoolProp access (cached: the solvers call these thousands of times with
# repeated arguments)
# ---------------------------------------------------------------------------

@functools.lru_cache(maxsize=8192)
def _sat_mass(prop, T, Q):
    """CoolProp saturation property `prop` at temperature T and quality Q
    (0 = saturated liquid, 1 = saturated vapour), mass-based SI units."""
    try:
        return PropsSI(prop, "T", T, "Q", Q, _FLUID)
    except Exception as exc:
        raise RuntimeError(
            f"CoolProp failed to evaluate '{prop}' for saturated N2O at "
            f"T = {T!r} K, Q = {Q}: {exc}") from exc


@functools.lru_cache(maxsize=8192)
def _T_from_P(P):
    try:
        return PropsSI("T", "P", P, "Q", 0, _FLUID)
    except Exception as exc:
        raise RuntimeError(
            f"CoolProp failed to evaluate the saturation temperature at "
            f"P = {P!r} Pa: {exc}") from exc


def _load_saturation_table():
    """
    Parse n2o_saturation_table.csv and return its tables as parallel lists,
    ready for interpolation. Result is cached in _table_cache.

    The CSV has four sections:
      [TABLE_A1] -- McGill/Perry: nu_v, h_l, h_v.            LEGACY (test cross-check only)
      [TABLE_A2] -- derivatives of Table A.1 (not parsed).   LEGACY
      [TABLE_A3] -- NIST: saturated vapour viscosity mu_v.   USED by the model
      [TABLE_A4] -- NIST: cp_l, mu_l, s_l, s_v.              mu_l USED; the rest
                    are LEGACY (test cross-check only)

    Returns
    -------
    dict
        "T_K"     : Table A.1 temperatures, K
        "nu_v"    : m^3/kmol       (legacy)
        "h_l"     : kJ/kmol        (legacy)
        "h_v"     : kJ/kmol        (legacy)
        "T_K_muv" : Table A.3 temperatures, K
        "mu_v"    : Pa.s (converted from uPa.s)
        "T_K_a4"  : Table A.4 temperatures, K
        "cp_l"    : J/(kg.K)       (legacy)
        "mu_l"    : Pa.s
        "s_l"     : J/(mol.K) == kJ/(kmol.K)   (legacy)
        "s_v"     : J/(mol.K) == kJ/(kmol.K)   (legacy)
    """
    global _table_cache
    if _table_cache is not None:
        return _table_cache

    T_list, nu_v_list, h_l_list, h_v_list = [], [], [], []
    T_muv_list, mu_v_list = [], []

    with open(_TABLE_PATH, "r") as f:
        lines = f.readlines()

    # --- Parse TABLE_A1 (legacy) ---
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
            # h_l and h_v are in kJ/mol in source; convert to kJ/kmol (x1000)
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
        "cp_l":    cp_l_list,
        "mu_l":    mu_l_list,
        "s_l":     s_l_list,
        "s_v":     s_v_list,
    }
    return _table_cache


def _interp(T, x_values, y_values):
    """
    Linear interpolation of y(T), given the table's x_values (temperatures,
    strictly ascending) and y_values (the property at each x_values[i]):
        y(T) = y_i + (y_(i+1) - y_i) / (x_(i+1) - x_i) * (T - x_i)
    for the bracketing pair x_i <= T <= x_(i+1). Implemented by hand to keep
    the viscosity path free of extra dependencies.

    T must already be range-checked by the caller (_check_range_a3 or
    _check_range_a4).
    """
    # Linear scan (not bisection): the tables have ~26 rows, so performance
    # is a non-issue and a plain scan is easier to verify by hand.
    for i in range(len(x_values) - 1):
        x_i, x_ip1 = x_values[i], x_values[i + 1]
        if x_i <= T <= x_ip1:
            y_i, y_ip1 = y_values[i], y_values[i + 1]
            fraction = (T - x_i) / (x_ip1 - x_i)
            return y_i + fraction * (y_ip1 - y_i)

    # Unreachable if the caller's range check passed -- fail loudly anyway.
    raise RuntimeError(
        f"T = {T:.2f} K fell outside the table's rows despite passing the "
        "module's range check -- this indicates a mismatch between the "
        "range constants and the table's actual first/last row."
    )


# ---------------------------------------------------------------------------
# Thermodynamic saturation properties (CoolProp)
# ---------------------------------------------------------------------------

def P_sat(T):
    """
    Saturation (vapour) pressure of N2O at temperature T (Lemmon & Span 2006
    equation of state, via CoolProp).

    Parameters
    ----------
    T : float
        Temperature in Kelvin, within [T_MIN, T_MAX].

    Returns
    -------
    float
        Saturation pressure in Pa.
    """
    _check_range(T)
    return _sat_mass("P", T, 0)


def dP_sat_dT(T):
    """
    Slope of the saturation curve, dP_sat/dT, from the Clausius-Clapeyron
    relation evaluated with equation-of-state quantities,

        dP_sat/dT = h_fg / (T * (1/rho_v - 1/rho_l))      [mass basis],

    which is an exact thermodynamic identity for the saturation curve (no
    finite differencing).

    Returns
    -------
    float
        dP/dT in Pa/K.
    """
    _check_range(T)
    h_fg_mass = _sat_mass("H", T, 1) - _sat_mass("H", T, 0)      # J/kg
    dv = 1.0 / _sat_mass("D", T, 1) - 1.0 / _sat_mass("D", T, 0)  # m^3/kg
    return h_fg_mass / (T * dv)


def T_sat(P):
    """
    Saturation temperature corresponding to a given pressure P -- the
    inverse of P_sat(T), obtained directly from the equation of state (no
    iteration in this module). This is the function used to compute the
    subcooling margin delta_T_sub = T_sat(P) - T at any point along the feed
    system, per Section 1.4 of the theory documents.

    (Before the CoolProp integration this was a damped Newton iteration on
    the Perry correlation; the extra parameters of that solver were removed.)

    Parameters
    ----------
    P : float
        Pressure in Pa, within [P_sat(T_MIN), P_sat(T_MAX)] (about 0.88 to
        72.4 bar).

    Returns
    -------
    float
        Saturation temperature in Kelvin.
    """
    P_low, P_high = P_sat(T_MIN), P_sat(T_MAX)
    if not (P_low <= P <= P_high):
        raise ValueError(
            f"P = {P/1e5:.3f} bar is outside the range covered by the "
            f"saturation properties: [{P_low/1e5:.3f}, {P_high/1e5:.3f}] bar "
            f"(T in [{T_MIN:.2f}, {T_MAX:.2f}] K)."
        )
    return _T_from_P(P)


def rho_liquid_sat(T):
    """
    Saturated liquid density of N2O at temperature T, kg/m^3.
    """
    _check_range(T)
    return _sat_mass("D", T, 0)


def nu_vapor_sat(T):
    """
    Saturated vapour molar volume of N2O at temperature T, m^3/kmol
    (= M_N2O / rho_v; the molar form is kept for API compatibility with the
    HEM/Dyer modules, which compute rho_v = M_N2O / nu_vapor_sat(T)).
    """
    _check_range(T)
    return M_N2O / _sat_mass("D", T, 1)


def h_liquid_sat(T):
    """
    Saturated liquid molar enthalpy, kJ/kmol. Only differences of enthalpy
    are physically meaningful here (see the module note on reference states).
    """
    _check_range(T)
    return _sat_mass("H", T, 0) * M_N2O / 1000.0


def h_vapor_sat(T):
    """
    Saturated vapour molar enthalpy, kJ/kmol (same reference as h_liquid_sat).
    """
    _check_range(T)
    return _sat_mass("H", T, 1) * M_N2O / 1000.0


def h_fg(T):
    """
    Latent heat of vaporization of N2O at temperature T,
        h_fg(T) = h_v(T) - h_l(T),  kJ/kmol,
    the energy required to vaporize a unit (molar) amount of liquid at that
    temperature (Sections 1.6 and 3.2 of the theory docs). It vanishes at the
    critical point.
    """
    return h_vapor_sat(T) - h_liquid_sat(T)


def cp_liquid_sat(T):
    """
    Isobaric heat capacity of saturated liquid N2O, J/(kg.K).

    Current use: none inside the model. (henry_fauske_critical_flow() turned
    out to need the liquid-entropy derivative, not cp_l.) Kept as a tested
    property accessor for future work (e.g. the tank thermal model,
    Priority 7). cp_l diverges at the critical point, so values very close to
    T_MAX should be used with caution.
    """
    _check_range(T)
    return _sat_mass("C", T, 0)


def s_liquid_sat(T):
    """
    Saturated liquid molar entropy, kJ/(kmol.K) (numerically equal to
    J/(mol.K)). Needed by the isentropic choking limit and by the
    Henry-Fauske model: the two-phase speed of sound is a constant-ENTROPY
    derivative, c^2 = (dP/drho)_s, whereas the real state at an orifice exit
    follows the isenthalpic path (see injector_two_phase.py).

    With the CoolProp backend the range is the full thermodynamic range
    [T_MIN, T_MAX]; the 307.33 K limit that applied when entropy came from
    NIST Table A.4 no longer exists.
    """
    _check_range(T)
    return _sat_mass("S", T, 0) * M_N2O / 1000.0


def s_vapor_sat(T):
    """Saturated vapour molar entropy, kJ/(kmol.K); see s_liquid_sat."""
    _check_range(T)
    return _sat_mass("S", T, 1) * M_N2O / 1000.0


def s_fg(T):
    """
    Entropy of vaporization, s_fg(T) = s_v(T) - s_l(T), kJ/(kmol.K): the
    entropy-domain analogue of h_fg(T), used as the denominator of the
    isentropic vapour-quality formula

        x_is(P2) = (s_up - s_l(T_sat(P2))) / s_fg(T_sat(P2)).

    Positive throughout the valid range, and -> 0 at the critical point.
    """
    return s_vapor_sat(T) - s_liquid_sat(T)


# ---------------------------------------------------------------------------
# Dynamic viscosity of saturated N2O (liquid and vapour) -- NIST tables
# ---------------------------------------------------------------------------

# MU_LIQUID_N2O: fallback constant, used only where T lies outside Table A.4's
# range (T > 307.33 K). Wherever T is inside that range, mu_liquid_sat(T) is
# used instead -- feed_line.py does so by default.
MU_LIQUID_N2O = 1.5e-4   # Pa.s, mid-range estimate (conservative)


def mu_liquid_sat(T):
    """
    Dynamic viscosity of saturated liquid N2O, Pa.s, interpolated from
    Table A.4 (NIST WebBook data). At design conditions (220-300 K):
    60-180 uPa.s; decreases strongly with temperature.

    Range: Table A.4 stops at 307.33 K (T_MAX_A4), narrower than the
    thermodynamic range -- see _check_range_a4.
    """
    _check_range_a4(T)
    tbl = _load_saturation_table()
    return _interp(T, tbl["T_K_a4"], tbl["mu_l"])


def mu_vapor_sat(T):
    """
    Dynamic viscosity of saturated N2O vapour at temperature T, Pa.s,
    interpolated linearly from Table A.3 (NIST WebBook data). NIST
    uncertainty: ~2 % at T > 150 K; higher near the critical point
    (T > 295 K).

    Range: 182.33 K to 307.33 K, enforced by _check_range_a3.

    Physical note: vapour viscosity INCREASES with T (unlike liquids). At
    design conditions (220-300 K): 11-18 uPa.s, roughly 10x below the liquid.
    In the two-phase feed-line model, mu_mix is dominated by the liquid
    fraction except at very high vapour quality.
    """
    _check_range_a3(T)
    tbl = _load_saturation_table()
    return _interp(T, tbl["T_K_muv"], tbl["mu_v"])


def mu_mixture(x, T=None, mu_l=None):
    """
    Dynamic viscosity of a two-phase liquid-vapour N2O mixture, Pa.s, using
    the McAdams mixing rule (linear in mass quality):
        mu_mix = (1 - x) * mu_l + x * mu_v

    Parameters
    ----------
    x : float
        Vapour quality (mass fraction of vapour), [0, 1].
    T : float or None
        Temperature, K. If provided, mu_l = mu_liquid_sat(T) (unless
        overridden) and mu_v = mu_vapor_sat(T); T must then lie within
        Tables A.3/A.4's range (<= 307.33 K). If None, uses the constant
        MU_LIQUID_N2O and a mid-range mu_v estimate of 13e-6 Pa.s (~250 K).
    mu_l : float or None
        Override for liquid viscosity, Pa.s.

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


def degree_of_subcooling(T, P):
    """
    Degree of subcooling, delta_T_sub = T_sat(P) - T, per Section 1.4.

    Positive: liquid is subcooled (has margin before saturation).
    Zero: liquid is exactly saturated.
    Negative: not physically meaningful for a pure liquid state -- signals
    that, at this pressure, the fluid at temperature T would already be at or
    past saturation (i.e. two-phase or superheated vapour).
    """
    return T_sat(P) - T


if __name__ == "__main__":
    # --- Self-check of the CoolProp backend ---
    print("N2O properties -- CoolProp (Lemmon & Span 2006) backend")
    print("-" * 70)
    print(f"Molar mass {M_N2O:.4f} kg/kmol | T_crit {T_CRIT:.2f} K | "
          f"P_crit {P_CRIT/1e5:.2f} bar | valid T: [{T_MIN:.2f}, {T_MAX:.2f}] K")
    print("-" * 70)

    # Literature reference values quoted in docs/01_n2o_thermodynamics.md
    for T_K, P_ref, label in [(273.15, 31.3e5, "0 degC"),
                              (293.15, 50.9e5, "20 degC"),
                              (T_MAX, P_CRIT, "T_max (0.02 K below T_crit)")]:
        P_calc = P_sat(T_K)
        print(f"  P_sat({label:28s}) = {P_calc/1e5:7.3f} bar "
              f"(reference {P_ref/1e5:.2f} bar, {100*(P_calc-P_ref)/P_ref:+.2f}%)")

    print("-" * 70)
    print("T_sat as the inverse of P_sat:")
    for T_K in (200.0, 273.15, 293.15, 307.0):
        T_back = T_sat(P_sat(T_K))
        print(f"  T = {T_K:7.2f} K -> recovered {T_back:.6f} K "
              f"(diff {T_back - T_K:+.1e} K)")

    print("-" * 70)
    T = 293.15
    print(f"At {T-273.15:.1f} degC: rho_l = {rho_liquid_sat(T):.1f} kg/m^3, "
          f"rho_v = {M_N2O/nu_vapor_sat(T):.1f} kg/m^3, "
          f"h_fg = {h_fg(T)/M_N2O:.1f} kJ/kg, s_fg = {s_fg(T)/M_N2O*1000:.1f} J/(kg.K)")
    print(f"h_fg trend approaching the critical point:")
    for T_check in (220.0, 260.0, 290.0, 305.0, T_MAX):
        print(f"  h_fg({T_check:.2f} K) = {h_fg(T_check)/1000:.2f} kJ/mol")

    print("-" * 70)
    print("Cross-check against the legacy tables kept in the CSV (informational):")
    tbl = _load_saturation_table()
    for T_check in (250.0, 280.0, 290.0, 300.0):
        nu_legacy = _interp(T_check, tbl["T_K"], tbl["nu_v"])
        hfg_legacy = (_interp(T_check, tbl["T_K"], tbl["h_v"])
                      - _interp(T_check, tbl["T_K"], tbl["h_l"]))
        print(f"  {T_check:.0f} K: nu_v {nu_vapor_sat(T_check):.5f} vs {nu_legacy:.5f} m^3/kmol | "
              f"h_fg {h_fg(T_check):.0f} vs {hfg_legacy:.0f} kJ/kmol")

    print("-" * 70)
    print(f"Out-of-range checks: viscosity tables stop at {T_MAX_A3} K, "
          f"thermodynamics at {T_MAX:.2f} K")
    for fname, f in (("s_liquid_sat", s_liquid_sat), ("mu_vapor_sat", mu_vapor_sat)):
        try:
            f(T_CRIT + 1.0)
            print(f"  UNEXPECTED: {fname} raised no error at T_crit + 1 K")
        except ValueError:
            print(f"  {fname}(T_crit + 1 K): correctly raised ValueError")
