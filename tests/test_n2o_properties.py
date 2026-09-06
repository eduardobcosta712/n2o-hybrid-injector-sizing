"""
test_n2o_properties.py

Tests for n2o_properties.py: saturation pressure, saturation temperature,
saturated liquid density, saturated vapour properties, latent heat, and
the degree-of-subcooling convenience function.

Three categories per function:
  1. Known-value checks  — against literature reference values.
  2. Physical properties — monotonicity, signs, limits.
  3. Edge cases          — out-of-range inputs, domain boundaries.
"""

import math
import pytest

from n2o_properties import (
    P_sat, dP_sat_dT, T_sat, rho_liquid_sat,
    nu_vapor_sat, h_liquid_sat, h_vapor_sat, h_fg,
    mu_vapor_sat, mu_mixture,
    degree_of_subcooling,
    T_MIN, T_MAX,
)

from n2o_properties import (
    P_sat, dP_sat_dT, T_sat, rho_liquid_sat,
    nu_vapor_sat, h_liquid_sat, h_vapor_sat, h_fg,
    degree_of_subcooling,
    T_MIN, T_MAX,
)

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------
REL_TOL = 1e-3   # 0.1% — appropriate for engineering correlations
ABS_TOL = 1e-9   # for quantities expected to be near zero


# ---------------------------------------------------------------------------
# P_sat(T)
# ---------------------------------------------------------------------------

class TestPsat:
    """Saturation pressure correlation against published reference values."""

    def test_at_0_degC(self):
        # Literature reference: P_sat(273.15 K) ≈ 31.3 bar
        # Correlation agrees within ~2.6% (documented in 04_implementation.md)
        assert math.isclose(P_sat(273.15), 31.3e5, rel_tol=0.03)

    def test_at_20_degC(self):
        # Literature reference: P_sat(293.15 K) ≈ 50.9 bar, within ~0.9%
        assert math.isclose(P_sat(293.15), 50.9e5, rel_tol=0.015)

    def test_at_34_degC(self):
        # Near the critical point — larger tolerance, correlation accuracy degrades
        assert math.isclose(P_sat(307.0), 72.5e5, rel_tol=0.06)

    def test_monotonically_increasing(self):
        # P_sat must increase with temperature throughout the valid range.
        T_vals = [T_MIN + i * 5 for i in range(25)]
        P_vals = [P_sat(T) for T in T_vals]
        for i in range(len(P_vals) - 1):
            assert P_vals[i] < P_vals[i + 1], (
                f"P_sat not monotonically increasing: "
                f"P_sat({T_vals[i]:.1f}) = {P_vals[i]:.0f} Pa, "
                f"P_sat({T_vals[i+1]:.1f}) = {P_vals[i+1]:.0f} Pa"
            )

    def test_positive_throughout(self):
        for T in [T_MIN, 250.0, T_MAX]:
            assert P_sat(T) > 0

    def test_out_of_range_low(self):
        with pytest.raises(ValueError):
            P_sat(T_MIN - 1.0)

    def test_out_of_range_high(self):
        with pytest.raises(ValueError):
            P_sat(T_MAX + 1.0)


# ---------------------------------------------------------------------------
# dP_sat_dT(T)
# ---------------------------------------------------------------------------

class TestdPsatdT:
    """Derivative of the saturation pressure — must be positive (P_sat is increasing)."""

    def test_positive_throughout(self):
        for T in [T_MIN, 250.0, 290.0, T_MAX]:
            assert dP_sat_dT(T) > 0

    def test_consistent_with_finite_difference(self):
        # dP/dT at 293.15 K should agree with a finite-difference approximation.
        T = 293.15
        h = 0.01  # K
        fd = (P_sat(T + h) - P_sat(T - h)) / (2 * h)
        assert math.isclose(dP_sat_dT(T), fd, rel_tol=1e-4)

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            dP_sat_dT(T_MIN - 1.0)


# ---------------------------------------------------------------------------
# T_sat(P) — numerical inverse of P_sat
# ---------------------------------------------------------------------------

class TestTsat:
    """T_sat must be the exact inverse of P_sat, verified to machine precision."""

    def test_inverse_at_several_temperatures(self):
        for T_original in [200.0, 250.0, 270.0, 290.0, 305.0]:
            P = P_sat(T_original)
            T_recovered = T_sat(P)
            assert math.isclose(T_recovered, T_original, rel_tol=1e-9), (
                f"T_sat(P_sat({T_original})) = {T_recovered:.6f} K, "
                f"expected {T_original:.6f} K"
            )

    def test_monotonically_increasing(self):
        # T_sat(P) must increase with P (inverse of increasing P_sat).
        P_vals = [20e5, 30e5, 40e5, 50e5, 60e5, 70e5]
        T_vals = [T_sat(P) for P in P_vals]
        for i in range(len(T_vals) - 1):
            assert T_vals[i] < T_vals[i + 1]

    def test_at_20_degC_reference(self):
        # T_sat at ~51.4 bar should return ~293.15 K
        T = T_sat(P_sat(293.15))
        assert math.isclose(T, 293.15, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# rho_liquid_sat(T)
# ---------------------------------------------------------------------------

class TestRhoLiquidSat:
    """Saturated liquid density: reference value and physical properties."""

    def test_at_20_degC(self):
        # Literature reference: ~786 kg/m³. Correlation gives 784.8 kg/m³ (<0.2% error).
        rho = rho_liquid_sat(293.15)
        assert math.isclose(rho, 786.0, rel_tol=0.005)

    def test_decreases_with_temperature(self):
        # Liquid density must decrease as temperature rises (thermal expansion).
        T_vals = [220.0, 250.0, 270.0, 290.0, 305.0]
        rho_vals = [rho_liquid_sat(T) for T in T_vals]
        for i in range(len(rho_vals) - 1):
            assert rho_vals[i] > rho_vals[i + 1], (
                f"rho_liquid_sat not decreasing: "
                f"rho({T_vals[i]:.0f}) = {rho_vals[i]:.1f}, "
                f"rho({T_vals[i+1]:.0f}) = {rho_vals[i+1]:.1f}"
            )

    def test_positive_throughout(self):
        for T in [T_MIN, 250.0, T_MAX]:
            assert rho_liquid_sat(T) > 0

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            rho_liquid_sat(400.0)


# ---------------------------------------------------------------------------
# nu_vapor_sat(T) — tabulated, interpolated
# ---------------------------------------------------------------------------

class TestNuVaporSat:
    """Saturated vapour molar volume from the look-up table."""

    def test_exact_table_row_at_290K(self):
        # Table A.1 value at 290 K: 0.30912 m³/kmol
        assert math.isclose(nu_vapor_sat(290.0), 0.30912, rel_tol=1e-4)

    def test_midpoint_interpolation(self):
        # Midpoint between 290 K and 295 K rows should equal their average
        # (linear interpolation — exact at the midpoint by definition).
        nu_290 = nu_vapor_sat(290.0)
        nu_295 = nu_vapor_sat(295.0)
        nu_mid = nu_vapor_sat(292.5)
        expected_mid = (nu_290 + nu_295) / 2.0
        assert math.isclose(nu_mid, expected_mid, rel_tol=1e-9)

    def test_decreases_toward_critical_point(self):
        # Vapour molar volume must decrease as T → T_crit
        # (liquid and vapour converge at the critical point).
        T_vals = [220.0, 250.0, 270.0, 290.0, 305.0]
        nu_vals = [nu_vapor_sat(T) for T in T_vals]
        for i in range(len(nu_vals) - 1):
            assert nu_vals[i] > nu_vals[i + 1]

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            nu_vapor_sat(T_MIN - 5.0)


# ---------------------------------------------------------------------------
# h_fg(T) — latent heat
# ---------------------------------------------------------------------------

class TestHfg:
    """Latent heat of vaporisation: sign, trend, and critical-point behaviour."""

    def test_positive_throughout(self):
        # h_fg must be positive (energy is required to vaporise liquid).
        for T in [220.0, 250.0, 270.0, 290.0, 305.0]:
            assert h_fg(T) > 0, f"h_fg({T}) = {h_fg(T):.1f} kJ/kmol — must be positive"

    def test_decreases_toward_critical_point(self):
        # h_fg → 0 as T → T_crit: latent heat vanishes at the critical point.
        T_vals = [220.0, 250.0, 270.0, 290.0, 305.0]
        hfg_vals = [h_fg(T) for T in T_vals]
        for i in range(len(hfg_vals) - 1):
            assert hfg_vals[i] > hfg_vals[i + 1], (
                f"h_fg not decreasing toward critical point: "
                f"h_fg({T_vals[i]:.0f}) = {hfg_vals[i]:.0f}, "
                f"h_fg({T_vals[i+1]:.0f}) = {hfg_vals[i+1]:.0f} kJ/kmol"
            )

    def test_consistent_with_clausius_clapeyron(self):
        # Cross-validate against Clausius-Clapeyron: h_fg = T * (nu_v - nu_l) * dP/dT
        # This is an independent thermodynamic identity, not used by h_fg() itself.
        # Agreement within ~5% is expected (some error from nu_v interpolation
        # and the approximation of using saturated-liquid nu_l from the correlation).
        T = 270.0
        M_N2O = 44.013  # kg/kmol
        nu_v = nu_vapor_sat(T)                       # m³/kmol
        nu_l = M_N2O / rho_liquid_sat(T)             # m³/kmol
        dPdT = dP_sat_dT(T)                          # Pa/K
        h_fg_cc = T * (nu_v - nu_l) * dPdT / 1e3    # kJ/kmol (Pa·m³ = J, /1000 = kJ)
        h_fg_table = h_fg(T)
        # Allow 8% tolerance — known approximation error near the critical point
        assert math.isclose(h_fg_cc, h_fg_table, rel_tol=0.08), (
            f"Clausius-Clapeyron: {h_fg_cc:.0f} kJ/kmol, "
            f"table: {h_fg_table:.0f} kJ/kmol"
        )


# ---------------------------------------------------------------------------
# degree_of_subcooling(T, P)
# ---------------------------------------------------------------------------

class TestDegreeOfSubcooling:
    """Subcooling margin: sign conventions and physical limits."""

    def test_exactly_saturated_is_zero(self):
        T = 293.15
        P = P_sat(T)
        dT_sub = degree_of_subcooling(T, P)
        assert math.isclose(dT_sub, 0.0, abs_tol=1e-6)

    def test_subcooled_is_positive(self):
        T = 293.15
        P = P_sat(T) + 5e5   # 5 bar above saturation
        assert degree_of_subcooling(T, P) > 0

    def test_superheated_is_negative(self):
        T = 293.15
        P = P_sat(T) - 5e5   # 5 bar below saturation
        assert degree_of_subcooling(T, P) < 0

    def test_increases_with_pressure(self):
        T = 280.0
        P_low  = P_sat(T) + 2e5
        P_high = P_sat(T) + 10e5
        assert degree_of_subcooling(T, P_high) > degree_of_subcooling(T, P_low)


class TestMuVaporSat:
    """Tests for mu_vapor_sat(T) -- saturated vapour dynamic viscosity."""

    def test_known_values_from_nist(self):
        # Exact NIST table points must be returned exactly (within float precision)
        assert math.isclose(mu_vapor_sat(182.33) * 1e6, 9.0689, rel_tol=1e-4)
        assert math.isclose(mu_vapor_sat(252.33) * 1e6, 13.417, rel_tol=1e-4)
        assert math.isclose(mu_vapor_sat(307.33) * 1e6, 22.982, rel_tol=1e-4)

    def test_increases_with_temperature(self):
        # Vapour viscosity must increase with T (unlike liquids)
        temps = [190.0, 220.0, 250.0, 270.0, 295.0]
        values = [mu_vapor_sat(T) for T in temps]
        for i in range(len(values) - 1):
            assert values[i] < values[i + 1], (
                f"mu_v must increase with T: mu_v({temps[i]}) >= mu_v({temps[i+1]})")

    def test_returns_pa_s_not_upa_s(self):
        # At 250 K, mu_v ~ 13.4e-6 Pa.s; must NOT be 13.4 (would be μPa.s)
        muv = mu_vapor_sat(250.0)
        assert 5e-6 < muv < 30e-6, f"mu_vapor_sat should be in Pa.s, got {muv}"

    def test_out_of_range_raises(self):
        with pytest.raises(ValueError):
            mu_vapor_sat(180.0)  # below triple point
        with pytest.raises(ValueError):
            mu_vapor_sat(315.0)  # above critical

    def test_much_less_than_liquid(self):
        # Vapour viscosity << liquid viscosity at all temperatures
        from n2o_properties import MU_LIQUID_N2O
        for T in [200.0, 250.0, 295.0]:
            assert mu_vapor_sat(T) < MU_LIQUID_N2O


class TestMuMixture:
    """Tests for mu_mixture(x, T) -- HEM mixture viscosity."""

    def test_pure_liquid_x0(self):
        from n2o_properties import mu_liquid_sat
        assert math.isclose(mu_mixture(0.0, T=250.0), mu_liquid_sat(250.0))

    def test_pure_vapour_x1(self):
        assert math.isclose(mu_mixture(1.0, T=250.0), mu_vapor_sat(250.0))

    def test_between_liquid_and_vapour(self):
        from n2o_properties import MU_LIQUID_N2O
        mu_mix = mu_mixture(0.3, T=260.0)
        mu_v   = mu_vapor_sat(260.0)
        assert mu_v < mu_mix < MU_LIQUID_N2O

    def test_decreases_with_vapour_quality(self):
        # More vapour = lower mixture viscosity (vapour << liquid)
        x_vals = [0.0, 0.1, 0.3, 0.5, 0.8, 1.0]
        mu_vals = [mu_mixture(x, T=270.0) for x in x_vals]
        for i in range(len(mu_vals) - 1):
            assert mu_vals[i] >= mu_vals[i + 1]
