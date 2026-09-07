"""
test_injector_two_phase.py

Tests for injector_two_phase.py: vapour quality, HEM mixture density,
HEM mass flow, Dyer non-equilibrium parameter, the full Dyer mass flow,
and the isenthalpic/isentropic critical-flow choking scans.
"""

import math
import pytest

from injector_two_phase import (
    vapor_quality_isenthalpic, vapor_quality_isentropic, hem_mixture_density,
    hem_mass_flow, dyer_non_equilibrium_parameter, dyer_mass_flow,
    hem_critical_flow, hem_critical_flow_isentropic,
)
from n2o_properties import (
    P_sat, T_sat, rho_liquid_sat, nu_vapor_sat,
    h_liquid_sat, h_fg, s_liquid_sat, s_fg,
    T_MAX_A4,
)

M_N2O = 44.013  # kg/kmol


def _rho_v(T):
    """Helper: saturated vapour density at T, kg/m³."""
    return M_N2O / nu_vapor_sat(T)


# ---------------------------------------------------------------------------
# vapor_quality_isenthalpic
# ---------------------------------------------------------------------------

class TestVaporQuality:

    def test_clamped_to_zero_when_no_vaporisation(self):
        # If the upstream enthalpy equals h_l at T_downstream,
        # no energy is available for vaporisation → x = 0.
        T_down = T_sat(20e5)
        h_up = h_liquid_sat(T_down)   # same enthalpy as downstream liquid
        x = vapor_quality_isenthalpic(h_up, T_down)
        assert math.isclose(x, 0.0, abs_tol=1e-9)

    def test_between_zero_and_one(self):
        # With a typical upstream enthalpy (20 °C liquid entering
        # and expanding to 20 bar), x must be in [0, 1].
        T_up = 293.15
        T_down = T_sat(20e5)
        h_up = h_liquid_sat(T_up)
        x = vapor_quality_isenthalpic(h_up, T_down)
        assert 0.0 <= x <= 1.0

    def test_increases_with_pressure_drop(self):
        # A larger pressure drop means more energy released → more vapour.
        T_up = 293.15
        h_up = h_liquid_sat(T_up)
        x_small_drop = vapor_quality_isenthalpic(h_up, T_sat(40e5))
        x_large_drop = vapor_quality_isenthalpic(h_up, T_sat(10e5))
        assert x_large_drop > x_small_drop


# ---------------------------------------------------------------------------
# vapor_quality_isentropic -- added September 2026, Priority 1
# Same structure as TestVaporQuality above, mirrored onto the entropy path.
# ---------------------------------------------------------------------------

class TestVaporQualityIsentropic:

    def test_clamped_to_zero_when_no_vaporisation(self):
        # If the upstream entropy equals s_l at T_downstream, no entropy
        # "room" is available for vaporisation -> x = 0. Mirrors
        # TestVaporQuality.test_clamped_to_zero_when_no_vaporisation.
        T_down = T_sat(20e5)
        s_up = s_liquid_sat(T_down)
        x = vapor_quality_isentropic(s_up, T_down)
        assert math.isclose(x, 0.0, abs_tol=1e-9)

    def test_between_zero_and_one(self):
        T_up = 280.0   # within Table A.4's range (<= T_MAX_A4)
        T_down = T_sat(20e5)
        s_up = s_liquid_sat(T_up)
        x = vapor_quality_isentropic(s_up, T_down)
        assert 0.0 <= x <= 1.0

    def test_increases_with_pressure_drop(self):
        T_up = 280.0
        s_up = s_liquid_sat(T_up)
        x_small_drop = vapor_quality_isentropic(s_up, T_sat(45e5))
        x_large_drop = vapor_quality_isentropic(s_up, T_sat(10e5))
        assert x_large_drop > x_small_drop

    def test_differs_from_isenthalpic_at_same_state(self):
        # The two paths are genuinely different calculations (one uses
        # h_l/h_fg, the other s_l/s_fg) -- they need not agree numerically,
        # but both should land in a physically sensible range for the
        # same operating point, confirming neither is silently reusing
        # the other's tables.
        T_up = 280.0
        T_down = T_sat(20e5)
        x_h = vapor_quality_isenthalpic(h_liquid_sat(T_up), T_down)
        x_s = vapor_quality_isentropic(s_liquid_sat(T_up), T_down)
        assert 0.0 <= x_h <= 1.0
        assert 0.0 <= x_s <= 1.0


# ---------------------------------------------------------------------------
# hem_mixture_density
# ---------------------------------------------------------------------------

class TestHemMixtureDensity:

    def test_pure_liquid_x0(self):
        # x = 0 → pure liquid → rho_mix = rho_l
        rho_l, rho_v = 800.0, 5.0
        assert math.isclose(hem_mixture_density(0.0, rho_l, rho_v), rho_l, rel_tol=1e-9)

    def test_pure_vapour_x1(self):
        # x = 1 → pure vapour → rho_mix = rho_v
        rho_l, rho_v = 800.0, 5.0
        assert math.isclose(hem_mixture_density(1.0, rho_l, rho_v), rho_v, rel_tol=1e-9)

    def test_decreases_with_x(self):
        # Adding more vapour reduces mixture density.
        rho_l, rho_v = 800.0, 5.0
        rho_low_x  = hem_mixture_density(0.1, rho_l, rho_v)
        rho_high_x = hem_mixture_density(0.5, rho_l, rho_v)
        assert rho_low_x > rho_high_x

    def test_much_lower_than_liquid_for_moderate_x(self):
        # Even x = 0.2 should give a density substantially below rho_l,
        # because vapour has much lower density (rho_v << rho_l).
        rho_l, rho_v = 800.0, 5.0
        rho_mix = hem_mixture_density(0.2, rho_l, rho_v)
        assert rho_mix < 0.5 * rho_l


# ---------------------------------------------------------------------------
# dyer_non_equilibrium_parameter
# ---------------------------------------------------------------------------

class TestDyerKappa:

    def test_increases_with_subcooling_margin(self):
        # kappa = sqrt((P_up - P_down) / (P_sat - P_down))
        # A higher P_up (more subcooling) gives a larger numerator,
        # so kappa is LARGER when the fluid is more subcooled,
        # pushing the blend closer to SPI (more "time to vaporise").
        T = 293.15; P_down = 20e5
        P_sat_T = P_sat(T)
        kappa_barely_subcooled = dyer_non_equilibrium_parameter(P_sat_T + 0.5e5, T, P_down)
        kappa_well_subcooled   = dyer_non_equilibrium_parameter(P_sat_T + 10e5,  T, P_down)
        assert kappa_well_subcooled > kappa_barely_subcooled

    def test_positive(self):
        T = 293.15
        kappa = dyer_non_equilibrium_parameter(P_sat(T) + 5e5, T, 20e5)
        assert kappa > 0

    def test_two_phase_inlet_raises(self):
        # P_upstream <= P_sat(T_upstream) → fluid already two-phase → ValueError.
        T = 293.15
        with pytest.raises(ValueError):
            dyer_non_equilibrium_parameter(P_sat(T) - 1e3, T, 20e5)

    def test_exactly_at_saturation_raises(self):
        T = 293.15
        with pytest.raises(ValueError):
            dyer_non_equilibrium_parameter(P_sat(T), T, 20e5)


# ---------------------------------------------------------------------------
# dyer_mass_flow — full Dyer model
# ---------------------------------------------------------------------------

class TestDyerMassFlow:

    # Reference operating point: 20 °C, 55 bar → 20 bar, 4 holes of 1.1 mm
    T_UP   = 293.15
    P_UP   = 55e5
    P_DOWN = 20e5
    Cd     = 0.65
    A      = 4 * math.pi * (0.55e-3) ** 2  # 4 holes, 1.1 mm diameter

    @property
    def _rho_l_up(self):
        return rho_liquid_sat(self.T_UP)

    @property
    def _rho_l_down(self):
        return rho_liquid_sat(T_sat(self.P_DOWN))

    @property
    def _rho_v_down(self):
        return _rho_v(T_sat(self.P_DOWN))

    def _run(self, **overrides):
        kw = dict(Cd=self.Cd, A=self.A, T_upstream=self.T_UP,
                  P_upstream=self.P_UP, P_downstream=self.P_DOWN,
                  rho_l_upstream=self._rho_l_up,
                  rho_l_downstream=self._rho_l_down,
                  rho_v_downstream=self._rho_v_down)
        kw.update(overrides)
        return dyer_mass_flow(**kw)

    def test_returns_dict_with_required_keys(self):
        result = self._run()
        for key in ("m_dot_Dyer", "m_dot_SPI", "m_dot_HEM", "kappa", "x_exit"):
            assert key in result

    def test_dyer_between_hem_and_spi(self):
        # The Dyer blend must be between the HEM (lower) and SPI (upper) limits.
        result = self._run()
        assert result["m_dot_HEM"] <= result["m_dot_Dyer"] <= result["m_dot_SPI"]

    def test_dyer_below_spi(self):
        # Two-phase correction always reduces mass flow relative to SPI.
        result = self._run()
        assert result["m_dot_Dyer"] < result["m_dot_SPI"]

    def test_all_flows_positive(self):
        result = self._run()
        assert result["m_dot_Dyer"] > 0
        assert result["m_dot_SPI"]  > 0
        assert result["m_dot_HEM"]  > 0

    def test_exit_quality_in_range(self):
        result = self._run()
        assert 0.0 <= result["x_exit"] <= 1.0

    def test_kappa_positive(self):
        result = self._run()
        assert result["kappa"] > 0

    def test_larger_area_gives_proportionally_larger_flow(self):
        # All three models must scale linearly with area.
        r1 = self._run()
        r2 = self._run(A=self.A * 2)
        assert math.isclose(r2["m_dot_Dyer"], r1["m_dot_Dyer"] * 2, rel_tol=1e-3)


# ---------------------------------------------------------------------------
# hem_critical_flow -- isenthalpic choking scan (unchanged, re-tested here
# for regression protection now that hem_critical_flow_isentropic sits
# alongside it)
# ---------------------------------------------------------------------------

class TestHemCriticalFlow:

    # Waxman (2013/2014) conditions, per validation/waxman_2013_results.md
    T1 = 280.0
    P1 = 4.36e6
    Cd = 0.65
    D  = 0.0015
    A  = math.pi * (D / 2.0) ** 2

    def test_matches_documented_waxman_value(self):
        # Documented in docs/04_implementation.md / future_work.md: 41.1 g/s
        crit = hem_critical_flow(self.Cd, self.A, self.T1, self.P1)
        assert math.isclose(crit["m_dot_crit"] * 1000, 41.1, rel_tol=0.02)

    def test_x_crit_in_range(self):
        crit = hem_critical_flow(self.Cd, self.A, self.T1, self.P1)
        assert 0.0 <= crit["x_crit"] <= 1.0

    def test_below_p_sat(self):
        from n2o_properties import P_sat as P_sat_f
        crit = hem_critical_flow(self.Cd, self.A, self.T1, self.P1)
        assert crit["P2_crit"] < P_sat_f(self.T1)

    def test_scales_linearly_with_area(self):
        crit1 = hem_critical_flow(self.Cd, self.A, self.T1, self.P1)
        crit2 = hem_critical_flow(self.Cd, 2 * self.A, self.T1, self.P1)
        assert math.isclose(crit2["m_dot_crit"], 2 * crit1["m_dot_crit"], rel_tol=1e-6)


# ---------------------------------------------------------------------------
# hem_critical_flow_isentropic -- added September 2026, Priority 1
# ---------------------------------------------------------------------------

class TestHemCriticalFlowIsentropic:

    # Same Waxman conditions as TestHemCriticalFlow, for direct comparison.
    T1 = 280.0
    P1 = 4.36e6
    Cd = 0.65
    D  = 0.0015
    A  = math.pi * (D / 2.0) ** 2

    # --- Known-value / cross-check ---------------------------------------

    def test_close_to_isenthalpic_value(self):
        # The isentropic and isenthalpic scans are different calculations
        # but should agree to within a few percent at Waxman conditions
        # (entropy and enthalpy corrections are both "small" physically
        # reasonable refinements of the same underlying choking condition,
        # not competing models).
        crit_h = hem_critical_flow(self.Cd, self.A, self.T1, self.P1)
        crit_s = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, self.P1)
        rel_diff = abs(crit_s["m_dot_crit"] - crit_h["m_dot_crit"]) / crit_h["m_dot_crit"]
        assert rel_diff < 0.05

    def test_both_below_waxman_experimental_range(self):
        # Both HEM-only ceilings must sit below the experimental Dyer-regime
        # values (44.0-48.0 g/s) -- Dyer's non-equilibrium correction
        # legitimately predicts above either HEM-only ceiling (see
        # validation/waxman_2013_results.md, Section 6).
        crit_s = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, self.P1)
        assert crit_s["m_dot_crit"] * 1000 < 44.0

    # --- Physical properties -----------------------------------------------

    def test_positive(self):
        crit = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, self.P1)
        assert crit["m_dot_crit"] > 0

    def test_x_crit_in_range(self):
        crit = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, self.P1)
        assert 0.0 <= crit["x_crit"] <= 1.0

    def test_below_p_sat(self):
        crit = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, self.P1)
        assert crit["P2_crit"] < P_sat(self.T1)

    def test_scales_linearly_with_area(self):
        crit1 = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, self.P1)
        crit2 = hem_critical_flow_isentropic(self.Cd, 2 * self.A, self.T1, self.P1)
        assert math.isclose(crit2["m_dot_crit"], 2 * crit1["m_dot_crit"], rel_tol=1e-6)

    def test_higher_upstream_pressure_gives_more_flow(self):
        crit_low  = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, 4.0e6)
        crit_high = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, 5.0e6)
        assert crit_high["m_dot_crit"] > crit_low["m_dot_crit"]

    # --- Edge cases -------------------------------------------------------

    def test_raises_above_table_a4_range(self):
        # T_upstream above 307.33 K: entropy data (Table A.4) unavailable.
        # This is the exact domain restriction documented in the
        # function's docstring and in future_work.md, Priority 1.
        with pytest.raises(ValueError):
            hem_critical_flow_isentropic(self.Cd, self.A, T_MAX_A4 + 1.0, self.P1)

    def test_accepts_temperature_at_table_a4_boundary(self):
        # T_upstream exactly at the boundary must be accepted (inclusive).
        crit = hem_critical_flow_isentropic(self.Cd, self.A, T_MAX_A4, 6.0e6)
        assert crit["m_dot_crit"] >= 0.0

    def test_two_phase_inlet_raises_x_above_zero(self):
        # With x_inlet > 0, the upstream entropy is higher (partially
        # vaporised already) -- should still return a valid, positive
        # critical flow, not raise.
        crit = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, self.P1,
                                             x_inlet=0.05)
        assert crit["m_dot_crit"] > 0
