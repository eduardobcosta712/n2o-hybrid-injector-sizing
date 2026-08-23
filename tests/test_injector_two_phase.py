"""
test_injector_two_phase.py

Tests for injector_two_phase.py: vapour quality, HEM mixture density,
HEM mass flow, Dyer non-equilibrium parameter, and the full Dyer mass flow.
"""

import math
import pytest

from injector_two_phase import (
    vapor_quality_isenthalpic, hem_mixture_density,
    hem_mass_flow, dyer_non_equilibrium_parameter, dyer_mass_flow,
)
from n2o_properties import (
    P_sat, T_sat, rho_liquid_sat, nu_vapor_sat,
    h_liquid_sat, h_fg,
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
