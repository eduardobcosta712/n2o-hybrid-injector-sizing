"""
test_injector_spi.py

Tests for injector_spi.py: SPI mass flow, area inversion,
and the SPI-sufficiency decision criterion.
"""

import math
import pytest

from injector_spi import spi_mass_flow, orifice_area_from_target_flow, spi_sufficient
from n2o_properties import P_sat, rho_liquid_sat


# ---------------------------------------------------------------------------
# spi_mass_flow
# ---------------------------------------------------------------------------

class TestSpiMassFlow:

    def test_known_value(self):
        # m_dot = Cd * A * sqrt(2 * rho * dP)
        # = 0.65 * 3.79e-6 * sqrt(2 * 785 * 30e5) = ?
        Cd, A, rho, dP = 0.65, 3.79e-6, 785.0, 30e5
        expected = Cd * A * math.sqrt(2 * rho * dP)
        result = spi_mass_flow(Cd, A, rho, dP)
        assert math.isclose(result, expected, rel_tol=1e-9)

    def test_scales_with_area(self):
        # Doubling the area should double the mass flow.
        rho, dP = 785.0, 20e5
        m1 = spi_mass_flow(0.65, 1e-6, rho, dP)
        m2 = spi_mass_flow(0.65, 2e-6, rho, dP)
        assert math.isclose(m2, 2 * m1, rel_tol=1e-9)

    def test_scales_with_sqrt_of_delta_p(self):
        # Quadrupling dP doubles mass flow (sqrt relationship).
        rho = 785.0
        m1 = spi_mass_flow(0.65, 3.79e-6, rho, 10e5)
        m2 = spi_mass_flow(0.65, 3.79e-6, rho, 40e5)
        assert math.isclose(m2, 2 * m1, rel_tol=1e-6)

    def test_positive(self):
        assert spi_mass_flow(0.65, 3.79e-6, 785.0, 20e5) > 0

    def test_negative_delta_p_raises(self):
        with pytest.raises(ValueError):
            spi_mass_flow(0.65, 3.79e-6, 785.0, -1e5)

    def test_zero_delta_p_gives_zero_flow(self):
        result = spi_mass_flow(0.65, 3.79e-6, 785.0, 0.0)
        assert result == 0.0


# ---------------------------------------------------------------------------
# orifice_area_from_target_flow — algebraic inverse of spi_mass_flow
# ---------------------------------------------------------------------------

class TestOrificeAreaFromTargetFlow:

    def test_round_trip_with_spi_mass_flow(self):
        # orifice_area_from_target_flow must be the exact inverse of spi_mass_flow:
        # A_required → m_dot → A_back should equal A_required.
        Cd, rho, dP = 0.65, 785.0, 25e5
        m_dot_target = 0.4  # kg/s
        A = orifice_area_from_target_flow(m_dot_target, Cd, rho, dP)
        m_dot_back = spi_mass_flow(Cd, A, rho, dP)
        assert math.isclose(m_dot_back, m_dot_target, rel_tol=1e-9)

    def test_larger_target_requires_larger_area(self):
        Cd, rho, dP = 0.65, 785.0, 20e5
        A1 = orifice_area_from_target_flow(0.3, Cd, rho, dP)
        A2 = orifice_area_from_target_flow(0.6, Cd, rho, dP)
        assert A2 > A1

    def test_positive(self):
        A = orifice_area_from_target_flow(0.4, 0.65, 785.0, 20e5)
        assert A > 0


# ---------------------------------------------------------------------------
# spi_sufficient
# ---------------------------------------------------------------------------

class TestSpiSufficient:

    def test_sufficient_when_downstream_above_psat(self):
        # P_downstream well above P_sat → SPI valid.
        T = 273.15   # 0 °C, P_sat ≈ 32 bar
        P_up   = 60e5
        P_down = 40e5   # still above P_sat(0°C) ≈ 32 bar
        assert spi_sufficient(P_up, T, P_down) is True

    def test_not_sufficient_when_downstream_below_psat(self):
        # P_downstream below P_sat → two-phase effects inside orifice.
        T = 293.15   # 20 °C, P_sat ≈ 51.4 bar
        P_up   = 60e5
        P_down = 20e5   # far below P_sat
        assert spi_sufficient(P_up, T, P_down) is False

    def test_borderline_case_exactly_at_psat(self):
        # P_downstream == P_sat → borderline, considered sufficient
        # (fluid just reaches saturation at the exit, not inside).
        T = 293.15
        Psat = P_sat(T)
        assert spi_sufficient(60e5, T, Psat) is True

    def test_returns_bool(self):
        result = spi_sufficient(60e5, 293.15, 20e5)
        assert isinstance(result, bool)
