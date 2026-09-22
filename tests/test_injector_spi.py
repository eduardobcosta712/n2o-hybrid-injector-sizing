"""TEST COPY of the repository's test_injector_spi.py (comments trimmed; unchanged by the audit)."""
import math
import pytest
from injector_spi import spi_mass_flow, orifice_area_from_target_flow, spi_sufficient
from n2o_properties import P_sat, rho_liquid_sat


class TestSpiMassFlow:
    def test_known_value(self):
        Cd, A, rho, dP = 0.65, 3.79e-6, 785.0, 30e5
        assert math.isclose(spi_mass_flow(Cd, A, rho, dP), Cd * A * math.sqrt(2 * rho * dP), rel_tol=1e-9)
    def test_scales_with_area(self):
        assert math.isclose(spi_mass_flow(0.65, 2e-6, 785.0, 20e5), 2 * spi_mass_flow(0.65, 1e-6, 785.0, 20e5), rel_tol=1e-9)
    def test_scales_with_sqrt_of_delta_p(self):
        assert math.isclose(spi_mass_flow(0.65, 3.79e-6, 785.0, 40e5), 2 * spi_mass_flow(0.65, 3.79e-6, 785.0, 10e5), rel_tol=1e-6)
    def test_positive(self):
        assert spi_mass_flow(0.65, 3.79e-6, 785.0, 20e5) > 0
    def test_negative_delta_p_raises(self):
        with pytest.raises(ValueError):
            spi_mass_flow(0.65, 3.79e-6, 785.0, -1e5)
    def test_zero_delta_p_gives_zero_flow(self):
        assert spi_mass_flow(0.65, 3.79e-6, 785.0, 0.0) == 0.0


class TestOrificeAreaFromTargetFlow:
    def test_round_trip_with_spi_mass_flow(self):
        A = orifice_area_from_target_flow(0.4, 0.65, 785.0, 25e5)
        assert math.isclose(spi_mass_flow(0.65, A, 785.0, 25e5), 0.4, rel_tol=1e-9)
    def test_larger_target_requires_larger_area(self):
        assert orifice_area_from_target_flow(0.6, 0.65, 785.0, 20e5) > orifice_area_from_target_flow(0.3, 0.65, 785.0, 20e5)
    def test_positive(self):
        assert orifice_area_from_target_flow(0.4, 0.65, 785.0, 20e5) > 0


class TestSpiSufficient:
    def test_sufficient_when_downstream_above_psat(self):
        assert spi_sufficient(60e5, 273.15, 40e5) is True
    def test_not_sufficient_when_downstream_below_psat(self):
        assert spi_sufficient(60e5, 293.15, 20e5) is False
    def test_borderline_case_exactly_at_psat(self):
        assert spi_sufficient(60e5, 293.15, P_sat(293.15)) is True
    def test_returns_bool(self):
        assert isinstance(spi_sufficient(60e5, 293.15, 20e5), bool)
