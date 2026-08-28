"""
test_full_system.py

Integration tests for full_system.py: the full tank → feed line →
injector path under three operating regimes:
  1. Flashing in the feed line (injector not evaluated)
  2. No flashing, SPI sufficient
  3. No flashing, SPI not sufficient (Dyer used)
"""

import math
import pytest

from full_system import evaluate_full_system
from n2o_properties import P_sat


# Shared feed line geometry (re-used across all test cases)
SEGMENTS = [
    {"type": "pipe",    "L": 1.0, "D": 0.008},
    {"type": "fitting", "D": 0.008, "K": 0.05},
    {"type": "pipe",    "L": 1.0, "D": 0.008},
    {"type": "fitting", "D": 0.008, "K": 0.30},
]
T_TANK = 293.15   # K, 20 °C
Cd     = 0.65
A      = 4 * math.pi * (0.55e-3) ** 2  # 4 holes, 1.1 mm


# ---------------------------------------------------------------------------
# Case 1: Flashing in the feed line
# ---------------------------------------------------------------------------

class TestFlashingInLine:
    """
    When flashing is detected in the feed line, the model now uses HEM with
    a two-phase inlet (hem_mass_flow_two_phase_inlet) instead of returning None.
    The fluid arrives partially vaporised; x_inlet > 0 is computed isentalpically.
    """

    def setup_method(self):
        # Tank exactly at saturation → zero initial margin → immediate flashing
        P_tank = P_sat(T_TANK)
        self.result = evaluate_full_system(
            0.5, T_TANK, P_tank, SEGMENTS, Cd, A, 20e5)

    def test_flashing_detected(self):
        assert self.result["feed_line_result"]["flashing_detected"] is True

    def test_x_inlet_positive(self):
        # With flashing, x_inlet must be > 0 (some vapour formed in the line).
        assert self.result["feed_line_result"]["x_inlet"] > 0.0

    def test_x_inlet_below_one(self):
        # x_inlet must be physically bounded [0, 1].
        assert self.result["feed_line_result"]["x_inlet"] < 1.0

    def test_m_dot_real_positive(self):
        # HEM with two-phase inlet must produce a positive mass flow.
        assert self.result["m_dot_real"] is not None
        assert self.result["m_dot_real"] > 0.0

    def test_regime_is_hem_two_phase(self):
        assert self.result["regime"] == "HEM_two_phase_inlet"

    def test_injector_result_has_x_exit(self):
        # The injector result dict must contain x_exit in [0, 1].
        ir = self.result["injector_result"]
        assert ir is not None
        assert 0.0 <= ir["x_exit"] <= 1.0

    def test_spi_sufficient_is_false(self):
        # Two-phase inlet: SPI is not sufficient (never valid here).
        assert self.result["spi_sufficient"] is False

    def test_p_injector_inlet_is_present(self):
        assert self.result["P_injector_inlet"] is not None


# ---------------------------------------------------------------------------
# Case 2: No flashing, SPI sufficient
# ---------------------------------------------------------------------------

class TestSpiSufficient:

    def setup_method(self):
        # 0 °C with plenty of subcooling, small chamber pressure drop
        # P_sat(0°C) ≈ 32 bar; tank at 60 bar → 28 bar subcooling
        # Chamber at 50 bar → P_downstream > P_sat(0°C) → SPI valid
        self.T = 273.15  # 0 °C
        self.P_tank = 60e5
        self.P_chamber = 50e5
        self.result = evaluate_full_system(
            0.3, self.T, self.P_tank, SEGMENTS, Cd, A, self.P_chamber)

    def test_no_flashing(self):
        assert not self.result["feed_line_result"]["flashing_detected"]

    def test_spi_sufficient(self):
        assert self.result["spi_sufficient"] is True

    def test_injector_result_is_none(self):
        # When SPI is sufficient, no two-phase result is computed
        assert self.result["injector_result"] is None

    def test_m_dot_real_positive(self):
        assert self.result["m_dot_real"] > 0

    def test_p_injector_inlet_below_p_tank(self):
        # Pressure must drop along the line
        assert self.result["P_injector_inlet"] < self.P_tank


# ---------------------------------------------------------------------------
# Case 3: No flashing, Dyer model used
# ---------------------------------------------------------------------------

class TestDyerUsed:

    def setup_method(self):
        # 20 °C, 55 bar tank, 20 bar chamber → crosses saturation in orifice
        self.P_tank = 55e5
        self.P_chamber = 20e5
        self.result = evaluate_full_system(
            0.5, T_TANK, self.P_tank, SEGMENTS, Cd, A, self.P_chamber)

    def test_no_flashing(self):
        assert not self.result["feed_line_result"]["flashing_detected"]

    def test_spi_not_sufficient(self):
        assert self.result["spi_sufficient"] is False

    def test_injector_result_present(self):
        assert self.result["injector_result"] is not None

    def test_m_dot_real_below_spi(self):
        # Dyer prediction must be below the (invalid) SPI prediction
        ir = self.result["injector_result"]
        assert self.result["m_dot_real"] < ir["m_dot_SPI"]

    def test_m_dot_real_positive(self):
        assert self.result["m_dot_real"] > 0

    def test_kappa_positive(self):
        assert self.result["injector_result"]["kappa"] > 0

    def test_exit_quality_in_range(self):
        x = self.result["injector_result"]["x_exit"]
        assert 0.0 <= x <= 1.0


# ---------------------------------------------------------------------------
# Physical consistency across cases
# ---------------------------------------------------------------------------

class TestPhysicalConsistency:

    def test_higher_tank_pressure_gives_more_flow(self):
        # More tank pressure → more driving dP across injector → more flow.
        r1 = evaluate_full_system(0.5, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5)
        r2 = evaluate_full_system(0.5, T_TANK, 65e5, SEGMENTS, Cd, A, 20e5)
        if r1["m_dot_real"] and r2["m_dot_real"]:
            assert r2["m_dot_real"] > r1["m_dot_real"]

    def test_result_has_all_required_keys(self):
        r = evaluate_full_system(0.5, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5)
        for key in ("feed_line_result", "P_injector_inlet",
                    "spi_sufficient", "m_dot_real", "injector_result"):
            assert key in r
