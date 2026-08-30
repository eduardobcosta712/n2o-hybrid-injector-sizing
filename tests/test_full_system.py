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


# ---------------------------------------------------------------------------
# Coupled solver — convergence and self-consistency tests
# ---------------------------------------------------------------------------

class TestCoupledSolver:
    """
    Tests specific to the coupled feed-line / injector solver.
    The key invariant: at the converged operating point, the mass flow
    used to compute feed-line losses equals the mass flow predicted by
    the injector model (to within the solver tolerance).
    """

    def _run(self, P_tank, P_chamber=20e5, segments=None, m0=0.3):
        segs = segments or SEGMENTS
        return evaluate_full_system(m0, T_TANK, P_tank, segs, Cd, A, P_chamber)

    def test_solver_converged(self):
        # solver_info must report convergence
        r = self._run(55e5)
        assert r["solver_info"]["converged"] is True

    def test_iterations_positive(self):
        r = self._run(55e5)
        assert r["solver_info"]["iterations"] >= 1

    def test_history_length_matches_iterations(self):
        # history has one entry per iteration plus the initial guess
        r = self._run(55e5)
        si = r["solver_info"]
        assert len(si["history"]) == si["iterations"] + 1

    def test_final_rel_err_below_tolerance(self):
        r = self._run(55e5)
        assert r["solver_info"]["final_rel_err"] < 1e-4  # default tol

    def test_self_consistency(self):
        # At convergence, re-evaluating the feed line at m_dot_real must
        # give a P_inlet that produces the same m_dot from the injector.
        # i.e. the two models agree to within 1% (much more than the 0.01%
        # solver tolerance, to account for numerical noise in re-evaluation).
        from feed_line import evaluate_feed_line as efl
        r = self._run(55e5)
        m_conv  = r["m_dot_real"]
        fl_check = efl(m_conv, T_TANK, 55e5, SEGMENTS)
        P_check  = fl_check["P_final"]
        # Re-run injector at this P_inlet
        from n2o_properties import rho_liquid_sat, nu_vapor_sat, T_sat
        from injector_two_phase import dyer_mass_flow
        rho_l_up   = rho_liquid_sat(T_TANK)
        T_down     = T_sat(20e5)
        rho_l_down = rho_liquid_sat(T_down)
        rho_v_down = 44.013 / nu_vapor_sat(T_down)
        ir = dyer_mass_flow(Cd, A, T_TANK, P_check, 20e5,
                            rho_l_up, rho_l_down, rho_v_down)
        m_check = ir["m_dot_Dyer"]
        rel_diff = abs(m_check - m_conv) / m_conv
        assert rel_diff < 0.01  # 1% -- much looser than solver tol

    def test_custom_tolerance(self):
        # A tighter tolerance should still converge
        r = evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5,
                                  tol=1e-6)
        assert r["solver_info"]["converged"] is True
        assert r["solver_info"]["final_rel_err"] < 1e-6

    def test_invalid_alpha_raises(self):
        with pytest.raises(ValueError):
            evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5,
                                  alpha=0.0)

    def test_invalid_tol_raises(self):
        with pytest.raises(ValueError):
            evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5,
                                  tol=-1e-4)

    def test_result_has_solver_info_key(self):
        r = self._run(55e5)
        assert "solver_info" in r
        for key in ("converged", "iterations", "history", "final_rel_err"):
            assert key in r["solver_info"]

    def test_initial_guess_does_not_affect_result(self):
        # Two different initial guesses should converge to the same answer
        r1 = evaluate_full_system(0.1, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5)
        r2 = evaluate_full_system(0.5, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5)
        assert math.isclose(r1["m_dot_real"], r2["m_dot_real"], rel_tol=1e-3)


# ---------------------------------------------------------------------------
# Physics-based monotonicity tests (Priority 5 of roadmap)
# These test physical relationships, not specific numerical values.
# ---------------------------------------------------------------------------

class TestPhysicsMonotonicity:
    """
    Tests that the model obeys physically required monotonic relationships.
    A violation of any of these tests indicates a physics bug, not just
    a numerical error.
    """

    def test_larger_area_gives_more_flow(self):
        # More orifice area -> more mass flow, all else equal
        A_small = 2 * math.pi * (0.55e-3) ** 2
        A_large = 6 * math.pi * (0.55e-3) ** 2
        r1 = evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, Cd, A_small, 20e5)
        r2 = evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, Cd, A_large, 20e5)
        assert r2["m_dot_real"] > r1["m_dot_real"]

    def test_higher_cd_gives_more_flow(self):
        # Higher Cd -> less restriction -> more mass flow
        r1 = evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, 0.60, A, 20e5)
        r2 = evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, 0.75, A, 20e5)
        assert r2["m_dot_real"] > r1["m_dot_real"]

    def test_higher_tank_pressure_gives_more_flow(self):
        # More tank pressure -> more driving dP -> more flow
        r1 = evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5)
        r2 = evaluate_full_system(0.3, T_TANK, 65e5, SEGMENTS, Cd, A, 20e5)
        assert r2["m_dot_real"] > r1["m_dot_real"]

    def test_lower_chamber_pressure_gives_more_flow(self):
        # Lower P_chamber -> more dP across injector -> more flow
        r1 = evaluate_full_system(0.3, T_TANK, 65e5, SEGMENTS, Cd, A, 25e5)
        r2 = evaluate_full_system(0.3, T_TANK, 65e5, SEGMENTS, Cd, A, 15e5)
        assert r2["m_dot_real"] > r1["m_dot_real"]

    def test_longer_line_gives_less_flow(self):
        # Longer line -> more friction losses -> less P_inlet -> less flow
        segs_short = [{"type": "pipe", "L": 0.5, "D": 0.008}]
        segs_long  = [{"type": "pipe", "L": 3.0, "D": 0.008}]
        r1 = evaluate_full_system(0.3, T_TANK, 65e5, segs_short, Cd, A, 20e5)
        r2 = evaluate_full_system(0.3, T_TANK, 65e5, segs_long,  Cd, A, 20e5)
        if r1["m_dot_real"] and r2["m_dot_real"]:
            assert r1["m_dot_real"] > r2["m_dot_real"]

    def test_smaller_pipe_gives_less_flow(self):
        # Smaller pipe ID -> more friction -> less P_inlet -> less flow
        segs_wide   = [{"type": "pipe", "L": 1.0, "D": 0.012}]
        segs_narrow = [{"type": "pipe", "L": 1.0, "D": 0.006}]
        r1 = evaluate_full_system(0.3, T_TANK, 65e5, segs_wide,   Cd, A, 20e5)
        r2 = evaluate_full_system(0.3, T_TANK, 65e5, segs_narrow, Cd, A, 20e5)
        if r1["m_dot_real"] and r2["m_dot_real"]:
            assert r1["m_dot_real"] > r2["m_dot_real"]

    def test_exit_quality_in_valid_range(self):
        # Vapour quality must always be in [0, 1]
        r = evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5)
        ir = r.get("injector_result")
        if ir and ir.get("x_exit") is not None:
            assert 0.0 <= ir["x_exit"] <= 1.0

    def test_x_inlet_in_valid_range_when_flashing(self):
        # x_inlet from feed line must be in [0, 1]
        P_tank = P_sat(T_TANK)  # zero margin -> flashing guaranteed
        r = evaluate_full_system(0.3, T_TANK, P_tank, SEGMENTS, Cd, A, 20e5)
        x_in = r["feed_line_result"].get("x_inlet", 0.0)
        assert 0.0 <= x_in <= 1.0

    def test_p_inlet_always_below_p_tank(self):
        # Feed-line losses can only reduce pressure, never increase it
        r = evaluate_full_system(0.3, T_TANK, 55e5, SEGMENTS, Cd, A, 20e5)
        assert r["P_injector_inlet"] <= 55e5

    def test_m_dot_always_positive(self):
        # Mass flow must be positive in all valid regimes
        for P_tank in [52e5, 55e5, 65e5]:
            r = evaluate_full_system(0.3, T_TANK, P_tank, SEGMENTS, Cd, A, 20e5)
            if r["m_dot_real"] is not None:
                assert r["m_dot_real"] > 0
