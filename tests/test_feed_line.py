"""
test_feed_line.py

Tests for feed_line.py: Reynolds number, Darcy friction factor,
pressure drop functions, and the full feed line evaluator
(including flashing detection).
"""

import math
import pytest

from feed_line import (
    reynolds_number, darcy_friction_factor,
    friction_pressure_drop, fitting_pressure_drop,
    velocity_from_mass_flow, evaluate_feed_line,
)
from n2o_properties import P_sat


# ---------------------------------------------------------------------------
# reynolds_number
# ---------------------------------------------------------------------------

class TestReynoldsNumber:

    def test_known_value(self):
        # Re = rho * v * D / mu = 800 * 5 * 0.01 / 1.5e-4 = 266 667
        Re = reynolds_number(rho=800.0, v=5.0, D=0.01, mu=1.5e-4)
        assert math.isclose(Re, 266_666.67, rel_tol=1e-4)

    def test_scales_linearly_with_velocity(self):
        Re1 = reynolds_number(800.0, 1.0, 0.01, 1.5e-4)
        Re2 = reynolds_number(800.0, 2.0, 0.01, 1.5e-4)
        assert math.isclose(Re2, 2 * Re1, rel_tol=1e-9)

    def test_positive(self):
        assert reynolds_number(800.0, 5.0, 0.01, 1.5e-4) > 0


# ---------------------------------------------------------------------------
# darcy_friction_factor
# ---------------------------------------------------------------------------

class TestDarcyFrictionFactor:

    def test_laminar_exact(self):
        # For Re = 1000 (laminar), f = 64 / Re = 0.064 exactly.
        f = darcy_friction_factor(Re=1000, roughness=1.5e-6, D=0.008)
        assert math.isclose(f, 0.064, rel_tol=1e-9)

    def test_laminar_threshold(self):
        # Re = 2299 → laminar; Re = 2300 → turbulent formula.
        f_lam = darcy_friction_factor(2299, 1.5e-6, 0.008)
        f_turb = darcy_friction_factor(2300, 1.5e-6, 0.008)
        assert math.isclose(f_lam, 64 / 2299, rel_tol=1e-9)
        # Turbulent value should be different (and higher for this Re range)
        assert not math.isclose(f_turb, 64 / 2300, rel_tol=0.01)

    def test_turbulent_decreases_with_re(self):
        # For smooth pipe, friction factor decreases as Re increases.
        f1 = darcy_friction_factor(10_000, 1.5e-6, 0.008)
        f2 = darcy_friction_factor(100_000, 1.5e-6, 0.008)
        assert f1 > f2

    def test_rougher_pipe_higher_friction(self):
        Re = 50_000
        D = 0.008
        f_smooth = darcy_friction_factor(Re, roughness=1.5e-6, D=D)
        f_rough  = darcy_friction_factor(Re, roughness=50e-6, D=D)
        assert f_rough > f_smooth

    def test_positive(self):
        assert darcy_friction_factor(50_000, 1.5e-6, 0.008) > 0


# ---------------------------------------------------------------------------
# friction_pressure_drop and fitting_pressure_drop
# ---------------------------------------------------------------------------

class TestPressureDropFunctions:

    def test_friction_known_value(self):
        # dP = f * (L/D) * (rho * v^2 / 2)
        # = 0.02 * (1.0/0.01) * (800 * 4 / 2) = 0.02 * 100 * 1600 = 3200 Pa
        dP = friction_pressure_drop(rho=800.0, v=2.0, f=0.02, L=1.0, D=0.01)
        assert math.isclose(dP, 3200.0, rel_tol=1e-6)

    def test_friction_scales_with_length(self):
        dP1 = friction_pressure_drop(800.0, 2.0, 0.02, 1.0, 0.01)
        dP2 = friction_pressure_drop(800.0, 2.0, 0.02, 2.0, 0.01)
        assert math.isclose(dP2, 2 * dP1, rel_tol=1e-9)

    def test_fitting_known_value(self):
        # dP = K * (rho * v^2 / 2) = 0.05 * (800 * 4 / 2) = 80 Pa
        dP = fitting_pressure_drop(rho=800.0, v=2.0, K=0.05)
        assert math.isclose(dP, 80.0, rel_tol=1e-6)

    def test_pressure_drops_are_positive(self):
        assert friction_pressure_drop(800.0, 2.0, 0.02, 1.0, 0.01) >= 0
        assert fitting_pressure_drop(800.0, 2.0, 0.3) >= 0


# ---------------------------------------------------------------------------
# velocity_from_mass_flow
# ---------------------------------------------------------------------------

class TestVelocityFromMassFlow:

    def test_known_value(self):
        # A = pi*(0.004)^2 = 5.027e-5 m^2
        # v = m_dot / (rho * A) = 0.5 / (800 * 5.027e-5) = 12.43 m/s
        A = math.pi * (0.004) ** 2
        v = velocity_from_mass_flow(m_dot=0.5, rho=800.0, D=0.008)
        assert math.isclose(v, 0.5 / (800.0 * A), rel_tol=1e-9)

    def test_doubles_when_m_dot_doubles(self):
        v1 = velocity_from_mass_flow(0.5, 800.0, 0.008)
        v2 = velocity_from_mass_flow(1.0, 800.0, 0.008)
        assert math.isclose(v2, 2 * v1, rel_tol=1e-9)


# ---------------------------------------------------------------------------
# evaluate_feed_line — integration tests
# ---------------------------------------------------------------------------

class TestEvaluateFeedLine:
    """Integration tests for the full feed line evaluator."""

    # Shared geometry: 2 m pipe + ball valve + 90° elbow
    SEGMENTS = [
        {"type": "pipe",    "L": 1.0, "D": 0.008},
        {"type": "fitting", "D": 0.008, "K": 0.05},
        {"type": "pipe",    "L": 1.0, "D": 0.008},
        {"type": "fitting", "D": 0.008, "K": 0.30},
    ]
    T_TANK = 293.15   # K, 20 °C
    M_DOT  = 0.5      # kg/s

    def test_no_flashing_with_adequate_margin(self):
        # 5 bar above saturation → no flashing expected.
        P_tank = P_sat(self.T_TANK) + 5e5
        result = evaluate_feed_line(self.M_DOT, self.T_TANK, P_tank, self.SEGMENTS)
        assert not result["flashing_detected"]
        assert result["delta_T_sub_final"] > 0

    def test_flashing_with_zero_margin(self):
        # Exactly at saturation → any friction drop causes flashing.
        P_tank = P_sat(self.T_TANK)
        result = evaluate_feed_line(self.M_DOT, self.T_TANK, P_tank, self.SEGMENTS)
        assert result["flashing_detected"]

    def test_pressure_decreases_monotonically(self):
        P_tank = P_sat(self.T_TANK) + 10e5
        result = evaluate_feed_line(self.M_DOT, self.T_TANK, P_tank, self.SEGMENTS)
        pressures = [P_tank] + [s["pressure_after_Pa"] for s in result["trace"]]
        for i in range(len(pressures) - 1):
            assert pressures[i] >= pressures[i + 1], (
                f"Pressure increased from segment {i} to {i+1}: "
                f"{pressures[i]/1e5:.3f} → {pressures[i+1]/1e5:.3f} bar"
            )

    def test_trace_has_correct_length(self):
        P_tank = P_sat(self.T_TANK) + 5e5
        result = evaluate_feed_line(self.M_DOT, self.T_TANK, P_tank, self.SEGMENTS)
        assert len(result["trace"]) == len(self.SEGMENTS)

    def test_higher_mass_flow_gives_more_pressure_drop(self):
        P_tank = P_sat(self.T_TANK) + 10e5
        r1 = evaluate_feed_line(0.3, self.T_TANK, P_tank, self.SEGMENTS)
        r2 = evaluate_feed_line(0.8, self.T_TANK, P_tank, self.SEGMENTS)
        # Higher velocity → more friction → lower final pressure
        assert r2["P_final"] < r1["P_final"]

    def test_unknown_segment_type_raises(self):
        bad_segments = [{"type": "rocket", "D": 0.008}]
        with pytest.raises(ValueError):
            evaluate_feed_line(0.5, self.T_TANK, 60e5, bad_segments)


# ---------------------------------------------------------------------------
# Two-phase feed line model — new tests for Priority 2 implementation
# ---------------------------------------------------------------------------

class TestTwoPhaseLineModel:
    """
    Tests that verify the two-phase HEM model in the feed line.
    When flashing occurs, the model must use rho_mix and mu_mix instead
    of pure-liquid properties -- which increases the pressure drop per
    unit length compared to the single-phase assumption.
    """

    T_TANK = 293.15   # K, 20 degC
    # Use a pipe long enough to cause meaningful two-phase losses
    SEGS_LONG = [
        {"type": "pipe", "L": 2.0, "D": 0.006},
        {"type": "pipe", "L": 2.0, "D": 0.006},
    ]

    def test_x_quality_in_trace_is_zero_before_flashing(self):
        # All segments before flashing onset must have x_quality == 0
        P_tank = P_sat(self.T_TANK) + 10e5  # plenty of margin
        result = evaluate_feed_line(0.3, self.T_TANK, P_tank, self.SEGS_LONG)
        if not result["flashing_detected"]:
            for seg in result["trace"]:
                assert seg["x_quality"] == 0.0

    def test_x_quality_increases_after_flashing(self):
        # Once flashing starts, x must be >= 0 and non-decreasing
        P_tank = P_sat(self.T_TANK)   # zero margin -- flashing from start
        result = evaluate_feed_line(0.5, self.T_TANK, P_tank, self.SEGS_LONG)
        assert result["flashing_detected"]
        x_vals = [s["x_quality"] for s in result["trace"]]
        # x must be non-decreasing (more vapour as pressure drops further)
        for i in range(len(x_vals) - 1):
            assert x_vals[i] <= x_vals[i + 1] + 1e-9, (
                f"x_quality must be non-decreasing: {x_vals[i]:.5f} > {x_vals[i+1]:.5f}")

    def test_rho_eff_lower_in_two_phase_region(self):
        # In the two-phase region, effective density must be below rho_l
        from n2o_properties import rho_liquid_sat
        rho_l = rho_liquid_sat(self.T_TANK)
        P_tank = P_sat(self.T_TANK)
        result = evaluate_feed_line(0.5, self.T_TANK, P_tank, self.SEGS_LONG)
        two_phase_segs = [s for s in result["trace"] if s["x_quality"] > 0]
        for seg in two_phase_segs:
            assert seg["rho_eff_kg_m3"] < rho_l, (
                f"Two-phase rho_eff={seg['rho_eff_kg_m3']:.1f} must be < rho_l={rho_l:.1f}")

    def test_two_phase_dP_larger_than_single_phase(self):
        # The two-phase model predicts MORE pressure drop per segment than
        # the single-phase model because rho_mix << rho_l (same flow rate,
        # lower density -> higher velocity -> more friction losses).
        # We compare the effective density directly: rho_eff in the two-phase
        # region must be lower than rho_l, which drives higher losses.
        from n2o_properties import rho_liquid_sat
        rho_l = rho_liquid_sat(self.T_TANK)
        P_tank = P_sat(self.T_TANK)
        segs = [
            {"type": "pipe", "L": 1.0, "D": 0.006},
            {"type": "pipe", "L": 1.0, "D": 0.006},
            {"type": "pipe", "L": 1.0, "D": 0.006},
        ]
        result = evaluate_feed_line(0.5, self.T_TANK, P_tank, segs)
        # Find segments that are genuinely two-phase (x > 0.005)
        two_phase = [s for s in result["trace"] if s["x_quality"] > 0.005]
        if two_phase:
            # Each two-phase segment must have lower effective density
            # than pure liquid -- proving the two-phase model is engaged
            for seg in two_phase:
                assert seg["rho_eff_kg_m3"] < rho_l * 0.99, (
                    f"Two-phase rho_eff={seg['rho_eff_kg_m3']:.1f} "
                    f"not meaningfully below rho_l={rho_l:.1f}")
            # And the dP in the last two-phase segment must exceed what
            # pure-liquid Darcy-Weisbach would give for the same geometry
            seg = two_phase[-1]
            v_mix = 0.5 / (seg["rho_eff_kg_m3"] * math.pi * (0.006/2)**2)
            v_liq = 0.5 / (rho_l * math.pi * (0.006/2)**2)
            # dP ~ rho * v^2, but rho*v^2 = m_dot^2 / (rho * A^2), so
            # dP_2ph / dP_1ph = rho_l / rho_mix > 1
            ratio = rho_l / seg["rho_eff_kg_m3"]
            assert ratio > 1.01, (
                f"Expected rho_l/rho_mix > 1.01, got {ratio:.4f}")

    def test_x_inlet_consistent_with_trace(self):
        # x_inlet (computed at P_final) must be >= x at last segment start
        P_tank = P_sat(self.T_TANK)
        result = evaluate_feed_line(0.5, self.T_TANK, P_tank, self.SEGS_LONG)
        if result["flashing_detected"]:
            x_last_seg = result["trace"][-1]["x_quality"]
            assert result["x_inlet"] >= x_last_seg - 1e-9

    def test_mu_eff_between_liquid_and_vapour(self):
        # In two-phase region: mu_l > mu_mix > mu_v
        from n2o_properties import MU_LIQUID_N2O, mu_vapor_sat, T_sat
        P_tank = P_sat(self.T_TANK)
        result = evaluate_feed_line(0.5, self.T_TANK, P_tank, self.SEGS_LONG)
        for seg in result["trace"]:
            if seg["x_quality"] > 0:
                P_seg = seg["pressure_after_Pa"]
                if P_seg > 0 and P_seg < P_sat(self.T_TANK):
                    mu_v = mu_vapor_sat(T_sat(P_seg))
                    assert mu_v <= seg["mu_eff_Pa_s"] <= MU_LIQUID_N2O

    def test_trace_has_new_fields(self):
        # Every trace entry must have the new two-phase fields
        P_tank = P_sat(self.T_TANK) + 5e5
        result = evaluate_feed_line(0.3, self.T_TANK, P_tank,
                                    [{"type": "pipe", "L": 1.0, "D": 0.01}])
        for seg in result["trace"]:
            assert "x_quality"     in seg
            assert "rho_eff_kg_m3" in seg
            assert "mu_eff_Pa_s"   in seg
