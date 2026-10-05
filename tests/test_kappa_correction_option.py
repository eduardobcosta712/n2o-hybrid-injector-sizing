"""
test_kappa_correction_option.py

Tests for the OPTIONAL, default-off `kappa_correction` argument of
full_system.evaluate_full_system() / design_injector_area()
(supercharge-gated correction of Dyer's kappa, exploratory --
docs/future_work.md, Priority 1).

Three groups, per project convention:
  1. The option is OFF by default and inert when it should be.
  2. Known / bracketing relations when it is ON.
  3. Edge cases (invalid arguments) and the continuity property that
     motivates the correction (corrected Dyer -> HEM as supercharge -> 0).
"""

import math
import pytest

from full_system import evaluate_full_system, design_injector_area
from injector_two_phase import dyer_mass_flow_corrected, hem_mass_flow_two_phase_inlet
from n2o_properties import (P_sat, rho_liquid_sat, nu_vapor_sat, T_sat, M_N2O)

PSI = 6894.757
KC = {"beta": 1.0, "supercharge_ref_Pa": 400.0 * PSI}

SEGMENTS = [
    {"type": "pipe",    "L": 1.0, "D": 0.008},
    {"type": "fitting", "D": 0.008, "K": 0.05},
    {"type": "pipe",    "L": 1.0, "D": 0.008},
    {"type": "fitting", "D": 0.008, "K": 0.30},
]
T_TANK = 293.15
Cd = 0.65
A = 4 * math.pi * (0.55e-3) ** 2
P_TANK = 55e5      # 4.5 bar above P_sat(20 degC): Dyer regime, low supercharge
P_CH = 20e5


def _run(kc=None):
    return evaluate_full_system(0.4, T_TANK, P_TANK, SEGMENTS, Cd, A, P_CH,
                                kappa_correction=kc)


# ---------------------------------------------------------------------------
# 1. Off by default / inert
# ---------------------------------------------------------------------------

class TestOptionIsOffByDefault:

    def test_none_equals_argument_omitted(self):
        r0 = evaluate_full_system(0.4, T_TANK, P_TANK, SEGMENTS, Cd, A, P_CH)
        r1 = _run(None)
        assert math.isclose(r0["m_dot_real"], r1["m_dot_real"], rel_tol=1e-12)
        assert "m_dot_Dyer_uncorrected" not in r1["injector_result"]

    def test_beta_zero_recovers_baseline(self):
        base = _run(None)
        r = _run({"beta": 0.0, "supercharge_ref_Pa": 400.0 * PSI})
        assert math.isclose(r["m_dot_real"], base["m_dot_real"], rel_tol=1e-9)

    def test_gate_closed_when_supercharge_above_reference(self):
        # reference far below the actual supercharge -> factor exactly 1
        base = _run(None)
        r = _run({"beta": 2.0, "supercharge_ref_Pa": 1.0})
        assert math.isclose(r["m_dot_real"], base["m_dot_real"], rel_tol=1e-9)

    def test_spi_regime_is_unaffected(self):
        # 0 degC, chamber above P_sat: single-phase through the orifice
        kw = dict(T_tank=273.15, P_tank=60e5, segments=SEGMENTS, Cd=Cd,
                  A_injector=A, P_chamber=50e5)
        r0 = evaluate_full_system(0.3, **kw)
        r1 = evaluate_full_system(0.3, kappa_correction=KC, **kw)
        assert r1["regime"] == "SPI"
        assert math.isclose(r0["m_dot_real"], r1["m_dot_real"], rel_tol=1e-12)


# ---------------------------------------------------------------------------
# 2. Option ON
# ---------------------------------------------------------------------------

class TestCorrectionActive:

    def setup_method(self):
        self.base = _run(None)
        self.corr = _run(KC)

    def test_still_dyer_regime_and_converged(self):
        assert self.corr["regime"] == "Dyer"
        assert self.corr["solver_info"]["converged"] is True

    def test_flow_lower_than_baseline_at_low_supercharge(self):
        # supercharge ~4.5 bar (~65 psi) << 400 psi: kappa is reduced, weight
        # moves towards HEM, flow drops.
        assert self.corr["m_dot_real"] < self.base["m_dot_real"]

    def test_extra_keys_present_and_consistent(self):
        ir = self.corr["injector_result"]
        for key in ("m_dot_Dyer_uncorrected", "kappa_corrected",
                    "kappa_correction_factor"):
            assert key in ir
        assert 0.0 < ir["kappa_correction_factor"] < 1.0
        assert ir["kappa_corrected"] < ir["kappa"]
        assert math.isclose(ir["kappa_corrected"],
                            ir["kappa"] * ir["kappa_correction_factor"], rel_tol=1e-9)

    def test_corrected_flow_bracketed_by_hem_and_uncorrected_dyer(self):
        ir = self.corr["injector_result"]
        m = self.corr["m_dot_real"]
        assert ir["m_dot_HEM"] * 0.999 <= m <= ir["m_dot_Dyer_uncorrected"] * 1.001

    def test_choked_flag_reevaluated_against_corrected_flow(self):
        ir = self.corr["injector_result"]
        if ir["m_dot_crit_HF"] is not None:
            assert ir["choked"] == (ir["m_dot_Dyer"] > ir["m_dot_crit_HF"])

    def test_more_tank_pressure_still_gives_more_flow(self):
        r1 = evaluate_full_system(0.4, T_TANK, 55e5, SEGMENTS, Cd, A, P_CH,
                                  kappa_correction=KC)
        r2 = evaluate_full_system(0.4, T_TANK, 65e5, SEGMENTS, Cd, A, P_CH,
                                  kappa_correction=KC)
        assert r2["m_dot_real"] > r1["m_dot_real"]


class TestDesignWithCorrection:

    T, Pin, Pc, m = 288.15, 59.75e5, 20e5, 0.5   # supercharge ~14.7 bar (~213 psi) < 400 psi

    def test_corrected_area_is_larger_and_hits_the_target(self):
        base = design_injector_area(self.m, Cd, self.T, self.Pin, self.Pc)
        corr = design_injector_area(self.m, Cd, self.T, self.Pin, self.Pc,
                                    kappa_correction=KC)
        assert corr["regime"] == "Dyer"
        assert math.isclose(corr["dyer_result"]["m_dot_Dyer"], self.m, rel_tol=1e-6)
        # lower flow per unit area -> more area needed for the same target
        assert corr["A_recommended"] > base["A_recommended"]

    def test_default_unchanged(self):
        d0 = design_injector_area(self.m, Cd, self.T, self.Pin, self.Pc)
        d1 = design_injector_area(self.m, Cd, self.T, self.Pin, self.Pc,
                                  kappa_correction=None)
        assert math.isclose(d0["A_recommended"], d1["A_recommended"], rel_tol=1e-12)


# ---------------------------------------------------------------------------
# 3. Edge cases and continuity
# ---------------------------------------------------------------------------

class TestInvalidArguments:

    @pytest.mark.parametrize("bad", [
        {"beta": -0.1, "supercharge_ref_Pa": 1e6},
        {"beta": 1.0, "supercharge_ref_Pa": 0.0},
        {"beta": 1.0, "supercharge_ref_Pa": -5.0},
        {"beta": 1.0},
        {"beta": 1.0, "supercharge_ref_Pa": 1e6, "extra": 1},
        [1.0, 1e6],
    ])
    def test_evaluate_full_system_rejects(self, bad):
        with pytest.raises(ValueError):
            evaluate_full_system(0.4, T_TANK, P_TANK, SEGMENTS, Cd, A, P_CH,
                                 kappa_correction=bad)

    def test_design_injector_area_rejects(self):
        with pytest.raises(ValueError):
            design_injector_area(0.5, Cd, 288.15, 59.75e5, 20e5,
                                 kappa_correction={"beta": -1.0,
                                                   "supercharge_ref_Pa": 1e6})


class TestContinuityAtTheFlashingThreshold:
    """The reason for the correction: with beta > 0 and the supercharge going
    to zero, kappa' -> 0, so the corrected Dyer flow tends to the HEM flow --
    which is exactly what the two-phase-inlet branch gives at x_inlet -> 0.
    (Uncorrected Dyer has kappa -> 1 there, i.e. (SPI + HEM)/2: a jump.)"""

    T, P2 = 280.0, 30e5

    def _flows(self, S_Pa):
        T_d = T_sat(self.P2)
        rl_d, rv_d = rho_liquid_sat(T_d), M_N2O / nu_vapor_sat(T_d)
        Pu = P_sat(self.T) + S_Pa
        d = dyer_mass_flow_corrected(Cd, A, self.T, Pu, self.P2,
                                     rho_liquid_sat(self.T), rl_d, rv_d,
                                     beta=1.0, supercharge_ref_Pa=400.0 * PSI)
        hem = hem_mass_flow_two_phase_inlet(Cd, A, self.T, 0.0, Pu, self.P2,
                                            rl_d, rv_d)["m_dot_HEM_2phase"]
        d0 = (d["kappa"] / (1 + d["kappa"])) * d["m_dot_SPI"] \
             + (1 / (1 + d["kappa"])) * d["m_dot_HEM"]
        return d["m_dot_Dyer_corrected"], hem, d0

    def test_corrected_dyer_joins_hem_as_supercharge_vanishes(self):
        corr, hem, _ = self._flows(1e3)          # 1 kPa above saturation
        assert math.isclose(corr, hem, rel_tol=5e-3)

    def test_uncorrected_dyer_does_not(self):
        _, hem, d0 = self._flows(1e3)
        assert d0 > 1.05 * hem                    # the discontinuity being fixed

    def test_gap_to_hem_shrinks_monotonically_with_supercharge(self):
        gaps = []
        for S in (1e3, 1e4, 1e5, 1e6):
            corr, hem, _ = self._flows(S)
            gaps.append(corr / hem - 1.0)
        assert all(g1 <= g2 + 1e-12 for g1, g2 in zip(gaps, gaps[1:]))
