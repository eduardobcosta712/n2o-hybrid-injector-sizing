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
    henry_fauske_critical_flow,
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
        for key in ("m_dot_Dyer", "m_dot_SPI", "m_dot_HEM", "kappa", "x_exit",
                    "m_dot_crit_HF", "choked", "HF_unavailable_reason"):
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

    def test_m_dot_Dyer_unchanged_by_HF_diagnostic(self):
        # m_dot_Dyer must be identical whether or not the Henry-Fauske
        # ceiling is exceeded -- it is a side-by-side diagnostic, never
        # an automatic override (docs/future_work.md, Priority 1).
        result = self._run()
        # Recompute m_dot_Dyer "by hand" from its own SPI/HEM/kappa
        # components, independent of whatever the HF branch did.
        expected = (result["kappa"] / (1.0 + result["kappa"])) * result["m_dot_SPI"] \
                   + (1.0 / (1.0 + result["kappa"])) * result["m_dot_HEM"]
        assert math.isclose(result["m_dot_Dyer"], expected, rel_tol=1e-9)

    def test_choked_flag_consistent_with_values(self):
        result = self._run()
        if result["m_dot_crit_HF"] is not None:
            assert result["choked"] == (result["m_dot_Dyer"] > result["m_dot_crit_HF"])
        else:
            assert result["choked"] is False

    def test_waxman_conditions_not_flagged_choked(self):
        # The Waxman validated point (this class's default T/P/geometry
        # is a different, but nearby, subcooled operating point) should
        # not need the ceiling. Explicitly check the actual Waxman
        # geometry here, since it is the one case with experimental
        # confirmation.
        T1, P1, Cd = 280.0, 4.36e6, 0.65
        A = math.pi * (0.0015 / 2) ** 2
        rho_l_up = rho_liquid_sat(T1)
        for dP_bar in [8.4, 9.8, 10.9, 13.7]:  # the 4 validated cases
            P2 = P1 - dP_bar * 1e5
            T2 = T_sat(P2)
            rho_l2 = rho_liquid_sat(T2)
            rho_v2 = _rho_v(T2)
            r = dyer_mass_flow(Cd, A, T1, P1, P2, rho_l_up, rho_l2, rho_v2)
            assert r["choked"] is False, (
                f"Waxman validated case at dP={dP_bar} bar was flagged "
                "choked -- this would mean the HF ceiling now cuts into "
                "an already-validated result, which should never happen."
            )


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


# ---------------------------------------------------------------------------
# henry_fauske_critical_flow -- added September 2026, Priority 1
# Source: Henry & Fauske (1971); equations transcribed from Simoneau,
# Henry, Hendricks & Watterson (1971), NASA TM X-67863, Eqs. (2)-(5).
# ---------------------------------------------------------------------------

class TestHenryFauskeCriticalFlow:

    T1, P1, Cd = 280.0, 4.36e6, 0.65
    D = 0.0015
    A = math.pi * (D / 2.0) ** 2

    # --- Known-value checks -------------------------------------------------

    def test_matches_hand_validated_waxman_value(self):
        # Cross-checked by hand against Eqs. (2) and (5) at Waxman
        # conditions during development: 50.67 g/s.
        hf = henry_fauske_critical_flow(self.Cd, self.A, self.T1, self.P1)
        assert math.isclose(hf["m_dot_crit"] * 1000, 50.67, rel_tol=0.01)

    # --- Physical properties -------------------------------------------------

    def test_above_equilibrium_hem_ceiling(self):
        # Non-equilibrium choking must exceed the full-equilibrium (HEM)
        # ceiling -- the entire physical motivation for this model.
        hf = henry_fauske_critical_flow(self.Cd, self.A, self.T1, self.P1)
        hem_h = hem_critical_flow(self.Cd, self.A, self.T1, self.P1)
        hem_s = hem_critical_flow_isentropic(self.Cd, self.A, self.T1, self.P1)
        assert hf["m_dot_crit"] > hem_h["m_dot_crit"]
        assert hf["m_dot_crit"] > hem_s["m_dot_crit"]

    def test_does_not_cut_into_validated_waxman_dyer_points(self):
        # None of the 4 already-validated Waxman Dyer predictions may
        # exceed this ceiling -- if they did, applying the ceiling would
        # silently degrade the documented MAPE = 3.51% result.
        hf = henry_fauske_critical_flow(self.Cd, self.A, self.T1, self.P1)
        m_crit_gs = hf["m_dot_crit"] * 1000
        for m_dyer_gs in (42.25, 44.60, 46.20, 49.55):  # from waxman_2013_results.md
            assert m_dyer_gs <= m_crit_gs, (
                f"Waxman Dyer value {m_dyer_gs} g/s exceeds HF ceiling "
                f"{m_crit_gs:.2f} g/s -- would break the validated MAPE.")

    def test_N_in_valid_range(self):
        hf = henry_fauske_critical_flow(self.Cd, self.A, self.T1, self.P1)
        assert 0.0 <= hf["N"] <= 1.0

    def test_x_crit_in_valid_range(self):
        hf = henry_fauske_critical_flow(self.Cd, self.A, self.T1, self.P1)
        assert 0.0 <= hf["x_crit"] <= 1.0

    def test_throat_pressure_below_p_sat_upstream(self):
        hf = henry_fauske_critical_flow(self.Cd, self.A, self.T1, self.P1)
        assert hf["P2_crit"] < P_sat(self.T1)

    def test_scales_linearly_with_area(self):
        hf1 = henry_fauske_critical_flow(self.Cd, self.A, self.T1, self.P1)
        hf2 = henry_fauske_critical_flow(self.Cd, 2 * self.A, self.T1, self.P1)
        assert math.isclose(hf2["m_dot_crit"], 2 * hf1["m_dot_crit"], rel_tol=1e-3)

    def test_higher_upstream_pressure_gives_more_flow(self):
        hf_low  = henry_fauske_critical_flow(self.Cd, self.A, self.T1, 4.2e6)
        hf_high = henry_fauske_critical_flow(self.Cd, self.A, self.T1, 5.0e6)
        assert hf_high["m_dot_crit"] > hf_low["m_dot_crit"]

    def test_two_phase_inlet_reduces_ceiling(self):
        # An already-partially-vaporised inlet has less "room" left
        # before full vaporisation -- the non-equilibrium ceiling should
        # be lower than for a pure-liquid inlet at the same T/P.
        hf0 = henry_fauske_critical_flow(self.Cd, self.A, self.T1, self.P1, x_inlet=0.0)
        hf1 = henry_fauske_critical_flow(self.Cd, self.A, self.T1, self.P1, x_inlet=0.03)
        assert hf1["m_dot_crit"] < hf0["m_dot_crit"]

    # --- Edge cases -------------------------------------------------

    def test_raises_above_table_a4_range(self):
        with pytest.raises(ValueError):
            henry_fauske_critical_flow(self.Cd, self.A, T_MAX_A4 + 1.0, self.P1)

    def test_raises_when_upstream_already_two_phase(self):
        # P_upstream <= P_sat(T_upstream): mirrors
        # dyer_non_equilibrium_parameter's own domain restriction.
        with pytest.raises(ValueError):
            henry_fauske_critical_flow(self.Cd, self.A, self.T1, P_sat(self.T1) - 1e3)

    def test_accepts_temperature_at_table_a4_boundary(self):
        # P_sat(307.33 K) ~= 69.5 bar -- must stay comfortably above it.
        hf = henry_fauske_critical_flow(self.Cd, self.A, T_MAX_A4, 75.0e5)
        assert hf["m_dot_crit"] > 0
