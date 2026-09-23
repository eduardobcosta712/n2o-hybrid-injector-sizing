import math
import pytest

from n2o_properties import (
    P_sat, dP_sat_dT, T_sat, rho_liquid_sat,
    nu_vapor_sat, h_liquid_sat, h_vapor_sat, h_fg,
    mu_vapor_sat, mu_mixture, mu_liquid_sat, cp_liquid_sat,
    s_liquid_sat, s_vapor_sat, s_fg,
    degree_of_subcooling,
    M_N2O, MU_LIQUID_N2O,
    T_MIN, T_MAX, T_MIN_A4, T_MAX_A4, T_MIN_A3, T_MAX_A3, T_CRIT, P_CRIT,
    _load_saturation_table, _interp,
)


def _perry_P_sat(T):
    c1, c2, c3, c4, c5 = 96.512, -4045.0, -12.277, 2.886e-5, 2.0
    return math.exp(c1 + c2 / T + c3 * math.log(T) + c4 * T ** c5)


def _perry_rho_l(T):
    c1, c2, c3, c4 = 2.781, 0.27244, 309.57, 0.2882
    return 44.013 / (c2 ** (1 + (1 - T / c3) ** c4) / c1)


def _mcgill(T, key):
    tbl = _load_saturation_table()
    return _interp(T, tbl["T_K"], tbl[key])


def _nist_a4(T, key):
    tbl = _load_saturation_table()
    return _interp(T, tbl["T_K_a4"], tbl[key])


class TestPsat:

    def test_at_0_degC(self):
        assert math.isclose(P_sat(273.15), 31.3e5, rel_tol=0.03)

    def test_at_20_degC(self):
        assert math.isclose(P_sat(293.15), 50.9e5, rel_tol=0.015)

    def test_critical_constants(self):
        assert math.isclose(T_CRIT, 309.52, rel_tol=1e-4)
        assert math.isclose(P_CRIT, 7.245e6, rel_tol=1e-3)

    def test_near_critical_point_reproduces_p_crit(self):
        assert math.isclose(P_sat(T_MAX), P_CRIT, rel_tol=0.005)

    def test_below_critical_pressure_just_under_critical_temperature(self):
        P = P_sat(307.0)
        assert P < P_CRIT
        assert P > 0.90 * P_CRIT

    def test_agrees_with_perry_correlation(self):
        for T in range(230, 306, 5):
            assert math.isclose(P_sat(float(T)), _perry_P_sat(float(T)), rel_tol=0.06), T

    def test_monotonically_increasing(self):
        T_vals = [T_MIN + i * 5 for i in range(25)]
        P_vals = [P_sat(T) for T in T_vals]
        for i in range(len(P_vals) - 1):
            assert P_vals[i] < P_vals[i + 1]

    def test_positive_throughout(self):
        for T in [T_MIN, 250.0, T_MAX]:
            assert P_sat(T) > 0

    def test_out_of_range_low(self):
        with pytest.raises(ValueError):
            P_sat(T_MIN - 1.0)

    def test_out_of_range_high(self):
        with pytest.raises(ValueError):
            P_sat(T_MAX + 1.0)


class TestdPsatdT:

    def test_positive_throughout(self):
        for T in [T_MIN, 250.0, 290.0, T_MAX]:
            assert dP_sat_dT(T) > 0

    def test_consistent_with_finite_difference(self):
        for T in (250.0, 293.15, 305.0):
            h = 0.05
            fd = (P_sat(T + h) - P_sat(T - h)) / (2 * h)
            assert math.isclose(dP_sat_dT(T), fd, rel_tol=0.05), T

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            dP_sat_dT(T_MIN - 1.0)


class TestTsat:

    def test_inverse_at_several_temperatures(self):
        for T_original in [200.0, 250.0, 270.0, 290.0, 305.0]:
            P = P_sat(T_original)
            T_recovered = T_sat(P)
            assert math.isclose(T_recovered, T_original, rel_tol=1e-7)

    def test_monotonically_increasing(self):
        P_vals = [20e5, 30e5, 40e5, 50e5, 60e5, 70e5]
        T_vals = [T_sat(P) for P in P_vals]
        for i in range(len(T_vals) - 1):
            assert T_vals[i] < T_vals[i + 1]

    def test_at_20_degC_reference(self):
        T = T_sat(P_sat(293.15))
        assert math.isclose(T, 293.15, rel_tol=1e-7)

    def test_pressure_above_critical_range_raises_clear_error(self):
        with pytest.raises(ValueError, match="outside the range"):
            T_sat(80e5)

    def test_pressure_below_triple_point_range_raises_clear_error(self):
        with pytest.raises(ValueError, match="outside the range"):
            T_sat(0.5e5)

    def test_interior_pressures_still_invert(self):
        for T in (200.0, 273.15, 307.0):
            assert math.isclose(T_sat(P_sat(T)), T, rel_tol=1e-7)


class TestRhoLiquidSat:

    def test_at_20_degC(self):
        assert math.isclose(rho_liquid_sat(293.15), 786.0, rel_tol=0.005)

    def test_agrees_with_perry_correlation(self):
        for T in range(230, 301, 10):
            assert math.isclose(rho_liquid_sat(float(T)), _perry_rho_l(float(T)), rel_tol=0.01), T

    def test_decreases_with_temperature(self):
        T_vals = [220.0, 250.0, 270.0, 290.0, 305.0]
        rho_vals = [rho_liquid_sat(T) for T in T_vals]
        for i in range(len(rho_vals) - 1):
            assert rho_vals[i] > rho_vals[i + 1]

    def test_positive_throughout(self):
        for T in [T_MIN, 250.0, T_MAX]:
            assert rho_liquid_sat(T) > 0

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            rho_liquid_sat(400.0)

    def test_molar_mass_constant_value(self):
        assert math.isclose(M_N2O, 44.013, rel_tol=1e-4)

    def test_project_modules_share_the_same_molar_mass(self):
        import feed_line, injector_two_phase, full_system
        assert feed_line.M_N2O is M_N2O
        assert injector_two_phase.M_N2O is M_N2O
        assert full_system.M_N2O is M_N2O


class TestNuVaporSat:

    def test_agrees_with_mcgill_table(self):
        for T in (230.0, 250.0, 270.0, 290.0, 300.0):
            assert math.isclose(nu_vapor_sat(T), _mcgill(T, "nu_v"), rel_tol=0.05), T

    def test_consistent_with_vapour_density(self):
        T = 285.0
        assert math.isclose(nu_vapor_sat(T) * (M_N2O / nu_vapor_sat(T)), M_N2O, rel_tol=1e-12)

    def test_decreases_toward_critical_point(self):
        T_vals = [220.0, 250.0, 270.0, 290.0, 305.0]
        nu_vals = [nu_vapor_sat(T) for T in T_vals]
        for i in range(len(nu_vals) - 1):
            assert nu_vals[i] > nu_vals[i + 1]

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            nu_vapor_sat(T_MIN - 5.0)


class TestHfg:

    def test_positive_throughout(self):
        for T in [220.0, 250.0, 270.0, 290.0, 305.0]:
            assert h_fg(T) > 0

    def test_decreases_toward_critical_point(self):
        T_vals = [220.0, 250.0, 270.0, 290.0, 305.0, T_MAX]
        hfg_vals = [h_fg(T) for T in T_vals]
        for i in range(len(hfg_vals) - 1):
            assert hfg_vals[i] > hfg_vals[i + 1]

    def test_nearly_vanishes_at_critical_point(self):
        assert h_fg(T_MAX) < 0.25 * h_fg(290.0)

    def test_consistent_with_clausius_clapeyron(self):
        for T in (250.0, 270.0, 290.0):
            nu_v = nu_vapor_sat(T)
            nu_l = M_N2O / rho_liquid_sat(T)
            h = 0.05
            dPdT = (P_sat(T + h) - P_sat(T - h)) / (2 * h)
            h_fg_cc = T * (nu_v - nu_l) * dPdT / 1e3
            assert math.isclose(h_fg_cc, h_fg(T), rel_tol=0.05), T

    def test_agrees_with_mcgill_table_within_known_error(self):
        for T in (230.0, 250.0, 270.0, 290.0, 300.0):
            legacy = _mcgill(T, "h_v") - _mcgill(T, "h_l")
            assert math.isclose(h_fg(T), legacy, rel_tol=0.06), T


class TestDegreeOfSubcooling:

    def test_exactly_saturated_is_zero(self):
        T = 293.15
        P = P_sat(T)
        dT_sub = degree_of_subcooling(T, P)
        assert math.isclose(dT_sub, 0.0, abs_tol=1e-5)

    def test_subcooled_is_positive(self):
        T = 293.15
        P = P_sat(T) + 5e5
        assert degree_of_subcooling(T, P) > 0

    def test_superheated_is_negative(self):
        T = 293.15
        P = P_sat(T) - 5e5
        assert degree_of_subcooling(T, P) < 0

    def test_increases_with_pressure(self):
        T = 280.0
        P_low  = P_sat(T) + 2e5
        P_high = P_sat(T) + 10e5
        assert degree_of_subcooling(T, P_high) > degree_of_subcooling(T, P_low)


class TestMuVaporSat:

    def test_known_values_from_nist(self):
        assert math.isclose(mu_vapor_sat(182.33) * 1e6, 9.0689, rel_tol=1e-4)
        assert math.isclose(mu_vapor_sat(252.33) * 1e6, 13.417, rel_tol=1e-4)
        assert math.isclose(mu_vapor_sat(307.33) * 1e6, 22.982, rel_tol=1e-4)

    def test_increases_with_temperature(self):
        temps = [190.0, 220.0, 250.0, 270.0, 295.0]
        values = [mu_vapor_sat(T) for T in temps]
        for i in range(len(values) - 1):
            assert values[i] < values[i + 1]

    def test_returns_pa_s_not_upa_s(self):
        muv = mu_vapor_sat(250.0)
        assert 5e-6 < muv < 30e-6

    def test_out_of_range_raises(self):
        with pytest.raises(ValueError):
            mu_vapor_sat(180.0)
        with pytest.raises(ValueError):
            mu_vapor_sat(315.0)

    def test_range_between_table_end_and_critical_point_raises_valueerror(self):
        for T in (T_MAX_A3 + 0.5, 308.0, 309.0, T_MAX):
            with pytest.raises(ValueError, match="Table A.3"):
                mu_vapor_sat(T)

    def test_table_a3_range_matches_table_a4_range(self):
        assert (T_MIN_A3, T_MAX_A3) == (T_MIN_A4, T_MAX_A4)

    def test_viscosity_tables_narrower_than_thermodynamic_range(self):
        assert T_MAX_A4 < T_MAX
        assert T_MIN_A4 <= T_MIN

    def test_much_less_than_liquid(self):
        for T in [200.0, 250.0, 295.0]:
            assert mu_vapor_sat(T) < MU_LIQUID_N2O


class TestMuMixture:

    def test_pure_liquid_x0(self):
        assert math.isclose(mu_mixture(0.0, T=250.0), mu_liquid_sat(250.0))

    def test_pure_vapour_x1(self):
        assert math.isclose(mu_mixture(1.0, T=250.0), mu_vapor_sat(250.0))

    def test_between_liquid_and_vapour(self):
        mu_mix = mu_mixture(0.3, T=260.0)
        mu_v   = mu_vapor_sat(260.0)
        assert mu_v < mu_mix < MU_LIQUID_N2O

    def test_decreases_with_vapour_quality(self):
        x_vals = [0.0, 0.1, 0.3, 0.5, 0.8, 1.0]
        mu_vals = [mu_mixture(x, T=270.0) for x in x_vals]
        for i in range(len(mu_vals) - 1):
            assert mu_vals[i] >= mu_vals[i + 1]

    def test_above_table_range_raises(self):
        with pytest.raises(ValueError):
            mu_mixture(0.1, T=308.5)

    def test_without_temperature_uses_constants(self):
        assert math.isclose(mu_mixture(0.0), MU_LIQUID_N2O)
        assert math.isclose(mu_mixture(1.0), 13e-6)


class TestMuLiquidSat:

    def test_exact_table_row(self):
        assert math.isclose(mu_liquid_sat(292.33), 6.918200e-05, rel_tol=1e-6)

    def test_decreases_with_temperature(self):
        T_vals = [200.0, 240.0, 270.0, 290.0, 305.0]
        mu_vals = [mu_liquid_sat(T) for T in T_vals]
        for i in range(len(mu_vals) - 1):
            assert mu_vals[i] > mu_vals[i + 1]

    def test_constant_fallback_is_above_room_temperature_value(self):
        assert mu_liquid_sat(293.15) < MU_LIQUID_N2O

    def test_out_of_range_raises(self):
        with pytest.raises(ValueError):
            mu_liquid_sat(T_MAX_A4 + 1.0)
        with pytest.raises(ValueError):
            mu_liquid_sat(T_MIN_A4 - 1.0)


class TestLiquidCp:

    def test_agrees_with_nist_table_rows(self):
        for T in (252.33, 272.33, 292.33):
            assert math.isclose(cp_liquid_sat(T), _nist_a4(T, "cp_l"), rel_tol=5e-3), T

    def test_increases_toward_critical_point(self):
        assert cp_liquid_sat(305.0) > cp_liquid_sat(290.0) > cp_liquid_sat(250.0)

    def test_out_of_range_raises(self):
        with pytest.raises(ValueError):
            cp_liquid_sat(T_MIN - 1.0)


class TestEntropyFunctions:

    def test_s_fg_agrees_with_nist_table_rows(self):
        for T in (252.33, 272.33, 292.33):
            nist = _nist_a4(T, "s_v") - _nist_a4(T, "s_l")
            assert math.isclose(s_fg(T), nist, rel_tol=5e-3), T

    def test_entropy_difference_agrees_with_nist_table(self):
        d_model = s_liquid_sat(292.33) - s_liquid_sat(252.33)
        d_nist = _nist_a4(292.33, "s_l") - _nist_a4(252.33, "s_l")
        assert math.isclose(d_model, d_nist, rel_tol=5e-3)

    def test_s_fg_equals_vapour_minus_liquid(self):
        assert math.isclose(s_fg(270.0), s_vapor_sat(270.0) - s_liquid_sat(270.0), rel_tol=1e-12)

    def test_s_fg_positive_throughout_range(self):
        for T in [190.0, 220.0, 260.0, 300.0, T_MAX]:
            assert s_fg(T) > 0

    def test_s_fg_decreases_toward_critical_point(self):
        T_vals = [200.0, 240.0, 270.0, 300.0, T_MAX]
        sfg_vals = [s_fg(T) for T in T_vals]
        for i in range(len(sfg_vals) - 1):
            assert sfg_vals[i] > sfg_vals[i + 1]

    def test_s_liquid_increases_with_temperature(self):
        T_vals = [190.0, 220.0, 260.0, 300.0]
        s_vals = [s_liquid_sat(T) for T in T_vals]
        for i in range(len(s_vals) - 1):
            assert s_vals[i] < s_vals[i + 1]

    def test_s_fg_units_are_kJ_per_kmolK(self):
        assert 10.0 < s_fg(260.0) < 100.0

    def test_raises_above_thermodynamic_range(self):
        s_liquid_sat(T_MAX_A4 + 1.0)
        with pytest.raises(ValueError):
            s_liquid_sat(T_MAX + 1.0)

    def test_raises_below_range(self):
        with pytest.raises(ValueError):
            s_vapor_sat(T_MIN - 1.0)

    def test_s_fg_raises_outside_range(self):
        with pytest.raises(ValueError):
            s_fg(T_MAX + 1.0)
