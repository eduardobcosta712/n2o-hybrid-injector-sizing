"""
test_grain_sizing.py

Tests for grain_sizing.py: Marxman regression rate, oxidiser mass flux,
fuel mass flow per port, initial port radius solver, and the top-level
size_grain() sizing function.

Three categories per function, per project convention:
  1. Known-value checks  — hand-computed or closed-form cross-checks.
  2. Physical properties — monotonicity, signs, special cases (n=0.5).
  3. Edge cases          — invalid inputs, unreachable targets.
"""

import math
import pytest

from grain_sizing import (
    regression_rate, oxidizer_mass_flux, fuel_mass_flow_rate_per_port,
    solve_initial_port_radius, size_grain, FUEL_PROPERTIES,
)


# ---------------------------------------------------------------------------
# regression_rate
# ---------------------------------------------------------------------------

class TestRegressionRate:

    def test_known_value(self):
        # r_dot = a * G_o^n = 1e-4 * 100^0.5 = 1e-4 * 10 = 1e-3
        assert math.isclose(regression_rate(1e-4, 0.5, 100.0), 1e-3, rel_tol=1e-9)

    def test_n_zero_gives_constant_rate(self):
        # n=0 -> r_dot = a regardless of G_o (degenerate, but must not error)
        assert math.isclose(regression_rate(5e-4, 0.0, 1.0), 5e-4, rel_tol=1e-9)
        assert math.isclose(regression_rate(5e-4, 0.0, 999.0), 5e-4, rel_tol=1e-9)

    def test_increases_with_G_o(self):
        r1 = regression_rate(1e-4, 0.5, 100.0)
        r2 = regression_rate(1e-4, 0.5, 400.0)
        assert r2 > r1

    def test_positive_for_positive_inputs(self):
        assert regression_rate(1e-4, 0.6, 250.0) > 0

    def test_raises_for_zero_or_negative_G_o(self):
        with pytest.raises(ValueError):
            regression_rate(1e-4, 0.5, 0.0)
        with pytest.raises(ValueError):
            regression_rate(1e-4, 0.5, -10.0)


# ---------------------------------------------------------------------------
# oxidizer_mass_flux
# ---------------------------------------------------------------------------

class TestOxidizerMassFlux:

    def test_known_value(self):
        # G_o = m_dot / (pi r^2) = 0.5 / (pi * 0.02^2) = 0.5 / 1.2566e-3 = 397.9
        G_o = oxidizer_mass_flux(0.5, 0.02)
        assert math.isclose(G_o, 0.5 / (math.pi * 0.02 ** 2), rel_tol=1e-9)

    def test_decreases_with_radius(self):
        G1 = oxidizer_mass_flux(0.5, 0.01)
        G2 = oxidizer_mass_flux(0.5, 0.03)
        assert G2 < G1

    def test_scales_linearly_with_mass_flow(self):
        G1 = oxidizer_mass_flux(0.5, 0.02)
        G2 = oxidizer_mass_flux(1.0, 0.02)
        assert math.isclose(G2, 2 * G1, rel_tol=1e-9)

    def test_raises_for_nonpositive_radius(self):
        with pytest.raises(ValueError):
            oxidizer_mass_flux(0.5, 0.0)
        with pytest.raises(ValueError):
            oxidizer_mass_flux(0.5, -0.01)


# ---------------------------------------------------------------------------
# fuel_mass_flow_rate_per_port
# ---------------------------------------------------------------------------

class TestFuelMassFlowRatePerPort:

    def test_known_value(self):
        # Hand-computed: G_o = 0.5/(pi*0.02^2), r_dot = a*G_o^n,
        # m_dot_fuel = rho*2*pi*r*L*r_dot
        a, n, rho, L, r, m_ox = 8e-5, 0.6, 900.0, 0.25, 0.02, 0.5
        G_o = oxidizer_mass_flux(m_ox, r)
        r_dot = regression_rate(a, n, G_o)
        expected = rho * 2 * math.pi * r * L * r_dot
        assert math.isclose(
            fuel_mass_flow_rate_per_port(a, n, rho, L, r, m_ox),
            expected, rel_tol=1e-9)

    def test_scales_linearly_with_length(self):
        kw = dict(a=8e-5, n=0.6, rho_fuel=900.0, r_port=0.02,
                  m_dot_ox_per_port=0.5)
        m1 = fuel_mass_flow_rate_per_port(L=0.2, **kw)
        m2 = fuel_mass_flow_rate_per_port(L=0.4, **kw)
        assert math.isclose(m2, 2 * m1, rel_tol=1e-9)

    def test_scales_linearly_with_density(self):
        kw = dict(a=8e-5, n=0.6, L=0.25, r_port=0.02,
                  m_dot_ox_per_port=0.5)
        m1 = fuel_mass_flow_rate_per_port(rho_fuel=800.0, **kw)
        m2 = fuel_mass_flow_rate_per_port(rho_fuel=1600.0, **kw)
        assert math.isclose(m2, 2 * m1, rel_tol=1e-9)

    def test_positive(self):
        assert fuel_mass_flow_rate_per_port(8e-5, 0.6, 900.0, 0.25, 0.02, 0.5) > 0


# ---------------------------------------------------------------------------
# solve_initial_port_radius
# ---------------------------------------------------------------------------

class TestSolveInitialPortRadius:

    # Shared realistic-ish parameters
    a, n, rho, L, m_ox_pp = 8.24e-5, 0.6, 900.0, 0.25, 0.5

    def test_round_trip_with_fuel_mass_flow_rate(self):
        # Solve for r_0 at a target, then re-evaluate the forward
        # function at that r_0 -- must recover the target.
        target = 0.08333  # kg/s
        r0 = solve_initial_port_radius(target, self.a, self.n, self.rho,
                                       self.L, self.m_ox_pp)
        recovered = fuel_mass_flow_rate_per_port(self.a, self.n, self.rho,
                                                  self.L, r0, self.m_ox_pp)
        assert math.isclose(recovered, target, rel_tol=1e-4)

    def test_matches_closed_form_for_n_not_half(self):
        # For n != 0.5: m_dot_fuel = K * r^(1-2n), K as derived in the
        # module docstring -- invert algebraically and compare.
        target = 0.08333
        K = (2 * math.pi * self.rho * self.L * self.a
             * math.pi ** (-self.n) * self.m_ox_pp ** self.n)
        r0_closed_form = (target / K) ** (1.0 / (1.0 - 2 * self.n))
        r0_bisection = solve_initial_port_radius(
            target, self.a, self.n, self.rho, self.L, self.m_ox_pp)
        assert math.isclose(r0_closed_form, r0_bisection, rel_tol=1e-4)

    def test_larger_target_requires_smaller_radius_when_n_above_half(self):
        # IMPORTANT PHYSICAL SUBTLETY (see solve_initial_port_radius's
        # docstring): m_dot_fuel(r) = K * r^(1-2n). For n > 0.5 (this
        # class's n=0.6), the exponent (1-2n) is NEGATIVE, so fuel flow
        # DECREASES as radius increases -- meaning a LARGER target flow
        # needs a SMALLER radius, the opposite of naive intuition.
        # Verified numerically before writing this test (an earlier,
        # wrong version of this test assumed the n<0.5 direction).
        r1 = solve_initial_port_radius(0.05, self.a, self.n, self.rho,
                                       self.L, self.m_ox_pp)
        r2 = solve_initial_port_radius(0.15, self.a, self.n, self.rho,
                                       self.L, self.m_ox_pp)
        assert r2 < r1

    def test_larger_target_requires_larger_radius_when_n_below_half(self):
        # Mirror check at n < 0.5: exponent (1-2n) is POSITIVE here, so
        # the naive "more flow needs a bigger hole" intuition holds.
        # Together with the test above, this pins down the full
        # direction-reversal behaviour around n=0.5 rather than assuming
        # either direction. Targets chosen well inside the achievable
        # range for these parameters (8.9-55.8 g/s across r in
        # [1e-4, 1] m) -- see test_raises_when_unbracketed... below for
        # what happens outside that range.
        n_low = 0.4
        r1 = solve_initial_port_radius(0.015, self.a, n_low, self.rho,
                                       self.L, self.m_ox_pp)
        r2 = solve_initial_port_radius(0.025, self.a, n_low, self.rho,
                                       self.L, self.m_ox_pp)
        assert r2 > r1

    def test_n_half_degenerate_case_unreachable_target_raises(self):
        # At n=0.5, fuel flow is independent of r -- a target that does
        # not match the unique achievable value has NO solution anywhere.
        n_half = 0.5
        achievable = fuel_mass_flow_rate_per_port(
            self.a, n_half, self.rho, self.L, 0.02, self.m_ox_pp)
        unreachable_target = achievable * 2.0  # deliberately different
        with pytest.raises(RuntimeError):
            solve_initial_port_radius(unreachable_target, self.a, n_half,
                                      self.rho, self.L, self.m_ox_pp)

    def test_raises_for_nonpositive_target(self):
        with pytest.raises(ValueError):
            solve_initial_port_radius(0.0, self.a, self.n, self.rho,
                                      self.L, self.m_ox_pp)
        with pytest.raises(ValueError):
            solve_initial_port_radius(-0.01, self.a, self.n, self.rho,
                                      self.L, self.m_ox_pp)

    def test_raises_when_unbracketed_outside_default_range(self):
        # A target far outside what [r_min, r_max] can reach must raise
        # RuntimeError with a clear message, not silently return a
        # meaningless endpoint.
        with pytest.raises(RuntimeError):
            solve_initial_port_radius(1e6, self.a, self.n, self.rho,
                                      self.L, self.m_ox_pp,
                                      r_min=1e-4, r_max=0.1)


# ---------------------------------------------------------------------------
# size_grain -- top-level sizing function
# ---------------------------------------------------------------------------

class TestSizeGrain:

    a, n, rho, L = 8.24e-5, 0.6, 900.0, 0.25

    def test_m_dot_fuel_matches_OF_definition(self):
        result = size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L)
        assert math.isclose(result["m_dot_fuel"], 0.5 / 6.0, rel_tol=1e-9)

    def test_higher_OF_gives_less_fuel_flow(self):
        r1 = size_grain(0.5, 4.0, self.a, self.n, self.rho, self.L)
        r2 = size_grain(0.5, 8.0, self.a, self.n, self.rho, self.L)
        assert r2["m_dot_fuel"] < r1["m_dot_fuel"]

    def test_multi_port_increases_total_perimeter_for_fixed_total_area(self):
        # Physical property highlighted in the module docstring: for a
        # FIXED total port area, splitting into N identical ports
        # increases total burning perimeter by a factor of sqrt(N).
        # This is a pure geometric fact about circles, independent of
        # the Marxman exponent n -- tested directly here (NOT through
        # solve_initial_port_radius, whose solved radius depends on n in
        # a way that does not hold total area fixed -- see
        # TestSolveInitialPortRadius for that n-dependent behaviour).
        A_total = 0.002  # m^2, arbitrary fixed total port area
        for N in (1, 4, 9):
            r_each = math.sqrt(A_total / (N * math.pi))
            total_perimeter = N * 2 * math.pi * r_each
            total_area = N * math.pi * r_each ** 2
            assert math.isclose(total_area, A_total, rel_tol=1e-9)
            expected_perimeter = 2 * math.sqrt(N * math.pi * A_total)
            assert math.isclose(total_perimeter, expected_perimeter, rel_tol=1e-9)
        # And the N=9 case must have 3x the perimeter of N=1 (sqrt(9)=3)
        r1 = math.sqrt(A_total / (1 * math.pi))
        r9 = math.sqrt(A_total / (9 * math.pi))
        perim1 = 1 * 2 * math.pi * r1
        perim9 = 9 * 2 * math.pi * r9
        assert math.isclose(perim9, 3 * perim1, rel_tol=1e-9)

    def test_more_ports_reduces_per_port_target_flow(self):
        # size_grain splits both m_dot_ox and the fuel-flow target
        # equally across N_ports -- this much is n-independent and
        # always true, unlike the resulting radius (see
        # TestSolveInitialPortRadius on why radius direction depends on
        # n and is therefore not asserted here).
        r1 = size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L,
                        N_ports=1)
        r4 = size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L,
                        N_ports=4)
        assert math.isclose(r4["m_dot_fuel_per_port"],
                            r1["m_dot_fuel_per_port"] / 4, rel_tol=1e-9)
        assert math.isclose(r4["m_dot_ox_per_port"],
                            r1["m_dot_ox_per_port"] / 4, rel_tol=1e-9)

    def test_result_has_required_keys(self):
        result = size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L)
        for key in ("m_dot_fuel", "m_dot_fuel_per_port", "m_dot_ox_per_port",
                    "r_0", "G_o_0", "r_dot_0", "N_ports", "L",
                    "r_final_estimate", "fuel_mass_consumed_estimate"):
            assert key in result

    def test_burn_time_none_gives_none_estimates(self):
        result = size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L,
                            burn_time=None)
        assert result["r_final_estimate"] is None
        assert result["fuel_mass_consumed_estimate"] is None

    def test_burn_time_gives_final_radius_larger_than_initial(self):
        result = size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L,
                            burn_time=8.0)
        assert result["r_final_estimate"] > result["r_0"]

    def test_burn_time_fuel_mass_consumed_positive(self):
        result = size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L,
                            burn_time=8.0)
        assert result["fuel_mass_consumed_estimate"] > 0

    def test_longer_burn_time_consumes_more_fuel(self):
        r_short = size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L,
                             burn_time=4.0)
        r_long = size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L,
                            burn_time=12.0)
        assert r_long["fuel_mass_consumed_estimate"] > r_short["fuel_mass_consumed_estimate"]

    def test_raises_for_nonpositive_m_dot_ox(self):
        with pytest.raises(ValueError):
            size_grain(0.0, 6.0, self.a, self.n, self.rho, self.L)

    def test_raises_for_nonpositive_OF(self):
        with pytest.raises(ValueError):
            size_grain(0.5, 0.0, self.a, self.n, self.rho, self.L)

    def test_raises_for_invalid_N_ports(self):
        with pytest.raises(ValueError):
            size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L, N_ports=0)

    def test_raises_for_nonpositive_burn_time(self):
        with pytest.raises(ValueError):
            size_grain(0.5, 6.0, self.a, self.n, self.rho, self.L,
                      burn_time=-1.0)


# ---------------------------------------------------------------------------
# FUEL_PROPERTIES -- structural checks (not a, n values, per module policy)
# ---------------------------------------------------------------------------

class TestFuelProperties:

    def test_expected_fuels_present(self):
        for fuel in ("Paraffin wax", "HTPB", "ABS", "PMMA"):
            assert fuel in FUEL_PROPERTIES

    def test_all_densities_positive_and_physically_plausible(self):
        # Sanity range for solid rocket-fuel-like polymers/waxes:
        # 700-1400 kg/m^3 comfortably covers all four fuels here.
        for fuel, props in FUEL_PROPERTIES.items():
            rho = props["rho_kg_m3"]
            assert 700.0 < rho < 1400.0, (
                f"{fuel} density {rho} kg/m^3 outside plausible range")

    def test_no_fuel_ships_a_or_n_as_a_usable_default(self):
        # Structural guard for the module's core policy: FUEL_PROPERTIES
        # must NOT contain a direct numeric "a" or "n" key that could be
        # mistaken for a usable default -- only descriptive reference
        # text under "a_n_reference_range".
        for fuel, props in FUEL_PROPERTIES.items():
            assert "a" not in props
            assert "n" not in props
            assert "a_n_reference_range" in props
