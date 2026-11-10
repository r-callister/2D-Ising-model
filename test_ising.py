"""
Tests for ising.py. Run with: pytest
"""
import numpy as np
import pytest

import ising
import ising_scan


def test_aligned_lattice_energy():
    """All spins aligned: every one of the 2N bonds contributes -J."""
    lattice = np.ones((8, 8), dtype=int)
    assert ising.total_energy(lattice) == -2 * 64 * ising.J


def test_energy_change_matches_full_recalculation():
    """energy_change must equal the difference in total energy after the flip."""
    rng = np.random.default_rng(1)
    lattice = ising.initial_lattice(6, rng)
    for i, j in [(0, 0), (2, 3), (5, 5), (0, 5)]:     # includes spins on the edges
        before = ising.total_energy(lattice)
        predicted = ising.energy_change(lattice, i, j)
        lattice[i, j] = -lattice[i, j]
        assert ising.total_energy(lattice) - before == predicted


def test_low_temperature_stays_ordered():
    """At T = 0.5 an aligned lattice stays almost fully aligned and near -2J per spin."""
    rng = np.random.default_rng(2)
    energies, magnetisations, _ = ising.simulate(16, 0.5, 100, rng, ordered=True)
    assert energies[-1] < -1.95
    assert magnetisations[-1] > 0.98


def test_high_temperature_is_disordered():
    """At T = 10 a lattice has almost no net magnetisation."""
    rng = np.random.default_rng(3)
    _, magnetisations, _ = ising.simulate(16, 10.0, 200, rng)
    assert np.mean(magnetisations[50:]) < 0.2


def test_energy_matches_onsager_at_t3():
    """Average energy at T = 3 agrees with Onsager's exact value of -0.8173."""
    rng = np.random.default_rng(4)
    energies, magnetisations, _ = ising.simulate(16, 3.0, 4000, rng)
    averages, _ = ising.calculate_averages(energies, magnetisations, 16, 3.0, 500)
    value, error = averages["energy"]
    assert abs(value + 0.8173) < 4 * error + 0.005


def test_drift_is_detected():
    """A steadily rising magnetisation is flagged, a steady one is not."""
    rng = np.random.default_rng(5)
    energies = np.full(2001, -1.7)
    steady = 0.9 + 0.03 * rng.standard_normal(2001)
    rising = np.linspace(0.2, 0.9, 2001) + 0.03 * rng.standard_normal(2001)
    _, drifted_steady = ising.calculate_averages(energies, steady, 16, 2.0, 0)
    _, drifted_rising = ising.calculate_averages(energies, rising, 16, 2.0, 0)
    assert not drifted_steady
    assert drifted_rising


def test_jackknife_error_for_independent_data():
    """For independent samples the error of the mean should be close to sigma / sqrt(n)."""
    rng = np.random.default_rng(6)
    n = 5000
    energies = rng.standard_normal(n + 1)     # sigma = 1; entry 0 is skipped as the start
    averages, _ = ising.calculate_averages(energies, np.zeros(n + 1), 8, 2.0, 0)
    expected = 1 / np.sqrt(n)
    assert 0.6 * expected < averages["energy"][1] < 1.4 * expected


def test_onsager_magnetisation():
    """Onsager's magnetisation is 0 above T_c and has known values below it."""
    m = ising_scan.exact_magnetisation([1.0, 2.0, 2.5, 3.0])
    assert m[0] == pytest.approx(0.99928, abs=1e-5)
    assert m[1] == pytest.approx(0.91132, abs=1e-5)
    assert m[2] == 0 and m[3] == 0


def test_too_few_sweeps_raises():
    """Too few sweeps after equilibration is reported as an error."""
    with pytest.raises(ValueError):
        ising.calculate_averages(np.zeros(10), np.zeros(10), 8, 2.0, 4)
