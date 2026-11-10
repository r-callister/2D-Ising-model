"""
2D Ising Model

Simulates a square lattice of spins (each +1 or -1) with the Metropolis Monte Carlo
algorithm. The user chooses a lattice size, a temperature and the number of sweeps. The
program shows how the lattice evolves towards equilibrium at that temperature, then
averages over the equilibrium part of the run to find the energy, magnetisation and
susceptibility.

Units: the coupling J and the Boltzmann constant k_B are both set to 1,
so energies are in units of J and temperatures are in units of J / k_B.
"""
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

J = 1.0     # coupling constant: J > 0 means neighbouring spins prefer to point the same way
N_BLOCKS = 10   # number of blocks the measurement sweeps are split into for the errors


# --- The lattice and its energy ---

def initial_lattice(size, rng, ordered=False):
    """
    Returns a size x size array of spins, each +1 or -1.
    ordered=True gives every spin +1 (all aligned, as at very low temperature).
    ordered=False gives random spins (as at very high temperature).
    """
    if ordered:
        return np.ones((size, size), dtype=int)
    return rng.choice([-1, 1], size=(size, size))


def total_energy(lattice):
    """
    Total energy E = -J * (sum over neighbouring pairs of s_i * s_j).

    Each spin has 4 neighbours (up, down, left, right), but each pair of neighbours is
    one bond, so we count only the right and down neighbours of every spin. That counts
    every bond exactly once.

    np.roll makes a shifted copy of the array, so that position [i, j] holds the
    spin that was at [i, j+1]. Spins on the edge wrap around to the other side
    (periodic boundary conditions), so there are no special edge spins and the lattice
    behaves like the surface of a torus.
    """
    right_neighbours = np.roll(lattice, -1, axis=1)
    down_neighbours = np.roll(lattice, -1, axis=0)
    return -J * np.sum(lattice * right_neighbours + lattice * down_neighbours)


def magnetisation_per_spin(lattice):
    """
    Average spin, M = (1/N) * sum of all spins. It is +1 when all spins point up,
    -1 when all point down, and close to 0 when the spins are random.
    """
    return np.mean(lattice)


# --- The Metropolis algorithm ---

def energy_change(lattice, i, j):
    """
    The change in total energy if the spin at row i, column j were flipped.
    Only the 4 neighbours of that spin share a bond with it, so
        delta_E = 2 * J * s_ij * (sum of the 4 neighbouring spins).
    The % size makes the indices wrap around the edges (periodic boundaries).
    """
    size = lattice.shape[0]
    neighbours = (lattice[(i + 1) % size, j] + lattice[(i - 1) % size, j]
                  + lattice[i, (j + 1) % size] + lattice[i, (j - 1) % size])
    return 2 * J * lattice[i, j] * neighbours


def metropolis_sweep(lattice, temperature, rng):
    """
    One sweep = N attempted spin flips (N = number of spins), at randomly chosen spins.
    A flip is accepted if it lowers the energy (delta_E <= 0), or otherwise with
    probability exp(-delta_E / T). This changes the lattice in place.
    """
    size = lattice.shape[0]
    n_spins = size * size
    rows = rng.integers(0, size, n_spins)
    columns = rng.integers(0, size, n_spins)
    randoms = rng.random(n_spins)
    for k in range(n_spins):
        i, j = rows[k], columns[k]
        delta_e = energy_change(lattice, i, j)
        if delta_e <= 0 or randoms[k] < math.exp(-delta_e / temperature):
            lattice[i, j] = -lattice[i, j]


# --- Running a simulation ---

def simulate(size, temperature, n_sweeps, rng, ordered=False):
    """
    Starts from a random (or all aligned) lattice and runs n_sweeps Metropolis sweeps.
    Returns the energy per spin and |magnetisation| after every sweep (entry 0 is the
    starting lattice), and copies of the lattice at three moments: the start, a tenth of
    the way through, and the end.
    """
    lattice = initial_lattice(size, rng, ordered)
    n_spins = size * size
    energies = np.zeros(n_sweeps + 1)
    magnetisations = np.zeros(n_sweeps + 1)
    energies[0] = total_energy(lattice) / n_spins
    magnetisations[0] = abs(magnetisation_per_spin(lattice))

    snapshot_sweeps = [0, max(1, n_sweeps // 10), n_sweeps]
    snapshots = {0: lattice.copy()}
    for sweep in range(1, n_sweeps + 1):
        metropolis_sweep(lattice, temperature, rng)
        energies[sweep] = total_energy(lattice) / n_spins
        magnetisations[sweep] = abs(magnetisation_per_spin(lattice))
        if sweep in snapshot_sweeps:
            snapshots[sweep] = lattice.copy()
    return energies, magnetisations, snapshots


def observables(e, m, n_spins, temperature):
    """
    The three averaged quantities from the energy per spin e and |magnetisation| m
    recorded at each sweep:
        energy          <e>
        magnetisation   <m>
        susceptibility  chi = N (<m^2> - <m>^2) / T          (per spin)
    where N is the number of spins and <...> is the average over the sweeps given.

    m here is |magnetisation|, not the signed value. On a finite lattice the whole lattice
    occasionally flips between mostly up and mostly down, so the signed <m> averages to
    zero even below T_c, and N <m^2> / T would grow with N instead of showing a peak.
    Using |m| is the standard choice for finite lattices.
    """
    return np.array([np.mean(e), np.mean(m), n_spins * np.var(m) / temperature])


def jackknife_errors(e_blocks, m_blocks, n_spins, temperature):
    """
    Standard errors of the three observables, using the blocked jackknife.
    Consecutive sweeps are correlated, so the sweeps cannot be treated as independent
    samples. Blocks much longer than the correlation time are close to independent.
    The observables are recalculated with each block left out in turn, and
        error^2 = (n - 1) / n * sum over k of (estimate_k - mean of the estimates)^2
    for n blocks, which works for the variance-based quantities as well as for means.
    """
    n_blocks = len(e_blocks)
    estimates = []
    for k in range(n_blocks):
        keep = np.arange(n_blocks) != k
        estimates.append(observables(e_blocks[keep].ravel(), m_blocks[keep].ravel(),
                                     n_spins, temperature))
    estimates = np.array(estimates)
    return np.sqrt((n_blocks - 1) / n_blocks * np.sum((estimates - estimates.mean(axis=0))**2,
                                                       axis=0))


def has_drifted(blocks, threshold=3.0):
    """
    True if the block means in the first half of the measurement window differ from those
    in the second half by more than threshold standard errors of that difference.
    A system in equilibrium has the same mean in both halves, so a clear difference means
    the run was still settling and its averages should not be trusted.
    """
    means = blocks.mean(axis=1)
    half = len(means) // 2
    first, second = means[:half], means[half:2 * half]
    spread = np.sqrt(np.var(first, ddof=1) / half + np.var(second, ddof=1) / half)
    return abs(first.mean() - second.mean()) > threshold * spread


def split_into_blocks(data):
    """
    Splits data into N_BLOCKS equal consecutive blocks, as a 2D array with one block per
    row. Any leftover sweeps at the end are dropped.
    """
    block_length = len(data) // N_BLOCKS
    return data[:block_length * N_BLOCKS].reshape(N_BLOCKS, block_length)


def calculate_averages(energies, magnetisations, size, temperature, n_equilibration):
    """
    Averages over the equilibrium part of a run. The first n_equilibration sweeps are
    discarded, because the lattice is still settling and is not yet typical of
    this temperature. Entry 0 of each array is the starting lattice, so the sweeps
    kept are entries n_equilibration + 1 onwards.

    Returns a dictionary mapping each of "energy", "magnetisation" and "susceptibility"
    to (value, error), and a flag that is True if the energy or the magnetisation
    drifted between the two halves of the measurement window.
    """
    if len(energies) - n_equilibration - 1 < 2 * N_BLOCKS:
        raise ValueError("Not enough sweeps left after equilibration to form blocks.")
    e_blocks = split_into_blocks(energies[n_equilibration + 1:])
    m_blocks = split_into_blocks(magnetisations[n_equilibration + 1:])

    n_spins = size * size
    values = observables(e_blocks.ravel(), m_blocks.ravel(), n_spins, temperature)
    errors = jackknife_errors(e_blocks, m_blocks, n_spins, temperature)
    names = ["energy", "magnetisation", "susceptibility"]
    results = dict(zip(names, zip(values, errors)))
    return results, has_drifted(e_blocks) or has_drifted(m_blocks)


# --- Input ---

def get_number(prompt, default, minimum, number_type):
    """
    Asks the user for a number and repeats until a valid one is given.
    Pressing Enter on its own uses the default.
    """
    while True:
        text = input(f"{prompt} [default {default}]: ").strip()
        if text == "":
            return default
        try:
            value = number_type(text)
        except ValueError:
            print("Enter a number.")
            continue
        if value < minimum:
            print(f"Enter a value of at least {minimum}.")
            continue
        return value


def get_start():
    """
    Asks whether to start from a random lattice or an all-aligned one.
    """
    while True:
        text = input("Starting lattice: (r)andom or (a)ligned [default r]: ").strip().lower()
        if text in ("", "r"):
            return False
        if text == "a":
            return True
        print("Enter r or a.")


# --- Plotting ---

def plot_run(snapshots, energies, magnetisations, temperature, n_equilibration):
    """
    Plots the lattice at the start, a tenth of the way through and the end,
    and how the energy per spin and |magnetisation| change with the number of sweeps.
    The dashed line marks the end of the equilibration sweeps that are discarded.
    """
    colours = ListedColormap(["tab:orange", "tab:blue"])   # spin -1 is orange, +1 is blue
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.8], hspace=0.3, wspace=0.3)

    for k, sweep in enumerate(sorted(snapshots)):
        ax = fig.add_subplot(gs[0, k])
        ax.imshow(snapshots[sweep], cmap=colours, vmin=-1, vmax=1, interpolation="nearest")
        ax.set_title(f"After {sweep} sweeps")
        ax.set_xticks([])
        ax.set_yticks([])

    sweeps = np.arange(len(energies))
    ax_energy = fig.add_subplot(gs[1, 0:2])
    ax_energy.plot(sweeps, energies, color="tab:blue")
    ax_energy.set_ylabel("Energy per spin (J)")
    ax_energy.set_title(f"Energy during the run at T = {temperature:g}")

    ax_magnetisation = fig.add_subplot(gs[1, 2])
    ax_magnetisation.plot(sweeps, magnetisations, color="tab:orange")
    ax_magnetisation.set_ylabel("|Magnetisation per spin|")
    ax_magnetisation.set_ylim(0, 1.05)
    ax_magnetisation.set_title("Magnetisation")

    for ax in (ax_energy, ax_magnetisation):
        ax.axvline(n_equilibration, color="black", linestyle="--", linewidth=1,
                   label="End of equilibration")
        ax.set_xlabel("Sweep")
    ax_energy.legend()

    fig.savefig("ising_run.png", dpi=200, bbox_inches="tight")
    plt.show()


def print_results(averages, drifted):
    """
    Prints each observable with its error, and a warning if the run drifted.
    """
    labels = {"energy": ("Average energy per spin", "J"),
              "magnetisation": ("Average |magnetisation| per spin", ""),
              "susceptibility": ("Susceptibility per spin", "")}
    for name, (label, unit) in labels.items():
        value, error = averages[name]
        print(f"{label}: {value:.3f} +/- {error:.3f} {unit}")
    if drifted:
        print("WARNING: the energy or magnetisation drifted during the measurement sweeps, "
              "so the lattice was not in equilibrium. Increase the equilibration sweeps.")


def main():
    """
    Asks for the simulation settings, runs the simulation and shows the results.
    """
    size = get_number("Lattice size (L, giving an L x L lattice)", 32, 2, int)
    temperature = get_number("Temperature (units of J/k_B; the transition is near 2.27)",
                             2.0, 0.1, float)
    n_equilibration = get_number("Equilibration sweeps (discarded)", 1000, 0, int)
    n_measurements = get_number("Measurement sweeps (averaged)", 2000, 2 * N_BLOCKS, int)
    seed = get_number("Random seed (0 gives a different run each time)", 0, 0, int)
    ordered = get_start()

    rng = np.random.default_rng(seed if seed != 0 else None)
    energies, magnetisations, snapshots = simulate(
        size, temperature, n_equilibration + n_measurements, rng, ordered)
    averages, drifted = calculate_averages(energies, magnetisations, size, temperature,
                                           n_equilibration)

    print_results(averages, drifted)
    plot_run(snapshots, energies, magnetisations, temperature, n_equilibration)


if __name__ == "__main__":
    main()
