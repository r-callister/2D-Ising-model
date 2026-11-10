"""
Temperature scan of the 2D Ising model

Runs the simulation from ising.py at a range of temperatures for several lattice sizes,
and plots the |magnetisation| and the susceptibility against temperature. The
magnetisation is compared with Onsager's exact result for an infinite lattice, and the
susceptibility peak is used to estimate the critical temperature.

The random seed is fixed, so running this again reproduces the same figure.
The full scan takes about 15 minutes.
"""
import numpy as np
import matplotlib.pyplot as plt

from ising import simulate, calculate_averages, J

SIZES = [8, 16, 32]
# Coarse steps far from the transition, fine steps of 0.05 around it
TEMPERATURES = np.round(np.concatenate([np.linspace(1.0, 1.75, 4),
                                        np.linspace(2.0, 2.75, 16),
                                        np.linspace(3.0, 4.0, 5)]), 2)
N_EQUILIBRATION = 500
N_MEASUREMENTS = 6000
SEED = 1
T_CRITICAL = 2 * J / np.log(1 + np.sqrt(2))      # Onsager's exact critical temperature
COLOURS = ["#2a78d6", "#eb6834", "#e34948"]      # one colour per lattice size


def exact_magnetisation(temperatures):
    """
    Onsager's exact |magnetisation| per spin of an infinite lattice,
        M = (1 - sinh(2J/T)^-4)^(1/8)   for T < T_c,   and 0 for T > T_c.
    """
    temperatures = np.asarray(temperatures, dtype=float)
    magnetisation = np.zeros_like(temperatures)
    below = temperatures < T_CRITICAL
    magnetisation[below] = (1 - np.sinh(2 * J / temperatures[below])**-4)**0.125
    return magnetisation


def scan_size(size, rng):
    """
    Runs every temperature for one lattice size. Every run starts with all spins aligned,
    which avoids the long-lived stripe states that random starts can fall into below T_c.
    Returns the |magnetisation| and susceptibility with their errors at each temperature,
    and the temperatures at which the run drifted (and so should not be trusted).
    """
    values = {"magnetisation": [], "susceptibility": []}
    errors = {"magnetisation": [], "susceptibility": []}
    drifted_at = []
    for temperature in TEMPERATURES:
        energies, magnetisations, _ = simulate(
            size, temperature, N_EQUILIBRATION + N_MEASUREMENTS, rng, ordered=True)
        averages, drifted = calculate_averages(
            energies, magnetisations, size, temperature, N_EQUILIBRATION)
        for name, series in values.items():
            series.append(averages[name][0])
            errors[name].append(averages[name][1])
        if drifted:
            drifted_at.append(temperature)
    return values, errors, drifted_at


def plot_scan(results):
    """
    Plots |magnetisation| and susceptibility against temperature for every lattice size.
    """
    fig, (ax_m, ax_chi) = plt.subplots(1, 2, figsize=(12, 5))
    smooth_temperatures = np.linspace(1.0, 4.0, 400)
    ax_m.plot(smooth_temperatures, exact_magnetisation(smooth_temperatures),
              color="black", linestyle="--", label="Onsager (infinite lattice)")
    for (size, (values, errors)), colour in zip(results.items(), COLOURS):
        ax_m.errorbar(TEMPERATURES, values["magnetisation"], yerr=errors["magnetisation"],
                      fmt="o", markersize=4, capsize=2, color=colour, label=f"L = {size}")
        ax_chi.errorbar(TEMPERATURES, values["susceptibility"], yerr=errors["susceptibility"],
                        fmt="o", markersize=4, capsize=2, color=colour, label=f"L = {size}")
    for ax in (ax_m, ax_chi):
        ax.axvline(T_CRITICAL, color="grey", linestyle=":", label=f"Exact $T_c$ = {T_CRITICAL:.3f}")
        ax.set_xlabel("Temperature (J / $k_B$)")
        ax.legend()
    ax_m.set_ylabel("|Magnetisation| per spin")
    ax_m.set_title("Magnetisation is lost at the critical temperature")
    ax_chi.set_ylabel("Susceptibility per spin")
    ax_chi.set_title("The susceptibility peaks at the transition")
    fig.tight_layout()
    fig.savefig("ising_scan.png", dpi=200)
    plt.show()


def main():
    """
    Scans all lattice sizes, prints where the susceptibility peaks for each, and plots.
    """
    rng = np.random.default_rng(SEED)
    results = {}
    for size in SIZES:
        print(f"Scanning L = {size} ...")
        values, errors, drifted_at = scan_size(size, rng)
        results[size] = (values, errors)
        peak = TEMPERATURES[np.argmax(values["susceptibility"])]
        print(f"  susceptibility peaks at T = {peak:.2f} (exact T_c = {T_CRITICAL:.3f})")
        if drifted_at:
            print(f"  warning: drift detected at T = {drifted_at}")
    plot_scan(results)


if __name__ == "__main__":
    main()
