# rejection.py
#
# Author: Seda den Boer
# Date: 04-02-2024
#
# Description:

import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
from classes.universe import Universe
from classes.simulation import Simulation
import copy


def plot_acceptance(chains: int, steps: int, lambd: int, total_time: int, initial_slice_size: int) -> float:
    """
    Determine the acceptance rate of the MCMC algorithm over a number of
    iterations.
    """
    x = np.arange(1, steps + 1, 1)
    all_acceptance_rates = []
    all_add_rates = []
    all_delete_rates = []
    all_flip_rates = []

    for _ in range(chains):
        universe = Universe(total_time=total_time, initial_slice_size=initial_slice_size)
        simulation = Simulation(universe, lambd=lambd, target_volume=0, epsilon=0)
        simulation.progress_universe(steps=steps, silence=False)
        all_acceptance_rates.append(simulation.acceptance_rates)
        all_add_rates.append(simulation.add_rates)
        all_delete_rates.append(simulation.delete_rates)
        all_flip_rates.append(simulation.flip_rates)

    mean_acceptance_rates, std_acceptance_rates = np.mean(np.array(all_acceptance_rates), axis=0), np.std(np.array(all_acceptance_rates), axis=0)
    mean_add_rates, std_add_rates = np.mean(np.array(all_add_rates), axis=0), np.std(np.array(all_add_rates), axis=0)
    mean_delete_rates, std_delete_rates = np.mean(np.array(all_delete_rates), axis=0), np.std(np.array(all_delete_rates), axis=0)
    mean_flip_rates, std_flip_rates = np.mean(np.array(all_flip_rates), axis=0), np.std(np.array(all_flip_rates), axis=0)

    plt.figure()
    plt.plot(x, mean_acceptance_rates)
    plt.fill_between(x, mean_acceptance_rates - std_acceptance_rates, mean_acceptance_rates + std_acceptance_rates, alpha=0.3)
    plt.xlim(0, steps)
    plt.title(f"Total acceptance rate of MCMC algorithm ($\lambda$={lambd:.2f}, $N_{{v,i}}$={total_time * initial_slice_size})")
    plt.xlabel("Iteration")
    plt.ylabel("Acceptance rate")
    plt.grid(True)
    plt.savefig(f"plots/acceptance_rate_lambda={lambd:.2f}_steps={steps}_size={total_time * initial_slice_size}.png", dpi=400, bbox_inches='tight')
    
    plt.figure()
    plt.plot(x, mean_add_rates, label='Add Rate')
    plt.plot(x, mean_delete_rates, label='Delete Rate')
    plt.plot(x, mean_flip_rates, label='Flip Rate')
    plt.fill_between(x, mean_add_rates - std_add_rates, mean_add_rates + std_add_rates, alpha=0.3)
    plt.fill_between(x, mean_delete_rates - std_delete_rates, mean_delete_rates + std_delete_rates, alpha=0.3)
    plt.fill_between(x, mean_flip_rates - std_flip_rates, mean_flip_rates + std_flip_rates, alpha=0.3)
    plt.xlim(0, steps)
    plt.title(f"Acceptance rates of MCMC algorithm ($\lambda$={lambd:.2f}, $N_{{v,i}}$={total_time * initial_slice_size})")
    plt.xlabel("Iteration")
    plt.ylabel("Acceptance Rate")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"plots/sep_acceptance_rates_lambda={lambd:.2f}_steps={steps}_size={total_time * initial_slice_size}.png", dpi=400, bbox_inches='tight')

def plot_size_ratio_lambda(chains: int, thermalisation_time: int, steps: int, total_time: int, initial_slice_size: int):
    """
    Plot the ratio of the final size of the universe to the initial size of the universe for different values of lambda.
    """
    lambdas = np.linspace(0.68, 0.78, 100)
    all_size_ratios = []

    for lambd in lambdas:
        print(f"\nLAMBDA: {lambd}")
        universe = Universe(total_time=total_time, initial_slice_size=initial_slice_size)
        simulation = Simulation(universe, lambd=lambd, target_volume=0, epsilon=0)
        simulation.progress_universe(steps=thermalisation_time, silence=False)
        size_ratios = []
        for _ in range(chains):
            simulation_copy = copy.deepcopy(simulation)
            simulation_copy.progress_universe(steps=steps, silence=False)
            size_ratios.append(simulation_copy.universe.get_total_size() / (total_time * initial_slice_size))

        all_size_ratios.append(size_ratios)

    size_ratios_matrix = np.array(all_size_ratios)
    mean_size_ratios = np.mean(size_ratios_matrix, axis=1)
    std_size_ratios = np.std(size_ratios_matrix, axis=1)

    plt.figure()
    plt.axvline(x=np.log(2), color='r', linestyle='--', label="$\lambda$=ln(2)")
    plt.plot(lambdas, mean_size_ratios, label="Mean size ratio")
    plt.fill_between(lambdas, mean_size_ratios - std_size_ratios, mean_size_ratios + std_size_ratios, alpha=0.3)
    plt.title("Ratio of final to initial universe size as function of $\lambda$")
    plt.xlabel("$\lambda$")
    plt.ylabel("$N_{{v,f}} / N_{{v,i}}$")
    plt.legend()
    plt.grid(True)
    plt.savefig("plots/size_ratio_lambda.png", dpi=400, bbox_inches='tight')
    plt.show()

def total_volume_change(chains: int, steps: int, lambd: int, total_time: int, initial_slice_size: int):
    """
    Plot the total volume change as a function over iterations.
    """
    x = np.arange(1, steps + 1, 1)
    all_volume_changes = []
    for _ in range(chains):
        universe = Universe(total_time=total_time, initial_slice_size=initial_slice_size)
        simulation = Simulation(universe, lambd=lambd, target_volume=0, epsilon=0)
        simulation.progress_universe(steps=steps, silence=False)
        all_volume_changes.append(simulation.volume_changes)

    all_volume_changes_array = np.array(all_volume_changes)
    mean_volume_changes = np.mean(all_volume_changes_array, axis=0)
    std_volume_changes = np.std(all_volume_changes_array, axis=0)

    plt.figure()
    plt.plot(x, mean_volume_changes)
    plt.fill_between(x, mean_volume_changes - std_volume_changes, mean_volume_changes + std_volume_changes, alpha=0.3)
    plt.xlim(0, steps)
    plt.title(f"Total volume of MCMC algorithm ($\lambda$={lambd:.2f}, $N_{{v,i}}$={total_time * initial_slice_size})")
    plt.xlabel("Iteration")
    plt.ylabel("Total volume")
    plt.grid(True)
    plt.savefig(f"plots/total_volume_lambda={lambd:.2f}_steps={steps}_size={total_time * initial_slice_size}.png", dpi=400, bbox_inches='tight')

if __name__ == "__main__":
    total_time = 50
    initial_slice_size = 40
    steps = 1000000
    thermalisation_time = 200000
    chains = 10

    # # TOTAL VOLUME CHANGE OVER ITERATIONS
    # total_volume_change(chains=chains, steps=steps, lambd=np.log(2), total_time=total_time, initial_slice_size=initial_slice_size)
    # total_volume_change(chains=chains, steps=steps, lambd=0.8, total_time=total_time, initial_slice_size=initial_slice_size)
    # total_volume_change(chains=chains, steps=steps, lambd=0.5, total_time=total_time, initial_slice_size=initial_slice_size)

    # # SIZE RATIO & LAMBDA
    # steps = 10000
    # thermalisation_time = 200000
    # plot_size_ratio_lambda(chains=chains,
    #                        thermalisation_time=thermalisation_time,
    #                        steps=steps,
    #                        total_time=total_time,
    #                        initial_slice_size=initial_slice_size
    #                     )

    # # ACCEPTANCE RATE
    # steps = 1000000
    # plot_acceptance(chains=chains, steps=steps, lambd=np.log(2), total_time=total_time, initial_slice_size=initial_slice_size)
    # plot_acceptance(chains=chains, steps=steps, lambd=0.5, total_time=total_time, initial_slice_size=initial_slice_size)
    # plot_acceptance(chains=chains, steps=steps, lambd=0.8, total_time=total_time, initial_slice_size=initial_slice_size)

    # plt.show()
