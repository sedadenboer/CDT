import sys
sys.path.append('..')
from typing import List
from classes.universe import Universe
from classes.simulation import Simulation
import multiprocessing as mp
import time
import os


def run_simulation(universe: Universe, chain: int, k0: float, k3: float):
    print(f'Running simulation with k0={k0} and k3={k3}')

    simulation = Simulation(
        universe=universe,
        seed=chain,
        k0=k0,
        k3=k3,
        tune_flag=True,
        thermal_sweeps=1000,
        sweeps=0,
        k_steps=300000,
        target_volume=3000, # Without tune does not do anything
        observables=['n_vertices', 'n_tetras', 'n_tetras_31', 'n_tetras_22', 'slice_sizes', 'slab_sizes', 'curvature'],
        include_mcmc_data=True,
        measuring_interval=1, # Measure every sweep
        measuring_thermal=True,
        save_thermal=True,
        saving_interval=2, # When to save geometry files
    )

    simulation.start(
        outfile=f'outfile_k0={simulation.k0}_tswps={simulation.thermal_sweeps}_swps={simulation.sweeps}_kstps={simulation.k_steps}_chain={chain}'
    )

    return simulation.observables

def run_parallel_simulations(universe: Universe, chain: int):
    k0_values: List[int] = [0, 1, 2, 3, 4, 5, 6, 7]
    k3_init_guesses: List[float] = [0.7, 0.9, 1.1, 1.3, 1.5, 1.7 ,1.9, 2.1]
    args = [(universe, chain, k0_values[i], k3_init_guesses[i]) for i in range(len(k0_values))]

    with mp.Pool(int(os.environ['SLURM_CPUS_PER_TASK'])) as pool:
        results = pool.starmap(run_simulation, args)

    return results

def run_sequential(universe: Universe, chain: int):
    k0_values: List[int] = [3, 4, 5, 6, 7]
    k3_init_guesses: List[float] = [1.3, 1.5, 1.7, 1.9, 2.1]

    for i in range(len(k0_values)):
        run_simulation(universe, chain, k0_values[i], k3_init_guesses[i])


if __name__ == "__main__":
    universe = Universe(geometry_infilename='classes/initial_universes/sample-g0-T3.cdt', strictness=3)
    run_parallel_simulations(universe=universe, chain=0)