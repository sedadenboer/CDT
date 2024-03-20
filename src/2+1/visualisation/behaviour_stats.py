import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from classes.vertex import Vertex
from classes.universe import Universe
from classes.simulation import Simulation
import copy
import matplotlib as mpl
import pandas as pd
# mpl.rcParams['font.size'] = 16


def critical_k3(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3, chains: int = 10):
    # Set parameters
    k0_values = [0, 1, 2, 3, 4, 5, 6, 7]
    k3_init_guess = [0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
    all_data = []
    sweeps = 0
    thermal_sweeps = 100
    k_steps = 1000000
    target_volume = 10000
    target2_volume = 0

    for chain in range(chains):
        critical_k3_values = []

        # Run simulations
        for i, k0 in enumerate(k0_values):
            # Create the universe and simulation
            universe = Universe(geometry_infilename=geometry_infilename, strictness=strictness)
            simulation = Simulation(universe)
            simulation.saving_interval = 100
            print(f"Running k0={k0}, k3_init_guess={k3_init_guess[i]}, chain={chain}")
            simulation.start(
                k0=k0,
                k3=k3_init_guess[i],
                sweeps=sweeps,
                thermal_sweeps=thermal_sweeps,
                k_steps=k_steps,
                target_volume=target_volume,
                target2_volume=target2_volume,
                volfix_switch=1,
                seed=chain,
                outfile=f"k3/swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtvol={target_volume}_k0={k0}_chain={chain}",
                validity_check=True,
                v1=1,
                v2=1,
                v3=1
            )

            # Save critical k3 value
            critical_k3_values.append(simulation.k3)
            print(f"Critical k3 value: {simulation.k3}\n")
        
        all_data.append(critical_k3_values)

    # Save the critical k3 values
    np.savetxt(f'saved_universes/k3/measurements/critical_k3_values.txt', all_data, delimiter=',')

def simulation_worker(geometry_infilename, strictness, k0, k3_init_guess, sweeps, thermal_sweeps, k_steps, target_volume, target2_volume, chain):
    universe = Universe(geometry_infilename=geometry_infilename, strictness=strictness)
    simulation = Simulation(universe)
    simulation.saving_interval = 100
    print(f"Running k0={k0}, k3_init_guess={k3_init_guess}, chain={chain}")
    simulation.start(
        k0=k0,
        k3=k3_init_guess,
        sweeps=sweeps,
        thermal_sweeps=thermal_sweeps,
        k_steps=k_steps,
        target_volume=target_volume,
        target2_volume=target2_volume,
        volfix_switch=1,
        seed=chain,
        outfile=f"k3/swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtvol={target_volume}_k0={k0}_chain={chain}",
        validity_check=True,
        v1=1,
        v2=1,
        v3=1
    )
    return simulation.k3

def critical_k3_parallel(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3, chains: int = 10):
    # Set parameters
    k0_values = [0, 1, 2, 3, 4, 5, 6, 7]
    k3_init_guess = [0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
    all_data = []
    sweeps = 0
    thermal_sweeps = 100
    k_steps = 1000000
    target_volume = 10000
    target2_volume = 0

    for chain in range(chains):
        with multiprocessing.Pool() as pool:
            # Prepare arguments
            items = [(geometry_infilename, strictness, k0_values[i], k3_init_guess[i], sweeps, thermal_sweeps, k_steps, target_volume, target2_volume, chain) for i, k0 in enumerate(k0_values)]
            critical_k3_values = pool.starmap(simulation_worker, items)
            all_data.append(critical_k3_values)

    # Save the critical k3 values as a dataframe
    df = pd.DataFrame(all_data, columns=k0_values)
    df.index.name = 'chain'
    df.columns.name = 'k0'
    print(df)
    df.to_csv('saved_universes/k3/measurements/critical_k3_values.csv', index=True)

if __name__ == "__main__":
    critical_k3_parallel()
    print("Done")

