import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
from classes.vertex import Vertex
from classes.universe import Universe
from classes.simulation import Simulation
import copy
import matplotlib as mpl
# mpl.rcParams['font.size'] = 16


def critical_k3(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3, chains: int = 10):
    # Create the universe and simulation
    universe = Universe(geometry_infilename=geometry_infilename, strictness=strictness)
    simulation = Simulation(universe)

    # Set parameters
    k0_values = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    k3_init_guess = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    all_critical_k3_values = []
    sweeps = 0
    thermal_sweeps = 10
    k_steps = 100000
    target_volume = 10000
    target2_volume = 0

    # Run simulations
    for i, k0 in enumerate(k0_values):
        critical_k3_values = []

        for j in range(chains):
            print(f"Running k0={k0}, k3_init_guess={k3_init_guess[i]}, chain={j}")
            simulation.start(
                k0=k0,
                k3=k3_init_guess[i],
                sweeps=sweeps,
                thermal_sweeps=thermal_sweeps,
                k_steps=k_steps,
                target_volume=target_volume,
                target2_volume=target2_volume,
                volfix_switch=1,
                seed=1,
                outfile=f"k3/swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtvol={target_volume}_k0={k0}_chain={j}",
                validity_check=True,
                v1=1,
                v2=1,
                v3=1
            )
    
            # Save critical k3 value
            critical_k3_values.append(simulation.k3)
            print(f"Critical k3 value: {simulation.k3}\n")

        all_critical_k3_values.append(critical_k3_values)

    # Save the critical k3 values
    np.savetxt(f'k3/critical_k3_values.txt', all_critical_k3_values)

def plot_volumeratio_k0(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3):
    # Create the universe and simulation
    universe = Universe(geometry_infilename=geometry_infilename, strictness=strictness)
    simulation = Simulation(universe)
    
    # Set parameters
    k0_values = np.linspace(5, 7, 10)
    k3 = 1.703
    sweeps = 1000
    thermal_sweeps = 1000
    k_steps = 1000000
    target_volume = 0
    target2_volume = 0

    # Run simulations
    ratios = []
    for k0 in k0_values:
        simulation.start(
            k0=k0,
            k3_i=k3,
            sweeps=sweeps,
            thermal_sweeps=thermal_sweeps,
            k_steps=k_steps,
            target_volume=target_volume,
            target2_volume=target2_volume,
            volfix_switch=0,
            seed=1,
            out_file=f"N/swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_k0={k0}",
            v1=1,
            v2=1,
            v3=1
        )

        # Save O2 value
        # ratios.append(simulation.O2_values[-1])



if __name__ == "__main__":
    critical_k3()
    print("Done")

