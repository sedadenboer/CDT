import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
from classes.vertex import Vertex
from classes.universe import Universe
from classes.simulation import Simulation
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import datetime
mpl.rcParams['font.size'] = 16


def critical_k3(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3, chains: int = 10):
    """
    Runs the critical k3 simulations sequentially.

    Args:
        geometry_infilename (str, optional): The filename of the geometry input file. Defaults to '../classes/initial_universes/sample-g0-T3.cdt'.
        strictness (int, optional): The strictness parameter. Defaults to 3.
        chains (int, optional): The number of chains. Defaults to 10.
    """
    # Set parameters
    k0_values = [0, 1, 2, 3, 4, 5, 6, 7]
    k3_init_guess = [0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
    all_data = []
    sweeps = 0
    thermal_sweeps = 50
    k_steps = 10000
    target_volume = 1000
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
        np.savetxt(f'saved_universes/k3/measurements/critical_k3_values_chain={chain}.txt', critical_k3_values, delimiter=',')

    # Save the critical k3 values
    df = pd.DataFrame(all_data, columns=k0_values)
    df.index.name = 'chain'
    df.columns.name = 'k0'
    df.to_csv('saved_universes/k3/measurements/critical_k3_values.csv', index=True)

def simulation_worker(geometry_infilename: str, strictness: int, k0: int, k3_init_guess: float, sweeps: int, thermal_sweeps: int, k_steps: int, target_volume: int, target2_volume: int, chain: int) -> float:
    """
    Runs a simulation for a given set of parameters and returns the critical k3 value.

    Args:
        geometry_infilename (str): The filename of the geometry input file.
        strictness (int): The strictness parameter.
        k0 (int): The k0 value.
        k3_init_guess (float): The initial guess for k3.
        sweeps (int): The number of sweeps.
        thermal_sweeps (int): The number of thermal sweeps.
        k_steps (int): The number of k steps.
        target_volume (int): The target volume.
        target2_volume (int): The second target volume.
        chain (int): The chain number.

    Returns:
        float: The critical k3 value.
    """
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

def critical_k3_parallel(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3, chains: int = 10) -> pd.DataFrame:
    """
    Runs the critical k3 simulations in parallel.

    Args:
        geometry_infilename (str, optional): _description_. Defaults to '../classes/initial_universes/sample-g0-T3.cdt'.
        strictness (int, optional): _description_. Defaults to 3.
        chains (int, optional): _description_. Defaults to 10.

    Returns:
        pd.DataFrame: Dataframe containing the critical k3 values for each chain and k0 value.
    """
    # Set parameters
    k0_values = [0, 1, 2, 3, 4, 5, 6, 7]
    k3_init_guess = [0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
    k0_values = [7]
    k3_init_guess = [2.1]
    all_data = []
    sweeps = 0
    thermal_sweeps = 50
    k_steps = 10000
    target_volume = 1000
    target2_volume = 0

    for chain in range(chains):
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            # Prepare arguments
            items = [(geometry_infilename, strictness, k0, k3_init_guess[i], sweeps, thermal_sweeps, k_steps, target_volume, target2_volume, chain) for i, k0 in enumerate(k0_values)]
            critical_k3_values = pool.starmap(simulation_worker, items)

            # Save critical k3 values
            timestamp = datetime.datetime.now().strftime("%H.%M")
            np.savetxt(f'saved_universes/k3/measurements/critical_k3_values_{timestamp}.txt', critical_k3_values, delimiter=',')
            all_data.append(critical_k3_values)

    # Save the critical k3 values as a dataframe
    df = pd.DataFrame(all_data, columns=k0_values)
    df.index.name = 'chain'
    df.columns.name = 'k0'
    df.to_csv('saved_universes/k3/measurements/critical_k3_values.csv', index=True)

    return df

def plot_critical_k3(load_df: bool = False, df_path: str = ""):
    if load_df:
        df = pd.read_csv(df_path, index_col=0)
    else:
        df = critical_k3_parallel()
        # test = {0: [0.75, 0.73, 0.7],
        #         1: [0.91, 0.89, 0.87],
        #         2: [1.1, 1.15, 1.2],
        #         3: [1.3, 1.35, 1.4],
        #         4: [1.5, 1.55, 1.53],
        #         5: [1.7, 1.78, 1.73],
        #         6: [1.9, 1.95, 1.98],
        #         7: [2.1, 2.05, 2.03]}
        
        # df = pd.DataFrame(test, columns=test.keys()) 
        # df.index.name = 'chain'
        # df.columns.name = 'k0'

    print(df)

    # Plot the critical k3 values
    fig, ax = plt.subplots(figsize=(8, 6))
    df.mean().plot(ax=ax, yerr=df.std(), fmt='.', capsize=5, color='royalblue')
    ax.set_title('Critical k3 values for different k0 values', fontsize=18)
    ax.set_xlabel('k0')
    ax.set_ylabel('k3')
    ax.grid(True, which="both", ls="-", alpha=0.6)
    ax.bbox_inches = 'tight'
    plt.savefig('plots/critical_k3_values.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    # critical_k3(chains=1)
    critical_k3_parallel(chains=1)

    # plot_critical_k3(load_df=False)
    print("Done")

