import sys
sys.path.append('..')
import matplotlib.pyplot as plt
import multiprocessing
from classes.universe import Universe
from classes.simulation import Simulation
import matplotlib as mpl
import pandas as pd
mpl.rcParams['font.size'] = 16


def critical_k3_simulation_worker(geometry_infilename: str, strictness: int, k0: int, k3_init_guess: float, sweeps: int, thermal_sweeps: int, k_steps: int, target_volume: int, target2_volume: int, chain: int) -> float:
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
    simulation.save_process = True
    simulation.save_final = True
    simulation.as_pickle = True

    print(f"Running k0={k0}, k3_init_guess={k3_init_guess}, chain={chain}")
    simulation.start(
        k0=k0,
        k3=k3_init_guess,
        sweeps=sweeps,
        thermal_sweeps=thermal_sweeps,
        k_steps=k_steps,
        target_volume=target_volume,
        target2_volume=target2_volume,
        volfix_switch=0,
        seed=chain,
        outfile=f"k3/swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtvol={target_volume}_k0={k0}_chain={chain}",
        validity_check=True,
        v1=1,
        v2=1,
        v3=1
    )
    print(f"Critical k3 value of k0 {k0}: {simulation.k3}\n")

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
    all_data = []
    sweeps = 0
    thermal_sweeps = 50
    k_steps = 100000
    target_volume = 10000
    target2_volume = 0

    for chain in range(chains):
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            # Prepare arguments
            items = [(geometry_infilename, strictness, k0, k3_init_guess[i], sweeps, thermal_sweeps, k_steps, target_volume, target2_volume, chain) for i, k0 in enumerate(k0_values)]
            critical_k3_values = pool.starmap(critical_k3_simulation_worker, items)

            # Save critical k3 values
            all_data.append(critical_k3_values)

        # Close the pool to prevent any more tasks from being submitted
        pool.close()
        # Wait for all tasks to finish
        pool.join()

    # Save the critical k3 values as a dataframe
    df = pd.DataFrame(all_data, columns=k0_values)
    df.index.name = 'chain'
    df.columns.name = 'k0'
    df.to_csv(f'saved_universes/k3/measurements/critical_k3_values_swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtv={target_volume}.csv', index=True)

    return df

def plot_critical_k3(load_df: bool = False, df_path: str = ""):
    """
    Plots the critical k3 values for different k0 values.

    Args:
        load_df (bool, optional): Whether to load the dataframe from a file. Defaults to False.
        df_path (str, optional): The path to the dataframe file. Defaults to "".
    """
    if load_df:
        df = pd.read_csv(df_path, index_col=0)
    else:
        df = critical_k3_parallel()

    print(df)

    # Plot the critical k3 values
    fig, ax = plt.subplots(figsize=(8, 6))
    df.mean().plot(ax=ax, yerr=df.std(), fmt='.', capsize=5, color='royalblue')
    ax.set_title('Critical $k_0$ values for different $k_3$ values', fontsize=18)
    ax.set_xlabel('$k_0$')
    ax.set_ylabel('$k_3$')
    ax.grid(True, which="both", ls="-", alpha=0.6)
    ax.bbox_inches = 'tight'
    plt.savefig('plots/critical_k3_values.png', dpi=400)
    plt.show()

def phase_transition_worker(geometry_infilename: str, strictness: int, k0: int, k3_init_guess: float, sweeps: int, thermal_sweeps: int, k_steps: int, target_volume: int, target2_volume: int, chain: int):
    """
    Runs a simulation for a given set of parameters and returns the N22/N31 ratio.

    Args:
        geometry_infilename (str): The filename of the geometry input file.
        strictness (int): The manifold strictness parameter.
        k0 (int): The k0 value.
        k3_init_guess (float): The initial guess for k3.
        sweeps (int): The number of sweeps.
        thermal_sweeps (int): The number of thermal sweeps.
        k_steps (int): The number of k steps.
        target_volume (int): The target volume.
        target2_volume (int): The second target volume.
        chain (int): The chain number.

    Returns:
        float: The N22/N31 ratio.
    """
    universe = Universe(geometry_infilename=geometry_infilename, strictness=strictness)
    simulation = Simulation(universe)
    simulation.save_process = True
    simulation.save_final = True
    simulation.as_pickle = True
    volfix_switch = 0
    print(f"Running k0={k0}, k3_init_guess={k3_init_guess}, chain={chain}")
    simulation.start(
        k0=k0,
        k3=k3_init_guess,
        sweeps=sweeps,
        thermal_sweeps=thermal_sweeps,
        k_steps=k_steps,
        target_volume=target_volume,
        target2_volume=target2_volume,
        volfix_switch=volfix_switch,
        seed=chain,
        outfile=f"N/measurements/swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtvol={target_volume}_k0={k0}_chain={chain}",
        validity_check=False,
        v1=1,
        v2=1,
        v3=1
    )

    # N3 = simulation.universe.tetrahedron_pool.get_number_occupied()
    # N31 = simulation.universe.tetras_31.get_number_occupied()
    # N22 = simulation.universe.tetras_22.get_number_occupied()
    # print(f"N3: {N3}, N31: {N31}, N13: {N3 - N31 - N22}, N22: {N22}")

    if volfix_switch == 1:
        result = simulation.expected_N22_N3
    else:
        result = simulation.expected_N22_N31

    return result

def phase_transition_parallel(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3, chains: int = 10) -> pd.DataFrame:
    """
    Runs the phase transition simulations in parallel.

    Args:
        geometry_infilename (str, optional): The filename of the geometry input file.
                                             Defaults to '../classes/initial_universes/sample-g0-T3.cdt'.
        strictness (int, optional): The manifold strictness parameter. Defaults to 3.
        chains (int, optional): The number of chains. Defaults to 10.

    Returns:
        pd.DataFrame: Dataframe containing the N22/N3 values for each chain and k0 value.
    """
    # Set parameters
    k0_values = [0, 1, 2, 3, 4, 5, 6, 7]
    k3_init_guess = [0.72, 0.92, 1.11, 1.31, 1.52, 1.75, 1.95, 2.15]
    all_data = []
    sweeps = 10
    thermal_sweeps = 50
    k_steps = 100000
    target_volume = 10000
    target2_volume = 0
    processes  = len(k0_values)
    assert processes <= multiprocessing.cpu_count(), "Number of processes exceeds the number of available CPUs."

    for chain in range(chains):
        with multiprocessing.Pool(processes=processes) as pool:
            # Prepare arguments
            items = [(geometry_infilename, strictness, k0, k3_init_guess[i], sweeps, thermal_sweeps, k_steps, target_volume, target2_volume, chain) for i, k0 in enumerate(k0_values)]
            phase_transition_values = pool.starmap(phase_transition_worker, items)

            # Save critical k3 values
            all_data.append(phase_transition_values)

        # Close the pool to prevent any more tasks from being submitted
        pool.close()
        # Wait for all tasks to finish
        pool.join()


    # Save the critical k3 values as a dataframe
    df = pd.DataFrame(all_data, columns=k0_values)
    df.index.name = 'chain'
    df.columns.name = 'k0'
    df.to_csv(f'saved_universes/N/phase_transition_values_swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtv={target_volume}.csv', index=True)

    print(df)

    return df

def plot_phase_transition(load_df: int = 0, fixvol_type: int = 1, df_path: str = ""):
    """
    Plots the N22/N31 values for different k0 values.

    Args:
        load_df (int): 
            0: Load the dataframe from a file.
            1: Run the phase transition simulations in parallel.
        fixvol_type (int):
            1: N22/N3 values.
            2: N22/N31 values.
        df_path (str, optional): The path to the dataframe file. Defaults to "".
    """
    if load_df == 0:
        df = pd.read_csv(df_path, index_col=0)
    else:
        df = phase_transition_parallel()

    print(df)

    # Plot the N22/N3 values
    fig, ax = plt.subplots(figsize=(8, 6))
    df.mean().plot(ax=ax, yerr=df.std(), fmt='.', capsize=5, color='royalblue')
    ax.set_xlabel('$k_0$')
    ax.grid(True, which="both", ls="-", alpha=0.6)
    ax.bbox_inches = 'tight'
    if fixvol_type == 1:
        ax.set_title('$N_{22}/N_{3}$ values for different $k_0$ values', fontsize=18)
        ax.set_ylabel('$N_{22}/N_{3}$')
        plt.savefig('plots/N22_N3_values.png', dpi=400)
    else:
        ax.set_title('$N_{22}/N_{31}$ values for different $k_0$ values', fontsize=18)
        ax.set_ylabel('$N_{22}/N_{31}$')
        plt.savefig('plots/N22_N31_values.png', dpi=400)

    plt.show()


if __name__ == "__main__":
    # critical_k3_parallel(chains=10)
    # plot_critical_k3(load_df=True, df_path='saved_universes/k3/measurements/critical_k3_values_swps=0_tswps=50_kstps=100000_trgtv=10000.csv')
    
    phase_transition_parallel(chains=1)
    # plot_phase_transition(load_df=1, df_path='saved_universes/N/phase_transition_values_swps=10_tswps=50_kstps=100000_trgtv=10000.csv')
    print("Done")

