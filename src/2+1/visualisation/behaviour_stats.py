import sys
sys.path.append('..')
import matplotlib.pyplot as plt
from typing import Tuple
import multiprocessing
from functools import partial
import random
from classes.universe import Universe
from classes.simulation import Simulation
import matplotlib as mpl
import pandas as pd
import numpy as np
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
        seed=np.random.randint(0, 1000000),
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
        with multiprocessing.Pool() as pool:
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
    ax.set_title('Critical $k_3$ values for different $k_0$ values', fontsize=18)
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
        seed=np.random.randint(0, 1000000),
        outfile=f"N/measurements/swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtvol={target_volume}_k0={k0}_chain={chain}",
        validity_check=False,
        v1=1,
        v2=1,
        v3=1
    )

    # Return the N22/N3 or N22/N31 ratio based on which tetra volume is fixed
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
    k3_init_guess = [0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
    all_data = []
    sweeps = 50
    thermal_sweeps = 50
    k_steps = 100000
    target_volume = 10000
    target2_volume = 0

    for chain in range(chains):
        with multiprocessing.Pool() as pool:
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
    df.to_csv(f'saved_universes/N/4filling_phase_transition_values_swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtv={target_volume}.csv', index=True)

    print(df)

    return df

def phase_transition_sequential(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3, chains: int = 10) -> pd.DataFrame:
    """
    Runs the phase transition simulations sequentially.

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
    k3_init_guess = [0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
    all_data = []
    sweeps = 50
    thermal_sweeps = 50
    k_steps = 100000
    target_volume = 10000
    target2_volume = 0

    for chain in range(chains):
        results = []

        # Prepare arguments
        for i, k0 in enumerate(k0_values):
            phase_transition_values = phase_transition_worker(geometry_infilename, strictness, k0, k3_init_guess[i], sweeps, thermal_sweeps, k_steps, target_volume, target2_volume, chain)
            results.append(phase_transition_values)

        all_data.append(results)

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
        df_full = pd.read_csv(df_path, index_col=0)
        df = df_full.drop(columns=['2.5', '3.2', '3.5'])
        df.columns = [0, 1, 2, 3, 4, 5, 6, 7]
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

def acceptance_worker(chain_num: int, geometry_infilename: str, strictness: int, k0: int, k3: float, sweeps: int, thermal_sweeps: int, k_steps: int, target_volume: int, target2_volume: int) -> Tuple[float]:
    """
    CDT simulation.
    
    Args:
        chain_num (int): The chain number.
        geometry_infilename (str): The filename of the geometry input file.
        strictness (int): The manifold strictness parameter.
        k0 (int): The k0 value.
        k3 (float): The k3 value.
        sweeps (int): The number of sweeps.
        thermal_sweeps (int): The number of thermal sweeps.
        k_steps (int): The number of k steps.
        target_volume (int): The target volume.
        target2_volume (int): The second target volume.
    
    Returns:
        Tuple[float]: The acceptance probabilities and success rates for different moves.
    """
    universe = Universe(geometry_infilename=geometry_infilename, strictness=strictness)
    simulation = Simulation(universe)
    simulation.save_process = False
    simulation.save_final = False
    simulation.as_pickle = True
    volfix_switch = 0
    print(f"Running k0={k0}, k3={k3}, chain={chain_num}")

    simulation.start(
        k0=k0,
        k3=k3,
        sweeps=sweeps,
        thermal_sweeps=thermal_sweeps,
        k_steps=k_steps,
        target_volume=target_volume,
        target2_volume=target2_volume,
        volfix_switch=volfix_switch,
        seed=np.random.randint(0, 1000000),
        outfile=f"acceptance/swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtvol={target_volume}_k0={k0}_chain={chain_num}",
        validity_check=False,
        v1=1,
        v2=1,
        v3=1
    )

    return simulation.add_ap, simulation.delete_ap, simulation.flip_ap, simulation.shift_ap, simulation.ishift_ap, simulation.succes_rates[1], simulation.succes_rates[2], simulation.succes_rates[3], simulation.succes_rates[4], simulation.succes_rates[5]

def plot_acceptance_parallel(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3, k0: int = 0, k3: float = 0.7, sweeps: int = 50, thermal_sweeps: int = 50, k_steps: int = 100000, target_volume: int = 10000, target2_volume: int = 0, chain: int = 0, savefig: bool = True):
    """
    Plots the acceptance probabilities and success rates for different moves in parallel.

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
    """
    assert chain <= multiprocessing.cpu_count(), "The number of chains must be less than or equal to the number of CPU cores."
    with multiprocessing.Pool(processes=chain) as pool:
        simulate_func = partial(acceptance_worker, geometry_infilename=geometry_infilename, strictness=strictness, k0=k0, k3=k3, sweeps=sweeps, thermal_sweeps=thermal_sweeps, k_steps=k_steps, target_volume=target_volume, target2_volume=target2_volume)
        results = pool.map(simulate_func, range(chain))
    
    pool.close()
    pool.join()

    add_ap_list, del_ap_list, flip_ap_list, shift_ap_list, ishift_ap_list, add_sr_list, del_sr_list, flip_sr_list, shift_sr_list, ishift_sr_list = zip(*results)

    all_add_ap_array = np.array(add_ap_list)
    all_del_ap_array = np.array(del_ap_list)
    all_flip_ap_array = np.array(flip_ap_list)
    all_shift_ap_array = np.array(shift_ap_list)
    all_ishift_ap_array = np.array(ishift_ap_list)

    all_add_sr_array = np.array(add_sr_list)
    all_del_sr_array = np.array(del_sr_list)
    all_flip_sr_array = np.array(flip_sr_list)
    all_shift_sr_array = np.array(shift_sr_list)
    all_ishift_sr_array = np.array(ishift_sr_list)

    # Calculate means and std
    add_ap_mean = np.mean(all_add_ap_array, axis=0)
    add_ap_std = np.std(all_add_ap_array, axis=0)
    del_ap_mean = np.mean(all_del_ap_array, axis=0)
    del_ap_std = np.std(all_del_ap_array, axis=0)
    flip_ap_mean = np.mean(all_flip_ap_array, axis=0)
    flip_ap_std = np.std(all_flip_ap_array, axis=0)
    shift_ap_mean = np.mean(all_shift_ap_array, axis=0)
    shift_ap_std = np.std(all_shift_ap_array, axis=0)
    ishift_ap_mean = np.mean(all_ishift_ap_array, axis=0)
    ishift_ap_std = np.std(all_ishift_ap_array, axis=0)

    add_sr_mean = np.mean(all_add_sr_array, axis=0)
    add_sr_std = np.std(all_add_sr_array, axis=0)
    del_sr_mean = np.mean(all_del_sr_array, axis=0)
    del_sr_std = np.std(all_del_sr_array, axis=0)
    flip_sr_mean = np.mean(all_flip_sr_array, axis=0)
    flip_sr_std = np.std(all_flip_sr_array, axis=0)
    shift_sr_mean = np.mean(all_shift_sr_array, axis=0)
    shift_sr_std = np.std(all_shift_sr_array, axis=0)
    ishift_sr_mean = np.mean(all_ishift_sr_array, axis=0)
    ishift_sr_std = np.std(all_ishift_sr_array, axis=0)

    # Plot the acceptance probabilities
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(add_ap_mean, label='Add', color='b')
    ax.fill_between(range(len(add_ap_mean)), add_ap_mean - add_ap_std, add_ap_mean + add_ap_std, alpha=0.3, color='b')
    ax.plot(del_ap_mean, label='Delete', color='r')
    ax.fill_between(range(len(del_ap_mean)), del_ap_mean - del_ap_std, del_ap_mean + del_ap_std, alpha=0.3, color='r')
    ax.plot(flip_ap_mean, label='Flip', color='g')
    ax.fill_between(range(len(flip_ap_mean)), flip_ap_mean - flip_ap_std, flip_ap_mean + flip_ap_std, alpha=0.3, color='g')
    ax.plot(shift_ap_mean, label='Shift', color='y')
    ax.fill_between(range(len(shift_ap_mean)), shift_ap_mean - shift_ap_std, shift_ap_mean + shift_ap_std, alpha=0.3, color='y')
    ax.plot(ishift_ap_mean, label='Inverse shift', color='purple')
    ax.fill_between(range(len(ishift_ap_mean)), ishift_ap_mean - ishift_ap_std, ishift_ap_mean + ishift_ap_std, alpha=0.3, color='purple')
    ax.set_title(f'Acceptance ratios of MCMC algorithm ($k_0 = {k0}$)', fontsize=18.5)
    ax.set_xlabel('Thermal sweep', fontsize=18)  
    ax.set_ylabel('Acceptance ratio', fontsize=18)
    ax.legend(fontsize=12, loc='upper left')
    ax.grid(True, which="both", ls="-", alpha=0.6)
    ax.bbox_inches = 'tight'
    step = 1
    ax.set_xticks(np.arange(0, len(add_ap_mean) + step, step))
    ax.set_xticklabels(np.arange(0, len(add_ap_mean) + step, step, dtype=int))
    if savefig:
        plt.savefig(f'plots/acceptance_ratios_k0={k0}.png', dpi=400)

    # Plot the success rates
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(add_sr_mean, label='Add', color='b')
    ax.fill_between(range(len(add_sr_mean)), add_sr_mean - add_sr_std, add_sr_mean + add_sr_std, alpha=0.3, color='b')
    ax.plot(del_sr_mean, label='Delete', color='r')
    ax.fill_between(range(len(del_sr_mean)), del_sr_mean - del_sr_std, del_sr_mean + del_sr_std, alpha=0.3, color='r')
    ax.plot(flip_sr_mean, label='Flip', color='g')
    ax.fill_between(range(len(flip_sr_mean)), flip_sr_mean - flip_sr_std, flip_sr_mean + flip_sr_std, alpha=0.3, color='g')
    ax.plot(shift_sr_mean, label='Shift', color='y')
    ax.fill_between(range(len(shift_sr_mean)), shift_sr_mean - shift_sr_std, shift_sr_mean + shift_sr_std, alpha=0.3, color='y')
    ax.plot(ishift_sr_mean, label='Inverse shift', color='purple')
    ax.fill_between(range(len(ishift_sr_mean)), ishift_sr_mean - ishift_sr_std, ishift_sr_mean + ishift_sr_std, alpha=0.3, color='purple')
    ax.set_title(f'Success rates of MCMC algorithm ($k_0 = {k0})$', fontsize=18.5)
    ax.set_xlabel('Thermal sweep', fontsize=18)
    ax.set_ylabel('Success rate', fontsize=18)
    ax.legend(fontsize=12, loc='upper left')
    ax.grid(True, which="both", ls="-", alpha=0.6)
    ax.bbox_inches = 'tight'
    step = 1
    ax.set_xticks(np.arange(0, len(add_sr_mean) + step, step))
    ax.set_xticklabels(np.arange(0, len(add_sr_mean) + step, step, dtype=int))
    if savefig:
        plt.savefig(f'plots/success_rates_k0={k0}.png', dpi=400)

    plt.show()

def plot_acceptance_sequential(geometry_infilename: str = '../classes/initial_universes/sample-g0-T3.cdt', strictness: int = 3, k0: int = 0, k3: float = 0.7, sweeps: int = 50, thermal_sweeps: int = 50, k_steps: int = 100000, target_volume: int = 10000, target2_volume: int = 0, chain: int = 0):
    """
    Plots the acceptance probabilities and success rates for different moves.

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
    """
    add_ap_list = []
    del_ap_list = []
    flip_ap_list = []
    shift_ap_list = []
    ishift_ap_list = []

    add_sr_list = []
    del_sr_list = []
    flip_sr_list = []
    shift_sr_list = []
    ishift_sr_list = []

    for _ in range(chain):
        universe = Universe(geometry_infilename=geometry_infilename, strictness=strictness)
        simulation = Simulation(universe)
        simulation.save_process = False
        simulation.save_final = False
        simulation.as_pickle = True
        volfix_switch = 0
        print(f"Running k0={k0}, k3={k3}, chain={chain}")

        simulation.start(
            k0=k0,
            k3=k3,
            sweeps=sweeps,
            thermal_sweeps=thermal_sweeps,
            k_steps=k_steps,
            target_volume=target_volume,
            target2_volume=target2_volume,
            volfix_switch=volfix_switch,
            seed=np.random.randint(0, 1000000),
            outfile=f"acceptance/swps={sweeps}_tswps={thermal_sweeps}_kstps={k_steps}_trgtvol={target_volume}_k0={k0}_chain={chain}",
            validity_check=False,
            v1=1,
            v2=1,
            v3=1
        )
        
        # Append to the general lists
        add_ap_list.append(simulation.add_ap)
        del_ap_list.append(simulation.delete_ap)
        flip_ap_list.append(simulation.flip_ap)
        shift_ap_list.append(simulation.shift_ap)
        ishift_ap_list.append(simulation.ishift_ap)
        
        add_sr_list.append(simulation.succes_rates[1])
        del_sr_list.append(simulation.succes_rates[2])
        flip_sr_list.append(simulation.succes_rates[3])
        shift_sr_list.append(simulation.succes_rates[4])
        ishift_sr_list.append(simulation.succes_rates[5])

    all_add_ap_array = np.array(add_ap_list)
    all_del_ap_array = np.array(del_ap_list)
    all_flip_ap_array = np.array(flip_ap_list)
    all_shift_ap_array = np.array(shift_ap_list)
    all_ishift_ap_array = np.array(ishift_ap_list)

    all_add_sr_array = np.array(add_sr_list)
    all_del_sr_array = np.array(del_sr_list)
    all_flip_sr_array = np.array(flip_sr_list)
    all_shift_sr_array = np.array(shift_sr_list)
    all_ishift_sr_array = np.array(ishift_sr_list)
    
    # Mean and std
    add_ap_mean = np.mean(all_add_ap_array, axis=0)
    add_ap_std = np.std(all_add_ap_array, axis=0)
    del_ap_mean = np.mean(all_del_ap_array, axis=0)
    del_ap_std = np.std(all_del_ap_array, axis=0)
    flip_ap_mean = np.mean(all_flip_ap_array, axis=0)
    flip_ap_std = np.std(all_flip_ap_array, axis=0)
    shift_ap_mean = np.mean(all_shift_ap_array, axis=0)
    shift_ap_std = np.std(all_shift_ap_array, axis=0)
    ishift_ap_mean = np.mean(all_ishift_ap_array, axis=0)
    ishift_ap_std = np.std(all_ishift_ap_array, axis=0)

    add_sr_mean = np.mean(all_add_sr_array, axis=0)
    add_sr_std = np.std(all_add_sr_array, axis=0)
    del_sr_mean = np.mean(all_del_sr_array, axis=0)
    del_sr_std = np.std(all_del_sr_array, axis=0)
    flip_sr_mean = np.mean(all_flip_sr_array, axis=0)
    flip_sr_std = np.std(all_flip_sr_array, axis=0)
    shift_sr_mean = np.mean(all_shift_sr_array, axis=0)
    shift_sr_std = np.std(all_shift_sr_array, axis=0)
    ishift_sr_mean = np.mean(all_ishift_sr_array, axis=0)
    ishift_sr_std = np.std(all_ishift_sr_array, axis=0)

    # Plot the acceptance probabilities
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(add_ap_mean, label='Add', color='b')
    ax.fill_between(range(len(add_ap_mean)), add_ap_mean - add_ap_std, add_ap_mean + add_ap_std, alpha=0.3, color='b')
    ax.plot(del_ap_mean, label='Delete', color='r')
    ax.fill_between(range(len(del_ap_mean)), del_ap_mean - del_ap_std, del_ap_mean + del_ap_std, alpha=0.3, color='r')
    ax.plot(flip_ap_mean, label='Flip', color='g')
    ax.fill_between(range(len(flip_ap_mean)), flip_ap_mean - flip_ap_std, flip_ap_mean + flip_ap_std, alpha=0.3, color='g')
    ax.plot(shift_ap_mean, label='Shift', color='y')
    ax.fill_between(range(len(shift_ap_mean)), shift_ap_mean - shift_ap_std, shift_ap_mean + shift_ap_std, alpha=0.3, color='y')
    ax.plot(ishift_ap_mean, label='Inverse shift', color='purple')
    ax.fill_between(range(len(ishift_ap_mean)), ishift_ap_mean - ishift_ap_std, ishift_ap_mean + ishift_ap_std, alpha=0.3, color='purple')
    ax.set_title(f'Acceptance probabilities of MCMC algorithm ($k_0 = {k0}$)', fontsize=18.5)
    ax.set_xlabel('Sweep', fontsize=18)  
    ax.set_ylabel('Acceptance ratio', fontsize=18)
    ax.legend(fontsize=14, loc='lower right')
    ax.grid(True, which="both", ls="-", alpha=0.6)
    ax.bbox_inches = 'tight'
    ax.set_xticks(np.arange(len(add_ap_mean)))
    ax.set_xticklabels(np.arange(len(add_ap_mean), dtype=int))
    # plt.savefig(f'plots/acceptance_ratios_k0={k0}.png', dpi=400)

    # Plot the success rates
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(add_sr_mean, label='Add', color='b')
    ax.fill_between(range(len(add_sr_mean)), add_sr_mean - add_sr_std, add_sr_mean + add_sr_std, alpha=0.3, color='b')
    ax.plot(del_sr_mean, label='Delete', color='r')
    ax.fill_between(range(len(del_sr_mean)), del_sr_mean - del_sr_std, del_sr_mean + del_sr_std, alpha=0.3, color='r')
    ax.plot(flip_sr_mean, label='Flip', color='g')
    ax.fill_between(range(len(flip_sr_mean)), flip_sr_mean - flip_sr_std, flip_sr_mean + flip_sr_std, alpha=0.3, color='g')
    ax.plot(shift_sr_mean, label='Shift', color='y')
    ax.fill_between(range(len(shift_sr_mean)), shift_sr_mean - shift_sr_std, shift_sr_mean + shift_sr_std, alpha=0.3, color='y')
    ax.plot(ishift_sr_mean, label='Inverse shift', color='purple')
    ax.fill_between(range(len(ishift_sr_mean)), ishift_sr_mean - ishift_sr_std, ishift_sr_mean + ishift_sr_std, alpha=0.3, color='purple')
    ax.set_title(f'Success rates of MCMC algorithm ($k_0 = {k0})$', fontsize=18.5)
    ax.set_xlabel('Sweep', fontsize=18)
    ax.set_ylabel('Success rate', fontsize=18)
    ax.legend(fontsize=14, loc='lower right')
    ax.grid(True, which="both", ls="-", alpha=0.6)
    ax.bbox_inches = 'tight'
    ax.set_xticks(np.arange(len(add_sr_mean)))
    ax.set_xticklabels(np.arange(len(add_sr_mean), dtype=int))
    # plt.savefig(f'plots/success_rates_k0={k0}.png', dpi=400)

    plt.show()


if __name__ == "__main__":
    # critical_k3_parallel(chains=10)
    # plot_critical_k3(load_df=True, df_path='saved_universes/k3/measurements/critical_k3_values_swps=0_tswps=50_kstps=100000_trgtv=10000.csv')
    
    # phase_transition_parallel(chains=5)
    # plot_phase_transition(load_df=0, fixvol_type=0, df_path='saved_universes/N/phase_transition_values_swps=50_tswps=50_kstps=100000_trgtv=10000.csv')

    # plot_acceptance_parallel(geometry_infilename='../classes/initial_universes/sample-g0-T3.cdt', strictness=3, k0=0, k3=0.7, sweeps=0, thermal_sweeps=50, k_steps=100000, target_volume=10000, target2_volume=0, chain=5, saevfig=True)
    # plot_acceptance_parallel(k0=3, k3=1.3, sweeps=0, thermal_sweeps=10, k_steps=100000, target_volume=10000, chain=5, savefig=True)
    plot_acceptance_sequential(k0=7, k3=2.1, sweeps=50, thermal_sweeps=50, k_steps=100000, target_volume=10000, chain=1)
    print("Done")

