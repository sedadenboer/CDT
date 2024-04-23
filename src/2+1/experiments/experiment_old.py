import sys
sys.path.append('..')
from typing import Dict, List
from classes.universe import Universe
from classes.simulation import Simulation
import pickle
import gzip
import multiprocessing as mp


def run_simulation(
        geometry_infilename: str,
        strictness: int,
        seed: int,
        k0: int,
        k3_init_guess: float,
        thermal_sweeps: int,
        sweeps: int,
        k_steps: int,
        target_volume: int = 10000, 
        target2_volume: int = 0,
        volfix_switch: int = 0,
        chain: int = 0,
        validity_check: bool = False,
        save_thermal: bool = False,
        save_main: bool = False,
        as_pickle: bool = False,
        measuring: bool = True,
        save_observables: bool = False,
        saving_interval: int = 1,
        measuring_interval: int = 1,
        ) -> Dict[str, any]:
    """
    Run a simulation with the given parameters and save the observables of the final state.

    Args:
        geometry_infilename (str): The path to the input file containing the initial universe.
        strictness (int): The strictness of the universe.
        seed (int): Seed for random generator.
        k0 (int): The value of the k0 parameter.
        k3_init_guess (float): The initial guess for the k3 parameter.
        thermal_sweeps (int): The number of thermal sweeps to perform.
        sweeps (int): The number of sweeps to perform.
        k_steps (int): The number of k steps to perform per (thermal) sweep.
        target_volume (int): The target (3,1)/(3)-volume. Defaults to 10000.
        target2_volume (int): The target triangle volume. Defaults to 0.
        volfix_switch (int): The volume fix switch, 0 for type (3,1), 1 for (3) tetrahedra. Defaults to 0.
        chain (int): The chain number. Defaults to 0.
        validity_check (bool): Whether to perform a validity check. Defaults to False.
        save_thermal (bool): Whether to save geometries from thermal sweeps. Defaults to False.
        save_main (bool): Whether to save geometries from main sweeps. Defaults to False.
        as_pickle (bool): Whether to save the final state as a pk.gz file. Defaults to False.
        measuring (bool): Whether to perform measurements during the main sweeps. Defaults to True.
        save_observables (bool): Whether to save the observables of the final state. Defaults to False.
        saving_interval (int): When to output geometries (every n sweeps). Defaults to 1.
        meauring_interval (int): When to measure observables (every n sweeps). Defaults to 1.

    Returns:
        Dict[str, any]: A dictionary containing the observables of the final state.
    """
    # Create universe and simulation
    universe = Universe(geometry_infilename=geometry_infilename, strictness=strictness)
    simulation = Simulation(universe)
    simulation.save_thermal = save_thermal
    simulation.save_main = save_main
    simulation.as_pickle = as_pickle
    simulation.measuring = measuring
    simulation.saving_interval = saving_interval
    simulation.measuring_interval = measuring_interval

    # Run simulation
    print(f"Running k0={k0}, k3_init_guess={k3_init_guess}, chain={chain}")
    simulation.start(
        k0=k0,
        k3=k3_init_guess,
        thermal_sweeps=thermal_sweeps,
        sweeps=sweeps,
        k_steps=k_steps,
        target_volume=target_volume,
        target2_volume=target2_volume,
        volfix_switch=volfix_switch,
        seed=seed,
        outfile=f"universe_k0={k0}_tswps={thermal_sweeps}_swps={sweeps}_kstps={k_steps}" +
                f"_trgtv={target_volume}_trgtv2={target2_volume}_fx={volfix_switch}_chn={chain}",
        validity_check=validity_check,
        v1=1,
        v2=1,
        v3=1
    )

    # k3 value
    k3: int = simulation.k3

    # Acceptance ratios and succeses (measured the whole simulation)
    acceptance_ratios: Dict[int, List[float]] = simulation.acceptance_ratios
    succes_rates: Dict[int, List[float]] = simulation.succes_rates

    if sweeps > 0 and measuring:
        # Save observables after mcmc is done
        n_slices: int = universe.n_slices
        expected_n22_n31: float = simulation.expected_N22_N31
        expected_n22_n3: float = simulation.expected_N22_N3
    
        # Observables measured over time (measured every main sweep)
        total_vertices: List[int] = simulation.total_vertices
        total_tetrahedra: List[int] = simulation.total_tetrahedra
        total_31_tetrahedra: List[int] = simulation.total_31_tetrahedra
        total_22_tetrahedra: List[int] = simulation.total_22_tetrahedra
        slice_sizes: List[Dict[int, int]] = simulation.measured_slice_sizes
        slab_sizes: List[Dict[int, int]] = simulation.measured_slab_sizes
        curvature_profiles: List[Dict[int, int]] = simulation.measured_curvature_profiles

        all_data = {
            'k3': k3,
            'n_slices': n_slices,
            'expected_n22_n31': expected_n22_n31,
            'expected_n22_n3': expected_n22_n3,
            'acceptance_ratios': acceptance_ratios,
            'succes_rates': succes_rates,
            'total_vertices': total_vertices,
            'total_tetrahedra': total_tetrahedra,
            'total_31_tetrahedra': total_31_tetrahedra,
            'total_22_tetrahedra': total_22_tetrahedra,
            'slice_sizes': slice_sizes,
            'slab_sizes': slab_sizes,
            'curvature_profiles': curvature_profiles,
        }

        # Save observables
        if save_observables:
            with gzip.open(f"measurements/simulation_k0={k0}_tswps={thermal_sweeps}_swps={sweeps}_kstps={k_steps}" +
                        f"_trgtv={target_volume}_trgtv2={target2_volume}_fx={volfix_switch}_chn={chain}.pkl.gz", 'wb', compresslevel=6) as f:
                pickle.dump(all_data, f)

        return all_data
    else:
        mcmc_data = {
            'k3': k3,
            'acceptance_ratios': acceptance_ratios,
            'succes_rates': succes_rates,
        }

        # Save observables
        if save_observables:
            with gzip.open(f"measurements/thermal_k0={k0}_tswps={thermal_sweeps}_swps={sweeps}_kstps={k_steps}" +
                        f"_trgtv={target_volume}_trgtv2={target2_volume}_fx={volfix_switch}_chn={chain}.pkl.gz", 'wb', compresslevel=6) as f:
                pickle.dump(mcmc_data, f)

        # Just return k3 as default
        return simulation.k3

def worker(geometry_infilename, strictness, seed, k0, k3_init_guess, thermal_sweeps, sweeps, k_steps, target_volume, target2_volume, volfix_switch, chain, validity_check, save_thermal, save_main, as_pickle, measuring, save_observables, saving_interval, measuring_interval):
    """
    Worker function to be used with multiprocessing.
    """
    return run_simulation(geometry_infilename, strictness, seed, k0, k3_init_guess, thermal_sweeps, sweeps, k_steps, target_volume, target2_volume, volfix_switch, chain, validity_check, save_thermal, save_main, as_pickle, measuring, save_observables, saving_interval, measuring_interval)

def critical_k3_parallel(
        geometry_infilename: str,
        output_filename: str,
        chain: int = 0,
        k0_values: List[float] = [0, 1, 2, 3, 4, 5, 6, 7],
        k3_init_guesses: List[float] = [0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
        ) -> Dict[int, List[float]]:
    """
    Runs n experiments with the given parameters for each k0 value and initial k3 guess, 
    and saves the final k3 values.

    Args:
        geometry_infilename (str): The path to the input file containing the initial universe.
    """
    # Prepare arguments for multiprocessing
    args_list = []
    for i, k0 in enumerate(k0_values):
        args_list.append((
            geometry_infilename,
            3, # strictness
            chain, # seed
            k0, 
            k3_init_guesses[i],
            1000, # thermal sweeps
            0,  # sweeps=0 to just measure k3
            1000000, # k_steps
            10000, # target volume
            0, # target 2 volume
            0, # volfix switch
            chain, 
            False, # validity check
            True,  # save_thermal
            False,  # save_main
            True,  # as_pickle
            False,  # save_observables
            100, # saving interval
            10 # measuring interval
        ))

    # Use multiprocessing Pool for parallel execution
    with mp.Pool() as pool:
        k3_final = pool.starmap(worker, args_list)

    # Save
    with gzip.open(output_filename + '.pkl.gz', 'wb', compresslevel=6) as f:
        pickle.dump(k3_final, f)

    print(k3_final)
    return k3_final

def generate_samples_parallel(
        chain: int = 0,
        k0_values: List[float] = [0, 1, 2, 3, 4, 5, 6, 7],
        k3_init_guesses: List[float] = [1.0417799999999955, 1.1760799999999827, 1.3237199999999882, 1.4782399999999793, 1.6480799999999842, 1.8284799999999846, 2.0621400000000256, 2.3117000000000187]
    ):
    # Prepare arguments for multiprocessing
    args_list = []
    for i, k0 in enumerate(k0_values):
        args_list.append((
            f'saved_universes/thermal_1000/universe_k0={k0}_tswps=1000_swps=0_kstps=1000000_trgtv=10000_trgtv2=0_fx=0_chn=0_thermal_1000.pkl.gz',
            3, # strictness
            chain, # seed
            k0, 
            k3_init_guesses[i],
            0, # thermal sweeps=0 since we take already generated thermalized samples
            1000,  # sweeps
            1000000, # k_steps
            10000, # target volume
            0, # target 2 volume
            0, # volfix switch
            chain, 
            False, # validity check
            False,  # save_thermal
            True,  # save_main
            True,  # as_pickle
            True, # measuring
            True,  # save_observables
            100, # saving interval
            10 # measuring interval
        ))

    # Use multiprocessing Pool for parallel execution
    with mp.Pool() as pool:
        pool.starmap(worker, args_list)

def main(geometry_infilename: str, seed: int):
    """
    Runs simulation with the given parameters for each k0 value and initial k3 guess.

    Args:
        geometry_infilename (str): The path to the input file containing the initial universe.
        seed (int): Seed for the random generator.
    """

    k0_values: List[int] = [0, 1, 2, 3, 4, 5, 6, 7]
    k3_init_guesses: List[float] = [1.0417799999999955, 1.1760799999999827, 1.3237199999999882, 1.4782399999999793, 1.6480799999999842, 1.8284799999999846, 2.0621400000000256, 2.3117000000000187]

    for i, k0 in enumerate(k0_values):
        run_simulation(
            geometry_infilename=geometry_infilename,
            strictness=3,
            seed=seed,
            k0=k0,
            k3_init_guess=k3_init_guesses[i],
            thermal_sweeps=10,
            sweeps=2, # Should be > 0, otherwise only k3 will be measured
            k_steps=1000000,
            target_volume=10000,
            target2_volume=0,
            volfix_switch=0,
            chain=seed,
            validity_check=False,
            save_thermal=False,
            save_main=True,
            as_pickle=True,
            save_observables=True,
            saving_interval=1
        )


if __name__ == "__main__":
    # main(geometry_infilename='../classes/initial_universes/sample-g0-T3.cdt', seed=0)

    # critical_k3_parallel(
    #     geometry_infilename='../classes/initial_universes/sample-g0-T3.cdt',
    #     output_filename='measurements/critical_k3_T10_trgtvN31=10000',
    #     chain=0,
    #     k0_values=[0, 1, 2, 3, 4, 5, 6, 7],
    #     k3_init_guesses=[1.0417799999999955, 1.1760799999999827, 1.3237199999999882, 1.4782399999999793, 1.6480799999999842, 1.8284799999999846, 2.0621400000000256, 2.3117000000000187]
    #     )
    
    generate_samples_parallel()
