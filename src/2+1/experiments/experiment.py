import sys
sys.path.append('..')
from typing import Tuple, Dict, List
import multiprocessing
from functools import partial
from classes.universe import Universe
from classes.simulation import Simulation
import numpy as np
import pickle
import gzip


def simulation_worker(
        geometry_infilename: str,
        strictness: int,
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
        save_process: bool = False,
        save_final: bool = False,
        as_pickle: bool = False,
        save_observables: bool = False
        ) -> Dict[str, any]:
    """
    Run a simulation with the given parameters and save the observables of the final state.

    Args:
        geometry_infilename (str): The path to the input file containing the initial universe.
        strictness (int): The strictness of the universe.
        k0 (int): The value of the k0 parameter.
        k3_init_guess (float): The initial guess for the k3 parameter.
        thermal_sweeps (int): The number of thermal sweeps to perform.
        sweeps (int): The number of sweeps to perform.
        k_steps (int): The number of k steps to perform per (thermal) sweep.
        target_volume (int, optional): The target (3,1)/(3)-volume. Defaults to 10000.
        target2_volume (int, optional): The target triangle volume. Defaults to 0.
        volfix_switch (int, optional): The volume fix switch, 0 for type (3,1), 1 for (3) tetrahedra. Defaults to 0.
        chain (int, optional): The chain number. Defaults to 0.
        validity_check (bool, optional): Whether to perform a validity check. Defaults to False.
        save_process (bool, optional): Whether to save the process. Defaults to False.
        save_final (bool, optional): Whether to save the final state. Defaults to False.
        as_pickle (bool, optional): Whether to save the final state as a pk.gz file. Defaults to False.
        save_observables (bool, optional): Whether to save the observables of the final state. Defaults to False.

    Returns:
        Dict[str, any]: A dictionary containing the observables of the final state.
    """
    # Create universe and simulation
    universe = Universe(geometry_infilename=geometry_infilename, strictness=strictness)
    simulation = Simulation(universe)
    simulation.save_process = save_process
    simulation.save_final = save_final
    simulation.as_pickle = as_pickle

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
        seed=np.random.randint(0, 1000000),
        outfile=f"universe_k0={k0}_tswps={thermal_sweeps}_swps={sweeps}_kstps={k_steps}" +
                f"_trgtv={target_volume}_trgtv2={target2_volume}_fx={volfix_switch}_chn={chain}",
        validity_check=validity_check,
        v1=1,
        v2=1,
        v3=1
    )

    # Save observables of final state
    k3: int = simulation.k3
    slice_sizes: Dict[int] = simulation.universe.slice_sizes
    slab_sizes: Dict[int] = simulation.universe.slab_sizes
    n_slices: int = simulation.universe.n_slices
    expected_n22_n31: float = simulation.expected_N22_N31
    expected_n22_n3: float = simulation.expected_N22_N3
    curvature_profile: Dict[int, List[int]] = simulation.universe.get_curvature_profile()

    # Observables measured over time (lists/dicts)
    acceptance_ratios: Dict[int, List[float]] = simulation.acceptance_ratios
    succes_rates: Dict[int, List[float]] = simulation.succes_rates
    total_vertices: List[int] = simulation.total_vertices
    total_tetrahedra: List[int] = simulation.total_tetrahedra
    total_31_tetrahedra: List[int] = simulation.total_31_tetrahedra
    total_22_tetrahedra: List[int] = simulation.total_22_tetrahedra

    all_data = {
        'k3': k3,
        'slice_sizes': slice_sizes,
        'slab_sizes': slab_sizes,
        'n_slices': n_slices,
        'expected_n22_n31': expected_n22_n31,
        'expected_n22_n3': expected_n22_n3,
        'curvature_profile': curvature_profile,
        'acceptance_ratios': acceptance_ratios,
        'succes_rates': succes_rates,
        'total_vertices': total_vertices,
        'total_tetrahedra': total_tetrahedra,
        'total_31_tetrahedra': total_31_tetrahedra,
        'total_22_tetrahedra': total_22_tetrahedra
    }

    # Save observables
    if save_observables:
        with gzip.open(f"measurements/simulation_k0={k0}_tswps={thermal_sweeps}_swps={sweeps}_kstps={k_steps}" +
                       f"_trgtv={target_volume}_trgtv2={target2_volume}_fx={volfix_switch}_chn={chain}.pkl.gz", 'wb', compresslevel=1) as f:
            pickle.dump(all_data, f)

    return all_data

def main(geometry_infilename: str):
    """
    Run n experiments with the given parameters for each k0 value and initial k3 guess.

    Args:
        geometry_infilename (str): The path to the input file containing the initial universe.
    """
    chains: int = 1
    k0_values: List[int] = [0, 1, 2, 3, 4, 5, 6, 7]
    k3_init_guesses: List[float] = [0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1]
    k0_values = [0]
    k3_init_guesses = [0.7]
    for chain in range(chains):
        for i, k0 in enumerate(k0_values):
            simulation_worker(
                geometry_infilename=geometry_infilename,
                strictness=3,
                k0=k0,
                k3_init_guess=k3_init_guesses[i],
                thermal_sweeps=10,
                sweeps=0,
                k_steps=100000,
                target_volume=10000,
                target2_volume=0,
                volfix_switch=0,
                chain=chain,
                validity_check=False,
                save_process=False,
                save_final=False,
                as_pickle=True,
                save_observables=True
        )


if __name__ == "__main__":
    main(geometry_infilename='../classes/initial_universes/output_g=0_T=10.txt')
