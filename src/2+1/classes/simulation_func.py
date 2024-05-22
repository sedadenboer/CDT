import random
import numpy as np
import gc
from universe_mt import Universe
from observable import Observable
from typing import TYPE_CHECKING, List, Tuple, Dict, Any
if TYPE_CHECKING:
    from universe import Universe
from helper_functions.helpers import total_size
import multiprocessing
from multiprocessing import Pool, shared_memory
from check_mt import *
import pickle
import ctypes

N_VERTICES_TETRA = 4
N_VERTICES_TRIANGLE = 3
N_MOVES = 5


def run(universe: Universe, seed: int,
        k0: int, k3: int, tune_flag: bool = True,
        thermal_sweeps: int = 10, sweeps: int = 10, k_steps: int = 1000,
        volfix_switch: int = 0, target_volume: int = 0, target2_volume: int = 0, epsilon: float = 0.00005,
        observables: List[Observable] = [], include_mcmc_data: bool = True,
        measuring_interval: int = 1, measuring_thermal: bool = False, measuring_main: bool = False,
        save_main: bool = False, save_thermal: bool = False, saving_interval: int = 1,
        validity_check: bool = False,
        n_proposals: int = 5, outfile: str = 'output') -> None:
    """
    Starts the MCMC CDT 2+1 simulation.
    """
    global rng
    global pool
    rng = random.Random(seed)

    # Initialize data structures for MCMC data and set the frequencies of the moves
    acceptance_ratios: np.darray = np.zeros(N_MOVES)
    observables: Dict[str, Observable] = {obs: Observable(obs, thermal_sweeps, sweeps, k0, measuring_interval) for obs in observables}
    
    # If MCMC data is included, save the success and fail counts, acceptance ratios and k3 values
    successes: List[np.ndarray] = []
    fails: List[np.darray] = []
    acceptance_ratios = []
    k3_values = []
    
    assert n_proposals <= multiprocessing.cpu_count(), "Number of proposals must be less than or equal to the number of CPUs."
    pool = Pool(n_proposals)

    for obs in observables.values():
        obs.clear_data()

    # Measure at the start 
    if measuring_thermal or measuring_main:
        # Measure observables
        for name, obs in observables.items():
            obs.measure(get_observable_data(name))

    # Thermal sweeps
    if thermal_sweeps > 0:
        print("========================================\n")
        print("THERMAL SWEEPS\n")
        print("----------------------------------------\n")
        print(f"k0 = {k0}, k3 = {k3}, epsilon = {epsilon}, thermal = {thermal_sweeps}, sweeps = {sweeps}, target = {target_volume}, target2d = {target2_volume}\n")
        print("----------------------------------------\n")

        for i in range(1, thermal_sweeps + 1):
            # Get the current state of the universe and print it
            n0 = universe.vertex_pool.get_number_occupied()
            n31 = universe.tetras_31.get_number_occupied()
            n3 = universe.tetrahedron_pool.get_number_occupied()
            n22 = universe.tetras_22.get_number_occupied()
            print(f"\nThermal i: {i} \t N0: {n0}, N3: {n3}, N31: {n31}, N13: {n3 - n31 - n22}, N22: {n22}, k0: {k0}, k3: {k3}")
        
            # Perform sweeps and save general mcmc move stats
            thermal_successes, thermal_fails = perform_sweep(k_steps)
            if include_mcmc_data:
                successes.append(thermal_successes)
                fails.append(thermal_fails)
                acceptance_ratios.append([get_acceptance_probability(i) for i in range(1, N_MOVES + 1)])
                k3_values.append(k3)
            
            # Tune the k3 parameter
            if tune_flag:
                tune()

            # Check if the universe geometry is valid
            if validity_check:
                # Update geometry
                universe.update_vertices()
                universe.log()
                universe.check_validity()
            
            # Measure observables
            if measuring_thermal and i % measuring_interval == 0:
                # Measure observables
                for name, obs in observables.items():
                    obs.measure(get_observable_data(name))

            # Save the universe (optional)
            if save_thermal and i % saving_interval == 0:
                universe.export_geometry(outfile + f"_thermal_{i}", k0=k0)

            # Print the sizes of the observables every 100 sweeps
            if i % 100 == 0:
                for obs, data in observables.items():
                    print(f"{obs} data size: {total_size(data.get_data()) / 1024 / 1024} MB")

            # Garbage collection every 10 sweeps
            if i % 10 == 0:
                gc.collect()

    if sweeps > 0:
        # Main sweeps
        print("========================================\n")
        print("MAIN SWEEPS\n")
        print("----------------------------------------\n")
        print(f"k0 = {k0}, k3 = {k3}, epsilon = {epsilon}\n")
        print("----------------------------------------\n")

        for i in range(1, sweeps + 1):
            # Get the current state of the universe and print it
            n0 = universe.vertex_pool.get_number_occupied()
            n31 = universe.tetras_31.get_number_occupied()
            n3 = universe.tetrahedron_pool.get_number_occupied()
            n22 = universe.tetras_22.get_number_occupied()
            print(f"Main i: {i} target: {target_volume} target2d: {target2_volume} k0: {k0} k3: {k3} \t CURRENT N0: {n0}, N3: {n3}, N31: {n31}, N13: {n3 - n31 - n22}, N22: {n22}\n")

            # Perform sweeps and save general mcmc move stats
            main_successes, main_fails = perform_sweep(k_steps)
            if include_mcmc_data:
                successes.append(main_successes)
                fails.append(main_fails)
                acceptance_ratios.append([get_acceptance_probability(i) for i in range(1, N_MOVES + 1)])

            # Ensure that the universe is at the target volume (if specified)
            if target_volume > 0 and i % measuring_interval == 0:
                # Flag to track if the target volume is reached
                compare = 0

                # Attempt moves until the target volume is reached
                while compare != target_volume:
                    move = choose_move()
                    get_move(move)
                    # Update the compare variable based on the volume switch
                    compare = universe.tetras_31.get_number_occupied() if volfix_switch == 0 else universe.tetrahedron_pool.get_number_occupied()
        
            # Check if there's a 2D target volume specified for the timeslices
            if target2_volume > 0 and i % measuring_interval == 0:
                # Flag to track if the 2D target volume is reached
                hit = False

                # Attempt moves until the 2D target volume is reached
                while not hit:
                    move = choose_move()
                    get_move(move)

                    # Check if any slice matches the 2D target volume
                    for s in universe.slice_sizes:
                        if s == target2_volume:
                            hit = True
                            break

            # Check if the universe geometry is valid
            if validity_check:
                # Update geometry 
                universe.update_vertices()
                universe.log()
                universe.check_validity()

            # Measure observables
            if measuring_main and i % measuring_interval == 0:
                # Measure observables             
                for name, obs in observables.items():
                    obs.measure(get_observable_data(name))

            # Save universe (optional)
            if save_main and i % saving_interval == 0:
                universe.export_geometry(outfile + f"_main_{i}", k0=k0)

            # Print the sizes of the observables every 100 sweeps
            if i % 100 == 0:
                for obs, data in observables.items():
                    print(f"{obs} data size: {total_size(data.get_data()) / 1024 / 1024} MB")

            # Garbage collection every 10 sweeps
            if i % 10 == 0:
                gc.collect()

    # Make success rates, acceptance ratios and k3 values observables
    if include_mcmc_data:
        observables['successes'] = Observable('successes', thermal_sweeps, sweeps, k0, measuring_interval)
        observables['successes'].data = successes
        observables['fails'] = Observable('fails', thermal_sweeps, sweeps, k0, measuring_interval)
        observables['fails'].data = fails
        observables['acceptance_ratios'] = Observable('acceptance_ratios', thermal_sweeps, sweeps,k0, measuring_interval)
        observables['acceptance_ratios'].data = acceptance_ratios
        observables['k3_values'] = Observable('k3_values', thermal_sweeps, sweeps, k0, measuring_interval)
        observables['k3_values'].data = k3_values

    # Save observables
    if measuring_thermal or measuring_main:
        for name, obs in observables.items():
            obs.save_data(outfile + f"_{name}")

    # Print the sizes of the observables final
    for obs, data in observables.items():
        print(f"{obs} data size: {total_size(data.get_data()) / 1024 / 1024} MB") 

    # Close the pool
    pool.close()
    pool.join()

def get_observable_data(name: str):
    if name == 'n_vertices':
        return universe.vertex_pool.get_number_occupied()
    elif name == 'n_tetras':
        return universe.tetrahedron_pool.get_number_occupied()
    elif name == 'n_tetras_31':
        return universe.tetras_31.get_number_occupied()
    elif name == 'n_tetras_22':
        return universe.tetras_22.get_number_occupied()
    elif name == 'slice_sizes':
        return universe.slice_sizes
    elif name == 'slab_sizes':
        return universe.slab_sizes
    elif name == 'curvature':
        return universe.get_curvature_profile()
    elif name == 'connections':
        universe.update_vertices()
        return universe.vertex_neighbours
    else:
        raise ValueError(f"Observable {name} not found.")

def get_acceptance_probability(move: int) -> float:
    """
    Calculates the acceptance probabilities for the different moves.

    Args:
        move (int): The move to calculate the acceptance probability for.
    """
    # Get relevant variables 
    n31 = universe.tetras_31.get_number_occupied()
    n3 = universe.tetrahedron_pool.get_number_occupied()

    # Add
    if move == 1:
        add_ap = (n31 / (n31 + 2)) * np.exp(k0 - 4 * k3)
        # If the target volume is specified, adjust AP according to the volume switch
        if target_volume > 0:
            if volfix_switch == 0:
                add_ap *= np.exp(4 * epsilon * (target_volume - n31 - 1))
            else:
                add_ap *= np.exp(8 * epsilon * (target_volume - n3 - 2))
        return add_ap
    # Delete
    elif move == 2:
        delete_ap = (n31 / (n31 - 2)) * np.exp(-k0 + 4 * k3)
        # If the target volume is specified, adjust AP according to the volume switch
        if target_volume > 0:
            if volfix_switch == 0:
                delete_ap *= np.exp(-4 * epsilon * (target_volume - n31 - 1))
            else:
                delete_ap *= np.exp(-8 * epsilon * (target_volume - n3 - 2))
        return delete_ap
    # Flip
    elif move == 3:
        return 1
    # Shift
    elif move == 4:
        shift_ap = np.exp(-k3)
        # If the target volume is specified, adjust AP according to the volume switch
        if target_volume > 0 and volfix_switch == 1:
            shift_ap *= np.exp(epsilon * (2 * target_volume - 2 * n3 - 1))
        return shift_ap
    # Inverse shift
    elif move == 5:
        ishift_ap = np.exp(k3)
        # If the target volume is specified, adjust AP according to the volume switch
        if target_volume > 0 and volfix_switch == 1:
            ishift_ap *= np.exp(-epsilon * (2 * target_volume - 2 * n3 - 1))
        return ishift_ap
    
def mcmc_check(acceptance_probability: float) -> bool:
    """
    Metropolis-Hastings check for acceptance of a move.

    Args:
        acceptance_probability (float): The acceptance probability of the move.

    Returns:
        bool: True if the move was accepted, False otherwise.
    """
    # Reject
    if acceptance_probability < 1:
        random_number = rng.random()
        if random_number > acceptance_probability:
            return False
        
    # Accept
    return True

def get_move_check(move: str, *args) -> int:
    """
    """
    if move == 'delete':
        return check_delete(*args)
    elif move == 'flip':
        return check_flip(*args)
    elif move == 'shift_u':
        return check_shift_u(*args)
    elif move == 'shift_d':
        return check_shift_d(*args)
    elif move == 'ishift_u':
        return check_ishift_u(*args)
    elif move == 'ishift_d':
        return check_ishift_d(*args)
    
def choose_move():
    # Choose a move with a weight proportional to the acceptance ratio
    moves = {'add': min(get_acceptance_probability(1), 1),
        'delete': min(get_acceptance_probability(2), 1),
        'flip': min(get_acceptance_probability(3), 1),
        'shift_u': min(get_acceptance_probability(4), 1) / 2,
        'shift_d': min(get_acceptance_probability(4), 1) / 2,
        'ishift_u': min(get_acceptance_probability(5), 1) / 2,
        'ishift_d': min(get_acceptance_probability(5), 1) / 2
    }
    
    return rng.choices(list(moves.keys()), weights=list(moves.values()))[0]

def add_task(func, *args):
    """
    Add a new task to the pool.

    Args:
        func (function): The function to execute.
        *args: The arguments to pass to the function.

    Returns:
        Any: The result of the function.
    """
    return pool.apply_async(func, args)

def spawn_move(move: str) -> List[int]:
    """
    Execute all tasks in the pool and collect results.
    """
    if move == 'delete':
        args = (universe.vertex_pool.elements, universe.vertex_pool.size)
    elif move in ['flip', 'shift_u', 'shift_d', 'ishift_u', 'ishift_d']:
        args = (universe.tetras_31.elements, universe.tetrahedron_pool.elements, universe.tetras_31.size)

    results = [add_task(get_move_check, move, *args) for _ in range(n_proposals)]
    output = []

    for result in results:
        try:
            output.append(result.get())
        except Exception as e:
            print(f"Error: {e}")
    
    print(f"Task completed: {output}")
    return output

def get_move(move: str) -> bool:
    """
    Get the move to perform and execute it.

    Args:
        move (str): The move to perform.

    Returns:
        bool: True if the move was executed, False otherwise.
    """
    if move == 'add':
        return move_add()
    elif move == 'delete':
        return move_delete()
    elif move == 'flip':
        return move_flip()
    elif move == 'shift_u':
        return move_shift_u()
    elif move == 'shift_d':
        return move_shift_d()
    elif move == 'ishift_u':
        return move_ishift_u()
    elif move == 'ishift_d':
        return move_ishift_d()

def perform_sweep(n: int) -> Tuple[List[int], List[int]]:
    """
    Perform a sweep of the simulation.

    Args:
        n (int): The number of moves to perform.
        which_sweep (str): The type of sweep ('thermal' or 'main').

    Returns:
        Tuple[List[int], List[int]]: The gathered counts and failed counts.
    """
    move_map = {'add': 0,'delete': 1,'flip': 2,'shift_u': 3,'shift_d': 3,'ishift_u': 4,'ishift_d': 4}
    successes = np.zeros(N_MOVES, dtype=int)
    fails = np.zeros(N_MOVES, dtype=int)

    # Perform n moves
    for i in range(n):
        move = choose_move()
        print(f"\n{i}, MOVE CHOSEN: {move}")
        passed = get_move(move)
        if passed == True:
            print(f"Move {move} passed.")
            successes[move_map[move]] += 1
        else:
            fails[move_map[move]] += 1

    print(f"Successes: {successes}, \nFails: {fails}")
    return successes, fails

def tune() -> float:
    """
    Tunes the k3 parameter of the simulation based on the difference between the
    target volume and the fixvolume.
    """
    global k3
    delta_k3 = 0.000001
    ratio = 100

    # Define different borders based on the target volume (adaptive tuning)
    border_far = target_volume * 0.5
    border_close = target_volume * 0.05
    border_vclose = target_volume * 0.002
    border_vvclose = target_volume * 0.0001

    # Calculate the fixvolume based on the volfix_switch
    fixvolume = 0
    if volfix_switch == 0:
        fixvolume = universe.tetras_31.get_number_occupied()
    else:
        fixvolume = universe.tetrahedron_pool.get_number_occupied()
        
    # Adjust k3 based on the difference between target volume and fixvolume
    if (target_volume - fixvolume) > border_far:
        k3 -= delta_k3 * ratio * 1000
    elif (target_volume - fixvolume) < -border_far:
        k3 += delta_k3 * ratio * 1000
    elif (target_volume - fixvolume) > border_close:
        k3 -= delta_k3 * 1000
    elif (target_volume - fixvolume) < -border_close:
        k3 += delta_k3 * 1000
    elif (target_volume - fixvolume) > border_vclose:
        k3 -= delta_k3 * 100
    elif (target_volume - fixvolume) < -border_vclose:
        k3 += delta_k3 * 100
    elif (target_volume - fixvolume) > border_vvclose:
        k3 -= delta_k3 * 20
    elif (target_volume - fixvolume) < -border_vvclose:
        k3 += delta_k3 * 20
        
def move_add() -> bool:
    """
    Metropolis-Hastings move to add a tetrahedron. This move is always possible.

    Returns:
        bool: True if the move was accepted, False otherwise.
    """
    # Perform MCMC check for acceptance
    if not mcmc_check(get_acceptance_probability(1)):
        return False
    
    # Perform the move
    tetra31_label = universe.tetras_31.pick()
    return universe.add(tetra31_id=tetra31_label, perform=True)

def move_delete() -> bool:
    """
    Metropolis-Hastings move to delete a tetrahedron.

    Returns:
        bool: True if the move was accepted, False otherwise.
    """
    # Perform MCMC check for acceptance
    if not mcmc_check(get_acceptance_probability(2)):
        return False
    
    # Do move in parallel
    output = spawn_move('delete')

    # Filter out the valid proposals, i.e. entries that are not -1
    valid_proposals = [tetra_label for tetra_label in output if tetra_label != -1]
    print(f"Output: {output}, valid proposals: {valid_proposals}")
    # Perform a random move from the valid proposals
    if valid_proposals:
        random_vertex_label = rng.choice(valid_proposals)
        print(f"Random vertex: {random_vertex_label}")
        return universe.delete(vertex_id=random_vertex_label, perform=True)
    
def move_flip() -> bool:
    """
    Metropolis-Hastings move to flip a tetrahedron.
    This move always has an acceptance ratio of 1.
    
    Returns:
        bool: True if the move was accepted, False otherwise.
    """
    # Do move in parallel
    output = spawn_move('flip')

    # Filter out the valid proposals, i.e. entries that are not -1
    valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
    print(f"Output: {output}, valid proposals: {valid_proposals}")
    # Perform a random move from the valid proposals
    if valid_proposals:
        random_tetra012_label, random_tetra230_label = rng.choice(valid_proposals)
        print(f"Random tetra012: {random_tetra012_label}, Random tetra230: {random_tetra230_label}")
        return universe.flip(tetra012_id=random_tetra012_label, tetra230_id=random_tetra230_label, perform=True)
    
    return False

def move_shift_u() -> bool:
    """
    Metropolis-Hastings move to perform a shift move with a (3,1)-
    tetrahedron.

    Returns:
        bool: True if the move was accepted, False otherwise.
    """
    # Perform MCMC check for acceptance
    if not mcmc_check(get_acceptance_probability(4)):
        return False

    # Do move in parallel
    output = spawn_move('shift_u')

    # Filter out the valid proposals, i.e. entries that are not -1
    valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
    print(f"Output: {output}, valid proposals: {valid_proposals}")
    # Perform a random move from the valid proposals
    if valid_proposals:
        random_tetra31_label, random_tetra22_label = rng.choice(valid_proposals)
        print(f"Random tetra31: {random_tetra31_label}, Random tetra22: {random_tetra22_label}")
        return universe.shift_u(tetra31_id=random_tetra31_label, tetra22_id=random_tetra22_label, perform=True)
    
    return False

def move_shift_d() -> bool:
    """
    Metropolis-Hastings move to perform a shift move with a (1,3)-
    tetrahedron.

    Returns:
        bool: True if the move was accepted, False otherwise.
    """
    # Perform MCMC check for acceptance
    if not mcmc_check(get_acceptance_probability(4)):
        return False
    
    # Do move in parallel
    output = spawn_move('shift_d')

    # Filter out the valid proposals, i.e. entries that are not -1
    valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
    print(f"Output: {output}, valid proposals: {valid_proposals}")
    # Perform a random move from the valid proposals
    if valid_proposals:
        random_tetra13_label, random_tetra22_label = rng.choice(valid_proposals)
        print(f"Random tetra13: {random_tetra13_label}, Random tetra22: {random_tetra22_label}")
        return universe.shift_d(tetra13_id=random_tetra13_label, tetra22_id=random_tetra22_label, perform=True)
    
    return False

def move_ishift_u() -> bool:
    """
    Metropolis-Hastings move to perform an inverse shift move with a
    (3,1)-tetrahedron.

    Returns:
        bool: True if the move was accepted, False otherwise.
    """
    # Perform MCMC check for acceptance
    if not mcmc_check(get_acceptance_probability(5)):
        return False
    
    # Do move in parallel
    output = spawn_move('ishift_u')

    # Filter out the valid proposals, i.e. entries that are not -1
    valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
    print(f"Output: {output}, valid proposals: {valid_proposals}")
    # Perform a random move from the valid proposals
    if valid_proposals:
        random_tetra31_label, random_tetra22l_label, random_tetra22r_label = rng.choice(valid_proposals)
        print(f"Random tetra31: {random_tetra31_label}, Random tetra22l: {random_tetra22l_label}, Random tetra22r: {random_tetra22r_label}")
        return universe.ishift_u(tetra31_id=random_tetra31_label, tetra22l_id=random_tetra22l_label, tetra22r_id=random_tetra22r_label, perform=True)

    return False
            
def move_ishift_d() -> bool:
    """
    Metropolis-Hastings move to perform an inverse shift move with a
    (1,3)-tetrahedron.

    Returns:
        bool: True if the move was accepted, False otherwise.
    """
    # Perform MCMC check for acceptance
    if not mcmc_check(get_acceptance_probability(5)):
        return False

    # Do move in parallel
    output = spawn_move('ishift_d')

    # Filter out the valid proposals, i.e. entries that are not -1
    valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
    print(f"Output: {output}, valid proposals: {valid_proposals}")
    # Perform a random move from the valid proposals
    if valid_proposals:
        random_tetra13_label, random_tetra22l_label, random_tetra22r_label = rng.choice(valid_proposals)
        print(f"Random tetra13: {random_tetra13_label}, Random tetra22l: {random_tetra22l_label}, Random tetra22r: {random_tetra22r_label}")
        return universe.ishift_d(tetra13_id=random_tetra13_label, tetra22l_id=random_tetra22l_label, tetra22r_id=random_tetra22r_label, perform=True)
    
    return False

def pick(array, size, rng) -> int:
    """
    Picks a random object from the array.
    """
    if len(array) == 0:
        raise Exception("Empty.")
    else:
        # Pick random element from pool
        random_element = array[rng.randint(0, size - 1)]
        while not random_element or random_element == -1:
            random_element = array[rng.randint(0, size - 1)]

        return random_element

def check_delete(vertex_pool, vertex_pool_size) -> int:
    """
    Helper function to check if a vertex can be deleted.
    """
    # Get a random vertex
    rng = random.Random()
    vertex = pick(vertex_pool, vertex_pool_size, rng)
    t01 = vertex.get_tetra()
    tv01 = t01.get_tetras()[3]
    
    # Get the vertex index in the tetrahedron
    vpos = np.where(t01.get_vertices() == vertex)[0][0]
    assert vpos >= 0

    # Get the base vertices of the two tetrahedra
    v0 = t01.get_vertices()[(vpos + 1) % 3]
    v1 = t01.get_vertices()[(vpos + 2) % 3]
    v2 = t01.get_vertex_opposite(v0)

    # Get all tetrahedra that need to be updated
    t12 = t01.get_tetra_opposite(v0)
    t20 = t01.get_tetra_opposite(v1)
    tv12 = tv01.get_tetra_opposite(v0)
    tv20 = tv01.get_tetra_opposite(v1)

    if (
        t01.is_31()
        and t12.is_31()
        and t20.is_31()
        and tv01.is_13()
        and tv12.is_13()
        and tv20.is_13()
        and v0.scnum >= 4
        and v1.scnum >= 4
        and v2.scnum >= 4
        and vertex.cnum == 6
        and vertex.scnum == 3
    ):
        # print("delete worked: ", vertex_label)
        return vertex.ID
    else:
        # print("delete failed: ", vertex_label)
        return -1
    
def check_flip(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a flip move is possible.
    If possible, returns the labels of the tetrahedra to flip.
    """
    # Pick a random (3,1)-tetrahedron
    rng = random.Random()
    t012_label = pick(tetras_31, tetras_31_size, rng)
    t012 = tetrahedron_pool[t012_label]
    
    # Get random neighbour of t012
    random_neighbour = rng.randint(0, 2)
    t230 = t012.get_tetras()[random_neighbour]
    
    # Get the opposite (1,3)-tetrahedra
    tv012 = t012.get_tetras()[3]
    tv230 = t230.get_tetras()[3]
    
    # # Get apex of the opposing tetrahedra
    # vt = t012.get_vertices()[3]
    # vb = tv012.get_vertices()[0]

    # Get the vertices of the base triangles that are going to be linked
    v1 = t012.get_vertex_opposite_tetra(t230)
    v3 = t230.get_vertex_opposite_tetra(t012)
    
    # Get the remaining base vertices
    v1pos = np.where(t012.get_vertices() == v1)[0][0]
    v2 = t012.get_vertices()[(v1pos + 1) % 3]
    v0 = t012.get_vertices()[(v1pos + 2) % 3]

    # Get opposite neighbouring tetrahedra
    ta01 = t012.get_tetra_opposite(v2)
    # ta12 = t012.get_tetra_opposite(v0)
    ta23 = t230.get_tetra_opposite(v0)
    # ta30 = t230.get_tetra_opposite(v2)
    tva01 = tv012.get_tetra_opposite(v2)
    tva12 = tv012.get_tetra_opposite(v0)
    tva23 = tv230.get_tetra_opposite(v0)
    # tva30 = tv230.get_tetra_opposite(v2)
    
    # Check if the tetrahedron is actually flippable (opposite tetras should also be neighbours)
    if (
        t230.is_31()
        and t012.get_tetras()[3].check_neighbours_tetra(t230.get_tetras()[3])
        and t012.get_tetras()[3].is_13()
        and t230.get_tetras()[3].is_13()
        and t012.is_31()
        and t230.is_31()
        and tv012.is_13()
        and tv230.is_13()
        and tv012.check_neighbours_tetra(tv230)
        and v1 != v3
        and v0.scnum >= 4
        and v2.scnum >= 4
        and ta01 != t230
        and ta23 != t012
        and tva01 != tv230
        and tva23 != tv012
        and not v1.check_vertex_neighbour(v3)
    ):
        return t012_label, t230.ID
    else:
        return -1
        
def check_shift_u(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a shift move is possible.
    If possible, returns the labels of the tetrahedra to shift.
    """
    # Pick a random (3,1)-tetrahedron
    rng = random.Random()
    t31_label = pick(tetras_31, tetras_31_size, rng)
    t31 = tetrahedron_pool[t31_label]

    # Get random neighbour of t31
    random_neighbour = rng.randint(0, 2)
    t22 = t31.get_tetras()[random_neighbour]

    # Get the vertices that will be linked
    v0 = t31.get_vertex_opposite_tetra(t22)
    v1 = t22.get_vertex_opposite_tetra(t31)
    
    # The remaining vertices
    # v3 = t31.get_vertices()[3]
    v0pos = np.where(t31.get_vertices() == v0)[0][0]
    v2 = t31.get_vertices()[(v0pos + 1) % 3]
    v4 = t31.get_vertices()[(v0pos + 2) % 3]

    # Get neighbouring tetrahedra that need to be updated after the move
    ta023 = t31.get_tetra_opposite(v4)
    ta034 = t31.get_tetra_opposite(v2)
    ta123 = t22.get_tetra_opposite(v4)
    # ta124 = t22.get_tetra_opposite(v3)
    ta134 = t22.get_tetra_opposite(v2)
    
    # Check if the move is valid
    if (
        t22.is_22()
        and not ta023.has_vertex(v1)
        and not ta123.has_vertex(v0)
        and not ta034.has_vertex(v1)
        and not ta134.has_vertex(v0)
        and not v0.check_vertex_neighbour(v1)
    ):
        return t31_label, t22.ID
    else:
        return -1
        
def check_shift_d(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a shift move is possible.
    If possible, returns the labels of the tetrahedra to shift.
    """
    # Pick a random (1,3)-tetrahedron
    rng = random.Random()  
    t31_label = pick(tetras_31, tetras_31_size, rng)
    t31 = tetrahedron_pool[t31_label]
    t13 = t31.get_tetras()[3]

    # Get random neighbour of t13
    random_neighbour = rng.randint(1, 3)
    t22 = t13.get_tetras()[random_neighbour]

    # Get the vertices that will be linked
    v0 = t13.get_vertex_opposite_tetra(t22)
    v1 = t22.get_vertex_opposite_tetra(t13)

    # The remaining vertices
    # v3 = t13.get_vertices()[0] # Top
    v0pos = np.where(t31.get_vertices() == v0)[0][0]
    v2 = t31.get_vertices()[(v0pos + 1) % 3]
    v4 = t31.get_vertices()[(v0pos + 2) % 3]

    # Get the neighbouring tetrahedra
    ta023 = t13.get_tetra_opposite(v4)
    ta034 = t13.get_tetra_opposite(v2)
    ta123 = t22.get_tetra_opposite(v4)
    # ta124 = t22.get_tetra_opposite(v3)
    ta134 = t22.get_tetra_opposite(v2)

    # Check if the tetrahedron is actually of type (2,2)
    if (
        t22.is_22()
        and not ta023.has_vertex(v1)
        and not ta123.has_vertex(v0)
        and not ta034.has_vertex(v1)
        and not ta134.has_vertex(v0)
        and not v0.check_vertex_neighbour(v1)
    ):
        return t13.ID, t22.ID
    else:
        return -1

def check_ishift_u(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int, int], int]:
    """
    Helper function to check if an inverse shift move is possible.
    If possible, returns the labels of the tetrahedra to inverse shift.
    """
    # Pick a random (3,1)-tetrahedron
    rng = random.Random()
    t31_label = pick(tetras_31, tetras_31_size, rng)
    t31 = tetrahedron_pool[t31_label]

    # Get random neighbour of t31
    random_neighbour = rng.randint(0, 2)
    t22l = t31.get_tetras()[random_neighbour]
    t22r = t31.get_tetras()[(random_neighbour + 2) % 3]

     # Get the vertices of the interior triangle
    v1 = t31.get_vertices()[3]
    v3 = t22l.get_vertex_opposite_tetra(t31)
    v4 = t31.get_vertex_opposite_tetra(t22l)

    # The remaining vertices
    v4pos = np.where(t31.get_vertices() == v4)[0][0]
    v0 = t31.get_vertices()[(v4pos + 1) % 3]
    v2 = t31.get_vertices()[(v4pos + 2) % 3]

    # Get neighbouring tetrahedra that need to be updated after the move
    ta023 = t22l.get_tetra_opposite(v1)
    ta034 = t22r.get_tetra_opposite(v1)
    ta123 = t22l.get_tetra_opposite(v0)
    ta124 = t31.get_tetra_opposite(v0)
    ta134 = t22r.get_tetra_opposite(v0)
    
    # Count the number of shared vertices between t22l and t22r
    shared_vertices = 0
    for i in range(N_VERTICES_TETRA):
        if t22r.has_vertex(t22l.get_vertices()[i]):
            shared_vertices += 1

    # Make sure the tetra is of type (2,2) and that they are neighbours and have 3 shared vertices
    if (
        t22l.is_22()
        and t22r.is_22()
        and t22l.check_neighbours_tetra(t22r)
        and shared_vertices == 3
        and ta023.has_vertex(v4)
        and ta123.has_vertex(v4)
        and ta034.has_vertex(v2)
        and ta124.has_vertex(v3)
        and ta134.has_vertex(v2)
    ):
        return t31_label, t22l.ID, t22r.ID
    else:
        return -1

def check_ishift_d(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int, int], int]:
    """
    Helper function to check if an inverse shift move is possible.
    If possible, returns the labels of the tetrahedra to inverse shift.
    """
    # Pick a random (1,3)-tetrahedron
    rng = random.Random()
    t31_label = pick(tetras_31, tetras_31_size, rng)
    t31 = tetrahedron_pool[t31_label]
    t13 = t31.get_tetras()[3]

    # Get random (2,2) neighbours of t13
    random_neighbour = rng.randint(0, 2)
    t22l = t13.get_tetras()[1 + random_neighbour]
    t22r = t13.get_tetras()[1 + (random_neighbour + 2) % 3]

    # Get the vertices of the inner triangle
    v1 = t13.get_vertices()[0]
    v3 = t22l.get_vertex_opposite_tetra(t13)
    v4 = t13.get_vertex_opposite_tetra(t22l)

    # Get the remaining vertices
    v4pos = np.where(t31.get_vertices() == v4)[0][0]
    v0 = t31.get_vertices()[(v4pos + 1) % 3]
    v2 = t31.get_vertices()[(v4pos + 2) % 3]

    # Get the neighbouring tetrahedra
    ta023 = t22l.get_tetra_opposite(v1)
    ta034 = t22r.get_tetra_opposite(v1)
    ta123 = t22l.get_tetra_opposite(v0)
    ta124 = t13.get_tetra_opposite(v0)
    ta134 = t22r.get_tetra_opposite(v0)
    
    # Count the number of shared vertices between t22l and t22r
    shared_vertices = 0
    for i in range(N_VERTICES_TETRA):
        if t22r.has_vertex(t22l.get_vertices()[i]):
            shared_vertices += 1

    # Make sure the tetra is of type (2,2) and that they are neighbours and have 3 shared vertices
    if (
        t22l.is_22()
        and t22r.is_22()
        and t22l.check_neighbours_tetra(t22r)
        and shared_vertices == 3
        and ta023.has_vertex(v4)
        and ta123.has_vertex(v4)
        and ta034.has_vertex(v2)
        and ta124.has_vertex(v3)
        and ta134.has_vertex(v2)
    ):
        return t13.ID, t22l.ID, t22r.ID
    else:
        return -1
    

if __name__ == "__main__":
    # global universe
    # global thermal_sweeps
    # global sweeps
    # global target_volume
    # global n_proposals
    # global volfix_switch
    # global epsilon
    # global k0

    universe = Universe(geometry_infilename='../classes/initial_universes/sample-g0-T3.cdt', strictness=3)
    observables = ['n_vertices', 'n_tetras', 'n_tetras_31', 'n_tetras_22', 'slice_sizes', 'slab_sizes', 'curvature', 'connections']

    seed = 0
    k0 = 0.0
    k3 = 0.8
    tune_flag = True
    thermal_sweeps = 1
    sweeps = 0
    k_steps = 50
    volfix_switch = 0
    target_volume = 3000
    epsilon = 0.00005
    measuring_interval = 1
    measuring_thermal = False
    measuring_main = False
    save_main = False
    save_thermal = False
    saving_interval = 100
    validity_check = False
    n_proposals = 8

    # Create observables
    run(
        universe=universe,
        seed=seed,
        k0=k0,
        k3=k3,
        tune_flag=tune_flag,
        thermal_sweeps=thermal_sweeps,
        sweeps=sweeps,
        k_steps=k_steps,
        volfix_switch=volfix_switch,
        target_volume=target_volume,
        epsilon=epsilon,
        measuring_interval=measuring_interval,
        measuring_thermal=measuring_thermal,
        measuring_main=measuring_main,
        save_main=save_main,
        save_thermal=save_thermal,
        saving_interval=saving_interval,
        validity_check=validity_check
    )

    universe.check_validity()