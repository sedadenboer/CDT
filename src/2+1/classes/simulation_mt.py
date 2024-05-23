# simulation_mt.py
# 
# Author: Seda den Boer
# Date: 15-05-2024
#
# Descriptiopn: Defines the multiple-try version
# of the original Simulation class.

import random
import numpy as np
import gc
from universe_mt import Universe
from observable import Observable
from typing import TYPE_CHECKING, List, Tuple, Dict, Any, Union
if TYPE_CHECKING:
    from universe import Universe
from helper_functions.helpers import total_size
import multiprocessing
from mpire import WorkerPool
# from check_mt import *
import pickle
import ctypes
N_VERTICES_TETRA = 4
N_VERTICES_TRIANGLE = 3
N_MOVES = 5


class Simulation:
    """
    The Simulation class is responsible for all procedures related
    to the actual Monte Carlo simulation. It proposes moves and computes
    the detailed-balance conditions. If it decides a move should be
    accepted, it calls the Universe class to carry out the move at a
    given location. It also triggers the measurement of observables.

    Args:
        universe (Universe): The universe object representing the system.
        seed (int): The seed for the random number generator.
        k0 (int): The number of k0 moves to perform.
        k3 (int): The number of k3 moves to perform.
        tune_flag (bool): Flag to tune the k3 parameter. Defaults to True.
        thermal_sweeps (int): The number of thermal sweeps to perform.
        sweeps (int): The number of sweeps to perform.
        k_steps (int): The number of k steps to perform.
        v1 (int, optional): The frequency of the v1 move. Defaults to 1.
        v2 (int, optional): The frequency of the v2 move. Defaults to 1.
        v3 (int, optional): The frequency of the v3 move. Defaults to 1.
        volfix_switch (int, optional): The volfix switch. Defaults to 0.
        target_volume (int, optional): The target volume. Defaults to 0.
        target2_volume (int, optional): The second target volume. Defaults to 0.
        epsilon (float, optional): The epsilon parameter. Defaults to 0.00005.
        observables (List[str], optional): The list of observables to measure. Defaults to []. 
                                           Options are: 'n_vertices', 'n_tetras', 'n_tetras_31', 'n_tetras_22',
                                           'slice_sizes', 'slab_sizes','curvature', 'connections'.
        include_mcmc_data (bool): If True, includes successes, fails, acceptance ratios and k3 values to observables. Defaults to True.
        measuring_interval (int, optional): The measuring interval. Defaults to 1.
        measuring_thermal (bool, optional): Flag to measure thermal data. Defaults to False.
        measuring_main (bool, optional): Flag to measure main data. Defaults to False.
        save_main (bool, optional): Flag to save main data. Defaults to False.
        save_thermal (bool, optional): Flag to save thermal data. Defaults to False.
        saving_interval (int, optional): The saving interval. Defaults to 1.
        validity_check (bool, optional): Flag to perform validity check. Defaults to False.
        n_proposals (int): The number of proposals to make in the multiple-try version. Defaults to 5.
    """
    def __init__(self,
                universe: Universe, seed: int,
                k0: int, k3: int, tune_flag: bool = True,
                thermal_sweeps: int = 10, sweeps: int = 10, k_steps: int = 1000,
                v1: int = 1, v2: int = 1, v3: int = 1,
                volfix_switch: int = 0, target_volume: int = 0, target2_volume: int = 0, epsilon: float = 0.00005,
                observables: List[Observable] = [], include_mcmc_data: bool = True,
                measuring_interval: int = 1, measuring_thermal: bool = False, measuring_main: bool = False,
                save_main: bool = False, save_thermal: bool = False, saving_interval: int = 1,
                validity_check: bool = False,
                n_proposals: int = 5
        ):
        self.universe: Universe = universe
        self.rng: random.Random = random.Random(seed)
        self.k0: int = k0
        self.k3: int = k3
        self.tune_flag: bool = tune_flag
        self.thermal_sweeps = thermal_sweeps
        self.sweeps: int = sweeps
        self.k_steps: int = k_steps
        self.volfix_switch: int = volfix_switch
        self.target_volume: int = target_volume
        self.target2_volume: int = target2_volume
        self.epsilon: float = epsilon
        self.validity_check: bool = validity_check
        self.save_main: bool = save_main
        self.save_thermal: bool = save_thermal
        self.saving_interval: int = saving_interval
        self.measuring_interval: int = measuring_interval
        self.measuring_thermal: bool = measuring_thermal
        self.measuring_main: bool = measuring_main
        self.include_mcmc_data: bool = include_mcmc_data
        self.n_proposals: int = n_proposals
        assert self.n_proposals <= multiprocessing.cpu_count(), "Number of proposals should be less than or equal to the number of CPU cores."
        
        # Initialize data structures for MCMC data and set the frequencies of the moves
        self.acceptance_ratios: np.darray = np.zeros(N_MOVES)
        self.move_freqs: Tuple[int, int, int] = (v1, v2, v3)
        self.observables: Dict[str, Observable] = {obs: Observable(obs, thermal_sweeps, sweeps, k0, measuring_interval) for obs in observables}
        
        # If MCMC data is included, save the success and fail counts, acceptance ratios and k3 values
        self.successes: List[np.ndarray] = []
        self.fails: List[np.darray] = []
        self.acceptance_ratios = []
        self.k3_values = []

        # Calculate cumulative frequencies
        self.cum_freqs = np.cumsum(self.move_freqs)
        self.freq_total = sum(self.move_freqs)

        # # Create a pool for parallel execution
        # shared_objects = {
        #     'vertex_pool': self.universe.vertex_pool.elements,
        #     'tetrahedron_pool': self.universe.tetrahedron_pool.elements,
        #     'tetras_31': self.universe.tetras_31.elements
        # }
        # self.pool = WorkerPool(n_jobs=self.n_proposals, shared_objects=shared_objects, start_method='fork')

    def start(self, outfile: str = 'output'):
        """
        Starts the MCMC CDT 2+1 simulation.
        """
        for obs in self.observables.values():
            obs.clear_data()

        # Measure at the start 
        if self.measuring_thermal or self.measuring_main:
            # Measure observables
            for name, obs in self.observables.items():
                obs.measure(self.get_observable_data(name))

        # Thermal sweeps
        if self.thermal_sweeps > 0:
            print("========================================\n")
            print("THERMAL SWEEPS\n")
            print("----------------------------------------\n")
            print(f"k0 = {self.k0}, k3 = {self.k3}, epsilon = {self.epsilon}, thermal = {self.thermal_sweeps}, sweeps = {self.sweeps}, target = {self.target_volume}, target2d = {self.target2_volume}\n")
            print("----------------------------------------\n")

            for i in range(1, self.thermal_sweeps + 1):
                # Get the current state of the universe and print it
                n0 = self.universe.vertex_pool.get_number_occupied()
                n31 = self.universe.tetras_31.get_number_occupied()
                n3 = self.universe.tetrahedron_pool.get_number_occupied()
                n22 = self.universe.tetras_22.get_number_occupied()
                print(f"\nThermal i: {i} \t N0: {n0}, N3: {n3}, N31: {n31}, N13: {n3 - n31 - n22}, N22: {n22}, k0: {self.k0}, k3: {self.k3}")
            
                # Perform sweeps and save general mcmc move stats
                thermal_successes, thermal_fails = self.perform_sweep(self.k_steps)
                if self.include_mcmc_data:
                    self.successes.append(thermal_successes)
                    self.fails.append(thermal_fails)
                    self.acceptance_ratios.append([self.get_acceptance_probability(i) for i in range(1, N_MOVES + 1)])
                    self.k3_values.append(self.k3)
                
                # Tune the k3 parameter
                if self.tune_flag:
                    self.tune()

                # Check if the universe geometry is valid
                if self.validity_check:
                    # Update geometry
                    self.prepare()
                    self.universe.log()
                    self.universe.check_validity()
                
                # Measure observables
                if self.measuring_thermal and i % self.measuring_interval == 0:
                    # Measure observables
                    for name, obs in self.observables.items():
                        obs.measure(self.get_observable_data(name))

                # Save the universe (optional)
                if self.save_thermal and i % self.saving_interval == 0:
                    self.universe.export_geometry(outfile + f"_thermal_{i}", k0=self.k0)

                # Print the sizes of the observables every 100 sweeps
                if i % 100 == 0:
                    for obs, data in self.observables.items():
                        print(f"{obs} data size: {total_size(data.get_data()) / 1024 / 1024} MB")

                # Garbage collection every 10 sweeps
                if i % 10 == 0:
                    gc.collect()

        if self.sweeps > 0:
            # Main sweeps
            print("========================================\n")
            print("MAIN SWEEPS\n")
            print("----------------------------------------\n")
            print(f"k0 = {self.k0}, k3 = {self.k3}, epsilon = {self.epsilon}\n")
            print("----------------------------------------\n")

            for i in range(1, self.sweeps + 1):
                # Get the current state of the universe and print it
                n0 = self.universe.vertex_pool.get_number_occupied()
                n31 = self.universe.tetras_31.get_number_occupied()
                n3 = self.universe.tetrahedron_pool.get_number_occupied()
                n22 = self.universe.tetras_22.get_number_occupied()
                print(f"Main i: {i} target: {self.target_volume} target2d: {self.target2_volume} k0: {self.k0} k3: {self.k3} \t CURRENT N0: {n0}, N3: {n3}, N31: {n31}, N13: {n3 - n31 - n22}, N22: {n22}\n")

                # Perform sweeps and save general mcmc move stats
                main_successes, main_fails = self.perform_sweep(self.k_steps)
                if self.include_mcmc_data:
                    self.successes.append(main_successes)
                    self.fails.append(main_fails)
                    self.acceptance_ratios.append([self.get_acceptance_probability(i) for i in range(1, N_MOVES + 1)])

                # Ensure that the universe is at the target volume (if specified)
                if self.target_volume > 0 and i % self.measuring_interval == 0:
                    # Flag to track if the target volume is reached
                    compare = 0

                    # Attempt moves until the target volume is reached
                    while compare != self.target_volume:
                        self.attempt_move()
                        # Update the compare variable based on the volume switch
                        compare = self.universe.tetras_31.get_number_occupied() if self.volfix_switch == 0 else self.universe.tetrahedron_pool.get_number_occupied()
            
                # Check if there's a 2D target volume specified for the timeslices
                if self.target2_volume > 0 and i % self.measuring_interval == 0:
                    # Flag to track if the 2D target volume is reached
                    hit = False

                    # Attempt moves until the 2D target volume is reached
                    while not hit:
                        self.attempt_move()

                        # Check if any slice matches the 2D target volume
                        for s in self.universe.slice_sizes:
                            if s == self.target2_volume:
                                hit = True
                                break

                # Check if the universe geometry is valid
                if self.validity_check:
                    # Update geometry 
                    self.prepare()
                    self.universe.log()
                    self.universe.check_validity()

                # Measure observables
                if self.measuring_main and i % self.measuring_interval == 0:
                    # Measure observables             
                    for name, obs in self.observables.items():
                        obs.measure(self.get_observable_data(name))

                # Save universe (optional)
                if self.save_main and i % self.saving_interval == 0:
                    self.universe.export_geometry(outfile + f"_main_{i}", k0=self.k0)

                # Print the sizes of the observables every 100 sweeps
                if i % 100 == 0:
                    for obs, data in self.observables.items():
                        print(f"{obs} data size: {total_size(data.get_data()) / 1024 / 1024} MB")

                # Garbage collection every 10 sweeps
                if i % 10 == 0:
                    gc.collect()

        # Make success rates, acceptance ratios and k3 values observables
        if self.include_mcmc_data:
            self.observables['successes'] = Observable('successes', self.thermal_sweeps, self.sweeps, self.k0, self.measuring_interval)
            self.observables['successes'].data = self.successes
            self.observables['fails'] = Observable('fails', self.thermal_sweeps, self.sweeps, self.k0, self.measuring_interval)
            self.observables['fails'].data = self.fails
            self.observables['acceptance_ratios'] = Observable('acceptance_ratios', self.thermal_sweeps, self.sweeps,self.k0, self.measuring_interval)
            self.observables['acceptance_ratios'].data = self.acceptance_ratios
            self.observables['k3_values'] = Observable('k3_values', self.thermal_sweeps, self.sweeps, self.k0, self.measuring_interval)
            self.observables['k3_values'].data = self.k3_values

        # Save observables
        if self.measuring_thermal or self.measuring_main:
            for name, obs in self.observables.items():
                obs.save_data(outfile + f"_{name}")

        # Print the sizes of the observables final
        for obs, data in self.observables.items():
            print(f"{obs} data size: {total_size(data.get_data()) / 1024 / 1024} MB") 

        # # Close the pool
        # self.pool.stop_and_join()

    def get_observable_data(self, name: str) -> Any:
        """
        Gets the observable object by name.
        Options are: 'n_vertices', 'n_tetras', 'n_tetras_31', 'n_tetras_22',
        'slice_sizes', 'slab_sizes','curvature', 'connections'.

        Args:
            name (str): The name of the observable.

        Returns:
            Any: The measured value of the observable.

        Raises:
            ValueError: If the observable is not found.
        """
        if name == 'n_vertices':
            return self.universe.vertex_pool.get_number_occupied()
        elif name == 'n_tetras':
            return self.universe.tetrahedron_pool.get_number_occupied()
        elif name == 'n_tetras_31':
            return self.universe.tetras_31.get_number_occupied()
        elif name == 'n_tetras_22':
            return self.universe.tetras_22.get_number_occupied()
        elif name == 'slice_sizes':
            return self.universe.slice_sizes
        elif name == 'slab_sizes':
            return self.universe.slab_sizes
        elif name == 'curvature':
            return self.universe.get_curvature_profile()
        elif name == 'connections':
            self.prepare()
            return self.universe.vertex_neighbours
        else:
            raise ValueError(f"Observable {name} not found.")
    
    def choose_move(self) -> int:
        """
        Chooses a move based on the acceptance probabilities.

        Returns:
            int: The move to perform.
        """
        # Choose a move with a weight proportional to the acceptance ratio
        moves = {'add': min(self.get_acceptance_probability(1), 1),
            'delete': min(self.get_acceptance_probability(2), 1),
            'flip': min(self.get_acceptance_probability(3), 1),
            'shift_u': min(self.get_acceptance_probability(4), 1) / 2,
            'shift_d': min(self.get_acceptance_probability(4), 1) / 2,
            'ishift_u': min(self.get_acceptance_probability(5), 1) / 2,
            'ishift_d': min(self.get_acceptance_probability(5), 1) / 2
        }
        
        return self.rng.choices(list(moves.keys()), weights=list(moves.values()))[0]
    
    def add_task(self, func, args):
        """
        Add a new task to the pool.

        Args:
            func (function): The function to execute.
            *args: The arguments to pass to the function.

        Returns:
            Any: The result of the function.
        """
        return self.pool.apply_async(func, args)

    def spawn_move(self, move: str) -> List[int]:
        """
        Execute all tasks in the pool and collect results.
        """
        # # Get the required elements and sizes based on the move type
        # if move == 'delete':
        #     # Generate the tasks
        #     results = [self.add_task(check_delete, self.universe.vertex_pool.size) for _ in range(self.n_proposals)]
        # elif move == 'flip':
        #     results = [self.add_task(check_flip, self.universe.tetras_31.size) for _ in range(self.n_proposals)]
        # elif move == 'shift_u':
        #     results = [self.add_task(check_shift_u, self.universe.tetras_31.size) for _ in range(self.n_proposals)]
        # elif move == 'shift_d':
        #     results = [self.add_task(check_shift_d, self.universe.tetras_31.size) for _ in range(self.n_proposals)]
        # elif move == 'ishift_u':
        #     results = [self.add_task(check_ishift_u, self.universe.tetras_31.size) for _ in range(self.n_proposals)]
        # elif move == 'ishift_d':
        #     results = [self.add_task(check_ishift_d, self.universe.tetras_31.size) for _ in range(self.n_proposals)]

        # Create a pool for parallel execution
        shared_objects = {
            'vertex_pool': self.universe.vertex_pool.elements,
            'tetrahedron_pool': self.universe.tetrahedron_pool.elements,
            'tetras_31': self.universe.tetras_31.elements
        }
        with WorkerPool(n_jobs=self.n_proposals, shared_objects=shared_objects, start_method='fork') as pool:
            if move == 'delete':
                # Generate the tasks
                results = pool.map(check_delete, [self.universe.vertex_pool.size for _ in range(self.n_proposals)])
            elif move == 'flip':
                results = pool.map(check_flip, [self.universe.tetras_31.size for _ in range(self.n_proposals)])
            elif move == 'shift_u':
                results = pool.map(check_shift_u, [self.universe.tetras_31.size for _ in range(self.n_proposals)])
            elif move == 'shift_d':
                results = pool.map(check_shift_d, [self.universe.tetras_31.size for _ in range(self.n_proposals)])
            elif move == 'ishift_u':
                results = pool.map(check_ishift_u, [self.universe.tetras_31.size for _ in range(self.n_proposals)])
            elif move == 'ishift_d':
                results = pool.map(check_ishift_d, [self.universe.tetras_31.size for _ in range(self.n_proposals)])
        
        # Generate the tasks
        # output = []

        # # Collect the results
        # for result in results:
        #     try:
        #         output.append(result.get())
        #     except Exception as e:
        #         print(f"Task generated an exception: {e}")
        
        # print(f"{move}: Task completed with result: {output}")
        
        return results
    
    def get_move(self, move: str) -> bool:
        """
        Get the move to perform and execute it.

        Args:
            move (str): The move to perform.

        Returns:
            bool: True if the move was executed, False otherwise.
        """
        if move == 'add':
            return self.move_add()
        elif move == 'delete':
            return self.move_delete()
        elif move == 'flip':
            return self.move_flip()
        elif move == 'shift_u':
            return self.move_shift_u()
        elif move == 'shift_d':
            return self.move_shift_d()
        elif move == 'ishift_u':
            return self.move_ishift_u()
        elif move == 'ishift_d':
            return self.move_ishift_d()

    def perform_sweep(self, n: int) -> Tuple[List[int], List[int]]:
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
            move = self.choose_move()
            # print(f"\n{i}, MOVE CHOSEN: {move}")
            passed = self.get_move(move)
            if passed == True:
                # print(f"Move {move} passed.")
                successes[move_map[move]] += 1
            else:
                fails[move_map[move]] += 1

        print(f"Successes: {successes}, \nFails: {fails}")
        return successes, fails

    def mcmc_check(self, acceptance_probability: float) -> bool:
        """
        Metropolis-Hastings check for acceptance of a move.

        Args:
            acceptance_probability (float): The acceptance probability of the move.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Reject
        if acceptance_probability < 1:
            random_number = self.rng.random()
            if random_number > acceptance_probability:
                return False
            
        # Accept
        return True
    
    def get_acceptance_probability(self, move: int):
        """
        Calculates the acceptance probabilities for the different moves.

        Args:
            move (int): The move to calculate the acceptance probability for.
        """
        # Get relevant variables 
        n31 = self.universe.tetras_31.get_number_occupied()
        n3 = self.universe.tetrahedron_pool.get_number_occupied()

        # Add
        if move == 1:
            add_ap = (n31 / (n31 + 2)) * np.exp(self.k0 - 4 * self.k3)
            # If the target volume is specified, adjust AP according to the volume switch
            if self.volfix_switch == 0:
                if self.target_volume > 0:
                    add_ap *= np.exp(4 * self.epsilon * (self.target_volume - n31 - 1))
            else:
                if self.target_volume > 0:
                    add_ap *= np.exp(8 * self.epsilon * (self.target_volume - n3 - 2))
                    
            return add_ap
        # Delete
        elif move == 2:
            delete_ap = (n31 / (n31 - 2)) * np.exp(-self.k0 + 4 * self.k3)

            # If the target volume is specified, adjust AP according to the volume switch
            if self.volfix_switch == 0:
                if self.target_volume > 0:
                    delete_ap *= np.exp(-4 * self.epsilon * (self.target_volume - n31 - 1))
            else:
                if self.target_volume > 0:
                    delete_ap *= np.exp(-8 * self.epsilon * (self.target_volume - n3 - 2))
    
            return delete_ap
        # Flip
        elif move == 3:
            return 1
        # Shift
        elif move == 4:
            shift_ap = np.exp(-self.k3)
            # If the target volume is specified, adjust AP according to the volume switch
            if self.volfix_switch == 1:
                if self.target_volume > 0:
                        shift_ap *= np.exp(self.epsilon * (2 * self.target_volume - 2 * n3 - 1))

            return shift_ap
        # Inverse shift
        elif move == 5:
            ishift_ap = np.exp(self.k3)
            # If the target volume is specified, adjust AP according to the volume switch
            if self.volfix_switch == 1:
                if self.target_volume > 0:
                    ishift_ap *= np.exp(-self.epsilon * (2 * self.target_volume - 2 * n3 - 1))

            return ishift_ap
        
    def move_add(self) -> bool:
        """
        Metropolis-Hastings move to add a tetrahedron. This move is always possible.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Perform MCMC check for acceptance
        if not self.mcmc_check(self.get_acceptance_probability(1)):
            return False
        
        # Perform the move
        tetra31_label = self.universe.tetras_31.pick()
        return self.universe.add(tetra31_id=tetra31_label, perform=True)

    def move_delete(self) -> bool:
        """
        Metropolis-Hastings move to delete a tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Perform MCMC check for acceptance
        if not self.mcmc_check(self.get_acceptance_probability(2)):
            return False
        
        # Do move in parallel
        output = self.spawn_move('delete')

        # Filter out the valid proposals, i.e. entries that are not -1
        valid_proposals = [tetra_label for tetra_label in output if tetra_label != -1]
        # print(f"Output: {output}, valid proposals: {valid_proposals}")
        # Perform a random move from the valid proposals
        if valid_proposals:
            random_vertex_label = self.rng.choice(valid_proposals)
            # print(f"Random vertex: {random_vertex_label}")
            return self.universe.delete(vertex_id=random_vertex_label, perform=True)
        
    def move_flip(self) -> bool:
        """
        Metropolis-Hastings move to flip a tetrahedron.
        This move always has an acceptance ratio of 1.
        
        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Do move in parallel
        output = self.spawn_move('flip')

        # Filter out the valid proposals, i.e. entries that are not -1
        valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
        # print(f"Output: {output}, valid proposals: {valid_proposals}")
        # Perform a random move from the valid proposals
        if valid_proposals:
            random_tetra012_label, random_tetra230_label = self.rng.choice(valid_proposals)
            # print(f"Random tetra012: {random_tetra012_label}, Random tetra230: {random_tetra230_label}")
            return self.universe.flip(tetra012_id=random_tetra012_label, tetra230_id=random_tetra230_label, perform=True)
        
        return False
    
    def move_shift_u(self) -> bool:
        """
        Metropolis-Hastings move to perform a shift move with a (3,1)-
        tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Perform MCMC check for acceptance
        if not self.mcmc_check(self.get_acceptance_probability(4)):
            return False

        # Do move in parallel
        output = self.spawn_move('shift_u')

        # Filter out the valid proposals, i.e. entries that are not -1
        valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
        # print(f"Output: {output}, valid proposals: {valid_proposals}")
        # Perform a random move from the valid proposals
        if valid_proposals:
            random_tetra31_label, random_tetra22_label = self.rng.choice(valid_proposals)
            # print(f"Random tetra31: {random_tetra31_label}, Random tetra22: {random_tetra22_label}")
            return self.universe.shift_u(tetra31_id=random_tetra31_label, tetra22_id=random_tetra22_label, perform=True)
        
        return False

    def move_shift_d(self) -> bool:
        """
        Metropolis-Hastings move to perform a shift move with a (1,3)-
        tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Perform MCMC check for acceptance
        if not self.mcmc_check(self.get_acceptance_probability(4)):
            return False
        
        # Do move in parallel
        output = self.spawn_move('shift_d')

        # Filter out the valid proposals, i.e. entries that are not -1
        valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
        # print(f"Output: {output}, valid proposals: {valid_proposals}")
        # Perform a random move from the valid proposals
        if valid_proposals:
            random_tetra13_label, random_tetra22_label = self.rng.choice(valid_proposals)
            # print(f"Random tetra13: {random_tetra13_label}, Random tetra22: {random_tetra22_label}")
            return self.universe.shift_d(tetra13_id=random_tetra13_label, tetra22_id=random_tetra22_label, perform=True)
        
        return False

    def move_ishift_u(self) -> bool:
        """
        Metropolis-Hastings move to perform an inverse shift move with a
        (3,1)-tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Perform MCMC check for acceptance
        if not self.mcmc_check(self.get_acceptance_probability(5)):
            return False
        
        # Do move in parallel
        output = self.spawn_move('ishift_u')

        # Filter out the valid proposals, i.e. entries that are not -1
        valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
        # print(f"Output: {output}, valid proposals: {valid_proposals}")
        # Perform a random move from the valid proposals
        if valid_proposals:
            random_tetra31_label, random_tetra22l_label, random_tetra22r_label = self.rng.choice(valid_proposals)
            # print(f"Random tetra31: {random_tetra31_label}, Random tetra22l: {random_tetra22l_label}, Random tetra22r: {random_tetra22r_label}")
            return self.universe.ishift_u(tetra31_id=random_tetra31_label, tetra22l_id=random_tetra22l_label, tetra22r_id=random_tetra22r_label, perform=True)

        return False
                
    def move_ishift_d(self) -> bool:
        """
        Metropolis-Hastings move to perform an inverse shift move with a
        (1,3)-tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Perform MCMC check for acceptance
        if not self.mcmc_check(self.get_acceptance_probability(5)):
            return False

        # Do move in parallel
        output = self.spawn_move('ishift_d')

        # Filter out the valid proposals, i.e. entries that are not -1
        valid_proposals = [tetra_labels for tetra_labels in output if tetra_labels != -1]
        # print(f"Output: {output}, valid proposals: {valid_proposals}")
        # Perform a random move from the valid proposals
        if valid_proposals:
            random_tetra13_label, random_tetra22l_label, random_tetra22r_label = self.rng.choice(valid_proposals)
            # print(f"Random tetra13: {random_tetra13_label}, Random tetra22l: {random_tetra22l_label}, Random tetra22r: {random_tetra22r_label}")
            return self.universe.ishift_d(tetra13_id=random_tetra13_label, tetra22l_id=random_tetra22l_label, tetra22r_id=random_tetra22r_label, perform=True)
        
        return False

    def prepare(self):
        """
        Prepares the universe for measurements in a sweep by updating the geometry.
        """
        self.universe.update_vertices()

    def tune(self):
        """
        Tunes the k3 parameter of the simulation based on the difference between the
        target volume and the fixvolume.
        """
        delta_k3 = 0.000001
        ratio = 100

        # Define different borders based on the target volume (adaptive tuning)
        border_far = self.target_volume * 0.5
        border_close = self.target_volume * 0.05
        border_vclose = self.target_volume * 0.002
        border_vvclose = self.target_volume * 0.0001

        # Calculate the fixvolume based on the self.volfix_switch
        fixvolume = 0
        if self.volfix_switch == 0:
            fixvolume = self.universe.tetras_31.get_number_occupied()
        else:
            fixvolume = self.universe.tetrahedron_pool.get_number_occupied()
            
        # Adjust k3 based on the difference between target volume and fixvolume
        if (self.target_volume - fixvolume) > border_far:
            self.k3 -= delta_k3 * ratio * 1000
        elif (self.target_volume - fixvolume) < -border_far:
            self.k3 += delta_k3 * ratio * 1000
        elif (self.target_volume - fixvolume) > border_close:
            self.k3 -= delta_k3 * 1000
        elif (self.target_volume - fixvolume) < -border_close:
            self.k3 += delta_k3 * 1000
        elif (self.target_volume - fixvolume) > border_vclose:
            self.k3 -= delta_k3 * 100
        elif (self.target_volume - fixvolume) < -border_vclose:
            self.k3 += delta_k3 * 100
        elif (self.target_volume - fixvolume) > border_vvclose:
            self.k3 -= delta_k3 * 20
        elif (self.target_volume - fixvolume) < -border_vvclose:
            self.k3 += delta_k3 * 20


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

def check_delete(shared_objects, vertex_pool_size) -> int:
    """
    Helper function to check if a vertex can be deleted.
    """
    vertex_pool = shared_objects['vertex_pool']
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
    
def check_flip(shared_objects, tetras_31_size) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a flip move is possible.
    If possible, returns the labels of the tetrahedra to flip.
    """
    tetrahedron_pool = shared_objects['tetrahedron_pool']
    tetras_31 = shared_objects['tetras_31']

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
        
def check_shift_u(shared_objects, tetras_31_size) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a shift move is possible.
    If possible, returns the labels of the tetrahedra to shift.
    """
    tetrahedron_pool = shared_objects['tetrahedron_pool']
    tetras_31 = shared_objects['tetras_31']

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
        
def check_shift_d(shared_objects, tetras_31_size) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a shift move is possible.
    If possible, returns the labels of the tetrahedra to shift.
    """
    tetrahedron_pool = shared_objects['tetrahedron_pool']
    tetras_31 = shared_objects['tetras_31']

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

def check_ishift_u(shared_objects, tetras_31_size) -> Union[Tuple[int, int, int], int]:
    """
    Helper function to check if an inverse shift move is possible.
    If possible, returns the labels of the tetrahedra to inverse shift.
    """
    tetrahedron_pool = shared_objects['tetrahedron_pool']
    tetras_31 = shared_objects['tetras_31']

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

def check_ishift_d(shared_objects, tetras_31_size) -> Union[Tuple[int, int, int], int]:
    """
    Helper function to check if an inverse shift move is possible.
    If possible, returns the labels of the tetrahedra to inverse shift.
    """
    tetrahedron_pool = shared_objects['tetrahedron_pool']
    tetras_31 = shared_objects['tetras_31']

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
    universe = Universe(geometry_infilename='../classes/initial_universes/sample-g0-T3.cdt', strictness=3)
    observables = ['n_vertices', 'n_tetras', 'n_tetras_31', 'n_tetras_22', 'slice_sizes', 'slab_sizes', 'curvature', 'connections']

    simulation = Simulation(
        universe=universe,
        seed=0,
        k0=0,
        k3=0.8,
        tune_flag=True,
        thermal_sweeps=10,
        sweeps=0,
        k_steps=300000,
        volfix_switch=0,
        target_volume=3000, # Without tune does not do anything
        observables=observables,
        include_mcmc_data=True,
        measuring_interval=1, # Measure every sweep
        measuring_thermal=False,
        measuring_main=False,
        save_main=False,
        save_thermal=False,
        saving_interval=100, # When to save geometry files
        validity_check=False,
        n_proposals=multiprocessing.cpu_count()
    )

    simulation.start(
        outfile=f'outfile_k0={simulation.k0}_tswps={simulation.thermal_sweeps}_swps={simulation.sweeps}_kstps={simulation.k_steps}_chain={0}'
    )

    simulation.universe.check_validity()
    # simulation.trial()
    
    # simulation.universe.check_validity()
    # observed = simulation.observables

    # for name, obs in observed.items():
    #     print(f"Observable: {name}")
    #     print(f"Thermal: {obs.data}\n")
    #     print()