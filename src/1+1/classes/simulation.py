# simulation.py
#
# Author: Seda den Boer
# Date: 02-01-2024
#
# Description: Defines the Simulation class,
# which is responsible for all procedures
# related to the actual Monte Carlo simulation.

import random
import numpy as np
import time
import os
from universe import Universe
from typing import Tuple


class Simulation:
    """
    The Simulation class is responsible for all procedures related
    to the actual Monte Carlo simulation. It proposes moves and computes
    the detailed-balance conditions. If it decides a move should be
    accepted, it calls the Universe class to carry out the move at a
    given location. It also triggers the measurement of observables.

    Args (Attributes):
        universe (Universe): Universe object.
        lambd (float): Lambda value.
        seed (int): Seed for random number generator.
        steps (int): Number of steps to progress the universe.
        weighted_moves (bool): Whether to use weighted moves.

    Attributes:
        rng (random.Random): Random number generator.
        move_freqs (List[int]): Frequencies of the different move types.
        volume_changes (List[int]): List of total sizes of the universe.
        ar_delete (List[float]): List of acceptance ratios for delete moves.
        ar_add (List[float]): List of acceptance ratios for add moves.
        ar_flip (List[float]): List of acceptance ratios for flip moves.
        count_add (List[int]): List of counts of add moves.
        count_delete (List[int]): List of counts of delete moves.
        count_flip (List[int]): List of counts of flip moves.
        failed_add (List[int]): List of counts of failed add moves.
        failed_delete (List[int]): List of counts of failed delete moves.
        failed_flip (List[int]): List of counts of failed flip moves.
    """
    N_MOVES = 3
    SLICE_SIZE_LIMIT = 3
    MAX_DELETE_TRIES = 1000

    def __init__(self, universe: Universe, lambd: float, seed: int = 0, steps: int = 1000, weighted_moves: bool = False):
        self.universe = universe
        self.lambd = lambd
        self.seed = seed
        self.rng = random.Random(seed)
        self.steps = steps
        self.weighted_moves = weighted_moves
        self.move_freqs = [1, 1]

        # Lists for statistics
        self.volume_changes = []
        self.ar_delete = []
        self.ar_add = []
        self.ar_flip = []
        self.count_add = [0]
        self.count_delete = [0]
        self.count_flip = [0]
        self.failed_add = [0]
        self.failed_delete = [0]
        self.failed_flip = [0]

    def save_data(self):
        """
        Save all the arrays to a numpy file as well as the geometries.

        Args:
            filename (str): Name of the file to save the data to.
        """
        # For using the lambda value in the filename
        if self.lambd == np.log(2):
            lambd_str = "ln2"
        else:
            lambd_str = str(self.lambd).replace(".", "_")

        # Create a directory for the measurements
        pathname = f'measurements/lambd={lambd_str}/'
        if not os.path.exists(pathname):
            os.makedirs(pathname)

        # Save the data to numpy files
        np.save(pathname + f"volume_changes_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.volume_changes))
        np.save(pathname + f"ar_delete_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.ar_delete))
        np.save(pathname + f"ar_add_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.ar_add))
        np.save(pathname + f"ar_flip_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.ar_flip))
        np.save(pathname + f"count_add_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.count_add))
        np.save(pathname + f"count_delete_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.count_delete))
        np.save(pathname + f"count_flip_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.count_flip))
        np.save(pathname + f"failed_add_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.failed_add))
        np.save(pathname + f"failed_delete_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.failed_delete))
        np.save(pathname + f"failed_flip_stps={self.steps}_lambd={lambd_str}_seed={self.seed}.npy", np.array(self.failed_flip))

        # Save the geometries
        self.universe.save_to_file(pathname + f"geometry_stps={self.steps}_lambd={lambd_str}_seed={self.seed}")

    def get_acceptance_ratio_and_object(self, move_type: int) -> Tuple[float, int]:
        """
        Compute the acceptance ratio for a given move type.

        Args:
            move_type (int): Type of move (1: add, 2: delete, 3: flip).

        Returns:
            float: Acceptance ratio for the given move type.
        """
        # Get the number of vertices and vertices of degree four
        n0 = self.universe.vertex_pool.get_number_occupied()
        n0_four = self.universe.four_vertices_bag.get_number_occupied()

        # Compute acceptance ratios for the different move types
        if move_type == 1:
            # Add
            triangle_id = self.universe.triangle_add_bag.pick()
            return (n0 / (n0_four + 1.0)) * np.exp(-2 * self.lambd), triangle_id
        elif move_type == 2:
            # Delete
            # If there are no vertices of degree four
            if n0_four == 0:
                return 0, 0
            
            vertex_id = self.universe.four_vertices_bag.pick()
            return ((n0_four + 1.0) / n0) * np.exp(2 * self.lambd), vertex_id
        elif move_type == 3:
            # Flip
            n3_flip = self.universe.triangle_flip_bag.get_number_occupied()

            # If there are no triangles to flip
            if n3_flip == 0:
                return 0, 0
            
            # Save the intial number of triangles
            new_n3_flip = n3_flip

            # Flip move needs a specific triangle to compute the acceptance ratio
            triangle_id = self.universe.triangle_flip_bag.pick()
            triangle = self.universe.triangle_pool.get(triangle_id)

            # Compute the new number of triangles that can be flipped
            if triangle.type == triangle.get_triangle_right().get_triangle_right().type:
                new_n3_flip += 1
            else:
                new_n3_flip -= 1

            return (n3_flip / new_n3_flip), triangle_id

        return -1, -1

    def mcmc_check(self, acceptance_ratio: float) -> bool:
        """
        Perform the MCMC check.

        Args:
            acceptance_ratio (float): Acceptance ratio for the move.

        Returns:
            bool: True if the move is accepted, False otherwise.
        """
        # Get the minimum of 1 and the acceptance ratio
        min_acceptance_ratio = min(1, acceptance_ratio)

        # Generate a random uniform number
        random_number = self.rng.random()

        # Check if the move is accepted
        if random_number < min_acceptance_ratio:
            return True
        
        return False
    
    def attempt_move(self) -> int:
        """
        Attempt a move.

        Returns:
            int: Move type (1: add, 2: delete, 3: flip).
        """
        # Get the acceptance ratio and the object for the move type
        add_ar, add_triangle_id = self.get_acceptance_ratio_and_object(1)
        flip_ar, flip_triangle_id = self.get_acceptance_ratio_and_object(3)
        delete_ar, del_vertex_id = self.get_acceptance_ratio_and_object(2)

        # Needed for checking if there are enough vertices to delete
        vertex_to_delete = self.universe.vertex_pool.get(del_vertex_id)
        slice_size = self.universe.slice_sizes[vertex_to_delete.time]
        tries = 0

        # If there are not enough vertices to delete, pick a new vertex
        while slice_size <= self.SLICE_SIZE_LIMIT:
            # If there are no vertices of degree four, return 0
            if tries > self.MAX_DELETE_TRIES:
                return 0
            
            # Get new acceptance ratio and object, and update tries
            delete_ar, del_vertex_id = self.get_acceptance_ratio_and_object(2)
            vertex_to_delete = self.universe.vertex_pool.get(del_vertex_id)
            slice_size = self.universe.slice_sizes[vertex_to_delete.time]
            tries += 1

        # Save the acceptance ratios for statistics
        self.ar_add.append(add_ar)
        self.ar_delete.append(delete_ar)
        self.ar_flip.append(flip_ar)

        # Choose move with probability proportional to acceptance ratios
        if self.weighted_moves:
            weighed_move_choice = self.rng.choices([1, 2, 3], weights=[add_ar, delete_ar, flip_ar])[0]

            if weighed_move_choice == 1:
                # Perform the MCMC check and if passed, insert the vertex
                if self.mcmc_check(add_ar):
                    self.universe.insert_vertex(add_triangle_id)
                    return 1
                else:
                    return -1
            elif weighed_move_choice == 2:
                # Perform the MCMC check and if passed, remove the vertex
                if self.mcmc_check(delete_ar):
                    self.universe.remove_vertex(del_vertex_id)
                    return 2
                else:
                    return -2
            elif weighed_move_choice == 3:
                # Perform the MCMC check and if passed, flip the edge
                if self.mcmc_check(flip_ar):
                    self.universe.flip_edge(flip_triangle_id)
                    return 3
                else:
                    return -3
        else:
            # Use cumulative frequencies to choose between add/delete and flip
            cum_freqs = [0, 0]
            tot_freq = 0
            prev_cum_freq = 0

            # Calculate cumulative frequencies for the move types based on their frequencies
            for i in range(len(self.move_freqs)):
                # Every time a move is attempted, the total frequency is updated
                tot_freq += self.move_freqs[i]
                cum_freqs[i] = prev_cum_freq + self.move_freqs[i]
                prev_cum_freq = cum_freqs[i]

            # Pick a random number and a random bin to determine the move type
            move = self.rng.randint(0, tot_freq)
        
            # Pick a random bin [0,1] to determine between add and delete
            bin_choice = self.rng.randint(0, 1)

            # Choose between add/delete or flip 
            if move < cum_freqs[0]:
                # Choose between add and delete based on bin_choice
                if bin_choice == 0:
                    # Perform the MCMC check and if passed, insert the vertex
                    if self.mcmc_check(add_ar):
                        self.universe.insert_vertex(add_triangle_id)
                        return 1
                    else: 
                        return -1
                else:
                    if self.mcmc_check(delete_ar):
                        # Perform the MCMC check and if passed, remove the vertex
                        self.universe.remove_vertex(del_vertex_id)
                        return 2
                    else:
                        return -2
            elif move >= cum_freqs[0]:
                if self.mcmc_check(flip_ar):
                    # Perform the MCMC check and if passed, flip the edge
                    self.universe.flip_edge(flip_triangle_id)
                    return 3
                else:
                    return -3
      
        return 0
    
    def progress_universe(self, silence: bool = True, save_data: bool = True):
        """
        Progress the universe by a given number of steps.

        Args:
            silence (bool, optional): Whether to print progress. Defaults to False.
            save_data (bool, optional): Whether to save the data. Defaults to True.
        """
        if not silence:
            print(f"Vertices: {self.universe.vertex_pool.get_number_occupied()}, Triangles: {self.universe.triangle_pool.get_number_occupied()}")

        add_count = 0
        delete_count = 0
        flip_count = 0
        add_failed_count = 0
        delete_failed_count = 0
        flip_failed_count = 0

        start = time.time()

        # Progress the universe by a given number of steps
        for _ in range(self.steps):
            # Attempt a move
            move_type = self.attempt_move()

            # Save the total size of the universe
            self.volume_changes.append(self.universe.get_total_size())

            # Update the counts
            if move_type == 1:
                add_count += 1
            elif move_type == 2:
                delete_count += 1
            elif move_type == 3:
                flip_count += 1
            elif move_type == -1:
                add_failed_count += 1
            elif move_type == -2:
                delete_failed_count += 1
            elif move_type == -3:
                flip_failed_count += 1
            
            # Save the counts
            self.count_add.append(add_count)
            self.count_delete.append(delete_count)
            self.count_flip.append(flip_count)
            self.failed_add.append(add_failed_count)
            self.failed_delete.append(delete_failed_count)
            self.failed_flip.append(flip_failed_count)

        end = time.time()

        # Check if the universe is still valid
        self.universe.check_validity()

        if not silence:
            print("...")
            print(f"Progressing the Universe (lambda = {self.lambd}) {self.steps} steps took {end-start} seconds")
            print(f"Add success count: {add_count}, delete success count: {delete_count}, flip success count: {flip_count}")
            print(f"Add failed count: {add_failed_count}, delete failed count: {delete_failed_count}, flip failed count: {flip_failed_count}")
            print(f"Total number of vertices: {self.universe.vertex_pool.get_number_occupied()}")
            print(f"Total number of triangles: {self.universe.triangle_pool.get_number_occupied()}")
            print(f"Ratio of order 4 vertices and normal vertices: {self.universe.four_vertices_bag.get_number_occupied() / self.universe.vertex_pool.get_number_occupied():.5f}")
            print()

        if save_data:
            self.save_data()


if __name__ == "__main__":
    lambd = 0.8
    steps = 10
    weighted_moves = False
    
    universe = Universe(total_time=50, initial_slice_size=40)
    simulation = Simulation(universe, lambd=lambd, seed=42, steps=steps, weighted_moves=weighted_moves)
    simulation.progress_universe(silence=False, save_data=False)
    simulation.universe.check_validity()