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
from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from triangle import Triangle


class Simulation:
    """
    The Simulation class is responsible for all procedures related
    to the actual Monte Carlo simulation. It proposes moves and computes
    the detailed-balance conditions. If it decides a move should be
    accepted, it calls the Universe class to carry out the move at a
    given location. It also triggers the measurement of observables.
    """
    class Constants:
        N_MOVES = 3
        SLICE_SIZE_LIMIT = 3

    def __init__(self, universe: Universe, lambd: float, seed: int = 0, steps: int = 1000, weighted_moves: bool = False):
        self.universe = universe
        self.lambd = lambd
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
        Save all the arrays to a numpy file.

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
        np.save(pathname + f"volume_changes_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.volume_changes))
        np.save(pathname + f"ar_delete_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.ar_delete))
        np.save(pathname + f"ar_add_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.ar_add))
        np.save(pathname + f"ar_flip_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.ar_flip))
        np.save(pathname + f"count_add_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.count_add))
        np.save(pathname + f"count_delete_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.count_delete))
        np.save(pathname + f"count_flip_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.count_flip))
        np.save(pathname + f"failed_add_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.failed_add))
        np.save(pathname + f"failed_delete_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.failed_delete))
        np.save(pathname + f"failed_flip_stps={self.steps}_lambd={lambd_str}.npy", np.array(self.failed_flip))

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
            triangle_id = self.universe.triangle_add_bag.pick()
            return (n0 / (n0_four + 1.0)) * np.exp(-2 * self.lambd), triangle_id
        elif move_type == 2:
            # If there are no vertices of degree four
            if n0_four == 0:
                return 0, 0
            
            vertex_id = self.universe.four_vertices_bag.pick()
            return ((n0_four + 1.0) / n0) * np.exp(2 * self.lambd), vertex_id
        elif move_type == 3:
            ntf = self.universe.triangle_flip_bag.get_number_occupied()

            # If there are no triangles to flip
            if ntf == 0:
                return 0, 0
            
            new_ntf = ntf

            # Flip move needs a specific triangle to compute the acceptance ratio
            triangle_id = self.universe.triangle_flip_bag.pick()
            triangle: Triangle = self.universe.triangle_pool.get(triangle_id)

            if triangle.type == triangle.get_triangle_right().get_triangle_right().type:
                new_ntf += 1
            else:
                new_ntf -= 1

            return (ntf / new_ntf), triangle_id

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

        vertex_to_delete = self.universe.vertex_pool.get(del_vertex_id)
        slice_size = self.universe.slice_sizes[vertex_to_delete.time]

        # If there are not enough vertices to delete, pick a new vertex
        while slice_size < self.Constants.SLICE_SIZE_LIMIT:
            delete_ar, del_vertex_id = self.get_acceptance_ratio_and_object(2)
            vertex_to_delete = self.universe.vertex_pool.get(del_vertex_id)
            slice_size = self.universe.slice_sizes[vertex_to_delete.time]

        # Save the acceptance ratios for statistics
        self.ar_add.append(add_ar)
        self.ar_delete.append(delete_ar)
        self.ar_flip.append(flip_ar)

        # Choose move based on acceptance ratios
        if self.weighted_moves:
            weighed_move_choice = self.rng.choices([1, 2, 3], weights=[add_ar, delete_ar, flip_ar])[0]
            if weighed_move_choice == 1:
                if self.mcmc_check(add_ar):
                    # Perform the add move
                    self.universe.insert_vertex(add_triangle_id)
                    return 1
                else:
                    return -1
            elif weighed_move_choice == 2:
                if self.mcmc_check(delete_ar):
                    # Perform the delete move
                    self.universe.remove_vertex(del_vertex_id)
                    return 2
                else:
                    return -2
            elif weighed_move_choice == 3:
                if self.mcmc_check(flip_ar):
                    # Perform the flip move
                    self.universe.flip_edge(flip_triangle_id)
                    return 3
                else:
                    return -3
        else:
            # Two bins for add/delete, and one for flip
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
                    if self.mcmc_check(add_ar):
                        # Perform the add move
                        self.universe.insert_vertex(add_triangle_id)
                        return 1
                    else: 
                        return -1
                else:
                    if self.mcmc_check(delete_ar):
                        # Perform the delete move
                        self.universe.remove_vertex(del_vertex_id)
                        return 2
                    else:
                        return -2
            elif move >= cum_freqs[0]:
                if self.mcmc_check(flip_ar):
                    # Perform the flip move
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

        for _ in range(self.steps):
            # Attempt a move
            move_type = self.attempt_move()

            # Save the total size of the universe
            self.volume_changes.append(self.universe.get_total_size())

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
            print(f"Progressing the Universe {self.steps} steps took {end-start} seconds")
            print(f"Add success count: {add_count}, delete success count: {delete_count}, flip success count: {flip_count}")
            print(f"Add failed count: {add_failed_count}, delete failed count: {delete_failed_count}, flip failed count: {flip_failed_count}")
            print(f"Total number of vertices: {self.universe.vertex_pool.get_number_occupied()}")
            print(f"Total number of triangles: {self.universe.triangle_pool.get_number_occupied()}")
            print(f"Ratio of order 4 vertices and normal vertices: {self.universe.four_vertices_bag.get_number_occupied() / self.universe.vertex_pool.get_number_occupied():.5f}")
            print()

        if save_data:
            self.save_data()


if __name__ == "__main__":
    # Set up the universe
    universe = Universe(total_time=20, initial_slice_size=20)
    simulation = Simulation(universe, lambd=np.log(2), seed=42, steps=1000, weighted_moves=True)

    # Progress the universe
    simulation.progress_universe(silence=False, save_data=True)
