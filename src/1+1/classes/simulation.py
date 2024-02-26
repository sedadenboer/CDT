# simulation.py
#
# Author: Seda den Boer
# Date: 02-01-2024
#
# Description:

import random
import numpy as np
import time
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
    def __init__(self, universe: Universe, lambd: float):
        self.universe = universe
        self.lambd = lambd
        self.add_count = 0
        self.delete_count = 0
        self.flip_count = 0
        self.selected_add = 0
        self.selected_delete = 0
        self.selected_flip = 0
        self.current_success = 0
        self.move_freqs = [1, 1]
        self.acceptance_rates = []
        self.volume_changes = []
        self.delete_rates = []
        self.add_rates = []
        self.flip_rates = []
        self.ar_delete = []
        self.ar_add = []
        self.ar_flip = []

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

        return 0

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
        random_number = np.random.uniform()

        # Check if the move is accepted
        if random_number < min_acceptance_ratio:
            return True
        
        return False
    
    def attempt_move(self) -> int:
        """
        Attempt a move.

        Returns:
            int: Move type (0: no move, 1: add, 2: delete, 3: flip).
        """
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
        move = np.random.default_rng().integers(0, tot_freq)
        # Pick a random bin [0,1] to determine between add and delete
        bin_choice = np.random.default_rng().integers(0, 2)

        # Get the acceptance ratio and the object for the move type
        add_ar, add_triangle_id = self.get_acceptance_ratio_and_object(1)
        flip_ar, flip_triangle_id = self.get_acceptance_ratio_and_object(3)
        delete_ar, del_vertex_id = self.get_acceptance_ratio_and_object(2)

        # If there are not enough vertices to delete, pick a new vertex
        vertex_to_delete = self.universe.vertex_pool.get(del_vertex_id)
        slice_size = self.universe.slice_sizes[vertex_to_delete.time]
        while slice_size < 3:
            delete_ar, del_vertex_id = self.get_acceptance_ratio_and_object(2)
            vertex_to_delete = self.universe.vertex_pool.get(del_vertex_id)
            slice_size = self.universe.slice_sizes[vertex_to_delete.time]

        # Save the acceptance ratios for statistics
        self.ar_add.append(add_ar)
        self.ar_delete.append(delete_ar)
        self.ar_flip.append(flip_ar)

        # Choose between add/delete or flip 
        if move < cum_freqs[0]:
            # Choose between add and delete based on bin_choice
            if bin_choice == 0:
                self.selected_add += 1
                if self.mcmc_check(add_ar):
                    # Perform the add move
                    # print(f"Added vertex {add_triangle_id}, at time {self.universe.triangle_pool.get(add_triangle_id).time}")
                    self.universe.insert_vertex(add_triangle_id)
                    # self.universe.print_state()
                    # print()
                    self.add_count += 1
                    return 1
            else:
                self.selected_delete += 1
                if self.mcmc_check(delete_ar):
                    # print(f"Deleted vertex {del_vertex_id}, at time {self.universe.vertex_pool.get(del_vertex_id).time}")
                    self.universe.remove_vertex(del_vertex_id)
                    self.delete_count += 1
                    # self.universe.print_state()
                    # print()
                    return 2
        elif move >= cum_freqs[0]:
            self.selected_flip += 1
            if self.mcmc_check(flip_ar):
                # Perform the flip move
                # print(f"Flipped triangle {flip_triangle_id}, at time {self.universe.triangle_pool.get(flip_triangle_id).time}")
                self.universe.flip_edge(flip_triangle_id)
                self.flip_count += 1
                # self.universe.print_state()
                # print()
                return 3
      
        return 0
    
    def get_total_moves(self) -> int:
        """
        Get the total number of performed moves.

        Returns:
            int: Total number of performed moves.
        """
        return self.add_count + self.delete_count + self.flip_count
    
    def progress_universe(self, steps: int, silence: bool = True):
        """
        Progress the universe by a given number of steps.

        Args:
            steps (int): Number of steps to progress the universe.
            silence (bool, optional): Whether to print progress. Defaults to False.
        """
        if not silence:
            print(f"Vertices: {self.universe.vertex_pool.get_number_occupied()}, Triangles: {self.universe.triangle_pool.get_number_occupied()}")

        start = time.time()

        for step in range(1, steps + 1):
            # Attempt a move
            move_type = self.attempt_move()

            # Check validity of the universe
            self.universe.check_validity()

            # Keep track of the number of successful moves
            if move_type != 0:
                self.current_success += 1
            
            # Save the total size of the universe
            self.volume_changes.append(self.universe.get_total_size())

            # Compute total acceptance rates
            self.acceptance_rates.append(self.current_success / step)

            # Compute individual acceptance rates
            if self.selected_delete != 0:
                self.delete_rates.append(self.delete_count / self.selected_delete)
            else:
                self.delete_rates.append(0)

            if self.selected_add != 0:
                self.add_rates.append(self.add_count / self.selected_add)
            else:
                self.add_rates.append(0)

            if self.selected_flip != 0:
                self.flip_rates.append(self.flip_count / self.selected_flip)
            else:
                self.flip_rates.append(0)

        end = time.time()

        if not silence:
            print("...")
            print(f"Progressing the Universe {steps} steps took {end-start} seconds")
            print(f"Add count: {self.add_count}, delete count: {self.delete_count}, flip count: {self.flip_count}")
            print(f"Ratio delete / add: {self.delete_count / self.add_count:.5f}. Ratio add + delete / flip: {(self.add_count + self.delete_count) / self.flip_count:.5f}")
            print(f"Total number of vertices: {self.universe.vertex_pool.get_number_occupied()}")
            print(f"Total number of triangles: {self.universe.triangle_pool.get_number_occupied()}")
            print(f"Ratio of order 4 vertices and normal vertices: {self.universe.four_vertices_bag.get_number_occupied() / self.universe.vertex_pool.get_number_occupied():.5f}")
            print()


if __name__ == "__main__":
    
    # Set up the universe
    universe = Universe(total_time=20, initial_slice_size=20)
    simulation = Simulation(universe, lambd=np.log(2))

    # Progress the universe
    simulation.progress_universe(10000, silence=False)
