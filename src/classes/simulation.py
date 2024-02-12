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
from typing import TYPE_CHECKING
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
    def __init__(self, universe: Universe, lambd: float, target_volume: int, epsilon: float):
        self.universe = universe
        self.lambd = lambd
        self.target_volume = target_volume
        self.epsilon = epsilon
        self.current_time = 0
        self.SWEEPSIZE = 100
        self.add_count = 0
        self.delete_count = 0
        self.flip_count = 0
        self.current_success = 0
        self.move_freqs = [1, 1]
        self.acceptance_rates = []
        self.volume_changes = []
        self.delete_rates = []
        self.add_rates = []
        self.flip_rates = []

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

        # Determine the move type
        if move < cum_freqs[0]:
            # Choose between add and delete based on bin_choice
            if bin_choice == 0:
                if self.add_move():
                    self.add_count += 1
                    return 1
            else:
                if self.delete_move():
                    self.delete_count += 1
                    return 2
        elif move >= cum_freqs[0]:
            # Chooste flip move
            if self.flip_move():
                self.flip_count += 1
                return 3

        # If no move is executed (due to rejection)
        return 0
    
    def get_total_moves(self) -> int:
        """
        Get the total number of performed moves.

        Returns:
            int: Total number of performed moves.
        """
        return self.add_count + self.delete_count + self.flip_count

    def add_move(self) -> bool:
        """
        Propose a move to add a vertex to the triangulation.

        Returns:
            bool: True if move was accepted, False otherwise.
        """        
        # Compute acceptance ratio
        n0 = self.universe.vertex_pool.get_number_occupied()
        n0_four = self.universe.four_vertices_bag.get_number_occupied()
        acceptance_ratio = (n0 / (n0_four + 1.0)) * np.exp(-2 * self.lambd)

        # Volume constraint
        if self.target_volume > 0:
            exp_epsilon = np.exp(2 * self.epsilon)
            acceptance_ratio *= self.universe.triangle_pool.get_number_occupied() if self.universe.triangle_pool.get_number_occupied() < self.target_volume else 1 / exp_epsilon

        # MCMC check
        if acceptance_ratio < 1.0:
            random_number = random.uniform(0.0, 1.0)
            if random_number > acceptance_ratio:
                return False

        # Pick a random vertex from the bag of triangles
        triangle_id = self.universe.triangle_add_bag.pick()

        # Add vertex to triangulation
        self.universe.insert_vertex(triangle_id)

        return True

    def delete_move(self) -> bool:
        """
        Propose a move to delete a vertex from the triangulation.

        Returns:
            bool: True if move was accepted, False otherwise.
        """
        # Get the number of vertices and vertices of degree four
        n0 = self.universe.vertex_pool.get_number_occupied()
        n0_four = self.universe.four_vertices_bag.get_number_occupied()

        # If there are no vertices of degree four, return False
        if n0_four == 0:
            return False
        
        # Compute acceptance ratio
        acceptance_ratio = ((n0_four + 1.0) / n0) * np.exp(2 * self.lambd)

        # Volume constraint
        if self.target_volume > 0:
            exp_epsilon = np.exp(2 * self.epsilon)
            acceptance_ratio *= 1 / exp_epsilon if self.universe.triangle_pool.get_number_occupied() < self.target_volume else exp_epsilon
        
        # MCMC check
        if acceptance_ratio < 1.0:
            random_number = random.uniform(0.0, 1.0)
            if random_number > acceptance_ratio:
                return False
        
        # Pick a random vertex from the bag of vertices of degree four
        vertex_id = self.universe.four_vertices_bag.pick()
        vertex = self.universe.vertex_pool.get(vertex_id)

        # Make sure that the slice size is at least 4
        if self.universe.slice_sizes[vertex.time] < 4:
            return False

        # Delete vertex from triangulation
        self.universe.remove_vertex(vertex_id)

        return True

    def flip_move(self) -> bool:
        """
        Propose a move to flip an edge in the triangulation.

        Returns:
            bool: True if move was accepted, False otherwise.
        """
        ntf = self.universe.triangle_flip_bag.get_number_occupied()

        if ntf == 0:
            return False
        
        # Pick a random triangle from the bag of triangles to flip
        triangle_id = self.universe.triangle_flip_bag.pick()
        triangle: Triangle = self.universe.triangle_pool.get(triangle_id)

        # Flip triangle
        tmp = ntf

        if triangle.type == triangle.get_triangle_right().get_triangle_right().type:
            tmp += 1
        else:
            tmp -= 1
        
        acceptance_ratio = 1.0 * ntf / tmp

        # MCMC check
        if acceptance_ratio < 1.0:
            random_number = random.uniform(0.0, 1.0)
            if random_number > acceptance_ratio:
                return False
            
        # Flip edge
        self.universe.flip_edge(triangle_id)

        return True

    def measure_observables(self):
        #TODO
        pass

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
            move_type = self.attempt_move()
            if move_type != 0:
                self.current_success += 1
            self.acceptance_rates.append(self.current_success / step)
            self.volume_changes.append(self.universe.get_total_size())
            self.delete_rates.append(self.delete_count / step)
            self.add_rates.append(self.add_count / step)
            self.flip_rates.append(self.flip_count / step)

        end = time.time()

        if not silence:
            print("...")
            print(f"Progressing the Universe {steps} steps took {end-start} seconds")
            print(f"Rejections: {steps - self.current_success}, Acceptance rate: {self.current_success / steps:.5f}")
            print(f"Add count: {self.add_count}, delete count: {self.delete_count}, flip count: {self.flip_count}")
            print(f"Ratio delete / add: {self.delete_count / self.add_count:.5f}. Ratio add + delete / flip: {(self.add_count + self.delete_count) / self.flip_count:.5f}")
            print(f"Total number of vertices: {self.universe.vertex_pool.get_number_occupied()}")
            print(f"Total number of triangles: {self.universe.triangle_pool.get_number_occupied()}")
            print(f"Ratio of order 4 vertices and normal vertices: {self.universe.four_vertices_bag.get_number_occupied() / self.universe.vertex_pool.get_number_occupied():.5f}")
            print()


if __name__ == "__main__":
    
    # Set up the universe
    universe = Universe(total_time=20, initial_slice_size=20)
    simulation = Simulation(universe, lambd=np.log(2), target_volume=0, epsilon=0.0)

    # Progress the universe
    simulation.progress_universe(10000, silence=False)
