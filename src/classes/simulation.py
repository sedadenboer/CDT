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
    def __init__(self, universe: Universe, time_step: int, lambd: float, target_volume: int, epsilon: float):
        self.universe = universe
        self.time_step = time_step
        self.lambd = lambd
        self.target_volume = target_volume
        self.epsilon = epsilon
        self.current_time = 0
        self.move_freqs = [1, 1]
        self.SWEEPSIZE = 100
        self.add_count = 0
        self.delete_count = 0
        self.flip_count = 0
        self.succes_count = 0

    def attempt_move(self) -> bool:
        """
        Attempt a move.

        Returns:
            bool: True if move was accepted, False otherwise.
        """
        # Initialisation
        cum_freqs = [0, 0]
        tot_freq = 0
        prev_cum_freq = 0

        # Calculate cumulative frequencies for the move types
        for i in range(len(self.move_freqs)):
            tot_freq += self.move_freqs[i]
            cum_freqs[i] = prev_cum_freq + self.move_freqs[i]
            prev_cum_freq = cum_freqs[i]

        # Pick a random number and a random bin to determine the move type
        move = np.random.default_rng().integers(0, tot_freq)
        bin_choice = np.random.default_rng().integers(0, 2)

        # Determine the move type
        if move < cum_freqs[0]:
            # Choose between add and delete based on bin_choice
            if bin_choice == 0:
                if self.add_move():
                    self.add_count += 1
                    self.succes_count += 1
                    return 1
                elif self.delete_move():
                    self.delete_count += 1
                    self.succes_count += 1
                    return 2
        elif move >= cum_freqs[0]:
            # Chooste flip move
            if self.flip_move():
                self.flip_count += 1
                self.succes_count += 1
                return 3
        
        # If no move is executed (due to rejection)
        return 0
            
    def sweep(self) -> None:
        """
        Perform a sweep of the simulation.
        """
        # Pick random uniform number
        uniform_int = random.uniform(0, 3)

        # Step 1: Perform moves
        moves = [0, 0, 0, 0]
        for _ in range(self.SWEEPSIZE * self.target_volume):
            moves[self.attempt_move()] += 1

        # Make sure that the triangle volume is equal to the target volume
        while self.universe.triangle_pool.get_number_occupied() != self.target_volume:
            self.attempt_move()
        
        # Step 2: Measure observables
        self.measure_observables()

    def add_move(self) -> bool:
        """
        Propose a move to add a vertex to the triangulation.

        Returns:
            bool: True if move was accepted, False otherwise.
        """        
        # Compute acceptance ratio
        n0 = self.universe.vertex_pool.get_number_occupied()
        n0_four = self.universe.four_vertices_bag.get_number_occupied()
        acceptance_ratio = n0 / (n0_four + 1.0) * np.exp(-2 * self.lambd)

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
        acceptance_ratio = n0_four / (n0 - 1.0) * np.exp(2 * self.lambd)

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
        pass

    def progress_universe(self, steps: int, silence: bool = False):
        """
        Progress the universe by a given number of steps.

        Args:
            steps (int): Number of steps to progress the universe.
            silence (bool, optional): Whether to print progress. Defaults to False.
        """
        if not silence:
            print(f"Initial number of vertices: {self.universe.vertex_pool.get_number_occupied()}")
            print(f"Initial number of triangles: {self.universe.triangle_pool.get_number_occupied()}")

        start = time.time()

        for _ in range(steps):
            self.attempt_move()

        end = time.time()

        if not silence:
            print("...")
            print(f"Progressing the Universe {steps} steps took {end-start} seconds")
            print(f"Rejected: {steps - self.succes_count}. Success rate: {self.succes_count / steps:.5f}")
            print(f"Add count: {self.add_count}, delete count: {self.delete_count}, flip count: {self.flip_count}")
            print(f"Ratio delete / add: {self.delete_count / self.add_count:.5f}. Ratio add + delete / flip: {(self.add_count + self.delete_count) / self.flip_count:.5f}")
            print(f"Total number of vertices: {self.universe.vertex_pool.get_number_occupied()}")
            print(f"Total number of triangles: {self.universe.triangle_pool.get_number_occupied()}")
            print(f"Number of order 4 vertices: {self.universe.four_vertices_bag.get_number_occupied()}")
            print(f"Ratio of order 4 vertices and normal vertices: {self.universe.four_vertices_bag.get_number_occupied() / self.universe.vertex_pool.get_number_occupied():.5f}")
            print()
            

if __name__ == "__main__":
    
    # Set up the universe
    universe = Universe(total_time=20, initial_slice_size=20)
    simulation = Simulation(universe, time_step=1, lambd=np.log(2), target_volume=0, epsilon=0.0)

    # Progress the universe
    simulation.progress_universe(100000, silence=False)
