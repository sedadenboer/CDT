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

    def attempt_move(self) -> bool:
        """
        Attempt a move.

        Returns:
            bool: True if move was accepted, False otherwise.
        """
        # Step 1: Initialisation
        cum_freqs = [0, 0]
        tot_freq = 0
        prev_cum_freq = 0

        # Step 2: Calculate cumulative frequencies
        for move_freq in self.move_freqs:
            total_freq += move_freq
            cum_freqs.append(prev_cum_freq + move_freq)
            prev_cum_freq = cum_freqs[-1]
        
        # Step 3: Pick a random number
        move = random.randint(0, tot_freq - 1)
        bin_choice = random.randint(0, 1)

        # Step 4: Determine the move type
        if move < cum_freqs[0]:
            # Choose between add and delete based on bin_choice
            move_type = 1 if bin_choice == 0 else 2
        else:
            # Chooste flip move
            move_type = 3

        # Step 5: Perform the move
        if move_type == 1:
            if self.add_move():
                return 1
        elif move_type == 2:
            if self.delete_move():
                return 2
        elif move_type == 3:
            if self.flip_move():
                return 3
        
        # If no move is executed (due to rejection)
        return 0
            
    def sweep(self) -> None:
        """
        Perform a sweep of the simulation.
        """
        uniform_int = random.randint(0, 3)

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
    
    def single_step(self):
        """
        Perform a single step of the simulation. Chances 
        between adding, deleting and flipping are equal.
        """
        # Pick a random move
        move = random.randint(0, 2)

        # Perform the move
        if move == 0:
            self.add_move()
            self.add_count += 1
        elif move == 1:
            self.delete_move()
            self.delete_count += 1
        elif move == 2:
            self.flip_move()
            self.flip_count += 1

    def progress_universe(self, steps: int):
        """
        Progress the universe by a given number of steps.

        Args:
            steps (int): Number of steps to progress the universe.
        """
        print("Initial number of vertices: %d"%(self.universe.vertex_pool.get_number_occupied()))
        print("Initial number of triangles: %d"%(self.universe.triangle_pool.get_number_occupied()))

        start = time.time()

        for _ in range(steps):
            self.single_step()

        end = time.time()

        print("...")
        print("Progressing the Universe %d steps took %d seconds"%(steps, end-start))
        print("Add count: %d, delete count: %d, flip count: %d"%(self.add_count, self.delete_count, self.flip_count))
        print("Total number of vertices: %d"%(self.universe.vertex_pool.get_number_occupied()))
        print("Total number of triangles: %d"%(self.universe.triangle_pool.get_number_occupied()))

if __name__ == "__main__":
    simulation = Simulation(
        universe=Universe(total_time=20, initial_slice_size=20),
        time_step=0,
        lambd=np.log(2),
        target_volume=0,
        epsilon=0.02
    )

    simulation.progress_universe(100)