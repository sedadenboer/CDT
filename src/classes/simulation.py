# simulation.py
#
# Author: Seda den Boer
# Date: 02-01-2024
#
# Description:

import random
import numpy as np
from universe import Universe

class Simulation:
    """
    The Simulation class is responsible for all procedures related
    to the actual Monte Carlo simulation. It proposes moves and computes
    the detailed-balance conditions. If it decides a move should be
    accepted, it calls the Universe class to carry out the move at a
    given location. It also triggers the measurement of observables.
    """
    def __init__(self, universe = Universe(total_time=3, initial_slice_size=3), time_step = 0, rescaled_cosmological_constant = np.log(2), target_volume = 0):
        self.universe = universe
        self.time_step = time_step
        self.rescaled_cosmological_constant = rescaled_cosmological_constant
        self.target_volume = target_volume
        self.current_time = 0
        self.add_move_count = 0
        self.delete_move_count = 0
        self.flip_move_count = 0

    def add_move(self) -> bool:
        """
        Propose a move to add a vertex to the triangulation.

        Returns:
            bool: True if move was accepted, False otherwise.
        """        
        # Compute acceptance ratio
        n0 = self.universe.vertex_bag.get_number_occupied()
        n0_four = self.universe.four_vertices_bag.get_number_occupied()
        acceptance_ratio = n0 / (n0_four + 1.0) * np.exp(-2 * self.rescaled_cosmological_constant)

        # Pick a random vertex from the bag of triangles
        triangle = self.universe.triangle_bag.pick()
        
        # Generate a random number between 0 and 1
        acceptance_threshold = random.uniform(0, 1)
        # If the acceptance criterion is greater than the random number,
        # accept the proposal state as the current state
        if acceptance_ratio > acceptance_threshold:
            # Add the vertex to the triangulation
            self.universe.insert_vertex(triangle)
            return True

        return False

    def delete_move(self):
        # Node requires only 1 future node and 1 past node (?)

        vertex = self.universe.get_random_node("delete")

    def flip_move(self):
        pass

    def make_a_move(self):
        """
        Propose a move and decide whether to accept it or not.
        """
        if random.random() < 0.5:
            if self.add_move():
                self.add_move_count += 1
        else:
            if self.delete_move():
                self.delete_move_count += 1

    def measure_observables(self):
        pass