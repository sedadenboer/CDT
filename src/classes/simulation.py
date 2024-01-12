# simulation.py
#
# Author: Seda den Boer
# Date: 02-01-2024
#
# Description:

import random
import numpy as np


class Simulation(object):
    """
    The Simulation class is responsible for all procedures related
    to the actual Monte Carlo simulation. It proposes moves and computes
    the detailed-balance conditions. If it decides a move should be
    accepted, it calls the Universe class to carry out the move at a
    given location. It also triggers the measurement of observables.
    """
    def __init__(self, universe, time_step):
        self.universe = universe
        self.time_step = time_step
        self.current_time = 0
        self.rescaled_cosmological_constant = np.log(2)
        self.add_move_count = 0
        self.delete_move_count = 0
        self.flip_move_count = 0
    
    def weighted_boltzmann_factor(self, rescaled_cosmological_constant: float) -> float:
        return np.exp(rescaled_cosmological_constant)
    
    def accept_move(self, p_move: float) -> bool:
        if random.random() > p_move:
            return False     

    def add_move(self) -> bool:
        """
        Let's consider an example: If split_f is 2, it means that the first two
        future nodes will be separated from the original node. The new number of
        future nodes (n_f_new) will be the original number of future nodes (n_f) minus
        the number of nodes being split (split_f), plus 1 for the new node created.
        The same logic applies to n_p_new for past nodes.

        Returns:
            bool: True if move was accepted, False otherwise.
        """
        vertex = self.universe.get_random_node("add")

        # Get the number of past and future nodes
        n_p = vertex.get_number_of_past_neighbours()
        n_f = vertex.get_number_of_future_neighbours()
        
        # Get the total number of spatial nodes
        N_v = self.universe.get_total_volume()

        # Choose a random split point
        split_f = random.randint(1, n_f)
        split_p = random.randint(1, n_p)
        
        # Compute the new number of past and future nodes
        n_f_new = n_f - split_f + 1
        n_p_new = n_p - split_p + 1

        # Compute the probability of the move
        P_move = (N_v / (N_v + 1)) * (
                 (n_p * n_f) / (n_p_new + n_f_new)) * (
                 self.weighted_boltzmann_factor(self.rescaled_cosmological_constant))

        # Accept or reject the move
        if self.accept_move(P_move):
            new_node = self.universe.create_node()
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