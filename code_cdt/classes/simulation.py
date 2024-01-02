# simulation.py
#
# Author: Seda den Boer
# Date: 02-01-2024
#
# Description:


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

    def compute_detailed_balance(self):
        pass
    
    def propose_move(self):
        pass

    def accept_move(self):
        pass
    
    def measure_observables(self):
        pass

    def update(self):
        pass