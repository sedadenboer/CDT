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
from universe import Universe
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from universe import Universe
    from observable import Observable


class Simulation:
    """
    The Simulation class is responsible for all procedures related
    to the actual Monte Carlo simulation. It proposes moves and computes
    the detailed-balance conditions. If it decides a move should be
    accepted, it calls the Universe class to carry out the move at a
    given location. It also triggers the measurement of observables.
    """
    class constants:
        N_VERTICES_TETRA = 4
        N_VERTICES_TRIANGLE = 3

    def __init__(self, universe: Universe):
        self.universe = universe 
        self.k0 = 0
        self.k3 = 0
        self.target_volume = 0
        self.target2_volume = 0
        self.rng = random.Random(0)
        self.epsilon = 0.00004
        self.measuring = False
        self.observables_3d = []
        self.observables_2d = []
        self.move_freqs = (1, 1, 1)
        self.k3_values = []
        self.total_vertices = []
        self.total_tetras = []

    def start(self, k0: int, k3: int,
                    sweeps: int, thermal_sweeps: int, k_steps: int,
                    target_volume: int, target2_volume: int,
                    seed: int, outfile: str,
                    v1: int, v2: int, v3: int):
        """
        Starts the MCMC CDT 2+1 simulation.

        Args:
            k0 (int): The k0 parameter.
            k3 (int): The k3 parameter.
            sweeps (int): The number of sweeps to perform.
            thermal_sweeps (int): The number of thermal sweeps to perform.
            k_steps (int): The number of steps to perform in each sweep.
            target_volume (int): The target volume of the simulation.
            target2_volume (int): The target 2D volume of the simulation.
            seed (int): The seed for the random number generator.
            outfile (str): The output file for the simulation.
            v1 (int): The first move frequency.
            v2 (int): The second move frequency.
            v3 (int): The third move frequency.
        """
        self.move_freqs = (v1, v2, v3)
        self.k0 = k0
        self.k3 = k3
        self.target_volume = target_volume
        self.target2_volume = target2_volume
        self.rng.seed(seed)
        self.measuring = True

        # Clear observables
        for o in self.observables_3d:
            o.clear()
        for o in self.observables_2d:
            o.clear()
        
        # Thermal sweeps
        print("========================================\n")
        print("THERMAL SWEEPS\n")
        print("----------------------------------------\n")
        print(f"k0 = {k0}, k3 = {k3}, epsilon = {self.epsilon}, thermal = {thermal_sweeps}, sweeps = {sweeps}, target = {target_volume}, target2d = {target2_volume}\n")
        print("----------------------------------------\n")

        for i in range(1, thermal_sweeps + 1):
            # Get the current state of the universe and print it
            total_2v = sum(self.universe.slice_sizes)
            n31 = self.universe.tetras_31.get_number_occupied()
            n3 = self.universe.tetrahedron_pool.get_number_occupied()
            print(f"Thermal: i: {i} \t tetra size: {n3} tetras31 size: {n31} k3: {self.k3} \n")

            # Perform sweeps and tune the k3 parameter
            self.perform_sweep(k_steps * 1000)
            self.tune()

            self.total_vertices.append(self.universe.vertex_pool.get_number_occupied())
            self.total_tetras.append(self.universe.tetrahedron_pool.get_number_occupied())

            # # Export the geometry every 10% of the thermal sweeps
            # if i % (thermal_sweeps / 10) == 0:
            #     self.universe.export_geometry(outfile)
            
            # # Update geometry and measure observables related to 3D structures
            # self.prepare()
            # for o in self.observables_3d:
            #     o.measure(self.universe)

        # Main sweeps
        print("========================================\n")
        print("MAIN SWEEPS\n")
        print("----------------------------------------\n")
        print(f"k0 = {k0}, k3 = {k3}, epsilon = {self.epsilon}\n")
        print("----------------------------------------\n")

        for i in range(1, sweeps + 1):
            # Get the current state of the universe and print it
            total_2v = sum(self.universe.slice_sizes)
            n31 = self.universe.tetras_31.get_number_occupied()
            n3 = self.universe.tetrahedron_pool.get_number_occupied()
            avg_2v = total_2v / self.universe.n_slices
            print(f"Main: i: {i} \t target: {self.target_volume} \t target2d: {self.target_2_volume} \t CURRENT n3: {n3} avgslice: {avg_2v}\n")

            # Perform sweeps and tune the k3 parameter
            self.perform_sweep(k_steps * 1000)

            self.total_vertices.append(self.universe.vertex_pool.get_number_occupied())
            self.total_tetras.append(self.universe.tetrahedron_pool.get_number_occupied())

            # Export the geometry every 10% of the main sweeps
            if i % (sweeps / 10) == 0:
                self.universe.export_geometry(outfile)
            
            # Check if there are any 3D observables to be measured
            if len(self.observables_3d) > 0: 
                # Remove if volume should fluctuate during measurements
                vol_switch = self.universe.volfix_switch

                # Ensure that the universe is at the target volume (if specified)
                if target_volume > 0:
                    # Flag to track if the target volume is reached
                    compare = 0

                    # Attempt moves until the target volume is reached
                    while compare != target_volume:
                        self.attempt_move()
                        # Update the compare variable based on the volume switch
                        compare = self.universe.tetras31.get_number_occupied() \
                            if vol_switch == 0 else self.universe.tetrahedron_pool.get_number_occupied()

                # # Update geometry and measure observables related to 3D structures
                # self.prepare()
                # for o in self.observables_3d:
                #     o.measure(self.universe)
            
            # Check if there's a 2D target volume specified for the timeslices
            if target2_volume > 0:
                # Flag to track if the 2D target volume is reached
                hit = False

                # Attempt moves until the 2D target volume is reached
                while not hit:
                    self.attempt_move()

                    # Check if any slice matches the 2D target volume
                    for s in self.universe.slice_sizes:
                         if s == target2_volume:
                             hit = True
                             break
                
                # # Update geometry and measure observables related to 2D structures
                # self.prepare()
                # for o in self.observables_2d:
                #     o.measure(self.universe)

    def attempt_move(self) -> int:
        """
        Attempt a move based on the move frequencies.

        Returns:
            int: The move number. 
            1 for add, 2 for delete, 3 for flip, 4 for shift, 5 for inverse shift.
            Negative numbers indicate failed moves.
        """
        cum_freqs = np.cumsum(self.move_freqs)
        freq_total = sum(self.move_freqs)

        # Generate a random number to select a move pair
        move = self.rng.randint(0, freq_total)

        if move < cum_freqs[0]:
            # Add or delete move
            if self.rng.randint(0, 1) == 0:
                return 1 if self.move_add() else -1
            else:
                return 2 if self.move_delete() else -2
        elif move < cum_freqs[1]:
            # Flip move
            return 3 if self.move_flip() else -3
        else:
            # Shift or inverse shift move
            if self.rng.randint(0, 1) == 0:
                # Shift (3,1) or (1,3) move
                if self.rng.randint(0, 1) == 0:
                    return 4 if self.move_shift_u() else -4
                else:
                    return 4 if self.move_shift_d() else -4
            else:
                # Inverse shift (3,1) or (1,3) move
                if self.rng.randint(0, 1) == 0:
                    return 5 if self.move_ishift_u() else -5
                else:
                    return 5 if self.move_ishift_d() else -5
                
    def perform_sweep(self, n: int) -> list[int]:
        """
        Perform a sweep of the simulation.

        Args:
            n (int): The number of moves to perform.

        Returns:
            list[int]: The number of moves performed for each move type.
        """
        moves = [0, 0, 0, 0, 0, 0]
        failed_moves = [0, 0, 0, 0, 0, 0]

        # Perform n moves
        for _ in range(n):
            move_num = self.attempt_move()
            move = abs(move_num)

            # Update the move counters
            moves[move] += 1
            if move_num < 0:
                failed_moves[move] += 1

        # Print the succes rate for each move
        print(f"Add: {moves[1]} \t Delete: {moves[2]} \t Flip: {moves[3]} \t Shift: {moves[4]} \t Inverse Shift: {moves[5]}")
        print(f"Failed Add: {failed_moves[1]} \t Failed Delete: {failed_moves[2]} \t Failed Flip: {failed_moves[3]} \t Failed Shift: {failed_moves[4]} \t Failed Inverse Shift: {failed_moves[5]}")

        return moves

    def mcmc_check(self, acceptance_probability: float) -> bool:
        """
        Metropolis-Hastings check for acceptance of a move.

        Args:
            acceptance_probability (float): The acceptance probability of the move.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Get the minimum between 1 and AP and generate random uniform number
        min_acceptance_ratio = min(1, acceptance_probability)
        random_number = self.rng.random()

        # Reject
        if random_number > min_acceptance_ratio:
            return False
        
        # Accept
        return True
    
    def move_add(self) -> bool:
        """
        Metropolis-Hastings move to add a tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Get relevant variables for calculating the add AP and calculate it
        n31 = self.universe.tetras_31.get_number_occupied()
        n3 = self.universe.tetrahedron_pool.get_number_occupied()
        n0 = self.universe.vertex_pool.get_number_occupied()
        vol_switch = self.universe.volfix_switch
        acceptance_probability = (n31 / (n0 + 1)) * np.exp(self.k0 - 4 * self.k3)

        # If the target volume is specified, adjust AP according to the volume switch
        if self.target_volume > 0:
            if vol_switch == 0:
                acceptance_probability *= np.exp(
                    4 * self.epsilon * (self.target_volume - n31 - 1))
            else:
                acceptance_probability *= np.exp(
                    8 * self.epsilon * (self.target_volume - n3 - 2))

        # Perform MCMC check for acceptance
        if not self.mcmc_check(acceptance_probability):
            return False
        
        # Perform the move
        tetra31_label = self.universe.tetras_31.pick()
        self.universe.add(tetra31_label)

        return True
    
    def move_delete(self) -> bool:
        """
        Metropolis-Hastings move to delete a tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Get relevant variables for calculating the delete AP and calculate it
        n31 = self.universe.tetras_31.get_number_occupied()
        n3 = self.universe.tetrahedron_pool.get_number_occupied()
        n0 = self.universe.vertex_pool.get_number_occupied()
        vol_switch = self.universe.volfix_switch
        acceptance_probability = ((n0 + 1) / n31) * np.exp(
            -self.k0 + 4 * self.k3)

        # If the target volume is specified, adjust AP according to the volume switch
        if self.target_volume > 0:
            if vol_switch == 0:
                acceptance_probability *= np.exp(
                    -4 * self.epsilon * (self.target_volume - n31 - 1))
            else:
                acceptance_probability *= np.exp(
                    -8 * self.epsilon * (self.target_volume - n3 - 2))

        # Perform MCMC check for acceptance
        if not self.mcmc_check(acceptance_probability):
            return False

        # Get a random vertex
        vertex_label = self.universe.vertex_pool.pick()
        vertex = self.universe.vertex_pool.get(vertex_label)

        # Check if the vertex is actually deletable
        if vertex.cnum != 6:
            return False
        if vertex.scnum != 3:
            return False
        
        # Perform the move
        self.universe.delete(vertex_label)

        # print(f"Deleted vertex: {vertex_label} \n")

        return True
    
    def move_flip(self) -> bool:
        """
        Metropolis-Hastings move to flip a tetrahedron.
        This move has an acceptance probability of 1.
        
        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        tetra012_label = self.universe.tetras_31.pick()
        tetra012 = self.universe.tetrahedron_pool.get(tetra012_label)
        
        # Get random neighbour of tetra012
        random_neighbour = self.rng.randint(0, 2)
        tetra230 = tetra012.get_tetras()[random_neighbour]

        assert tetra012.check_neighbours_tetra(tetra230), "SIM: Tetra012 and tetra230 are not neighbours"

        # Check if the tetrahedron is actually flippable (opposite tetras should also be neighbours)
        if not tetra230.is_31():
            return False
        if not tetra012.get_tetras()[3].check_neighbours_tetra(tetra230.get_tetras()[3]):
            return False
        
        # Try the move
        return self.universe.flip(tetra012_label, tetra230.ID)

    def move_shift_u(self) -> bool:
        """
        Metropolis-Hastings move to perform a shift move with a (3,1)-
        tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        acceptance_probability = np.exp(-self.k3)
        n3 = self.universe.tetrahedron_pool.get_number_occupied()
        vol_switch = self.universe.volfix_switch

        if vol_switch == 1 and self.target_volume > 0:
            acceptance_probability *= np.exp(
                self.epsilon * (2 * self.target_volume - 2 * n3 -1))

        # Perform MCMC check for acceptance
        if not self.mcmc_check(acceptance_probability):
            return False
        
        # Pick a random (3,1)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra31 = self.universe.tetrahedron_pool.get(tetra31_label)

        # Get random neighbour of tetra31
        random_neighbour = self.rng.randint(0, 2)
        tetra22 = tetra31.get_tetras()[random_neighbour]
        
        assert tetra22.check_neighbours_tetra(tetra31), "SIM: Tetra22 and tetra31 are not neighbours"

        # Check if the tetrahedron is actually of type (2,2)
        if not tetra22.is_22():
            return False
        
        # Try the move
        return self.universe.shift_u(tetra31_label, tetra22.ID)
    
    def move_shift_d(self) -> bool:
        """
        Metropolis-Hastings move to perform a shift move with a (1,3)-
        tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        acceptance_probability = np.exp(-self.k3)
        n3 = self.universe.tetrahedron_pool.get_number_occupied()
        vol_switch = self.universe.volfix_switch

        if vol_switch == 1 and self.target_volume > 0:
            acceptance_probability *= np.exp(
                self.epsilon * (2 * self.target_volume - 2 * n3 -1))

        # Perform MCMC check for acceptance
        if not self.mcmc_check(acceptance_probability):
            return False
        
        # Pick a random (1,3)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra13 = self.universe.tetrahedron_pool.get(tetra31_label).get_tetras()[3]

        # Get random neighbour of tetra13
        random_neighbour = self.rng.randint(1, 3)
        tetra22 = tetra13.get_tetras()[random_neighbour]
        
        # Check if the tetrahedron is actually of type (2,2)
        if not tetra22.is_22() or not tetra13.is_13():
            return False
        
        assert tetra22.check_neighbours_tetra(tetra13), "SIM: Tetra22 and tetra13 are not neighbours"
        
        # Try the move
        return self.universe.shift_d(tetra13.ID, tetra22.ID)

    def move_ishift_u(self) -> bool:
        """
        Metropolis-Hastings move to perform an inverse shift move with a
        (3,1)-tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        acceptance_probability = np.exp(self.k3)
        n3 = self.universe.tetrahedron_pool.get_number_occupied()
        vol_switch = self.universe.volfix_switch

        if vol_switch == 1 and self.target_volume > 0:
            acceptance_probability *= np.exp(
                -self.epsilon * (2 * self.target_volume - 2 * n3 -1))

        # Perform MCMC check for acceptance
        if not self.mcmc_check(acceptance_probability):
            return False
        
        # Pick a random (3,1)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra31 = self.universe.tetrahedron_pool.get(tetra31_label)

        # Get random neighbour of tetra31
        random_neighbour = self.rng.randint(0, 2)
        tetra22l = tetra31.get_tetras()[random_neighbour]
        tetra22r = tetra31.get_tetras()[(random_neighbour + 2) % 3]

        assert tetra31.check_neighbours_tetra(tetra22l) and tetra31.check_neighbours_tetra(tetra22r), \
            "SIM: Tetra31 and tetra22l or tetra22r are not neighbours"

        # Make sure the tetra is actually of type (2,2) and that the (2,2) tetras are neighbours
        if not tetra22l.is_22() or not tetra22r.is_22() or not tetra22l.check_neighbours_tetra(tetra22r):
            return False

        # Count the number of shared vertices between tetra22l and tetra22r
        shared_vertices = 0
        for i in range(self.constants.N_VERTICES_TETRA):
            if tetra22r.has_vertex(tetra22l.get_vertices()[i]):
                shared_vertices += 1

        # If the number of shared vertices is not 3, the move is not valid
        if shared_vertices != 3:
            return False
        
        print(f"ishift_u: tetra31: {tetra31_label} tetra22l: {tetra22l.ID} tetra22r: {tetra22r.ID}")
        
        # Try the move
        return self.universe.ishift_u(tetra31_label, tetra22l.ID, tetra22r.ID)

    def move_ishift_d(self) -> bool:
        """
        Metropolis-Hastings move to perform an inverse shift move with a
        (1,3)-tetrahedron.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        acceptance_probability = np.exp(self.k3)
        n3 = self.universe.tetrahedron_pool.get_number_occupied()
        vol_switch = self.universe.volfix_switch

        if vol_switch == 1 and self.target_volume > 0:
            acceptance_probability *= np.exp(
                -self.epsilon * (2 * self.target_volume - 2 * n3 -1))
            
        # Perform MCMC check for acceptance
        if not self.mcmc_check(acceptance_probability):
            return False
        
        # Pick a random (1,3)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra13 = self.universe.tetrahedron_pool.get(tetra31_label).get_tetras()[3]

        # Get random (2,2) neighbours of tetra13
        random_neighbour = self.rng.randint(0, 2)
        tetra22l = tetra13.get_tetras()[1 + random_neighbour]
        tetra22r = tetra13.get_tetras()[1 + (random_neighbour + 2) % 3]

        assert tetra13.check_neighbours_tetra(tetra22l) and tetra13.check_neighbours_tetra(tetra22r), \
            "SIM: Tetra13 and tetra22l or tetra22r are not neighbours"

        # Make sure the tetra is actually of type (2,2) and that the (2,2) tetras are neighbours
        if not tetra22l.is_22() or not tetra22r.is_22() or not tetra22l.check_neighbours_tetra(tetra22r):
            return False
        
        # Count the number of shared vertices between tetra22l and tetra22r
        shared_vertices = 0
        for i in range(self.constants.N_VERTICES_TETRA):
            if tetra22r.has_vertex(tetra22l.get_vertices()[i]):
                shared_vertices += 1

        # If the number of shared vertices is not 3, the move is not valid
        if shared_vertices != 3:
            return False
     
        # Try the move
        return self.universe.ishift_d(tetra13.ID, tetra22l.ID, tetra22r.ID)

    def prepare(self):
        """
        Prepares the universe for measurements in a sweep by updating the geometry.
        """
        self.universe.update_geometry()

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

        # Calculate the fixvolume based on the vol_switch
        vol_switch = self.universe.volfix_switch
        fixvolume = 0
        if vol_switch == 0:
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

        # Append the k3 value to the list
        self.k3_values.append(self.k3)

    def attempt_simple(self):
        # Generate a random number to select a move pair
        move = self.rng.randint(0, 1)

        if move == 0:
            # Add or delete move
            if self.rng.randint(0, 1) == 0:
                return 1 if self.move_add() else -1
            else:
                return 2 if self.move_delete() else -2
        else:
            # Flip move
            return 3 if self.move_flip() else -3
        
    def trial(self, k0, k3, seed, N):
        self.k0 = k0
        self.k3 = k3
        self.rng.seed(seed)
        
        test = {1: 0, -1: 0, 2: 0, -2: 0, 3: 0, -3: 0, 4: 0, -4: 0, 5: 0, -5: 0}
        for i in range(N):
            print(f"Trial: {i} \n")
            move_num = self.attempt_move()
            test[move_num] += 1
            print(f"Move: {move_num} \n")
            print(f"Test: {test} \n")
            self.universe.update_geometry()
            self.universe.check_validity()

if __name__ == "__main__":
    universe = Universe(geometry_infilename='initial_universes/test.dat', strictness=3, volfix_switch=0)
    # Start simulation
    simulation = Simulation(universe)

    # simulation.start(
    #     k0=1, k3=0, sweeps=10, thermal_sweeps=100, k_steps=100,
    #     target_volume=500, target2_volume=0, seed=0, outfile="output_test.txt",
    #     v1=1, v2=1, v3=1
    # )

    seed = random.randint(0, 1000000)
    # seed = 0
    simulation.trial(k0=1, k3=1, seed=seed, N=10000)