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
from universe import Universe
from typing import TYPE_CHECKING, List, Tuple
if TYPE_CHECKING:
    from universe import Universe


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
        self.volfix_switch = 0
        self.target_volume = 0
        self.target2_volume = 0
        self.rng = random.Random(0)
        self.epsilon = 0.00005
        self.measuring = False
        self.observables_3d = []
        self.observables_2d = []
        self.move_freqs = (1, 1, 1)
        self.k3_values = []
        self.N31_sizes = []
        self.N22_sizes = []
        self.N22_N31_ratio = []
        self.N22_N3_ratio = []
        self.total_tetrahedra = []
        self.total_vertices = []
        self.save_process = False
        self.save_final = False
        self.saving_interval = 10
        self.as_pickle = True

        # To save acceptance probabilities
        self.add_ap = []
        self.delete_ap = []
        self.flip_ap = []
        self.shift_ap = []
        self.ishift_ap = []

        # To save acceptance probabilities and rates
        self.succes_rates = {1: [], 2: [], 3: [], 4: [], 5: []}

    def start(self, k0: int, k3: int,
                    sweeps: int, thermal_sweeps: int, k_steps: int,
                    volfix_switch: int, target_volume: int, target2_volume: int,
                    seed: int, outfile: str, validity_check: bool,
                    v1: int, v2: int, v3: int):
        """
        Starts the MCMC CDT 2+1 simulation.

        Args:
            k0 (int): The k0 parameter.
            k3 (int): The k3 parameter.
            sweeps (int): The number of sweeps to perform.
            thermal_sweeps (int): The number of thermal sweeps to perform.
            k_steps (int): The number of steps to perform in each sweep.
            volfix_switch (int): The volume fix switch. 0 for (3,1)-tetrahedra, 1 for all tetrahedra.
            target_volume (int): The target volume of the simulation.
            target2_volume (int): The target 2D volume of the simulation.
            seed (int): The seed for the random number generator.
            outfile (str): The output file for the simulation.
            validity_check (bool): True if the validity of the universe should be checked, False otherwise.
            v1 (int): The first move frequency.
            v2 (int): The second move frequency.
            v3 (int): The third move frequency.
        """
        self.move_freqs = (v1, v2, v3)
        self.k0 = k0
        self.k3 = k3
        self.volfix_switch = volfix_switch
        self.target_volume = target_volume
        self.target2_volume = target2_volume
        self.rng.seed(seed)
        self.measuring = True
        self.saving_interval = thermal_sweeps

        # # Clear observables
        # for o in self.observables_3d:
        #     o.clear()
        # for o in self.observables_2d:
        #     o.clear()
        
        # Thermal sweeps
        print("========================================\n")
        print("THERMAL SWEEPS\n")
        print("----------------------------------------\n")
        print(f"k0 = {self.k0}, k3 = {self.k3}, epsilon = {self.epsilon}, thermal = {thermal_sweeps}, sweeps = {sweeps}, target = {target_volume}, target2d = {target2_volume}\n")
        print("----------------------------------------\n")

        for i in range(1, thermal_sweeps + 1):
            # Get the current state of the universe and print it
            total_2v = sum(self.universe.slice_sizes)
            n0 = self.universe.vertex_pool.get_number_occupied()
            n31 = self.universe.tetras_31.get_number_occupied()
            n3 = self.universe.tetrahedron_pool.get_number_occupied()
            n22 = self.universe.tetras_22.get_number_occupied()
            print(f"\nThermal i: {i} \t N0: {n0}, N3: {n3}, N31: {n31}, N13: {n3 - n31 - n22}, N22: {n22}, k0: {self.k0}, k3: {self.k3}")

            # Save acceptance probabilities
            add_ap_vals, delete_ap_vals, flip_ap_vals, shift_ap_vals, ishift_ap_vals = self.get_acceptance_probabilities()
            self.add_ap.append(add_ap_vals)
            self.delete_ap.append(delete_ap_vals)
            self.flip_ap.append(flip_ap_vals)
            self.shift_ap.append(shift_ap_vals)
            self.ishift_ap.append(ishift_ap_vals)
        
            # Perform sweeps and tune the k3 parameter
            self.perform_sweep(k_steps)
            self.tune()

            self.total_vertices.append(self.universe.vertex_pool.get_number_occupied())
            self.total_tetrahedra.append(n3)

            # Export the geometry every 10% of the thermal sweeps
            if self.save_process:
                if i % self.saving_interval == 0:
                    self.universe.export_geometry(outfile + f"_thermal_{i}", as_pickle=self.as_pickle)
                
            # # Update geometry and measure observables related to 3D structures
            # self.prepare()
            # for o in self.observables_3d:
            #     o.measure(self.universe)

        if validity_check:
            # self.universe.log()
            self.universe.check_validity()

        if sweeps > 0:
            # Main sweeps
            print("========================================\n")
            print("MAIN SWEEPS\n")
            print("----------------------------------------\n")
            print(f"k0 = {self.k0}, k3 = {self.k3}, epsilon = {self.epsilon}\n")
            print("----------------------------------------\n")

            for i in range(1, sweeps + 1):
                # Get the current state of the universe and print it
                total_2v = sum(self.universe.slice_sizes)
                n0 = self.universe.vertex_pool.get_number_occupied()
                n31 = self.universe.tetras_31.get_number_occupied()
                n3 = self.universe.tetrahedron_pool.get_number_occupied()
                n22 = self.universe.tetras_22.get_number_occupied()
                avg_2v = total_2v / self.universe.n_slices
                print(f"Main i: {i} target: {self.target_volume} target2d: {self.target2_volume} k0: {self.k0} k3: {self.k3} \t CURRENT N0: {n0}, N3: {n3}, N31: {n31}, N13: {n3 - n31 - n22}, N22: {n22}\n")

                # Save acceptance probabilities
                add_ap_vals, delete_ap_vals, flip_ap_vals, shift_ap_vals, ishift_ap_vals = self.get_acceptance_probabilities()
                self.add_ap.append(add_ap_vals)
                self.delete_ap.append(delete_ap_vals)
                self.flip_ap.append(flip_ap_vals)
                self.shift_ap.append(shift_ap_vals)
                self.ishift_ap.append(ishift_ap_vals)

                # Perform sweeps
                self.perform_sweep(k_steps)
                
                # Save sizes
                self.total_vertices.append(self.universe.vertex_pool.get_number_occupied())
                self.total_tetrahedra.append(n3)

                # Export the geometry every 10% of the main sweeps
                if self.save_process:
                    if i % self.saving_interval == 0:
                        self.universe.export_geometry(outfile + f"_main_{i}", as_pickle=self.as_pickle)
                
                # Check if there are any 3D observables to be measured
                if len(self.observables_3d) > 0: 
                    # Remove if volume should fluctuate during measurements
                    vol_switch = self.volfix_switch

                    # Ensure that the universe is at the target volume (if specified)
                    if target_volume > 0:
                        # Flag to track if the target volume is reached
                        compare = 0

                        # Attempt moves until the target volume is reached
                        while compare != target_volume:
                            self.attempt_move()
                            # Update the compare variable based on the volume switch
                            compare = self.universe.tetras31.get_number_occupied() if vol_switch == 0 else self.universe.tetrahedron_pool.get_number_occupied()

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
            
            if validity_check:
                self.universe.log()
                self.universe.check_validity()

        # Compute <N22/N3(1)> variable
        if sweeps > 0:
            # Remove thermal phase
            self.expected_N22_N31 = sum(self.N22_N31_ratio[(thermal_sweeps * k_steps):]) / (sweeps * k_steps)
            self.expected_N22_N3 = sum(self.N22_N3_ratio[(thermal_sweeps * k_steps):]) / (sweeps * k_steps)
        else:
            self.expected_N22_N31 = sum(self.N22_N31_ratio) / (thermal_sweeps * k_steps)
            self.expected_N22_N3 = sum(self.N22_N3_ratio) / (thermal_sweeps * k_steps)

        # Export the final geometry
        if self.save_final:
            self.universe.export_geometry(outfile + "_final")

    def attempt_move(self) -> int:
        """
        Attempt a move based on the move frequencies.

        Returns:
            int: The move number. 
            1 for add, 2 for delete, 3 for flip, 4 for shift, 5 for inverse shift.
            Negative numbers indicate failed moves.
        """        
        # Calculate cumulative frequencies
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
                
    def perform_sweep(self, n: int):
        """
        Perform a sweep of the simulation.

        Args:
            n (int): The number of moves to perform.
        """
        gathered_counts = [0, 0, 0, 0, 0]
        gathered_failed_counts = [0, 0, 0, 0, 0]
        
        # Perform n moves
        for _ in range(n):
            move_num = self.attempt_move()
            move = abs(move_num)

            # Update the move counters
            if move_num > 0:
                gathered_counts[move - 1] += 1
            else:
                gathered_failed_counts[move - 1] += 1

            # Append the sizes of N31, N3 and N22
            n31 = self.universe.tetras_31.get_number_occupied()
            n3 = self.universe.tetrahedron_pool.get_number_occupied()
            n22 = self.universe.tetras_22.get_number_occupied()
            self.N31_sizes.append(n31)
            self.N22_sizes.append(n22)
            self.N22_N31_ratio.append(n22 / n31)
            self.N22_N3_ratio.append(n22 / n3)

        for i, count in enumerate(gathered_counts):
            self.succes_rates[i + 1].append(count / (gathered_failed_counts[i] + count))
    
    def mcmc_check(self, acceptance_probability: float) -> bool:
        """
        Metropolis-Hastings check for acceptance of a move.

        Args:
            acceptance_probability (float): The acceptance probability of the move.

        Returns:
            bool: True if the move was accepted, False otherwise.
        """
        # Reject
        if acceptance_probability < 1:
            random_number = self.rng.random()
            if random_number > acceptance_probability:
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
        vol_switch = self.volfix_switch
        acceptance_probability = (n31 / (n0 + 1)) * np.exp(self.k0 - 4 * self.k3)

        # If the target volume is specified, adjust AP according to the volume switch
        if vol_switch == 0:
            if self.target_volume > 0:
                acceptance_probability *= np.exp(
                    4 * self.epsilon * (self.target_volume - n31 - 1))
        else:
            if self.target_volume > 0:
                acceptance_probability *= np.exp(
                    8 * self.epsilon * (self.target_volume - n3 - 2))

        # Perform MCMC check for acceptance
        if not self.mcmc_check(acceptance_probability):
            return False
        
        # Perform the move
        tetra31_label = self.universe.tetras_31.pick()
        return self.universe.add(tetra31_label)

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
        vol_switch = self.volfix_switch
        acceptance_probability = ((n0 + 1) / n31) * np.exp(-self.k0 + 4 * self.k3)
     
        # If the target volume is specified, adjust AP according to the volume switch
        if vol_switch == 0:
            if self.target_volume > 0:
                acceptance_probability *= np.exp(
                    -4 * self.epsilon * (self.target_volume - n31 - 1))
        else:
            if self.target_volume > 0:
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
        return self.universe.delete(vertex_label)

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

        # Check if the tetrahedron is actually flippable (opposite tetras should also be neighbours)
        if not tetra230.is_31():
            return False
        if not tetra012.get_tetras()[3].check_neighbours_tetra(tetra230.get_tetras()[3]):
            return False
        if not tetra012.get_tetras()[3].is_13():
            return False
        if not tetra230.get_tetras()[3].is_13():
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
        vol_switch = self.volfix_switch
     
        # If the target volume is specified, adjust AP according to the volume switch
        if vol_switch == 1:
            if self.target_volume > 0:
                acceptance_probability *= np.exp(self.epsilon * (2 * self.target_volume - 2 * n3 -1))

        # Perform MCMC check for acceptance
        if not self.mcmc_check(acceptance_probability):
            return False
        
        # Pick a random (3,1)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra31 = self.universe.tetrahedron_pool.get(tetra31_label)

        # Get random neighbour of tetra31
        random_neighbour = self.rng.randint(0, 2)
        tetra22 = tetra31.get_tetras()[random_neighbour]
        
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
        vol_switch = self.volfix_switch

        # If the target volume is specified, adjust AP according to the volume switch
        if vol_switch == 1:
            if self.target_volume > 0:
                acceptance_probability *= np.exp(self.epsilon * (2 * self.target_volume - 2 * n3 -1))

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
        if not tetra22.is_22():
            return False
        
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
        vol_switch = self.volfix_switch
       
        # If the target volume is specified, adjust AP according to the volume switch
        if vol_switch == 1:
            if self.target_volume > 0:
                acceptance_probability *= np.exp(-self.epsilon * (2 * self.target_volume - 2 * n3 -1))

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

        # Make sure the tetra is actually of type (2,2) and that the (2,2) tetras are neighbours
        if not tetra22l.is_22():
            return False
        if not tetra22r.is_22():
            return False
        if not tetra22l.check_neighbours_tetra(tetra22r):
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
        vol_switch = self.volfix_switch
       
        # If the target volume is specified, adjust AP according to the volume switch
        if vol_switch == 1:
            if self.target_volume > 0:
                acceptance_probability *= np.exp(-self.epsilon * (2 * self.target_volume - 2 * n3 -1))

        # Perform MCMC check for acceptance
        if not self.mcmc_check(acceptance_probability):
            return False
        
        # Pick a random (1,3)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra31 = self.universe.tetrahedron_pool.get(tetra31_label)
        tetra13 = tetra31.get_tetras()[3]

        # Get random (2,2) neighbours of tetra13
        random_neighbour = self.rng.randint(0, 2)
        tetra22l = tetra13.get_tetras()[1 + random_neighbour]
        tetra22r = tetra13.get_tetras()[1 + (random_neighbour + 2) % 3]

        # Make sure the tetra is actually of type (2,2) and that the (2,2) tetras are neighbours
        if not tetra22l.is_22():
            return False
        if not tetra22r.is_22():
            return False
        if not tetra22l.check_neighbours_tetra(tetra22r):
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
        vol_switch = self.volfix_switch
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

        # print(f"New k3: {self.k3}\n")
        # Append the k3 value to the list
        self.k3_values.append(self.k3)

    def get_acceptance_probabilities(self):
        """
        Calculates the acceptance probabilities for the different moves.
        """
        # Get relevant variables 
        n31 = self.universe.tetras_31.get_number_occupied()
        n3 = self.universe.tetrahedron_pool.get_number_occupied()
        n0 = self.universe.vertex_pool.get_number_occupied()
        vol_switch = self.volfix_switch
        add_ap = (n31 / (n0 + 1)) * np.exp(self.k0 - 4 * self.k3)
        delete_ap = ((n0 + 1) / n31) * np.exp(-self.k0 + 4 * self.k3)
        flip_ap = 1
        shift_ap = np.exp(-self.k3)
        ishift_ap = np.exp(self.k3)

        # If the target volume is specified, adjust AP according to the volume switch
        if vol_switch == 0:
            if self.target_volume > 0:
                add_ap *= np.exp(
                    4 * self.epsilon * (self.target_volume - n31 - 1))
                delete_ap *= np.exp(
                    -4 * self.epsilon * (self.target_volume - n31 - 1))
        else:
            if self.target_volume > 0:
                add_ap *= np.exp(
                    8 * self.epsilon * (self.target_volume - n3 - 2))
                delete_ap *= np.exp(
                    -8 * self.epsilon * (self.target_volume - n3 - 2))
                shift_ap *= np.exp(self.epsilon * (2 * self.target_volume - 2 * n3 -1))
                ishift_ap *= np.exp(-self.epsilon * (2 * self.target_volume - 2 * n3 -1))
        
        return add_ap, delete_ap, flip_ap, shift_ap, ishift_ap
        
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
        start = time.time()

        for i in range(N):
            # print(f"Trial: {i}")
            # move_num = self.attempt_move()
            move_num = self.attempt_move()
            test[move_num] += 1
            # print(f"Move: {move_num}, k3: {self.k3}")
            # print(f"Test: {test}")
            print(f"Total n31: {self.universe.tetras_31.get_number_occupied()}, Total n3: {self.universe.tetrahedron_pool.get_number_occupied()}, Total n0: {self.universe.vertex_pool.get_number_occupied()}\n")
            # self.universe.update_geometry()

            # if move_num > 0:
            #     self.universe.check_validity()

        print(f"{test}")
        self.universe.check_validity()

        end = time.time()
        print(f"Time: {end - start} \n")
        # self.universe.export_geometry(f"geometry_k0={k0}_k3={k3}_seed={seed}_N={N}")


if __name__ == "__main__":
    universe = Universe(geometry_infilename='initial_universes/sample-g0-T3.cdt', strictness=3)
    # universe_T32 = Universe(geometry_infilename='initial_universes/output_g=0_T=32.txt', strictness=3)
    
    simulation = Simulation(universe)
    simulation.save_process = False
    simulation.save_final = False
    simulation.as_pickle = True
    simulation.start( 
        k0=0, k3=0.72, sweeps=2, thermal_sweeps=5, k_steps=100000,
        volfix_switch=0, target_volume=10000, target2_volume=0,
        seed=0, outfile="saved_universes/test_run/output", validity_check=True,
        v1=1, v2=1, v3=1
    )
    print(f"Expected N22/N31: {simulation.expected_N22_N31_with_thermal}")
    print(f"Expected N22/N3: {simulation.expected_N22_N3_with_thermal}")

    # # seed = random.randint(0, 1000000)
    # seed = 1
    # simulation.trial(k0=7, k3=2.1, seed=seed, N=100000)