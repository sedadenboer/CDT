from typing import Tuple, Union
import random
from multiprocessing import Pool, shared_memory

class Constants:
    N_VERTICES_TETRA = 4
    N_VERTICES_TRIANGLE = 3
    N_MOVES = 5

class CheckMT:
    """
    Shared memory class to store the universe dictionary.
    """
    def __init__(self, universe, n_proposals):
        self.universe = universe
        self.n_proposals = n_proposals

        # Share memory between processes
        self.universe_shm = shared_memory.SharedMemory(create=True, size=self.universe.size)
        self.universe_shm.buf[:self.universe.size] = self.universe.buf

    def try_move(self, move: int) -> Union[Tuple[int], int]:
        """
        Helper function to try a move.
        """
        # Try out n times a move in parellel with shared memory resources and save the valid ones
        with Pool() as pool:
            results = pool.map(self.try_move_parallel, [move] * self.n_proposals)

        # Check if any of the results are valid
        for result in results:
            if result != -1:
                return result
            
        return -1
    
    def try_move_parallel(self, move: int) -> Union[Tuple[int], int]:
        """
        Helper function to try a move in parallel.
        """
        if move == 2:
            return self.check_delete()
        elif move == 3:
            return self.check_flip()
        elif move == 4:
            return self.check_shift_u()
        elif move == 5:
            return self.check_shift_d()
        elif move == 6:
            return self.check_ishift_u()
        elif move == 7:
            return self.check_ishift_d()
        else:
            raise ValueError("Invalid move type.")
            
    def check_delete(self) -> int:
        """
        Helper function to check if a vertex can be deleted.
        """
        # Get a random vertex
        vertex_label = self.universe.vertex_pool.pick()
        vertex = self.universe.vertex_pool.get(vertex_label)

        # Check if the vertex is actually deletable
        if (
            vertex.cnum == 6
            and vertex.scnum == 3
            and self.universe.delete(vertex_id=vertex_label, perform=False
            )
        ):
            # print("delete worked: ", vertex_label)
            return vertex_label
        else:
            # print("delete failed: ", vertex_label)
            return -1
        
    def check_flip(self) -> Union[Tuple[int, int], int]:
        """
        Helper function to check if a flip move is possible.
        If possible, returns the labels of the tetrahedra to flip.
        """
        rng = random.Random()  # Create a new RNG instance for each process
        tetra012_label = self.universe.tetras_31.pick()
        tetra012 = self.universe.tetrahedron_pool.get(tetra012_label)
        
        # Get random neighbour of tetra012
        random_neighbour = rng.randint(0, 2)
        tetra230 = tetra012.get_tetras()[random_neighbour]

        # Check if the tetrahedron is actually flippable (opposite tetras should also be neighbours)
        if (
            tetra230.is_31()
            and tetra012.get_tetras()[3].check_neighbours_tetra(tetra230.get_tetras()[3])
            and tetra012.get_tetras()[3].is_13()
            and tetra230.get_tetras()[3].is_13()
            and self.universe.flip(
                tetra012_id=tetra012_label,
                tetra230_id=tetra230.ID,
                perform=False
            )
        ):
            return tetra012_label, tetra230.ID
        else:
            return -1
            
    def check_shift_u(self) -> Union[Tuple[int, int], int]:
        """
        Helper function to check if a shift move is possible.
        If possible, returns the labels of the tetrahedra to shift.
        """
        rng = random.Random()  # Create a new RNG instance for each process
        # Pick a random (3,1)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra31 = self.universe.tetrahedron_pool.get(tetra31_label)

        # Get random neighbour of tetra31
        random_neighbour = rng.randint(0, 2)
        tetra22 = tetra31.get_tetras()[random_neighbour]

        # Check if the tetrahedron is actually of type (2,2)
        if (
            tetra22.is_22()
            and self.universe.shift_u(
                tetra31_id=tetra31_label,
                tetra22_id=tetra22.ID,
                perform=False
            )
        ):
            return tetra31_label, tetra22.ID
        else:
            return -1
            
    def check_shift_d(self) -> Union[Tuple[int, int], int]:
        """
        Helper function to check if a shift move is possible.
        If possible, returns the labels of the tetrahedra to shift.
        """
        rng = random.Random()  # Create a new RNG instance for each process
        # Pick a random (1,3)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra13 = self.universe.tetrahedron_pool.get(tetra31_label).get_tetras()[3]

        # Get random neighbour of tetra13
        random_neighbour = rng.randint(1, 3)
        tetra22 = tetra13.get_tetras()[random_neighbour]

        # Check if the tetrahedron is actually of type (2,2)
        if (
            tetra22.is_22()
            and self.universe.shift_d(
                tetra13_id=tetra13.ID,
                tetra22_id=tetra22.ID,
                perform=False
            )
        ):
            return tetra13.ID, tetra22.ID
        else:
            return -1

    def check_ishift_u(self) -> Union[Tuple[int, int, int], int]:
        """
        Helper function to check if an inverse shift move is possible.
        If possible, returns the labels of the tetrahedra to inverse shift.
        """
        rng = random.Random()  # Create a new RNG instance for each process
        # Pick a random (3,1)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra31 = self.universe.tetrahedron_pool.get(tetra31_label)

        # Get random neighbour of tetra31
        random_neighbour = rng.randint(0, 2)
        tetra22l = tetra31.get_tetras()[random_neighbour]
        tetra22r = tetra31.get_tetras()[(random_neighbour + 2) % 3]

        # Count the number of shared vertices between tetra22l and tetra22r
        shared_vertices = 0
        for i in range(Constants.N_VERTICES_TETRA):
            if tetra22r.has_vertex(tetra22l.get_vertices()[i]):
                shared_vertices += 1

        # Make sure the tetra is of type (2,2) and that they are neighbours and have 3 shared vertices
        if (
            tetra22l.is_22()
            and tetra22r.is_22()
            and tetra22l.check_neighbours_tetra(tetra22r)
            and shared_vertices == 3
            and self.universe.ishift_u(
                tetra31_id=tetra31_label,
                tetra22l_id=tetra22l.ID,
                tetra22r_id=tetra22r.ID,
                perform=False
            )
        ):
            return tetra31_label, tetra22l.ID, tetra22r.ID
        else:
            return -1

    def check_ishift_d(self) -> Union[Tuple[int, int, int], int]:
        """
        Helper function to check if an inverse shift move is possible.
        If possible, returns the labels of the tetrahedra to inverse shift.
        """
        rng = random.Random()  # Create a new RNG instance for each process
        # Pick a random (1,3)-tetrahedron
        tetra31_label = self.universe.tetras_31.pick()
        tetra13 = self.universe.tetrahedron_pool.get(tetra31_label).get_tetras()[3]

        # Get random (2,2) neighbours of tetra13
        random_neighbour = rng.randint(0, 2)
        tetra22l = tetra13.get_tetras()[1 + random_neighbour]
        tetra22r = tetra13.get_tetras()[1 + (random_neighbour + 2) % 3]

        # Count the number of shared vertices between tetra22l and tetra22r
        shared_vertices = 0
        for i in range(Constants.N_VERTICES_TETRA):
            if tetra22r.has_vertex(tetra22l.get_vertices()[i]):
                shared_vertices += 1

        # Make sure the tetra is of type (2,2) and that they are neighbours and have 3 shared vertices
        if (
            tetra22l.is_22()
            and tetra22r.is_22()
            and tetra22l.check_neighbours_tetra(tetra22r)
            and shared_vertices == 3
            and self.universe.ishift_d(
                tetra13_id=tetra13.ID,
                tetra22l_id=tetra22l.ID,
                tetra22r_id=tetra22r.ID,
                perform=False)
        ):
            return tetra13.ID, tetra22l.ID, tetra22r.ID
        else:
            return -1
        
    def __del__(self):
        self.universe_shm.close()
        self.universe_shm.unlink()