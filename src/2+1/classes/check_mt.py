from typing import Tuple, Union

class Constants:
        N_VERTICES_TETRA = 4
        N_VERTICES_TRIANGLE = 3
        N_MOVES = 5

def check_delete(universe) -> int:
    """
    Helper function to check if a vertex can be deleted.
    """
    # Get a random vertex
    vertex_label = universe.vertex_pool.pick()
    vertex = universe.vertex_pool.get(vertex_label)

    # Check if the vertex is actually deletable
    if (
        vertex.cnum == 6
        and vertex.scnum == 3
        and universe(vertex_id=vertex_label, perform=False
        )
    ):
        return vertex_label
    else:
        return -1
    
def check_flip(universe, rng) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a flip move is possible.
    If possible, returns the labels of the tetrahedra to flip.
    """
    tetra012_label = universe.tetras_31.pick()
    tetra012 = universe.tetrahedron_pool.get(tetra012_label)
    
    # Get random neighbour of tetra012
    random_neighbour = rng.randint(0, 2)
    tetra230 = tetra012.get_tetras()[random_neighbour]

    # Check if the tetrahedron is actually flippable (opposite tetras should also be neighbours)
    if (
        tetra230.is_31()
        and tetra012.get_tetras()[3].check_neighbours_tetra(tetra230.get_tetras()[3])
        and tetra012.get_tetras()[3].is_13()
        and tetra230.get_tetras()[3].is_13()
        and universe.flip(
            tetra012_label=tetra012_label,
            tetra230_label=tetra230.ID,
            perform=False
        )
    ):
        return tetra012_label, tetra230.ID
    else:
        return -1

def check_shift_u(universe, rng) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a shift move is possible.
    If possible, returns the labels of the tetrahedra to shift.
    """
    # Pick a random (3,1)-tetrahedron
    tetra31_label = universe.tetras_31.pick()
    tetra31 = universe.tetrahedron_pool.get(tetra31_label)

    # Get random neighbour of tetra31
    random_neighbour = rng.randint(0, 2)
    tetra22 = tetra31.get_tetras()[random_neighbour]

    # Check if the tetrahedron is actually of type (2,2)
    if (
        tetra22.is_22()
        and universe.shift_u(
            tetra31_label=tetra31_label,
            tetra22_label=tetra22.ID,
            perform=False
        )
    ):
        return tetra31_label, tetra22.ID
    else:
        return -1
        
def check_shift_d(universe, rng) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a shift move is possible.
    If possible, returns the labels of the tetrahedra to shift.
    """
    # Pick a random (1,3)-tetrahedron
    tetra31_label = universe.tetras_31.pick()
    tetra13 = universe.tetrahedron_pool.get(tetra31_label).get_tetras()[3]

    # Get random neighbour of tetra13
    random_neighbour = rng.randint(1, 3)
    tetra22 = tetra13.get_tetras()[random_neighbour]

    # Check if the tetrahedron is actually of type (2,2)
    if (
        tetra22.is_22()
        and universe.shift_d(
            tetra13_label=tetra13.ID,
            tetra22_label=tetra22.ID,
            perform=False
        )
    ):
        return tetra13.ID, tetra22.ID
    else:
        return -1

def check_ishift_u(universe, rng) -> Union[Tuple[int, int, int], int]:
    """
    Helper function to check if an inverse shift move is possible.
    If possible, returns the labels of the tetrahedra to inverse shift.
    """
    # Pick a random (3,1)-tetrahedron
    tetra31_label = universe.tetras_31.pick()
    tetra31 = universe.tetrahedron_pool.get(tetra31_label)

    # Get random neighbour of tetra31
    random_neighbour = rng.randint(0, 2)
    tetra22l = tetra31.get_tetras()[random_neighbour]
    tetra22r = tetra31.get_tetras()[(random_neighbour + 2) % 3]

    # Count the number of shared vertices between tetra22l and tetra22r
    shared_vertices = 0
    for i in range(Constants.N_VERTICES_TETRA):
        if tetra22r.has_vertex(tetra22l.get_vertices()[i]):
            shared_vertices += 1

    # Make sure the tetra is of type (2,2) and that the (2,2) tetras are neighbours and have 3 shared vertices
    if (
        tetra22l.is_22()
        and tetra22r.is_22()
        and tetra22l.check_neighbours_tetra(tetra22r)
        and shared_vertices == 3
    ):
        return tetra31_label, tetra22l.ID, tetra22r.ID
    else:
        return -1

def check_ishift_d(universe, rng) -> Union[Tuple[int, int, int], int]:
    """
    Helper function to check if an inverse shift move is possible.
    If possible, returns the labels of the tetrahedra to inverse shift.
    """
    # Pick a random (1,3)-tetrahedron
    tetra31_label = universe.tetras_31.pick()
    tetra31 = universe.tetrahedron_pool.get(tetra31_label)
    tetra13 = tetra31.get_tetras()[3]

    # Get random (2,2) neighbours of tetra13
    random_neighbour = rng.randint(0, 2)
    tetra22l = tetra13.get_tetras()[1 + random_neighbour]
    tetra22r = tetra13.get_tetras()[1 + (random_neighbour + 2) % 3]

    # Count the number of shared vertices between tetra22l and tetra22r
    shared_vertices = 0
    for i in range(Constants.N_VERTICES_TETRA):
        if tetra22r.has_vertex(tetra22l.get_vertices()[i]):
            shared_vertices += 1

    # Make sure the tetra is of type (2,2) and that the (2,2) tetras are neighbours and have 3 shared vertices
    if (
        tetra22l.is_22()
        and tetra22r.is_22()
        and tetra22l.check_neighbours_tetra(tetra22r)
        and shared_vertices == 3
    ):
        return tetra13.ID, tetra22l.ID, tetra22r.ID
    else:
        return -1