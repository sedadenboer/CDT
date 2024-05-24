# check_mt_ray.py
#
# Author: Seda den Boer
# Date: 15-05-2024
#
# Descriptiopn: This file contains the functions that check if a move is possible 
# before it is executed. The functions are are similar to the ones in the
# Universe class but don't have the executement part where the move is performed, 
# whereas these functions have a check part that returns a the labels of the vertex/
# tetrahedra if the move is possible, or -1 if the move is not possible.

import ray
import random
import numpy as np
from typing import Tuple, Union

N_VERTICES_TETRA = 4
N_VERTICES_TRIANGLE = 3
N_MOVES = 5


def pick(array: np.ndarray, size: int, rng: random.Random) -> int:
    """
    Picks a random object from the array.

    Args:
        array (np.darray): The array to pick from.
        size (int): The size of the array.
        rng (random.Random): The random number generator.

    Returns:
        int: The index of the picked object.
    """
    if len(array) == 0:
        raise Exception("Empty.")
    else:
        # Pick random element from pool
        random_element = array[rng.randint(0, size - 1)]
        while not random_element or random_element == -1:
            random_element = array[rng.randint(0, size - 1)]

        return random_element

@ray.remote
def check_delete(vertex_pool, vertex_pool_size) -> int:
    """
    Helper function to check if a vertex can be deleted.
    """
    # Get a random vertex
    rng = random.Random()
    vertex = pick(vertex_pool, vertex_pool_size, rng)
    t01 = vertex.get_tetra()
    tv01 = t01.get_tetras()[3]
    
    # Get the vertex index in the tetrahedron
    vpos = np.where(t01.get_vertices() == vertex)[0][0]
    assert vpos >= 0

    # Get the base vertices of the two tetrahedra
    v0 = t01.get_vertices()[(vpos + 1) % 3]
    v1 = t01.get_vertices()[(vpos + 2) % 3]
    v2 = t01.get_vertex_opposite(v0)

    # Get all tetrahedra that need to be updated
    t12 = t01.get_tetra_opposite(v0)
    t20 = t01.get_tetra_opposite(v1)
    tv12 = tv01.get_tetra_opposite(v0)
    tv20 = tv01.get_tetra_opposite(v1)

    if (
        t01.is_31()
        and t12.is_31()
        and t20.is_31()
        and tv01.is_13()
        and tv12.is_13()
        and tv20.is_13()
        and v0.scnum >= 4
        and v1.scnum >= 4
        and v2.scnum >= 4
        and vertex.cnum == 6
        and vertex.scnum == 3
    ):
        return vertex.ID
    else:
        return -1

@ray.remote
def check_flip(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a flip move is possible.
    If possible, returns the labels of the tetrahedra to flip.
    """
    # Pick a random (3,1)-tetrahedron
    rng = random.Random()
    t012_label = pick(tetras_31, tetras_31_size, rng)
    t012 = tetrahedron_pool[t012_label]
    
    # Get random neighbour of t012
    random_neighbour = rng.randint(0, 2)
    t230 = t012.get_tetras()[random_neighbour]
    
    # Get the opposite (1,3)-tetrahedra
    tv012 = t012.get_tetras()[3]
    tv230 = t230.get_tetras()[3]

    # Get the vertices of the base triangles that are going to be linked
    v1 = t012.get_vertex_opposite_tetra(t230)
    v3 = t230.get_vertex_opposite_tetra(t012)
    
    # Get the remaining base vertices
    v1pos = np.where(t012.get_vertices() == v1)[0][0]
    v2 = t012.get_vertices()[(v1pos + 1) % 3]
    v0 = t012.get_vertices()[(v1pos + 2) % 3]

    # Get opposite neighbouring tetrahedra
    ta01 = t012.get_tetra_opposite(v2)
    ta23 = t230.get_tetra_opposite(v0)
    tva01 = tv012.get_tetra_opposite(v2)
    tva23 = tv230.get_tetra_opposite(v0)
    
    # Check if the tetrahedron is actually flippable (opposite tetras should also be neighbours)
    if (
        t230.is_31()
        and t012.get_tetras()[3].check_neighbours_tetra(t230.get_tetras()[3])
        and t012.get_tetras()[3].is_13()
        and t230.get_tetras()[3].is_13()
        and t012.is_31()
        and t230.is_31()
        and tv012.is_13()
        and tv230.is_13()
        and tv012.check_neighbours_tetra(tv230)
        and v1 != v3
        and v0.scnum >= 4
        and v2.scnum >= 4
        and ta01 != t230
        and ta23 != t012
        and tva01 != tv230
        and tva23 != tv012
        and not v1.check_vertex_neighbour(v3)
    ):
        return t012_label, t230.ID
    else:
        return -1

@ray.remote 
def check_shift_u(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a shift move is possible.
    If possible, returns the labels of the tetrahedra to shift.
    """
    # Pick a random (3,1)-tetrahedron
    rng = random.Random()
    t31_label = pick(tetras_31, tetras_31_size, rng)
    t31 = tetrahedron_pool[t31_label]

    # Get random neighbour of t31
    random_neighbour = rng.randint(0, 2)
    t22 = t31.get_tetras()[random_neighbour]

    # Get the vertices that will be linked
    v0 = t31.get_vertex_opposite_tetra(t22)
    v1 = t22.get_vertex_opposite_tetra(t31)
    
    # The remaining vertices
    v0pos = np.where(t31.get_vertices() == v0)[0][0]
    v2 = t31.get_vertices()[(v0pos + 1) % 3]
    v4 = t31.get_vertices()[(v0pos + 2) % 3]

    # Get neighbouring tetrahedra that need to be updated after the move
    ta023 = t31.get_tetra_opposite(v4)
    ta034 = t31.get_tetra_opposite(v2)
    ta123 = t22.get_tetra_opposite(v4)
    ta134 = t22.get_tetra_opposite(v2)
    
    # Check if the move is valid
    if (
        t22.is_22()
        and not ta023.has_vertex(v1)
        and not ta123.has_vertex(v0)
        and not ta034.has_vertex(v1)
        and not ta134.has_vertex(v0)
        and not v0.check_vertex_neighbour(v1)
    ):
        return t31_label, t22.ID
    else:
        return -1

@ray.remote      
def check_shift_d(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int], int]:
    """
    Helper function to check if a shift move is possible.
    If possible, returns the labels of the tetrahedra to shift.
    """
    # Pick a random (1,3)-tetrahedron
    rng = random.Random()  
    t31_label = pick(tetras_31, tetras_31_size, rng)
    t31 = tetrahedron_pool[t31_label]
    t13 = t31.get_tetras()[3]

    # Get random neighbour of t13
    random_neighbour = rng.randint(1, 3)
    t22 = t13.get_tetras()[random_neighbour]

    # Get the vertices that will be linked
    v0 = t13.get_vertex_opposite_tetra(t22)
    v1 = t22.get_vertex_opposite_tetra(t13)

    # The remaining vertices
    v0pos = np.where(t31.get_vertices() == v0)[0][0]
    v2 = t31.get_vertices()[(v0pos + 1) % 3]
    v4 = t31.get_vertices()[(v0pos + 2) % 3]

    # Get the neighbouring tetrahedra
    ta023 = t13.get_tetra_opposite(v4)
    ta034 = t13.get_tetra_opposite(v2)
    ta123 = t22.get_tetra_opposite(v4)
    ta134 = t22.get_tetra_opposite(v2)

    # Check if the tetrahedron is actually of type (2,2)
    if (
        t22.is_22()
        and not ta023.has_vertex(v1)
        and not ta123.has_vertex(v0)
        and not ta034.has_vertex(v1)
        and not ta134.has_vertex(v0)
        and not v0.check_vertex_neighbour(v1)
    ):
        return t13.ID, t22.ID
    else:
        return -1

@ray.remote
def check_ishift_u(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int, int], int]:
    """
    Helper function to check if an inverse shift move is possible.
    If possible, returns the labels of the tetrahedra to inverse shift.
    """
    # Pick a random (3,1)-tetrahedron
    rng = random.Random()
    t31_label = pick(tetras_31, tetras_31_size, rng)
    t31 = tetrahedron_pool[t31_label]

    # Get random neighbour of t31
    random_neighbour = rng.randint(0, 2)
    t22l = t31.get_tetras()[random_neighbour]
    t22r = t31.get_tetras()[(random_neighbour + 2) % 3]

     # Get the vertices of the interior triangle
    v1 = t31.get_vertices()[3]
    v3 = t22l.get_vertex_opposite_tetra(t31)
    v4 = t31.get_vertex_opposite_tetra(t22l)

    # The remaining vertices
    v4pos = np.where(t31.get_vertices() == v4)[0][0]
    v0 = t31.get_vertices()[(v4pos + 1) % 3]
    v2 = t31.get_vertices()[(v4pos + 2) % 3]

    # Get neighbouring tetrahedra that need to be updated after the move
    ta023 = t22l.get_tetra_opposite(v1)
    ta034 = t22r.get_tetra_opposite(v1)
    ta123 = t22l.get_tetra_opposite(v0)
    ta124 = t31.get_tetra_opposite(v0)
    ta134 = t22r.get_tetra_opposite(v0)
    
    # Count the number of shared vertices between t22l and t22r
    shared_vertices = 0
    for i in range(N_VERTICES_TETRA):
        if t22r.has_vertex(t22l.get_vertices()[i]):
            shared_vertices += 1

    # Make sure the tetra is of type (2,2) and that they are neighbours and have 3 shared vertices
    if (
        t22l.is_22()
        and t22r.is_22()
        and t22l.check_neighbours_tetra(t22r)
        and shared_vertices == 3
        and ta023.has_vertex(v4)
        and ta123.has_vertex(v4)
        and ta034.has_vertex(v2)
        and ta124.has_vertex(v3)
        and ta134.has_vertex(v2)
    ):
        return t31_label, t22l.ID, t22r.ID
    else:
        return -1

@ray.remote
def check_ishift_d(tetras_31, tetrahedron_pool, tetras_31_size) -> Union[Tuple[int, int, int], int]:
    """
    Helper function to check if an inverse shift move is possible.
    If possible, returns the labels of the tetrahedra to inverse shift.
    """
    # Pick a random (1,3)-tetrahedron
    rng = random.Random()
    t31_label = pick(tetras_31, tetras_31_size, rng)
    t31 = tetrahedron_pool[t31_label]
    t13 = t31.get_tetras()[3]

    # Get random (2,2) neighbours of t13
    random_neighbour = rng.randint(0, 2)
    t22l = t13.get_tetras()[1 + random_neighbour]
    t22r = t13.get_tetras()[1 + (random_neighbour + 2) % 3]

    # Get the vertices of the inner triangle
    v1 = t13.get_vertices()[0]
    v3 = t22l.get_vertex_opposite_tetra(t13)
    v4 = t13.get_vertex_opposite_tetra(t22l)

    # Get the remaining vertices
    v4pos = np.where(t31.get_vertices() == v4)[0][0]
    v0 = t31.get_vertices()[(v4pos + 1) % 3]
    v2 = t31.get_vertices()[(v4pos + 2) % 3]

    # Get the neighbouring tetrahedra
    ta023 = t22l.get_tetra_opposite(v1)
    ta034 = t22r.get_tetra_opposite(v1)
    ta123 = t22l.get_tetra_opposite(v0)
    ta124 = t13.get_tetra_opposite(v0)
    ta134 = t22r.get_tetra_opposite(v0)
    
    # Count the number of shared vertices between t22l and t22r
    shared_vertices = 0
    for i in range(N_VERTICES_TETRA):
        if t22r.has_vertex(t22l.get_vertices()[i]):
            shared_vertices += 1

    # Make sure the tetra is of type (2,2) and that they are neighbours and have 3 shared vertices
    if (
        t22l.is_22()
        and t22r.is_22()
        and t22l.check_neighbours_tetra(t22r)
        and shared_vertices == 3
        and ta023.has_vertex(v4)
        and ta123.has_vertex(v4)
        and ta034.has_vertex(v2)
        and ta124.has_vertex(v3)
        and ta134.has_vertex(v2)
    ):
        return t13.ID, t22l.ID, t22r.ID
    else:
        return -1