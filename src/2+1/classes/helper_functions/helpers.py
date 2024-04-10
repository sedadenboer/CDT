# Load packages
import sys
sys.path.append('..')
from typing import Dict, List
from classes.universe import Universe


def get_vertices_in_slice(universe: Universe) -> Dict[int, List[int]]:
    """
    Generates a dictionary with time t as keys and lists of vertex IDs with time t as values.

    Args:
        universe (Universe): The universe sample.
    """
    # Setup
    T = universe.n_slices
    vertices_in_slice = {t: [] for t in range(T)}

    # Get time of all vertices and add to dictionary under right key
    for vertex in universe.vertex_pool.get_objects():
        vertices_in_slice[vertex.time].append(vertex.ID)
    
    return vertices_in_slice

def get_spatial_neighbours(universe: Universe) -> Dict[int, List[int]]:
    """
    Generates a dictionary with vertex IDs as keys and a list with their spatial neighbours as values.

    Args:
        universe (Universe): The universe sample.
    """
    spatial_neighbours = {}

    # Filter the vertex_neighbours dict on spatial neighbours
    for vertex_id, neighbour_ids in universe.vertex_neighbours.items():
        vertex = universe.vertex_pool.get(vertex_id)

        for neighbour_id in neighbour_ids:
            neighbour = universe.vertex_pool.get(neighbour_id)

            if neighbour.time == vertex.time:
                if spatial_neighbours.get(vertex_id):
                    spatial_neighbours[vertex_id].append(neighbour_id)
                else:
                    spatial_neighbours[vertex_id] = [neighbour_id]
    
    return spatial_neighbours