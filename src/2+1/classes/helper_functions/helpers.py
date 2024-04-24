# Load packages
from __future__ import print_function
import sys
sys.path.append('..')
from typing import Dict, List
from classes.universe import Universe
from sys import getsizeof, stderr
from itertools import chain
from collections import deque
try:
    from reprlib import repr
except ImportError:
    pass

def total_size(o, handlers={}, verbose=False):
    """ Returns the approximate memory footprint an object and all of its contents.

    Automatically finds the contents of the following builtin containers and
    their subclasses:  tuple, list, deque, dict, set and frozenset.
    To search other containers, add handlers to iterate over their contents:

        handlers = {SomeContainerClass: iter,
                    OtherContainerClass: OtherContainerClass.get_elements}

    """
    dict_handler = lambda d: chain.from_iterable(d.items())
    all_handlers = {tuple: iter,
                    list: iter,
                    deque: iter,
                    dict: dict_handler,
                    set: iter,
                    frozenset: iter,
                   }
    all_handlers.update(handlers)     # user handlers take precedence
    seen = set()                      # track which object id's have already been seen
    default_size = getsizeof(0)       # estimate sizeof object without __sizeof__

    def sizeof(o):
        if id(o) in seen:       # do not double count the same object
            return 0
        seen.add(id(o))
        s = getsizeof(o, default_size)

        if verbose:
            print(s, type(o), repr(o), file=stderr)

        for typ, handler in all_handlers.items():
            if isinstance(o, typ):
                s += sum(map(sizeof, handler(o)))
                break
        return s

    return sizeof(o)

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