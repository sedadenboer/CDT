# universe.py
#
# Author: Seda den Boer
# Date: 02-01-2024
# 
# Description:

from __future__ import annotations
from typing import cast
from vertex import Vertex
from triangle import Triangle
from pool import Pool
from bag import Bag
import pickle
import resource
import sys

max_rec = 0x100000

# May segfault without this line. 0x100 is a guess at the size of each stack frame.
resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
sys.setrecursionlimit(max_rec)


class Universe:
    """
    The Universe class represents the current state of the triangulation
    and stores properties of the geometry in a convenient matter. It also
    provides member functions that carry out changes on the geometry.

    Attributes:
        total_time (int): Total number of time slices.
        initial_slice_size (int): Initial size of the time slices.
        VERTEX_CAPACITY (int): Maximum number of vertices in the triangulation.
        TRIANGLE_CAPACITY (int): Maximum number of triangles in the triangulation.
        vertex_pool (Pool): Pool of vertices.
        triangle_pool (Pool): Pool of triangles.
        triangle_add_bag (Bag): Bag of triangles that can be added.
        four_vertices_bag (Bag): Bag of vertices that have degree 4.
        triangle_flip_bag (Bag): Bag of triangles that can be flipped.
        n_vertices (int): Total number of vertices in the triangulation.
        n_triangles (int): Total number of triangles in the triangulation.
        slice_sizes (dict): Dictionary with the size of each time slice.
        triangle_up_count (int): Number of triangles with an upwards orientation.
        triangle_down_count (int): Number of triangles with a downwards orientation.
    """

    def __init__(self, total_time: int, initial_slice_size: int, VERTEX_CAPACITY: int = 1000000):
        if total_time < 3:
            raise ValueError("Total time must be greater than 3.")
        if initial_slice_size < 3:
            raise ValueError("Initial slice size must be greater than 3.")
        if VERTEX_CAPACITY < 9:
            raise ValueError("Vertex capacity must be greater than 9.")
        
        self.total_time = total_time
        self.initial_slice_size = initial_slice_size

        VERTEX_CAPACITY = VERTEX_CAPACITY
        TRIANGLE_CAPACITY = 2 * VERTEX_CAPACITY

        # Create pools for vertices and triangles
        self.vertex_pool = Pool(capacity=VERTEX_CAPACITY)
        self.triangle_pool = Pool(capacity=TRIANGLE_CAPACITY)

        # Create bags
        self.triangle_add_bag = Bag(pool_capacity=TRIANGLE_CAPACITY)
        self.four_vertices_bag = Bag(pool_capacity=VERTEX_CAPACITY)
        self.triangle_flip_bag = Bag(pool_capacity=TRIANGLE_CAPACITY)

        # Total size of the triangulation   
        self.n_vertices = total_time * initial_slice_size
        self.n_triangles = (total_time - 1) * (initial_slice_size - 1) * 2

        # Create a dictionary to store the size of each time slice
        self.slice_sizes = {t: initial_slice_size for t in range(total_time)}

        # Number of triangles with a certain orientation
        self.triangle_up_count = 0
        self.triangle_down_count = 0

        # Initialise the triangulation, store the vertices and triangles
        self.initialise_triangulation()

    def initialise_triangulation(self):
        """
        Initialise the grid with vertices and triangles. Creates
        vertices and triangles and set their connectivity.
        """
        total_time = self.total_time
        width = self.initial_slice_size
        
    def insert_vertex(self, triangle_id: int) -> tuple[Vertex, Triangle, Triangle]:
        """
        Insert a vertex into the triangulation.

        Args:
            triangle (Triangle): Triangle to insert vertex into.

        Returns:
            tuple[Vertex, Triangle, Triangle]: The inserted vertex and the two triangles
            that are created with it.
        """
        pass
    
    def remove_vertex(self, vertex_id: int) -> tuple[Vertex, Triangle, Triangle]:
        """
        Remove a vertex from the triangulation.

        Args:
            vertex (Vertex): Vertex to remove.

        Returns:
            tuple[Vertex, Triangle, Triangle]: The removed vertex and the two triangles
            that are removed with it.
        """
        pass

    def flip_edge(self, triangle_id : int) -> tuple[Triangle, Triangle]:
        """
        Flip an edge in the triangulation.

        Args:
            triangle (Triangle): Triangle to flip edge in.

        Returns:
            tuple[Triangle, Triangle]: The two triangles that are flipped (t, tr).
        """
        pass
    
    def is_four_vertex(self, vertex: Vertex) -> bool:
        """
        Checks if a vertex is of degree 4.

        Args:
            vertex (Vertex): Vertex to be checked.

        Returns:
            bool: True if vertex is of degree 4, otherwise False.
        """
        return (vertex.get_triangle_left().get_triangle_right() 
                == vertex.get_triangle_right()) and (
                vertex.get_triangle_left().get_triangle_center().get_triangle_right() 
                == vertex.get_triangle_right().get_triangle_center())

    def get_total_size(self) -> int:
        """
        Get the total size of the triangulation.

        Returns:
            int: Total size of the triangulation.
        """
        return self.vertex_pool.get_number_occupied()
    
    def sort_vertices_periodic(self, vertices):
        pass

    def get_triangulation_state(self):
        """
        Get the current state of the triangulation.
        """
        pass

    def print_state(self):
        """
        Print the current state of the triangulation.
        """
        pass

    def check_validity(self):
        """
        Check the validity of the triangulation.
        """
        pass

    def save_to_file(self, filename):
        """
        Save the state of the Universe to a file using pickle.

        Args:
            filename (str): The name of the file to save the state to.
        """
        with open(filename, 'wb') as file:
            pickle.dump(self.__dict__, file)

    def load_from_file(self, filename):
        """
        Load the state of the Universe from a file using pickle.

        Args:
            filename (str): The name of the file to load the state from.
        """
        with open(filename, 'rb') as file:
            state = pickle.load(file)
        self.__dict__.update(state)