# vertex.py
#
# Author: Seda den Boer
# Date: 17-02-2024
# 
# Description: Defines a vertex in the triangulation.

from __future__ import annotations
from typing import TYPE_CHECKING, List
if TYPE_CHECKING:
    from tetra import Tetrahedron
import numpy as np
from numba import jit


class Vertex:
    """
    Represents a vertex in the triangulation.

    Attributes:
        ID: Unique identifier for the vertex.
        time: Slice number.
        tetra: Tetrahedron that contains this vertex.
        cnum: Number of tetrahedra that share this vertex.
        scnum: Spatial coordination number.
    """
    def __init__(self, time: int):
        self.time = time
        self.ID: int = -1
        self.tetra: Tetrahedron = None
        self.cnum = 0
        self.scnum = 0

    def set_tetra(self, t: Tetrahedron) -> None:
        """
        Sets an (3,1)-tetrahedron that contains this vertex in its base.

        Args:
            t (Tetrahedron): A (3,1)-tetrahedron that contains this vertex in its base.
        """
        # Make sure that the vertex is in the base of the tetrahedron
        assert t.has_vertex(self)
        assert np.array_equal(t.get_vertices(), self) != 3
        self.tetra = t

    def get_tetra(self) -> int:
        """
        Returns an (3,1)-tetrahedron that contains this vertex in its base.

        Returns:
            int: A (3,1)-tetrahedron that contains this vertex in its base.
        """
        return self.tetra

    def check_vertex_neighbour(self, v: Vertex) -> bool:
        """
        Checks if the given vertex `v` is a neighbor of this vertex.

        Args:
            v (Vertex): The vertex to check for neighbor relationship.

        Returns:
            bool: True if `v` is a neighbor of this vertex, False otherwise.
        """
        if v == self:
            return False
        
        # Set to keep track of visited tetrahedra
        visited: set[Tetrahedron] = set()
        # Queue for breadth-first search
        queue: List[Tetrahedron] = [self.tetra]

        # Breadth-first search to find neighboring vertices
        while queue:
            current_tetra = queue.pop(0)

            # Skip if the tetrahedron has already been visited
            if current_tetra in visited:
                continue

            # Check if the current tetrahedron contains the vertex
            if current_tetra.has_vertex(v):
                return True

            # Mark the current tetrahedron as visited
            visited.add(current_tetra)

            # Check the neighboring tetrahedra
            for neighbour_tetra in current_tetra.get_tetras():
                # Skip if the tetrahedron has already been visited
                if neighbour_tetra not in visited:
                    # Add the neighboring tetrahedron to the queue
                    queue.append(neighbour_tetra)

        return False

    def log(self) -> str:
        """
        Returns a string representation of the vertex.

        Returns:
            str: A string representation of the vertex.
        """
        return f"Vertex {self.ID} @ time {self.time} with cnum {self.cnum} and scnum {self.scnum}."