# vertex.py
#
# Author: Seda den Boer
# Date: 17-02-2024
# 
# Description: Defines a vertex in the triangulation.
# Contains methods to get and set set the tetrahedron that contains the vertex, as well as
# methods to check if a vertex is a neighbour of this vertex and to clear all
# references to other objects.

from __future__ import annotations
from typing import TYPE_CHECKING, List, Dict
if TYPE_CHECKING:
    from tetra import Tetrahedron
import numpy as np


class Vertex:
    """
    Represents a vertex in the triangulation.

    Args (Attributes):
        time: Slice number.

    Attributes:
        ID: Unique identifier for the vertex.
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

    def set_tetra(self, t: Tetrahedron):
        """
        Sets an (3,1)-tetrahedron that contains this vertex in its base.

        Raises:
            AssertionError: If the vertex is not part of the tetrahedron
                            or not in the base of the tetrahedron.
        
        Args:
            t (Tetrahedron): A (3,1)-tetrahedron that contains this vertex in its base.
        """
        # Make sure that the vertex is in the base of the tetrahedron
        assert t.has_vertex(self), f"Vertex {self.ID} is not part of the tetrahedron {t.ID}."
        assert np.array_equal(t.get_vertices(), self) != 3, f"Vertex {self.ID} is not in the base of the tetrahedron."
        self.tetra = t

    def get_tetra(self) -> Tetrahedron:
        """
        Returns an (3,1)-tetrahedron that contains this vertex in its base.

        Returns:
            Tetrahedron: A (3,1)-tetrahedron that contains this vertex in its base.
        """
        return self.tetra

    def check_vertex_neighbour(self, v: Vertex) -> bool:
        """
        Checks if the given vertex `v` is a neighbor of this vertex.
        Uses a breadth-first search to do so.

        Args:
            v (Vertex): The vertex to check for neighbor relationship.

        Returns:
            bool: True if `v` is a neighbor of this vertex, False otherwise.
        """
        if v == self:
            return False
        
        t = self.tetra

        # Dictionary to keep track of triangles that have been checked
        done: Dict[Vertex, bool] = {}
        
        # List of triangles to check
        current: List[Tetrahedron] = [t]
        # List of triangles to check next
        next: List[Tetrahedron] = []

        # Breadth-first search to find neighboring vertices
        while current:
            # For each tetrahedron in the current list
            for current_tetra in current:
                # For each tetrahedron neighbor of the current tetrahedron
                for neighbour_current_tetra in current_tetra.get_tetras():
                    # If the current tetrahedron does not contain this vertex, skip
                    if not neighbour_current_tetra.has_vertex(self):
                        continue

                    # If the tetrahedron contains the vertex, return True
                    if neighbour_current_tetra not in done:
                        if neighbour_current_tetra.has_vertex(v):
                            return True
                        
                        done[neighbour_current_tetra] = True
                        next.append(neighbour_current_tetra)

            current = next
            next = []
    
        return False

    def clear_references(self):
        """
        Clears all references to other objects.
        """
        self.tetra = None
        
    def log(self) -> str:
        """
        Returns a string representation of the vertex.

        Returns:
            str: A string representation of the vertex.
        """
        tetra = self.tetra
        return f"Vertex {self.ID} @ time {self.time} with cnum {self.cnum} and scnum {self.scnum}. Tetrahedron: {tetra}"