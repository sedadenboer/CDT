# tetra.py
#
# Author: Seda den Boer
# Date: 17-02-2024
# 
# Description: Defines a tetrahedron in the triangulation.
# Contains functions to set the vertices and tetrahedron neighbours, as well as
# functions to get and set the vertices and tetrahedron neighbours of the tetrahedron, 
# check the type of the tetrahedron, and remove all references to other objects.

from __future__ import annotations
from typing import TYPE_CHECKING, Union
if TYPE_CHECKING:
    from vertex import Vertex
    from tetra import Tetrahedron
import numpy as np


class Tetrahedron:
    """
    Represents a tetrahedron in the triangulation.

    Attributes:
        ID: Unique identifier for the tetrahedron.
        time: The time of the tetrahedron.'
        type: The type of the tetrahedron.
        tnbr: The tetrahedron neighbours of the tetrahedron.
        vs: The vertices of the tetrahedron.
    """
    class Type:
        THREEONE = '31'
        ONETHREE = '13'
        TWOTWO = '22'

    def __init__(self):
        self.ID: int = -1
        self.time: int = -1
        self.type: Union[Tetrahedron.Type, None] = None
        self.tnbr: np.ndarray = np.empty(4, dtype=object)
        self.vs: np.ndarray = np.empty(4, dtype=object)

    def to_string(self, t: Type) -> str:
        """
        Returns a string representation of the tetrahedron.

        Args:
            t (Type): The type of tetrahedron.

        Returns:
            str: String representation of the tetrahedron.
        """
        return {
            Tetrahedron.Type.THREEONE: "31",
            Tetrahedron.Type.ONETHREE: "13",
            Tetrahedron.Type.TWOTWO: "22"
        }[t]
    
    def set_vertices(self, v0: Vertex, v1: Vertex, v2: Vertex, v3: Vertex):
        """
        Sets the vertices of the tetrahedron.

        Args:
            v0 (Vertex): First vertex.
            v1 (Vertex): Second vertex.
            v2 (Vertex): Third vertex.
            v3 (Vertex): Fourth vertex.
        """
        if v0.time == v1.time == v2.time:
            self.type = Tetrahedron.Type.THREEONE
        elif v1.time == v2.time == v3.time:
            self.type = Tetrahedron.Type.ONETHREE
        elif v0.time == v1.time and v2.time == v3.time:
            self.type = Tetrahedron.Type.TWOTWO

        assert v0.time != v3.time
        
        self.vs = np.array([v0, v1, v2, v3], dtype=object)
        self.time = v0.time
    
    def get_vertices(self) -> np.ndarray[Vertex]:
        """
        Returns the vertices of the tetrahedron. v1, v2, and v0 make up the base 
        of the tetrahedron, and v3 is the opposite vertex.

        Returns:
            np.ndarray[Vertex]: The vertices of the tetrahedron.
        """
        return self.vs
    
    def set_tetras(self, t0: Tetrahedron, t1: Tetrahedron, t2: Tetrahedron, t3: Tetrahedron):
        """
        Sets the tetrahedron neighbours of the tetrahedron. The order of the tetrahedra is
        important. t0 and t2 are (2,2)- oriented, and t1 is a(3,1) oriented. The last
        tetrahedron is opposite to the last vertex, and has orientation (1,3).

        Args:
            t0 (Tetrahedron): First tetrahedron.
            t1 (Tetrahedron): Second tetrahedron.
            t2 (Tetrahedron): Third tetrahedron.
            t3 (Tetrahedron): Fourth tetrahedron.
        """
        self.tnbr = np.array([t0, t1, t2, t3], dtype=object)

    def get_tetras(self) -> np.ndarray[Tetrahedron]:
        """
        Returns the tetrahedron neighbours of the tetrahedron.

        Returns:
            np.ndarray[Tetrahedron]: The tetrahedron neighbours of the tetrahedron.
        """
        return self.tnbr
    
    def is_31(self) -> bool:
        """
        Returns whether the tetrahedron is of type 31.

        Returns:
            bool: True if the tetrahedron is of type 31, False otherwise.
        """
        return self.type == Tetrahedron.Type.THREEONE
    
    def is_13(self) -> bool:
        """
        Returns whether the tetrahedron is of type 13.

        Returns:
            bool: True if the tetrahedron is of type 13, False otherwise.
        """
        return self.type == Tetrahedron.Type.ONETHREE
    
    def is_22(self) -> bool:
        """
        Returns whether the tetrahedron is of type 22.

        Returns:
            bool: True if the tetrahedron is of type 22, False otherwise.
        """
        return self.type == Tetrahedron.Type.TWOTWO
    
    def has_vertex(self, v: Vertex) -> bool:
        """
        Returns whether the tetrahedron has the given vertex.

        Args:
            v (Vertex): The vertex to check for.

        Returns:
            bool: True if the tetrahedron has the given vertex, False otherwise.
        """
        set_vs = set(self.vs)
        return v in set_vs
    
    def check_neighbours_tetra(self, t: Tetrahedron) -> bool:
        """
        Returns whether the given tetrahedron is a neighbour of the tetrahedron.

        Args:
            t (Tetrahedron): The tetrahedron to check for.

        Returns:
            bool: True if the given tetrahedron is a neighbour of the tetrahedron, False otherwise.
        """
        set_tnbr = set(self.tnbr)
        return t in set_tnbr
    
    def get_tetra_opposite(self, v: Vertex) -> Tetrahedron:
        """
        Returns the tetrahedron opposite to the given vertex.

        Args:
            v (Vertex): The vertex to find the opposite tetrahedron for.

        Returns:
            Tetrahedron: The tetrahedron opposite to the given vertex.

        Raises:
            ValueError: If the vertex is not in the tetrahedron.
        """
        vertex_index = np.where(self.vs == v)[0]

        # Check if the vertex is in the tetrahedron (index not empty)
        if len(vertex_index) > 0:
            return self.tnbr[vertex_index[0]]
            
        raise ValueError(f"Vertex {v.ID} is not in tetrahedron {self.ID}")

    def get_vertex_opposite(self, v: Vertex) -> Vertex:
        """
        Returns the vertex opposite to the given vertex in the opposite tetrahedron.

        Args:
            v (Vertex): The vertex to find the opposite vertex for.

        Returns:
            Vertex: The vertex opposite to the given vertex.
        """
        tn = self.get_tetra_opposite(v)

        face = set(self.vs)
        face.remove(v)

        # Find the vertex that is not in the face
        for tnv in tn.vs:
            if tnv not in face:
                return tnv

        raise ValueError(f"Vertex {v.ID} is not in tetrahedron {self.ID}")

    def get_vertex_opposite_tetra(self, tn: Tetrahedron) -> Vertex:
        """
        Returns the vertex opposite to the given tetrahedron.

        Args:
            tn (Tetrahedron): The tetrahedron to find the opposite vertex for.

        Returns:
            Vertex: The vertex opposite to the given tetrahedron.

        Raises:
            ValueError: If the given tetrahedron is not a neighbour of the tetrahedron.
        """
        tn_index = np.where(self.tnbr == tn)[0]

        # Check if the tetrahedron is a neighbour of the tetrahedron (index not empty)
        if len(tn_index) > 0:
            return self.vs[tn_index[0]]
            
        raise ValueError(f"Tetrahedron {tn.ID} is not a neighbour of tetrahedron {self.ID}")

    def exchange_tetra_opposite(self, v: Vertex, tn: Tetrahedron):
        """
        Exchanges the tetrahedron opposite to the given vertex with the given tetrahedron.

        Args:
            v (Vertex): The vertex to exchange the tetrahedron opposite for.
            tn (Tetrahedron): The tetrahedron to exchange the opposite tetrahedron for.

        Raises:
            ValueError: If the vertex is not in the tetrahedron.
        """
        index = np.where(self.vs == v)[0]

        # Check if the vertex is in the tetrahedron (index not empty)
        if len(index) > 0:
            self.tnbr[index[0]] = tn
        else:
            raise ValueError(f"Vertex {v.ID} is not in tetrahedron {self.ID}")

    def clear_references(self):
        """
        Clears all the tetrahedron neighbours and vertices.
        """
        self.tnbr = np.empty(4, dtype=object)
        self.vs = np.empty(4, dtype=object)
        
    def log(self):
        """
        Prints information about the tetrahedron.
        """
        print(f"Tetrahedron {self.ID} @ {self.time} with type {self.to_string(self.type)}")
        print(f"Vertices: {[v.ID for v in self.vs]}, time: {[v.time for v in self.vs]}")
        print(f"Neighbours: {[t.ID for t in self.tnbr]}, time: {[t.time for t in self.tnbr]}")
        print()