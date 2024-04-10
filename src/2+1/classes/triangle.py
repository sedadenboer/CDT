from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from vertex import Vertex
    from halfedge import HalfEdge
    from triangle import Triangle


class Triangle:
    """
    Represents a triangle in the triangulation.

    Attributes:
        ID: Unique identifier for the triangle.
    """

    def __init__(self):
        self.ID: int = -1
        self.vs: np.ndarray = np.empty(3, dtype=object)
        self.hes: np.ndarray = np.empty(3, dtype=object)
        self.trnbr: np.ndarray = np.empty(3, dtype=object)
    
    def set_vertices(self, v0: Vertex, v1: Vertex, v2: Vertex):
        """
        Sets the vertices of the triangle.

        Args:
            v0 (Vertex): First vertex.
            v1 (Vertex): Second vertex.
            v2 (Vertex): Third vertex.
        """
        self.vs = np.array([v0, v1, v2], dtype=object)
        assert v0.time == v1.time == v2.time
        self.time = v0.time
    
    def set_half_edges(self, h0: HalfEdge, h1: HalfEdge, h2: HalfEdge):
        """
        Sets the halfedges of the triangle.

        Args:
            h0 (HalfEdge): First halfedge.
            h1 (HalfEdge): Second halfedge.
            h2 (HalfEdge): Third halfedge.
        """
        self.hes = np.array([h0, h1, h2], dtype=object)
    
    def set_triangle_neighbours(self, t0: Triangle, t1: Triangle, t2: Triangle):
        """
        Sets the triangle neighbours of the triangle.

        Args:
            t0 (Triangle): First triangle.
            t1 (Triangle): Second triangle.
            t2 (Triangle): Third triangle.
        """
        self.trnbr = np.array([t0, t1, t2], dtype=object)
    
    def get_vertices(self) -> np.ndarray:
        """
        Returns the vertices of the triangle.

        Returns:
            np.ndarray: The vertices of the triangle.
        """
        return self.vs
    
    def get_half_edges(self) -> np.ndarray:
        """
        Returns the halfedges of the triangle.

        Returns:
            np.ndarray: The halfedges of the triangle.
        """
        return self.hes
    
    def get_triangle_neighbours(self) -> np.ndarray:
        """
        Returns the triangle neighbours of the triangle.

        Returns:
            np.ndarray: The triangle neighbours of the triangle.
        """
        return self.trnbr
    
    def has_vertex(self, v: Vertex) -> bool:
        """
        Checks if the given vertex `v` is a vertex of this triangle.

        Args:
            v (Vertex): The vertex to check for.

        Returns:
            bool: True if `v` is a vertex of this triangle, False otherwise.
        """
        return np.any(self.vs == v.ID)
