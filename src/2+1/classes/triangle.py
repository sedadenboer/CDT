# triangle.py
#
# Author: Seda den Boer
# Date: 17-02-2024
# 
# Description: Defines a triangle in the triangulation.

from __future__ import annotations
from typing import TYPE_CHECKING, Tuple
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
        self.vs: Tuple[int] = ()
        self.hes: Tuple[int] = ()
        self.trnbr: Tuple[int] = ()
    
    def set_vertices(self, v0: Vertex, v1: Vertex, v2: Vertex):
        """
        Sets the vertices of the triangle.

        Args:
            v0 (Vertex): First vertex.
            v1 (Vertex): Second vertex.
            v2 (Vertex): Third vertex.
        """
        self.vs = (v0, v1, v2)
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
        self.hes = (h0, h1, h2)
    
    def set_triangle_neighbours(self, t0: Triangle, t1: Triangle, t2: Triangle):
        """
        Sets the triangle neighbours of the triangle.

        Args:
            t0 (Triangle): First triangle.
            t1 (Triangle): Second triangle.
            t2 (Triangle): Third triangle.
        """
        self.trnbr = (t0, t1, t2)
    
    def get_vertices(self) -> list[int, int, int]:
        """
        Returns the vertices of the triangle.

        Returns:
            list[int, int, int]: The vertices of the triangle.
        """
        return self.vs
    
    def get_half_edges(self) -> list[int, int, int]:
        """
        Returns the halfedges of the triangle.

        Returns:
            list[int, int, int]: The halfedges of the triangle.
        """
        return self.hes
    
    def get_triangle_neighbours(self) -> list[int, int, int]:
        """
        Returns the triangle neighbours of the triangle.

        Returns:
            list[int, int, int]: The triangle neighbours of the triangle.
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
        return v.ID in self.vs



