# halfedge.py
#
# Author: Seda den Boer
# Date: 17-02-2024
#
# Description: Defines a halfedge in the triangulation.


from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from vertex import Vertex
    from triangle import Triangle
    from tetra import Tetrahedron


class HalfEdge:
    """
    Represents a halfedge in the triangulation.

    Attributes:
        vs (np.ndarray): Array of vertices.
        adj (HalfEdge): Adjacent halfedge.
        next (HalfEdge): Next halfedge.
        prev (HalfEdge): Previous halfedge.
        tetra (Tetrahedron): Tetrahedron that contains this halfedge.
        triangle (Triangle): Triangle that contains this halfedge.
    """

    def __init__(self):
        """
        Initializes a new instance of the HalfEdge class.
        """
        self.vs: np.ndarray = np.empty(2, dtype=object)
        self.adj: HalfEdge = None
        self.next: HalfEdge = None
        self.prev: HalfEdge = None
        self.tetra: Tetrahedron = None
        self.triangle: Triangle = None

    def set_vertices(self, ve: Vertex, vf: Vertex):
        """
        Set the vertices of the halfedge.

        Args:
            ve (Vertex): First vertex
            vf (Vertex): Second vertex
        """
        self.vs[0] = ve
        self.vs[1] = vf
    
    def get_vertices(self) -> np.ndarray:
        """
        Get the vertices of the halfedge.

        Returns:
            np.ndarray: Array of vertices
        """
        return self.vs

    def get_adjacent(self) -> HalfEdge:
        """
        Get the adjacent halfedge, i.e. edges that
        share the same vertices.

        Returns:
            HalfEdge: Adjacent halfedge
        """
        return self.adj

    def set_adjacent(self, he: HalfEdge):
        """
        Set the adjacent halfedge, i.e. edges that
        share the same vertices.

        Args:
            he (HalfEdge): Adjacent halfedge
        """
        he.adj = self
        self.adj = he
    
    def set_next(self, he: HalfEdge):
        """
        Set the next halfedge.

        Args:
            he (HalfEdge): Next halfedge
        """
        self.next = he
    
    def get_next(self) -> HalfEdge:
        """
        Get the next halfedge.

        Returns:
            HalfEdge: Next halfedge
        """
        return self.next
    
    def set_previous(self, he: HalfEdge):
        """
        Set the previous halfedge.

        Args:
            he (HalfEdge): Previous halfedge
        """
        self.prev = he
    
    def set_tetra(self, t: Tetrahedron):
        """
        Set the tetrahedron that contains this halfedge.

        Args:
            t (Tetrahedron): Tetrahedron
        """
        self.tetra = t
    
    def get_tetra(self) -> Tetrahedron:
        """
        Get the tetrahedron that contains this halfedge.

        Returns:
            Tetrahedron: Tetrahedron
        """
        return self.tetra
    
    def set_triangle(self, t: Triangle):
        """
        Set the triangle that contains this halfedge.

        Args:
            t (Triangle): Triangle 
        """
        self.triangle = t
    
    def get_triangle(self) -> Triangle:
        """
        Get the triangle that contains this halfedge.

        Returns:
            Triangle: Triangle
        """
        return self.triangle
