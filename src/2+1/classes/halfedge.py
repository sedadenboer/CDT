# halfedge.py
#
# Author: Seda den Boer
# Date: 17-02-2024
#
# Description: Defines a halfedge in the triangulation.

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from vertex import Vertex
    from triangle import Triangle
    from tetra import Tetrahedron


class HalfEdge:
    """
    Represents a halfedge in the triangulation.

    Attributes:
        vs (list[Vertex]): List of vertices.
        adj (HalfEdge): Adjacent halfedge.
        next (HalfEdge): Next halfedge.
        prev (HalfEdge): Previous halfedge.
        tetra (Tetrahedron): Tetrahedron that contains this halfedge.
        triangle (Triangle): Triangle that contains this halfedge.
    """

    def __init__(self) -> None:
        self.vs: list[Vertex] = [None, None]
        self.adj: HalfEdge
        self.next: HalfEdge
        self.prev: HalfEdge
        self.tetra: Tetrahedron
        self.triangle: Triangle

    def set_vertices(self, ve: Vertex, vf: Vertex):
        self.vs[0] = ve
        self.vs[1] = vf
    
    def get_vertices(self) -> list[Vertex]:
        return self.vs

    def get_adjacent(self) -> HalfEdge:
        return self.adj

    def set_adj(self, he: HalfEdge):
        he.adj = self
        self.adj = he
    
    def get_adj(self) -> HalfEdge:
        return self.adj
    
    def set_next(self, he: HalfEdge):
        self.next = he
    
    def get_next(self) -> HalfEdge:
        return self.next
    
    def set_prev(self, he: HalfEdge):
        self.prev = he
    
    def set_tetra(self, t: Tetrahedron):
        self.tetra = t
    
    def get_tetra(self) -> Tetrahedron:
        return self.tetra
    
    def set_triangle(self, t: Triangle):
        self.triangle = t
    
    def get_triangle(self) -> Triangle:
        return self.triangle
