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
    """

    def __init__(self) -> None:
        self.vs: list[Vertex] = [None, None]
        self.adj: HalfEdge
        self.next: HalfEdge
        self.prev: HalfEdge
        self.tetra: Tetrahedron
        self.triangle: Triangle

    def set_vertices(self, ve: Vertex, vf: Vertex) -> None:
        self.vs[0] = ve
        self.vs[1] = vf
    
    def get_vertices(self) -> list[Vertex]:
        return self.vs

    def get_adjacent(self) -> HalfEdge:
        return self.adj

    def set_adjacent(self, he: HalfEdge) -> None:
        he.adj = self
        self.adj = he
