# tetra.py
#
# Author: Seda den Boer
# Date: 17-02-2024
#
# Description: A tetrahedron in the triangulation.

from __future__ import annotations
from typing import TYPE_CHECKING, Tuple, Union
if TYPE_CHECKING:
    from triangle import Triangle
    from vertex import Vertex
    from tetra import Tetrahedron
    from halfedge import HalfEdge


class Tetrahedron:
    """
    Represents a tetrahedron in the triangulation.

    Attributes:
        ID: Unique identifier for the tetrahedron.
        time: The time of the tetrahedron.'
        type: The type of the tetrahedron.
        tnbr: The tetrahedron neighbours of the tetrahedron.
        vs: The vertices of the tetrahedron.
        hes: The halfedges of the tetrahedron.
    """

    class Type:
        THREEONE = '31'
        ONETHREE = '13'
        TWOTWO = '22'

    def __init__(self) -> None:
        self.ID: int
        self.time: int
        self.type: Union[Tetrahedron.Type, None] = None
        self.tnbr: Tuple[Tetrahedron] = ()
        self.vs: Tuple[Vertex] = ()
        self.hes: Tuple[HalfEdge] = ()

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
    
    def set_vertices(self, v0: Vertex, v1: Vertex, v2: Vertex, v3: Vertex) -> None:
        """
        Sets the vertices of the tetrahedron. v1, v2, and v0 make up the base 
        of the tetrahedron, and v3 is the opposite vertex.

        Args:
            v0 (Vertex): First vertex.
            v1 (Vertex): Second vertex.
            v2 (Vertex): Third vertex.
            v3 (Vertex): Fourth vertex.
        """
        if v0.time == v1.time and v0.time == v2.time:
            self.type = Tetrahedron.Type.THREEONE
        if v1.time == v2.time and v1.time == v3.time:
            self.type = Tetrahedron.Type.ONETHREE
        if v0.time == v1.time and v2.time == v3.time:
            self.type = Tetrahedron.Type.TWOTWO

        assert v0.time != v3.time
        
        self.vs = (v0, v1, v2, v3)
        self.time = v0.time
    
    def get_vertices(self) -> list[Vertex]:
        """
        Returns the vertices of the tetrahedron. v1, v2, and v0 make up the base 
        of the tetrahedron, and v3 is the opposite vertex.

        Returns:
            list[Vertex]: The vertices of the tetrahedron.
        """
        return self.vs
    
    def set_half_edges(self, h0: HalfEdge, h1: HalfEdge, h2: HalfEdge):
        """
        Sets the halfedges of the tetrahedron.

        Args:
            h0 (HalfEdge): First halfedge.
            h1 (HalfEdge): Second halfedge.
            h2 (HalfEdge): Third halfedge.
        """
        self.hes = (h0, h1, h2)
    
    def get_half_edges(self) -> list[HalfEdge]:
        """
        Returns the halfedges of the tetrahedron.

        Returns:
            list[HalfEdge]: The halfedges of the tetrahedron.
        """
        return self.hes

    def get_half_edge_from(self, v: Vertex) -> HalfEdge:
        """
        Returns the halfedge of the tetrahedron that starts at the given vertex.

        Args:
            v (Vertex): The vertex to start the halfedge from.

        Returns:
            HalfEdge: The halfedge of the tetrahedron that starts at the given vertex.
        """
        for he in self.hes:
            if he.vs[0] == v:
                return he
            
        return None
    
    def get_half_edge_to(self, v: Vertex) -> HalfEdge:
        """
        Returns the halfedge of the tetrahedron that ends at the given vertex.

        Args:
            v (Vertex): The vertex to end the halfedge at.

        Returns:
            HalfEdge: The halfedge of the tetrahedron that ends at the given vertex.
        """
        for he in self.hes:
            if he.vs[1] == v:
                return he
            
        return None
    
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
        self.tnbr = (t0, t1, t2, t3)

    def get_tetras(self) -> list[Tetrahedron]:
        """
        Returns the tetrahedron neighbours of the tetrahedron.

        Returns:
            list[Tetrahedron]: The tetrahedron neighbours of the tetrahedron.
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
        return v in self.vs
    
    def check_neighbours_tetra(self, t: Tetrahedron) -> bool:
        """
        Returns whether the given tetrahedron is a neighbour of the tetrahedron.

        Args:
            t (Tetrahedron): The tetrahedron to check for.

        Returns:
            bool: True if the given tetrahedron is a neighbour of the tetrahedron, False otherwise.
        """
        return t in self.tnbr
    
    def get_tetra_opposite(self, v: Vertex) -> Tetrahedron:
        """
        Returns the tetrahedron opposite to the given vertex.

        Args:
            v (Vertex): The vertex to find the opposite tetrahedron for.

        Returns:
            Tetrahedron: The tetrahedron opposite to the given vertex.
        """
        assert self.has_vertex(v)

        for i, v in enumerate(self.vs):
            if self.vs[i] == v:
                return self.tnbr[i]
            
        assert False

    def get_vertex_opposite(self, v: Vertex) -> Vertex:
        """
        Returns the vertex opposite to the given vertex.

        Args:
            v (Vertex): The vertex to find the opposite vertex for.

        Returns:
            Vertex: The vertex opposite to the given vertex.
        """
        opposite_tetra = self.get_tetra_opposite(v)
        face_vertices = [v_other for v_other in opposite_tetra.vs if v_other != v]
        for v_opposite in opposite_tetra.vs:
            if v_opposite not in face_vertices:
                return v_opposite
            
        assert False

    def get_vertex_opposite_tetra(self, tn: Tetrahedron) -> Vertex:
        """
        Returns the vertex opposite to the given tetrahedron.

        Args:
            tn (Tetrahedron): The tetrahedron to find the opposite vertex for.

        Returns:
            Vertex: The vertex opposite to the given tetrahedron.
        """
        for i, t in enumerate(self.tnbr):
            if t == tn:
                return self.vs[i]
            
        assert False
    
    def exchange_tetra_opposite(self, v: Vertex, tn: Tetrahedron) -> None:
        """
        Exchanges the tetrahedron opposite to the given vertex with the given tetrahedron.

        Args:
            v (Vertex): The vertex to exchange the tetrahedron opposite for.
            tn (Tetrahedron): The tetrahedron to exchange the opposite tetrahedron for.
        """
        for i, v_i in enumerate(self.vs):
            if v_i == v:
                self.tnbr[i] = tn
    
    def log(self):
        """
        Prints information about the tetrahedron.
        """
        print(f"Tetrahedron {self.ID} @ {self.time} with type {self.to_string(self.type)}")
        for i, v in enumerate(self.vs):
            print(f"Vertex {i} @ {v.time}")
        for i, t in enumerate(self.tnbr):
            print(f"Tetra {i} @ {t.time}")
        for i, h in enumerate(self.hes):
            print(f"Halfedge {i} @ {h.time}")
        print()