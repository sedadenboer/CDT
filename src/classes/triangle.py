# triangle.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:

from __future__ import annotations
from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from vertex import Vertex


class Triangle:
    """
    Represents a triangle in the triangulation.

    Attributes:
        ID: Unique identifier for the triangle.
        type: Type of the triangle (upwards or downwards).
        tl: Triangle to the left.
        tr: Triangle to the right.
        tc: Triangle to the center (apex).
        vl: Vertex to the left.
        vr: Vertex to the right.
        vc: Vertex to the center.
    """

    def __init__(self) -> None:
        self.ID: int
        self.type: str

        # The three triangles sharing and edge with triangle t
        self.tl_: Triangle = None # left
        self.tr_: Triangle = None # right
        self.tc_: Triangle = None # vertical

        # The three vertices comprising the triangle t
        self.vl_: Vertex = None # left
        self.vr_: Vertex = None # right
        self.vc_: Vertex = None # apex

        # Time can only be determined after initialisation
        self.time: int = None

    def get_triangle_left(self) -> Triangle:
        """
        Get the triangle object to the left of this triangle.

        Raises:
            ValueError: If there is no left triangle.

        Returns:
            Triangle: Triangle to the left.
        """
        if not self.tl_:
            raise ValueError("No triangle left.")
        
        return self.tl_

    def get_triangle_right(self) -> Triangle:
        """
        Get the triangle object to the right of this triangle.

        Raises:
            ValueError: If there is no right triangle.

        Returns:
            Triangle: Triangle to the right.
        """
        if not self.tr_:
            raise ValueError("No triangle right.")
        
        return self.tr_

    def get_triangle_center(self) -> Triangle:
        """
        Get the triangle object to the center of this triangle.

        Raises:
            ValueError: If there is no center triangle.

        Returns:
            Triangle: Triangle to the center.
        """        
        if not self.tc_:
            raise ValueError("No triangle center.")
        
        return self.tc_
    
    def get_triangles(self) -> Tuple[Triangle, Triangle, Triangle]:
        """
        Get the triangles to the left, right and center of this triangle.

        Raises:
            ValueError: If there are no triangles.

        Returns:
            Tuple[Triangle, Triangle, Triangle]: Triangles to the left, right and center.
        """
        if not self.tl_:
            raise ValueError("No triangle left.")
        if not self.tr_:
            raise ValueError("No triangle right.")
        if not self.tc_:
            raise ValueError("No triangle center.")
        
        return self.tl_, self.tr_, self.tc_

    def set_triangle_left(self, t: Triangle) -> None:
        """
        Set the triangle to the left of this triangle. Also ensures vice versa
        connections.

        Args:
            t (Triangle): Triangle to the left.
        """
        self.tl_ = t 

        if t:
            t.tr_ = self

    def set_triangle_right(self, t: Triangle) -> None:
        """
        Set the triangle to the right of this triangle. Also ensures vice versa
        connections.

        Args:
            t (Triangle): Triangle to the right.
        """
        self.tr_ = t

        if t:
            t.tl_ = self

    def set_triangle_center(self, t: Triangle) -> None:
        """
        Set the triangle to the center of this triangle. Also ensures vice versa
        connections.

        Args:
            t (Triangle): Triangle to the center.
        """
        self.tc_ = t

        if t:
            t.tc_ = self

    def set_triangles(self, tl: Triangle, tr : Triangle, tc: Triangle) -> None:
        """
        Set the triangles to the left, right and center of this triangle. Also
        ensures vice versa connections.

        Args:
            tl (Triangle): Triangle to the left.
            tr (Triangle): Triangle to the right.
            tc (Triangle): Triangle to the center.
        """
        self.tl_ = tl
        self.tr_ = tr
        self.tc_ = tc

        # Ensure vice versa connections
        if tl:
            tl.tr_ = self
        if tr:
            tr.tl_ = self
        if tc:
            tc.tc_ = self

    def get_vertex_left(self) -> Vertex:
        """
        Get the left vertex of this triangle.

        Raises:
            ValueError: If there is no left vertex.

        Returns:
            Vertex: Left vertex.
        """
        if not self.vl_:
            raise ValueError("No vertex left.")
        
        return self.vl_

    def get_vertex_right(self) -> Vertex:
        """
        Get the right vertex of this triangle.

        Raises:
            ValueError: If there is no right vertex.

        Returns:
            Vertex: Right vertex.
        """
        if not self.vr_:
            raise ValueError("No vertex right.")
        
        return self.vr_

    def get_vertex_center(self) -> Vertex:
        """
        Get the center vertex of this triangle.

        Raises:
            ValueError: If there is no center vertex.

        Returns:
            Vertex: Center vertex.
        """
        if not self.vc_:
            raise ValueError("No vertex center.")
        
        return self.vc_

    def get_vertices(self) -> Tuple[Vertex, Vertex, Vertex]:
        """
        Get the vertices of this triangle.

        Raises:
            ValueError: If there are no vertices.

        Returns:
            Tuple[Vertex, Vertex, Vertex]: Left, right and center vertex.
        """  
        return self.vl_, self.vr_, self.vc_

    def set_vertex_right(self, v: Vertex) -> None:
        """
        Set the right vertex of this triangle. Also updates the time of the triangle,
        and sets the triangle to the left of the right vertex.

        Args:
            v (Vertex): Right vertex.
        """
        self.vr_ = v
        self.time = v.time
        self.update_type()

        # Update vertex connections in space
        v.set_neighbour_left(self.vl_)

        if self.is_upwards():
            v.set_triangle_left(self)

    def set_vertices(self, vl: Vertex, vr: Vertex, vc: Vertex) -> None:
        """
        Set the vertices of this triangle. Also updates the time of the triangle,
        and sets the triangle to the right of the left vertex and the triangle to
        the left of the right vertex.

        Args:
            vl (Vertex): Left vertex.
            vr (Vertex): Right vertex.
            vc (Vertex): Center vertex.
        """
        self.vl_ = vl
        self.vr_ = vr
        self.vc_ = vc
        self.time = vl.time
        self.update_type()

        # Save vertex-vertex connections in space
        vl.set_neighbour_right(vr)
        vr.set_neighbour_left(vl)

        if self.is_upwards():
            # Save vertex-triangle connections for upwards triangles
            vl.set_triangle_right(self)
            vr.set_triangle_left(self)

            # Save vertex-vertex connections in time for upwards triangles
            vl.add_future_neighbour(vc)
            vr.add_future_neighbour(vc)
        else:
            # Save vertex-vertex connections in time for downwards triangles
            vl.add_past_neighbour(vc)
            vr.add_past_neighbour(vc)

    def is_upwards(self) -> bool:
        """
        Return True if the triangle is upwards, False otherwise.
        """
        return self.type == "UP"

    def is_downwards(self) -> bool:
        """
        Return True if the triangle is downwards, False otherwise.
        """
        return self.type == "DOWN"

    def update_type(self) -> None:
        """
        Update the type of the triangle based on the times of the vertices.
        """
        if self.vl_.time < self.vc_.time:
            self.type = "UP"
        else:
            self.type = "DOWN"

        # Correct for boundary conditions
        if self.vc_.time == 0 and self.vl_.time > 1:
            self.type = "UP"
        
        if self.vl_.time == 0 and self.vc_.time > 1:
            self.type = "DOWN"
