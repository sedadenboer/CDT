# triangle.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:

from __future__ import annotations
from typing import TYPE_CHECKING
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
        self.tl_: Triangle # left
        self.tr_: Triangle # right
        self.tc_: Triangle # vertical

        # The three vertices comprising the triangle t
        self.vl_: Vertex # left
        self.vr_: Vertex # right
        self.vc_: Vertex # apex

        # Time can only be determined after initialisation
        self.time: int

    def get_triangle_left(self) -> Triangle:
        """
        Get the triangle object to the left of this triangle.

        Returns:
            Triangle: Triangle to the left.
        """
        return self.tl_

    def get_triangle_right(self) -> Triangle:
        """
        Get the triangle object to the right of this triangle.

        Returns:
            Triangle: Triangle to the right.
        """
        return self.tr_

    def get_triangle_center(self) -> Triangle:
        """
        Get the triangle object to the center of this triangle.

        Returns:
            Triangle: Triangle to the center.
        """        
        return self.tc_

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

    def set_vertex_left(self, v: Vertex) -> None:
        """
        Set the left vertex of this triangle. Also updates the time of the triangle,
        and sets the triangle to the right of the left vertex.

        Args:
            v (Vertex): Left vertex.
        """
        self.vl_ = v
        self.time = v.time
        if self.is_upwards():
            v.set_triangle_right(self)

    def set_vertex_right(self, v: Vertex) -> None:
        """
        Set the right vertex of this triangle. Also updates the time of the triangle,
        and sets the triangle to the left of the right vertex.

        Args:
            v (Vertex): Right vertex.
        """
        self.vr_ = v
        self.time = v.time
        if self.is_upwards():
            v.set_triangle_left(self)

    def set_vertex_center(self, v: Vertex) -> None:
        """
        Set the center vertex of this triangle.

        Args:
            v (Vertex): Center vertex.
        """
        self.vc_ = v
    
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

        if self.is_upwards():
            vl.set_triangle_right(self)
            vr.set_triangle_left(self)

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
