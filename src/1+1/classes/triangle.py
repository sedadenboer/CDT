# triangle.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description: Represents a triangle in the triangulation.
# Contains methods to get and set the neighbouring triangles and vertices of the triangle,
# as well as methods to determine the type of the triangle and to clear all
# references to other objects.

from __future__ import annotations
from typing import TYPE_CHECKING, Tuple
if TYPE_CHECKING:
    from vertex import Vertex


class Triangle:
    """
    Represents a triangle in the triangulation.

    Attributes:
        ID (int): Unique identifier for the triangle.
        type (str): Type of the triangle (upwards or downwards).
        tl (Triangle): Triangle to the left.
        tr (Triangle): Triangle to the right.
        tc (Triangle): Triangle to the center (apex).
        vl (Vertex): Vertex to the left.
        vr (Vertex): Vertex to the right.
        vc (Vertex): Vertex to the center.
        time (int): Time of the triangle (time of the base).
    """

    def __init__(self):
        self.ID: int
        self.type: str

        # The three triangles sharing and edge with triangle t
        self.tl: Triangle = None 
        self.tr: Triangle = None
        self.tc: Triangle = None 

        # The three vertices comprising the triangle t
        self.vl: Vertex = None
        self.vr: Vertex = None
        self.vc: Vertex = None

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
        if not self.tl:
            raise ValueError("No triangle left.")
        
        return self.tl

    def get_triangle_right(self) -> Triangle:
        """
        Get the triangle object to the right of this triangle.

        Raises:
            ValueError: If there is no right triangle.

        Returns:
            Triangle: Triangle to the right.
        """
        if not self.tr:
            raise ValueError("No triangle right.")
        
        return self.tr

    def get_triangle_center(self) -> Triangle:
        """
        Get the triangle object to the center of this triangle.

        Raises:
            ValueError: If there is no center triangle.

        Returns:
            Triangle: Triangle to the center.
        """        
        if not self.tc:
            raise ValueError("No triangle center.")
        
        return self.tc
    
    def get_triangles(self) -> Tuple[Triangle, Triangle, Triangle]:
        """
        Get the triangles to the left, right and center of this triangle.

        Raises:
            ValueError: If there are no triangles.

        Returns:
            Tuple[Triangle, Triangle, Triangle]: Triangles to the left, right and center.
        """
        if not self.tl:
            raise ValueError("No triangle left.")
        if not self.tr:
            raise ValueError("No triangle right.")
        if not self.tc:
            raise ValueError("No triangle center.")
        
        return self.tl, self.tr, self.tc

    def set_triangle_left(self, t: Triangle):
        """
        Set the triangle to the left of this triangle. Also ensures vice versa
        connections.

        Args:
            t (Triangle): Triangle to the left.
        """
        self.tl = t 

        # Ensure vice versa connection
        if t:
            t.tr = self

    def set_triangle_right(self, t: Triangle):
        """
        Set the triangle to the right of this triangle. Also ensures vice versa
        connections.

        Args:
            t (Triangle): Triangle to the right.
        """
        self.tr = t

        # Ensure vice versa connection
        if t:
            t.tl_ = self

    def set_triangle_center(self, t: Triangle):
        """
        Set the triangle to the center of this triangle. Also ensures vice versa
        connections.

        Args:
            t (Triangle): Triangle to the center.
        """
        self.tc = t

        # Ensure vice versa connection
        if t:
            t.tc = self

    def set_triangles(self, tl: Triangle, tr : Triangle, tc: Triangle):
        """
        Set the triangles to the left, right and center of this triangle. Also
        ensures vice versa connections.

        Args:
            tl (Triangle): Triangle to the left.
            tr (Triangle): Triangle to the right.
            tc (Triangle): Triangle to the center.
        """
        self.tl = tl
        self.tr = tr
        self.tc = tc

        # Ensure vice versa connections
        if tl:
            tl.tr = self
        if tr:
            tr.tl_ = self
        if tc:
            tc.tc = self

    def get_vertex_left(self) -> Vertex:
        """
        Get the left vertex of this triangle.

        Raises:
            ValueError: If there is no left vertex.

        Returns:
            Vertex: Left vertex.
        """
        if not self.vl:
            raise ValueError("No vertex left.")
        
        return self.vl

    def get_vertex_right(self) -> Vertex:
        """
        Get the right vertex of this triangle.

        Raises:
            ValueError: If there is no right vertex.

        Returns:
            Vertex: Right vertex.
        """
        if not self.vr:
            raise ValueError("No vertex right.")
        
        return self.vr

    def get_vertex_center(self) -> Vertex:
        """
        Get the center vertex of this triangle.

        Raises:
            ValueError: If there is no center vertex.

        Returns:
            Vertex: Center vertex.
        """
        if not self.vc:
            raise ValueError("No vertex center.")
        
        return self.vc

    def get_vertices(self) -> Tuple[Vertex, Vertex, Vertex]:
        """
        Get the vertices of this triangle.

        Raises:
            ValueError: If there are no vertices.

        Returns:
            Tuple[Vertex, Vertex, Vertex]: Left, right and center vertex.
        """  
        return self.vl, self.vr, self.vc

    def set_vertex_right(self, v: Vertex):
        """
        Set the right vertex of this triangle. Also updates the time of the triangle,
        and sets the triangle to the left of the right vertex.

        Args:
            v (Vertex): Right vertex.
        """
        self.vr = v
        self.time = v.time
        self.update_type()

        # Update vertex connections in space
        v.set_neighbour_left(self.vl)

        # If the triangle is upwards, update vertex-triangle connections
        if self.is_upwards():
            v.set_triangle_left(self)

    def set_vertices(self, vl: Vertex, vr: Vertex, vc: Vertex):
        """
        Set the vertices of this triangle. Also updates the time of the triangle,
        and sets the triangle to the right of the left vertex and the triangle to
        the left of the right vertex.

        Args:
            vl (Vertex): Left vertex.
            vr (Vertex): Right vertex.
            vc (Vertex): Center vertex.
        """
        self.vl = vl
        self.vr = vr
        self.vc = vc
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

    def update_type(self):
        """
        Update the type of the triangle based on the times of the vertices.
        """
        if self.vl.time < self.vc.time:
            self.type = "UP"
        else:
            self.type = "DOWN"

        # Correct for boundary conditions
        if self.vc.time == 0 and self.vl.time > 1:
            self.type = "UP"
        
        if self.vl.time == 0 and self.vc.time > 1:
            self.type = "DOWN"

    def clear_references(self):
        """
        Clear all references to other objects.
        """
        self.tl = None
        self.tr = None
        self.tc = None
        self.vl = None
        self.vr = None
        self.vc = None
        self.time = None
        self.type = None