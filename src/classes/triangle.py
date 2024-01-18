# triangle.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:

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
    capacity = 100

    def __init__(self, ID: int) -> None:
        self.ID = ID
        self.type = None

        # The three triangles sharing and edge with triangle t
        self.tl = None # left
        self.tr = None # right
        self.tc = None # vertical

        # The three vertices comprising the triangle t
        self.vl = None # left
        self.vr = None # right
        self.vc = None # apex

        # Set the time of the triangle to the time of the base of the triangle
        self.time = self.vl.time if self.vl else None

    def get_triangle_left(self) -> "Triangle":
        """
        Get the triangle object to the left of this triangle.

        Raises:
            ValueError: If there is no triangle to the left.

        Returns:
            Triangle: Triangle to the left.
        """
        if not self.tl:
            raise ValueError("No triangle left.")
        
        return self.tl

    def get_triangle_right(self) -> "Triangle":
        """
        Get the triangle object to the right of this triangle.

        Raises:
            ValueError: If there is no triangle to the right.

        Returns:
            Triangle: Triangle to the right.
        """
        if not self.tr:
            raise ValueError("No triangle right.")
        
        return self.tr

    def get_triangle_center(self) -> "Triangle":
        """
        Get the triangle object to the center of this triangle.
        
        Raises:
            ValueError: If there is no triangle to the center.

        Returns:
            Triangle: Triangle to the center.
        """
        if not self.tc:
            raise ValueError("No triangle center.")
        
        return self.tc

    def set_triangle_left(self, t: "Triangle") -> None:
        """
        Set the triangle to the left of this triangle.

        Args:
            t (Triangle): Triangle to the left.
        """
        self.tl = t 

    def set_triangle_right(self, t: "Triangle") -> None:
        """
        Set the triangle to the right of this triangle.

        Args:
            t (Triangle): Triangle to the right.
        """
        self.tr = t

    def set_triangle_center(self, t: "Triangle") -> None:
        """
        Set the triangle to the center of this triangle.

        Args:
            t (Triangle): Triangle to the center.
        """
        self.tc = t

    def set_triangles(self, tl: "Triangle", tr : "Triangle", tc: "Triangle") -> None:
        """
        Set the triangles to the left, right and center of this triangle.

        Args:
            tl (Triangle): Triangle to the left.
            tr (Triangle): Triangle to the right.
            tc (Triangle): Triangle to the center.
        """
        self.tl = tl
        self.tr = tr
        self.tc = tc

    def get_vertex_left(self) -> "Vertex":
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

    def get_vertex_right(self) -> "Vertex":
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

    def get_vertex_center(self) -> "Vertex":
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

    def set_vertex_left(self, v: "Vertex") -> None:
        """
        Set the left vertex of this triangle.

        Args:
            v (Vertex): Left vertex.
        """
        self.vl = v
        self.time = v.time
        self.update_type()

    def set_vertex_right(self, v: "Vertex") -> None:
        """
        Set the right vertex of this triangle.

        Args:
            v (Vertex): Right vertex.
        """
        self.vr = v
        self.time = v.time
        self.update_type()

    def set_vertex_center(self, v: "Vertex") -> None:
        """
        Set the center vertex of this triangle.

        Args:
            v (Vertex): Center vertex.
        """
        self.vc = v
        self.update_type()
    
    def set_vertices(self, vl: "Vertex", vr: "Vertex", vc: "Vertex") -> None:
        """
        Set the vertices of this triangle.

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

    def is_upwards(self):
        """
        Return True if the triangle is upwards, False otherwise.
        """
        return self.type == "UP"

    def is_downwards(self):
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

        



# triangles = Triangle()
# for _ in range(5):
#     triangles.create()

# triangles.destroy(0)
# triangles.destroy(1)
# triangles.destroy(4)
