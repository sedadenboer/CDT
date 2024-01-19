# triangle.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:


class Triangle(object):
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
    capacity = 30

    def __init__(self) -> None:
        self.ID = None
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

    def get_triangle_left(self) -> object:
        """
        Get the triangle object to the left of this triangle.

        Returns:
            object: Triangle to the left.
        """
        return self.tl

    def get_triangle_right(self) -> object:
        """
        Get the triangle object to the right of this triangle.

        Returns:
            object: Triangle to the right.
        """
        return self.tr

    def get_triangle_center(self) -> object:
        """
        Get the triangle object to the center of this triangle.

        Returns:
            object: Triangle to the center.
        """        
        return self.tc

    def set_triangle_left(self, t: object) -> None:
        """
        Set the triangle to the left of this triangle.

        Args:
            t (object): Triangle to the left.
        """
        self.tl = t 

    def set_triangle_right(self, t: object) -> None:
        """
        Set the triangle to the right of this triangle.

        Args:
            t (object): Triangle to the right.
        """
        self.tr = t

    def set_triangle_center(self, t: object) -> None:
        """
        Set the triangle to the center of this triangle.

        Args:
            t (object): Triangle to the center.
        """
        self.tc = t

    def set_triangles(self, tl: object, tr : object, tc: object) -> None:
        """
        Set the triangles to the left, right and center of this triangle.

        Args:
            tl (object): Triangle to the left.
            tr (object): Triangle to the right.
            tc (object): Triangle to the center.
        """
        self.tl = tl
        self.tr = tr
        self.tc = tc

    def get_vertex_left(self) -> object:
        """
        Get the left vertex of this triangle.

        Raises:
            ValueError: If there is no left vertex.

        Returns:
            object: Left vertex.
        """
        if not self.vl:
            raise ValueError("No vertex left.")
        
        return self.vl

    def get_vertex_right(self) -> object:
        """
        Get the right vertex of this triangle.

        Raises:
            ValueError: If there is no right vertex.

        Returns:
            object: Right vertex.
        """
        if not self.vr:
            raise ValueError("No vertex right.")
        
        return self.vr

    def get_vertex_center(self) -> object:
        """
        Get the center vertex of this triangle.

        Raises:
            ValueError: If there is no center vertex.

        Returns:
            object: Center vertex.
        """
        if not self.vc:
            raise ValueError("No vertex center.")
        
        return self.vc

    def set_vertex_left(self, v: object) -> None:
        """
        Set the left vertex of this triangle.

        Args:
            v (object): Left vertex.
        """
        self.vl = v
        self.time = v.time
        self.update_type()

    def set_vertex_right(self, v: object) -> None:
        """
        Set the right vertex of this triangle.

        Args:
            v (object): Right vertex.
        """
        self.vr = v
        self.time = v.time
        self.update_type()

    def set_vertex_center(self, v: object) -> None:
        """
        Set the center vertex of this triangle.

        Args:
            v (object): Center vertex.
        """
        self.vc = v
        self.update_type()
    
    def set_vertices(self, vl: object, vr: object, vc: object) -> None:
        """
        Set the vertices of this triangle.

        Args:
            vl (object): Left vertex.
            vr (object): Right vertex.
            vc (object): Center vertex.
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

        
# Test
if __name__ == "__main__":
    triangles = Triangle()
    for _ in range(5):
        triangles.create()

    triangles.destroy(0)
    triangles.destroy(1)
    triangles.destroy(4)
