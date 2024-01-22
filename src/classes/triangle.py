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

    def __init__(self) -> None:
        self.ID = None
        self.type = None

        # The three triangles sharing and edge with triangle t
        self.tl_ = None # left
        self.tr_ = None # right
        self.tc_ = None # vertical

        # The three vertices comprising the triangle t
        self.vl_ = None # left
        self.vr_ = None # right
        self.vc_ = None # apex

        # Set the time of the triangle to the time of the base of the triangle
        self.time = self.vl.time if self.vl else None

    def get_triangle_left(self) -> object:
        """
        Get the triangle object to the left of this triangle.

        Returns:
            object: Triangle to the left.
        """
        return self.tl_

    def get_triangle_right(self) -> object:
        """
        Get the triangle object to the right of this triangle.

        Returns:
            object: Triangle to the right.
        """
        return self.tr_

    def get_triangle_center(self) -> object:
        """
        Get the triangle object to the center of this triangle.

        Returns:
            object: Triangle to the center.
        """        
        return self.tc_

    def set_triangle_left(self, t: object) -> None:
        """
        Set the triangle to the left of this triangle.

        Args:
            t (object): Triangle to the left.
        """
        self.tl_ = t 

    def set_triangle_right(self, t: object) -> None:
        """
        Set the triangle to the right of this triangle.

        Args:
            t (object): Triangle to the right.
        """
        self.tr_ = t

    def set_triangle_center(self, t: object) -> None:
        """
        Set the triangle to the center of this triangle.

        Args:
            t (object): Triangle to the center.
        """
        self.tc_ = t

    def set_triangles(self, tl: object, tr : object, tc: object) -> None:
        """
        Set the triangles to the left, right and center of this triangle.

        Args:
            tl (object): Triangle to the left.
            tr (object): Triangle to the right.
            tc (object): Triangle to the center.
        """
        self.tl_ = tl
        self.tr_ = tr
        self.tc_ = tc

    def get_vertex_left(self) -> object:
        """
        Get the left vertex of this triangle.

        Raises:
            ValueError: If there is no left vertex.

        Returns:
            object: Left vertex.
        """
        if not self.vl_:
            raise ValueError("No vertex left.")
        
        return self.vl_

    def get_vertex_right(self) -> object:
        """
        Get the right vertex of this triangle.

        Raises:
            ValueError: If there is no right vertex.

        Returns:
            object: Right vertex.
        """
        if not self.vr_:
            raise ValueError("No vertex right.")
        
        return self.vr_

    def get_vertex_center(self) -> object:
        """
        Get the center vertex of this triangle.

        Raises:
            ValueError: If there is no center vertex.

        Returns:
            object: Center vertex.
        """
        if not self.vc_:
            raise ValueError("No vertex center.")
        
        return self.vc_

    def set_vertex_left(self, v: object) -> None:
        """
        Set the left vertex of this triangle.

        Args:
            v (object): Left vertex.
        """
        self.vl_ = v
        self.time = v.time
        self.update_type()

    def set_vertex_right(self, v: object) -> None:
        """
        Set the right vertex of this triangle.

        Args:
            v (object): Right vertex.
        """
        self.vr_ = v
        self.time = v.time
        self.update_type()

    def set_vertex_center(self, v: object) -> None:
        """
        Set the center vertex of this triangle.

        Args:
            v (object): Center vertex.
        """
        self.vc_ = v
        self.update_type()
    
    def set_vertices(self, vl: object, vr: object, vc: object) -> None:
        """
        Set the vertices of this triangle.

        Args:
            vl (object): Left vertex.
            vr (object): Right vertex.
            vc (object): Center vertex.
        """
        self.vl_ = vl
        self.vr_ = vr
        self.vc_ = vc
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
        if self.vl_.time < self.vc_.time:
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
