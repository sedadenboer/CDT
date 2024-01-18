# vertex.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:

from triangle import Triangle


class Vertex:
    capacity = 100

    def __init__(self, ID, time):
        self.ID = ID
        self.time = time
        self.triangles = []
        # The two (2,1)-triangles sharing an edge with vertex v
        self.tr = None
        self.tl = None

    def set_triangle_left(self, t: "Triangle") -> None:
        """
        Set the (2,1)-triangle to the left of this vertex.

        Args:
            t (Triangle): (2,1)-triangle to set to the left of this vertex.
        """
        self.tl = t
    
    def set_triangle_right(self, t: "Triangle") -> None:
        """
        Set the (2,1)-triangle to the right of this vertex.

        Args:
            t (Triangle): (2,1)-triangle to set to the right of this vertex.
        """
        self.tr = t
    
    def get_triangle_left(self) -> "Triangle":
        """
        Get the (2,1)-triangle to the left of this vertex.

        Raises:
            ValueError: If there is no (2,1)-triangle to the left.

        Returns:
            Triangle: (2,1)-triangle to the left.
        """
        if not self.tl:
            raise ValueError("No (2,1)-triangle left.")
        
        return self.tl
    
    def get_triangle_right(self) -> "Triangle":
        """
        Get the (2,1)-triangle to the right of this vertex.

        Raises:
            ValueError: If there is no (2,1)-triangle to the right.

        Returns:
            Triangle: (2,1)-triangle to the right.
        """
        if not self.tr:
            raise ValueError("No (2,1)-triangle right.")
        
        return self.tr
