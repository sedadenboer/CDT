# vertex.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:


class Vertex(object):
    capacity = 30

    def __init__(self, time: int):
        self.ID = None
        self.time = time
        self.triangles = []
        # The two (2,1)-triangles sharing an edge with vertex v
        self.tr = None
        self.tl = None

    def set_triangle_left(self, t: object) -> None:
        """
        Set the (2,1)-triangle to the left of this vertex.

        Args:
            t (object): (2,1)-triangle to set to the left of this vertex.
        """
        self.tl = t
    
    def set_triangle_right(self, t: object) -> None:
        """
        Set the (2,1)-triangle to the right of this vertex.

        Args:
            t (object): (2,1)-triangle to set to the right of this vertex.
        """
        self.tr = t
    
    def get_triangle_left(self) -> object:
        """
        Get the (2,1)-triangle to the left of this vertex.

        Returns:
            object: (2,1)-triangle to the left.
        """
        return self.tl
    
    def get_triangle_right(self) -> object:
        """
        Get the (2,1)-triangle to the right of this vertex.

        Returns:
            object: (2,1)-triangle to the right.
        """      
        return self.tr
