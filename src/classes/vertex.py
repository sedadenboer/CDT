# vertex.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:


class Vertex(object):

    def __init__(self, time: int):
        self.ID = None
        self.time = time
        self.tr_iangles = []
        # The two (2,1)-triangles sharing an edge with vertex v
        self.tr_ = None
        self.tl_ = None

    def set_triangle_left(self, t: object) -> None:
        """
        Set the (2,1)-triangle to the left of this vertex.

        Args:
            t (object): (2,1)-triangle to set to the left of this vertex.
        """
        self.tl_ = t
    
    def set_triangle_right(self, t: object) -> None:
        """
        Set the (2,1)-triangle to the right of this vertex.

        Args:
            t (object): (2,1)-triangle to set to the right of this vertex.
        """
        self.tr_ = t
    
    def get_triangle_left(self) -> object:
        """
        Get the (2,1)-triangle to the left of this vertex.

        Returns:
            object: (2,1)-triangle to the left.
        """
        return self.tl_
    
    def get_triangle_right(self) -> object:
        """
        Get the (2,1)-triangle to the right of this vertex.

        Returns:
            object: (2,1)-triangle to the right.
        """      
        return self.tr_
