# vertex.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from triangle import Triangle


class Vertex(object):
    """
    Represents a vertex in the triangulation.

    Attributes:
        ID: Unique identifier for the vertex.
        time: Time of the vertex.
        tl: Triangle to the left.
        tr: Triangle to the right.
    """
    def __init__(self, time: int) -> None:
        self.time = time
        self.ID: int
        # The two (2,1)-triangles sharing an edge with vertex v
        self.tr_: Triangle
        self.tl_: Triangle

    def set_triangle_left(self, t: Triangle) -> None:
        """
        Set the (2,1)-triangle to the left of this vertex.

        Args:
            t (Triangle): (2,1)-triangle to set to the left of this vertex.
        """
        assert t.is_upwards() is not None, "Triangle is not upwards."
        self.tl_ = t
    
    def set_triangle_right(self, t: Triangle) -> None:
        """
        Set the (2,1)-triangle to the right of this vertex.

        Args:
            t (Triangle): (2,1)-triangle to set to the right of this vertex.
        """
        assert t.is_upwards() is not None, "Triangle is not upwards."
        self.tr_ = t
    
    def get_triangle_left(self) -> Triangle:
        """
        Get the (2,1)-triangle to the left of this vertex.

        Returns:
            Triangle: (2,1)-triangle to the left.
        """
        return self.tl_
    
    def get_triangle_right(self) -> Triangle:
        """
        Get the (2,1)-triangle to the right of this vertex.

        Returns:
            Triangle: (2,1)-triangle to the right.
        """      
        return self.tr_
