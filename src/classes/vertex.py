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
        self.tr_: Triangle = None
        self.tl_: Triangle = None
        self.vr_: Vertex = None
        self.vl_: Vertex = None
        self.past_neighbours = []
        self.future_neighbours = []

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

        Raises:
            ValueError: If there is no left triangle.
            
        Returns:
            Triangle: (2,1)-triangle to the left.
        """
        if not self.tl_:
            raise ValueError("No triangle left.")
        
        return self.tl_
    
    def get_triangle_right(self) -> Triangle:
        """
        Get the (2,1)-triangle to the right of this vertex.

        Raises:
            ValueError: If there is no right triangle.

        Returns:
            Triangle: (2,1)-triangle to the right.
        """
        if not self.tr_:
            raise ValueError("No triangle right.")
        
        return self.tr_

    def set_neighbour_left(self, v: Vertex) -> None:
        """
        Set the vertex to the left of this vertex.

        Args:
            v (Vertex): Vertex to set to the left of this vertex.
        """
        self.vl_ = v
        v.vr_ = self
    
    def set_neighbour_right(self, v: Vertex) -> None:
        """
        Set the vertex to the right of this vertex.

        Args:
            v (Vertex): Vertex to set to the right of this vertex.
        """
        self.vr_ = v
        v.vl_ = self
    
    def get_neighbour_left(self) -> Vertex:
        """
        Get the vertex to the left of this vertex.

        Raises:
            ValueError: If there is no left vertex.

        Returns:
            Vertex: Vertex to the left.
        """
        if not self.vl_:
            raise ValueError("No vertex left.")
        
        return self.vl_
    
    def get_neighbour_right(self) -> Vertex:
        """
        Get the vertex to the right of this vertex.

        Raises:
            ValueError: If there is no right vertex.

        Returns:
            Vertex: Vertex to the right.
        """
        if not self.vr_:
            raise ValueError("No vertex right.")
        
        return self.vr_

    def get_space_neighbours(self) -> list[Vertex]:
        """
        Get the space neighbours of the vertex.

        Returns:
            list[Vertex]: Space neighbours of the vertex (vl, vr).
        """
        return self.get_neighbour_left(), self.get_neighbour_right()
    
    def add_future_neighbour(self, v: Vertex) -> None:
        """
        Add a future neighbour to the vertex and vice versa.

        Args:
            v (Vertex): Vertex to add as future neighbour.
        """
        if not v in self.future_neighbours:
            self.future_neighbours.append(v)
            v.past_neighbours.append(self)

    def add_past_neighbour(self, v: Vertex) -> None:
        """
        Add a past neighbour to the vertex and vice versa.

        Args:
            v (Vertex): Vertex to add as past neighbour.
        """
        if self.time != 0:
            assert v.time == self.time - 1, "Past neighbour has wrong time."
        
        if not v in self.past_neighbours:
            self.past_neighbours.append(v)
            v.future_neighbours.append(self)
    
    def delete_future_neighbour(self, v: Vertex) -> None:
        """
        Delete a future neighbour from the vertex and vice versa.

        Args:
            v (Vertex): Vertex to delete as future neighbour.
        """
        if v in self.future_neighbours:
            self.future_neighbours.remove(v)
            v.past_neighbours.remove(self)
    
    def delete_past_neighbour(self, v: Vertex) -> None:
        """
        Delete a past neighbour from the vertex and vice versa.

        Args:
            v (Vertex): Vertex to delete as past neighbour.
        """
        if v in self.past_neighbours:
            self.past_neighbours.remove(v)
            v.future_neighbours.remove(self)

    def get_future_neighbours(self) -> list[Vertex]:
        """
        Get the future neighbours of the vertex.

        Returns:
            list[Vertex]: Future neighbours of the vertex.
        """
        return self.future_neighbours
    
    def get_past_neighbours(self) -> list[Vertex]:
        """
        Get the past neighbours of the vertex.

        Returns:
            list[Vertex]: Past neighbours of the vertex.
        """
        return self.past_neighbours