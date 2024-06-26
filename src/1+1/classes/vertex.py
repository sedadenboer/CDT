# vertex.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description: Represents a vertex in the triangulation.
# Contains methods to get and set the triangles and vertices of the vertex,
# as well as methods to add and remove neighbours and to clear all references
# to other objects.

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from triangle import Triangle


class Vertex:
    """
    Represents a vertex in the triangulation.

    Args (Attributes):
        time (int): Time of the vertex.

    Attributes:
        ID (int): Unique identifier for the vertex.
        tr (Triangle): Triangle to the right.
        tl (Triangle): Triangle to the left.
        past_neighbours (List[Vertex]): List of past neighbours.
        future_neighbours (List[Vertex]): List of future neighbours.
    """
    def __init__(self, time: int):
        self.time = time
        self.ID: int

        # The two (2,1)-triangles sharing an edge with vertex v
        self.tr: Triangle = None
        self.tl: Triangle = None

        # The two neighbouring vertices of vertex v
        self.vr: Vertex = None
        self.vl: Vertex = None
        
        self.past_neighbours = []
        self.future_neighbours = []

    def set_triangle_left(self, t: Triangle):
        """
        Set the (2,1)-triangle to the left of this vertex.

        Raises:
            AssertionError: If the triangle is not upwards.

        Args:
            t (Triangle): (2,1)-triangle to set to the left of this vertex.
        """
        assert t.is_upwards() is not None, "Triangle is not upwards."
        self.tl = t
    
    def set_triangle_right(self, t: Triangle):
        """
        Set the (2,1)-triangle to the right of this vertex.

        Raises:
            AssertionError: If the triangle is not upwards.

        Args:
            t (Triangle): (2,1)-triangle to set to the right of this vertex.
        """
        assert t.is_upwards() is not None, "Triangle is not upwards."
        self.tr = t
    
    def get_triangle_left(self) -> Triangle:
        """
        Get the (2,1)-triangle to the left of this vertex.

        Raises:
            ValueError: If there is no left triangle.
            
        Returns:
            Triangle: (2,1)-triangle to the left.
        """
        if not self.tl:
            raise ValueError("No triangle left.")
        
        return self.tl
    
    def get_triangle_right(self) -> Triangle:
        """
        Get the (2,1)-triangle to the right of this vertex.

        Raises:
            ValueError: If there is no right triangle.

        Returns:
            Triangle: (2,1)-triangle to the right.
        """
        if not self.tr:
            raise ValueError("No triangle right.")
        
        return self.tr

    def set_neighbour_left(self, v: Vertex):
        """
        Set the vertex to the left of this vertex.

        Args:
            v (Vertex): Vertex to set to the left of this vertex.
        """
        self.vl = v
        v.vr = self
    
    def set_neighbour_right(self, v: Vertex):
        """
        Set the vertex to the right of this vertex.

        Args:
            v (Vertex): Vertex to set to the right of this vertex.
        """
        self.vr = v
        v.vl = self
    
    def get_neighbour_left(self) -> Vertex:
        """
        Get the vertex to the left of this vertex.

        Raises:
            ValueError: If there is no left vertex.

        Returns:
            Vertex: Vertex to the left.
        """
        if not self.vl:
            raise ValueError("No vertex left.")
        
        return self.vl
    
    def get_neighbour_right(self) -> Vertex:
        """
        Get the vertex to the right of this vertex.

        Raises:
            ValueError: If there is no right vertex.

        Returns:
            Vertex: Vertex to the right.
        """
        if not self.vr:
            raise ValueError("No vertex right.")
        
        return self.vr

    def get_space_neighbours(self) -> list[Vertex]:
        """
        Get the space neighbours of the vertex.

        Returns:
            list[Vertex]: Space neighbours of the vertex (vl, vr).
        """
        return self.get_neighbour_left(), self.get_neighbour_right()
    
    def add_future_neighbour(self, v: Vertex):
        """
        Add a future neighbour to the vertex and vice versa.

        Args:
            v (Vertex): Vertex to add as future neighbour.
        """
        if not v in self.future_neighbours:
            self.future_neighbours.append(v)
            v.past_neighbours.append(self)

    def add_past_neighbour(self, v: Vertex):
        """
        Add a past neighbour to the vertex and vice versa.

        Raises:
            AssertionError: If the past neighbour has the wrong time.

        Args:
            v (Vertex): Vertex to add as past neighbour.
        """
        # Check if past neighbour has the right time, considering it is not at the boundary
        if self.time != 0:
            assert v.time == self.time - 1, "Past neighbour has wrong time."
        
        if not v in self.past_neighbours:
            self.past_neighbours.append(v)
            v.future_neighbours.append(self)
    
    def delete_future_neighbour(self, v: Vertex):
        """
        Delete a future neighbour from the vertex and vice versa.

        Args:
            v (Vertex): Vertex to delete as future neighbour.
        """
        if v in self.future_neighbours:
            self.future_neighbours.remove(v)
            v.past_neighbours.remove(self)
    
    def delete_past_neighbour(self, v: Vertex):
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
    
    def clear_references(self):
        """
        Clears all the references of the vertex.
        """
        self.tl = None
        self.tr = None
        self.vr = None
        self.vl = None
        self.past_neighbours = []
        self.future_neighbours = []