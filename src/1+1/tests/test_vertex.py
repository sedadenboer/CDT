# test_vertex.py
#
# Author: Seda den Boer
# Date: 01-02-2024
#
# Description: Test file for the Vertex class.

import pytest
from classes.vertex import Vertex
from classes.triangle import Triangle


def test_vertex_initialization():
    """
    Test if a vertex is initialized correctly.
    """
    time = 0
    vertex = Vertex(time)
    assert vertex.time == time
    assert not vertex.tr_ 
    assert not vertex.tl_ 

def test_set_triangle_left():
    """
    Test if a triangle is set correctly to the left of a vertex.
    """
    vertex = Vertex(time=0)
    triangle_left = Triangle()
    triangle_left.type = "UP"
    vertex.set_triangle_left(triangle_left)
    assert vertex.tl_ == triangle_left

def test_set_triangle_right():
    """
    Test if a triangle is set correctly to the right of a vertex.
    """
    vertex = Vertex(time=0)
    triangle_right = Triangle()
    triangle_right.type = "UP"
    vertex.set_triangle_right(triangle_right)
    assert vertex.tr_ == triangle_right

def test_get_triangle_left():
    """
    Test if a triangle is set correctly to the left of a vertex
    and if it is returned correctly.
    """
    vertex = Vertex(time=0)
    triangle_left = Triangle()
    triangle_left.type = "UP"
    vertex.set_triangle_left(triangle_left)
    assert vertex.get_triangle_left() == triangle_left

def test_get_triangle_right():
    """
    Test if a triangle is set correctly to the right of a vertex
    and if it is returned correctly.
    """
    vertex = Vertex(time=0)
    triangle_right = Triangle()
    triangle_right.type = "UP"
    vertex.set_triangle_right(triangle_right)
    assert vertex.get_triangle_right() == triangle_right

def test_value_errors():
    """
    Test if the correct value errors are raised.
    """
    vertex = Vertex(time=0)
    
    with pytest.raises(ValueError, match="No triangle left."):
        vertex.get_triangle_left()
    
    with pytest.raises(ValueError, match="No triangle right."):
        vertex.get_triangle_right()

def test_space_neighbours():
    """
    Test if the space neighbours are set correctly.
    """
    v1 = Vertex(1)
    v2 = Vertex(2)
    v3 = Vertex(3)

    v1.set_neighbour_left(v2)
    v1.set_neighbour_right(v3)

    assert v1.get_neighbour_left() == v2
    assert v1.get_neighbour_right() == v3
    assert v2.get_neighbour_right() == v1
    assert v3.get_neighbour_left() == v1

def test_future_past_neighbours():
    """
    Test if the future and past neighbours are set and deleted correctly.
    """
    v1 = Vertex(1)
    v2 = Vertex(2)
    v3 = Vertex(3)

    v1.add_future_neighbour(v2)
    v1.add_past_neighbour(v3)

    assert v1.get_future_neighbours() == [v2]
    assert v1.get_past_neighbours() == [v3]
    assert v2.get_past_neighbours() == [v1]
    assert v3.get_future_neighbours() == [v1]

    v1.delete_future_neighbour(v2)
    v3.delete_future_neighbour(v1)

    assert v1.get_future_neighbours() == []
    assert v1.get_past_neighbours() == []
    assert v2.get_past_neighbours() == []
    assert v3.get_future_neighbours() == []
