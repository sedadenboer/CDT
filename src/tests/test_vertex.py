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