# test_triangle.py
#
# Author: Seda den Boer
# Date: 01-02-2024
#
# Description: Test file for the Triangle class.

import pytest
from classes.triangle import Triangle
from classes.vertex import Vertex


def test_triangle_initialization():
    """
    Test if a triangle is initialized correctly.
    """
    triangle = Triangle()
    assert not triangle.tl_
    assert not triangle.tr_ 
    assert not triangle.tc_ 
    assert not triangle.vl_ 
    assert not triangle.vr_ 
    assert not triangle.vc_ 
    assert not triangle.time 

def test_set_triangle_left():
    """
    Test if a triangle is set correctly to the left of another triangle.
    """
    triangle1 = Triangle()
    triangle2 = Triangle()
    triangle1.set_triangle_left(triangle2)
    assert triangle1.tl_ == triangle2
    assert triangle2.tr_ == triangle1

def test_set_triangle_right():
    """
    Test if a triangle is set correctly to the right of another triangle.
    """
    triangle1 = Triangle()
    triangle2 = Triangle()
    triangle1.set_triangle_right(triangle2)
    assert triangle1.tr_ == triangle2
    assert triangle2.tl_ == triangle1

def test_set_triangle_center():
    """
    Test if a triangle is set correctly to the center of another triangle.
    """
    triangle1 = Triangle()
    triangle2 = Triangle()
    triangle1.set_triangle_center(triangle2)
    assert triangle1.tc_ == triangle2
    assert triangle2.tc_ == triangle1

def test_set_triangles():
    """
    Test if a triangle is set correctly to the left, right and center of another triangle.
    """
    triangle1 = Triangle()
    triangle2 = Triangle()
    triangle3 = Triangle()
    triangle1.set_triangles(triangle2, triangle3, None)
    assert triangle1.tl_ == triangle2
    assert triangle1.tr_ == triangle3
    assert triangle2.tr_ == triangle1
    assert triangle3.tl_ == triangle1

def test_set_vertex_left():
    """
    Test if a vertex is set correctly to the left of a triangle.
    """
    triangle = Triangle()
    triangle.type = "UP"
    vertex = Vertex(time=0)
    triangle.set_vertex_left(vertex)
    assert triangle.vl_ == vertex
    assert triangle.time == vertex.time

def test_set_vertex_right():
    """
    Test if a vertex is set correctly to the right of a triangle.
    """
    triangle = Triangle()
    triangle.type = "UP"
    vertex = Vertex(time=0)
    triangle.set_vertex_right(vertex)
    assert triangle.vr_ == vertex
    assert triangle.time == vertex.time

def test_set_vertex_center():
    """
    Test if a vertex is set correctly to the center of a triangle.
    """
    triangle = Triangle()
    triangle.type = "UP"
    vertex = Vertex(time=0)
    triangle.set_vertex_center(vertex)
    assert triangle.vc_ == vertex

def test_set_vertices():
    """
    Test if a vertex is set correctly to the left, right and center of a triangle.
    """
    triangle = Triangle()
    triangle.type = "UP"
    vertex_left = Vertex(time=0)
    vertex_right = Vertex(time=0)
    vertex_center = Vertex(time=1)
    triangle.set_vertices(vertex_left, vertex_right, vertex_center)
    assert triangle.vl_ == vertex_left
    assert triangle.vr_ == vertex_right
    assert triangle.vc_ == vertex_center
    assert triangle.time == vertex_left.time
    assert vertex_left.tr_ == triangle
    assert vertex_right.tl_ == triangle

def test_update_type():
    """
    Test if the type of a triangle is updated correctly.
    """
    triangle = Triangle()
    triangle.set_vertices(Vertex(time=0), Vertex(time=0), Vertex(time=1))
    triangle.update_type()
    assert triangle.is_downwards() == False
    assert triangle.is_upwards() == True
    triangle.vl_.time = 1
    triangle.vr_.time = 1
    triangle.vc_.time = 0
    triangle.update_type()
    assert triangle.is_upwards() == False
    assert triangle.is_downwards() == True

def test_value_errors():
    """
    Test if the correct errors are raised when a triangle or vertex is not set.
    """
    triangle = Triangle()

    with pytest.raises(ValueError, match="No triangle left."):
        triangle.get_triangle_left()

    with pytest.raises(ValueError, match="No triangle right."):
        triangle.get_triangle_right()

    with pytest.raises(ValueError, match="No triangle center."):
        triangle.get_triangle_center()

    with pytest.raises(ValueError, match="No triangle left."):
        triangle.get_triangles()

    with pytest.raises(ValueError, match="No vertex left."):
        triangle.get_vertex_left()

    with pytest.raises(ValueError, match="No vertex right."):
        triangle.get_vertex_right()

    with pytest.raises(ValueError, match="No vertex center."):
        triangle.get_vertex_center()