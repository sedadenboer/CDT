# test_universe.py
#
# Author: Seda den Boer
# Date: 31-01-2024
#
# Description: Test file for the Universe class.

import pytest
import copy
import numpy as np
from classes.universe import Universe


@pytest.fixture
def universe():
    return Universe(total_time=4, initial_slice_size=4, VERTEX_CAPACITY=20)

def test_initialization(universe):
    """
    Test the initialization of the Universe class by checking
    the number of triangles and vertices created, the orientation
    of the triangles and if the triangles are correctly added to
    the vertex objects.
    """
    # Check if the number of triangles and vertices is correct
    assert universe.triangle_pool.get_number_occupied() == 32
    assert universe.vertex_pool.get_number_occupied() == 16

    # Check the number of triangles in each bag
    assert universe.triangle_add_bag.get_number_occupied() == 32
    assert universe.triangle_flip_bag.get_number_occupied() == 32
    assert universe.four_vertices_bag.get_number_occupied() == 0

    # Check the orientation of triangles
    for id in np.arange(0, 32, 2):
        triangle = universe.triangle_pool.get(index=id)
        assert triangle.is_upwards()
    
    for id in np.arange(1, 33, 2):
        triangle = universe.triangle_pool.get(index=id)
        assert triangle.is_downwards()

    # Check if triangles are correctly added to vertex objects
    v_0 = universe.vertex_pool.get(index=0)
    v_1 = universe.vertex_pool.get(index=1)
    v_4 = universe.vertex_pool.get(index=4)
    v_5 = universe.vertex_pool.get(index=5)
    v_15 = universe.vertex_pool.get(index=15)

    assert v_0.get_triangle_left() == universe.triangle_pool.get(index=6)
    assert v_0.get_triangle_right() == universe.triangle_pool.get(index=0)
    assert v_1.get_triangle_left() == universe.triangle_pool.get(index=0)
    assert v_1.get_triangle_right() == universe.triangle_pool.get(index=2)
    assert v_4.get_triangle_left() == universe.triangle_pool.get(index=14)
    assert v_4.get_triangle_right() == universe.triangle_pool.get(index=8)
    assert v_5.get_triangle_left() == universe.triangle_pool.get(index=8)
    assert v_5.get_triangle_right() == universe.triangle_pool.get(index=10)
    assert v_15.get_triangle_left() == universe.triangle_pool.get(index=28)
    assert v_15.get_triangle_right() == universe.triangle_pool.get(index=30)

def test_insert_vertex(universe):
    """
    Test insert_vertex() operation in the Universe class by
    checking if after the insertion the vertices and triangles
    are correctly updated, if the correct number of triangles
    and vertices are created and if the correct number of
    triangles are added to the triangle_add_bag and triangle_flip_bag.
    """
    # Get triangle with id 8
    t = universe.triangle_pool.get(index=8)

    tl, tr, tc = t.get_triangles()
    vl, vr, vc = t.get_vertices()
    tc_tl, tc_tr, tc_tc = tc.get_triangles()
    tc_vl, tc_vr, tc_vc = tc.get_vertices()

    # Perform insert_vertex operation
    inserted_vertex, new_t1, new_t2 = universe.insert_vertex(triangle_id=8)

    # Check if vertices are updated correctly
    assert t.get_vertex_right() == inserted_vertex
    assert tc.get_vertex_right() == inserted_vertex
    assert new_t1.get_vertex_right() == vr
    assert new_t1.get_vertex_center() == vc
    assert new_t1.get_vertex_left() == inserted_vertex
    assert new_t2.get_vertex_left() == inserted_vertex
    assert new_t2.get_vertex_right() == vr
    assert new_t2.get_vertex_center() == tc_vc

    # Check if neighbors are updated correctly
    assert t.get_triangle_right() == new_t1
    assert tc.get_triangle_right() == new_t2
    assert new_t1.get_triangle_left() == t
    assert new_t2.get_triangle_left() == tc

    # Check if pools are updated correctly
    assert universe.vertex_pool.get_number_occupied() == 17
    assert universe.triangle_pool.get_number_occupied() == 34
    assert universe.vertex_pool.contains(inserted_vertex.ID)
    assert universe.triangle_pool.contains(new_t1.ID)
    assert universe.triangle_pool.contains(new_t2.ID)

    # Check if bags are updated correctly
    assert universe.triangle_add_bag.contains(new_t1.ID)
    assert universe.triangle_add_bag.contains(new_t2.ID)
    assert universe.triangle_flip_bag.contains(new_t1.ID)
    assert universe.triangle_flip_bag.contains(new_t2.ID)
    assert universe.four_vertices_bag.contains(inserted_vertex.ID)

def test_remove_vertex(universe):
    """
    Test remove_vertex() in the Universe class by
    checking if after the removal the vertices and triangles
    are correctly updated, if the correct number of triangles
    and vertices are created and if the correct (number of)
    triangles are added to the triangle_add_bag and triangle_flip_bag.
    """
    # Save the initial state for comparison
    initial_triangle_add_bag = copy.deepcopy(universe.triangle_add_bag.used_indices)
    initial_triangle_flip_bag = copy.deepcopy(universe.triangle_flip_bag.used_indices)
    initial_four_vertices_bag = copy.deepcopy(universe.four_vertices_bag.used_indices)

    # Get triangle with id 8
    t = universe.triangle_pool.get(index=8)

    tl, tr, tc = t.get_triangles()
    vl, vr, vc = t.get_vertices()
    tc_tl, tc_tr, tc_tc = tc.get_triangles()
    tc_vl, tc_vr, tc_vc = tc.get_vertices()

    # First perform insert_vertex operation
    inserted_vertex, new_t1, new_t2 = universe.insert_vertex(triangle_id=8)

    # Then perform remove_vertex operation
    removed_vertex, removed_t1, removed_t2 = universe.remove_vertex(vertex_id=inserted_vertex.ID)

    t = universe.triangle_pool.get(index=8)
    tc = t.get_triangle_center()

    # Check if vertices are updated correctly
    assert t.get_vertex_left() == vl
    assert t.get_vertex_right() == vr 
    assert t.get_vertex_center() == vc
    assert tc.get_vertex_left() == tc_vl
    assert tc.get_vertex_right() == tc_vr
    assert tc.get_vertex_center() == tc_vc

    # Check if neighbors are updated correctly
    assert t.get_triangle_right() == tr
    assert t.get_triangle_center() == tc
    assert t.get_triangle_left() == tl
    assert tc.get_triangle_right() == tc_tr
    assert tc.get_triangle_left() == tc_tl
    assert tc.get_triangle_center() == tc_tc
    assert inserted_vertex == removed_vertex
    assert new_t1 == removed_t1
    assert new_t2 == removed_t2

    # Check if pools are updated correctly
    assert universe.vertex_pool.get_number_occupied() == 16
    assert universe.triangle_pool.get_number_occupied() == 32
    assert not universe.vertex_pool.contains(removed_vertex.ID)
    assert not universe.triangle_pool.contains(removed_t1.ID)
    assert not universe.triangle_pool.contains(removed_t2.ID)

    # Check if bags are updated correctly
    assert universe.triangle_add_bag.used_indices == initial_triangle_add_bag
    assert universe.triangle_flip_bag.used_indices == initial_triangle_flip_bag
    assert universe.four_vertices_bag.used_indices == initial_four_vertices_bag

def test_flip_edge(universe):
    """
    Test flip_edge() in the Universe class by checking if after 
    two times flipping the same edge the triangles and vertices
    are correctly updated.
    """
    # Get triangle with id 8
    t = universe.triangle_pool.get(index=8)

    tl, tr, tc = t.get_triangles()
    vl, vr, vc = t.get_vertices()
    tr_tl, tr_tr, tr_tc = tr.get_triangles()
    tr_vl, tr_vr, tr_vc = tr.get_vertices()

    # Perform flip_edge operation
    flipped_t1, flipped_t2 = universe.flip_edge(triangle_id=8)

    # Check if vertices are updated correctly
    assert t.get_vertex_right() == tr_vr
    assert t.get_vertex_left() == vc
    assert t.get_vertex_center() == vl
    assert tr.get_vertex_right() == vr
    assert tr.get_vertex_left() == vl
    assert tr.get_vertex_center() == tr_vr

    # Check if neighbors are updated correctly
    assert t.get_triangle_left() == tl
    assert t.get_triangle_center() == tr_tc
    assert tr.get_triangle_center() == tc
    assert t.get_triangle_right() == flipped_t2
    assert tr.get_triangle_left() == flipped_t1
    
    # Flip back
    flipped_t3, flipped_t4 = universe.flip_edge(triangle_id=8)

    assert flipped_t3 == t
    assert flipped_t4 == tr
    assert flipped_t3.get_triangles() == (tl, tr, tc)
    assert flipped_t4.get_triangles() == (t, tr_tr, tr_tc)

def test_flip_four_vertices(universe):
    """
    Test if after the flip move the four vertices bag is updated correctly.
    """
    centered_vertex = universe.vertex_pool.get(index=9)

    t_16, t_17 = universe.flip_edge(triangle_id=16)
    assert not universe.four_vertices_bag.contains(t_17.get_vertex_right().ID)

    t_10, t_11 = universe.flip_edge(triangle_id=10)

    # Check the IDs of the vertices of the triangles
    assert t_10.get_vertex_left().ID == 9
    assert t_10.get_vertex_right().ID == 10
    assert t_10.get_vertex_center().ID == 5
    assert t_11.get_vertex_left().ID == 5
    assert t_11.get_vertex_right().ID == 6
    assert t_11.get_vertex_center().ID == 10
    assert t_16.get_vertex_left().ID == 12
    assert t_16.get_vertex_right().ID == 13
    assert t_16.get_vertex_center().ID == 8
    assert t_17.get_vertex_left().ID == 8
    assert t_17.get_vertex_right().ID == 9
    assert t_17.get_vertex_center().ID == 13

    # Check the triangles saved in vertex 9 are correct
    assert centered_vertex.get_triangle_left() == t_17
    assert centered_vertex.get_triangle_right() == universe.triangle_pool.get(index=18)

    # Check if bags are updated correctly
    assert t_10.get_vertex_left() == t_17.get_vertex_right()
    assert universe.is_four_vertex(t_10.get_vertex_left())
    assert universe.four_vertices_bag.contains(9)

    t_4, t_5 = universe.flip_edge(triangle_id=4)

    # Check the IDs of the vertices of the triangles
    assert t_4.get_vertex_left().ID == 6
    assert t_4.get_vertex_right().ID == 7
    assert t_4.get_vertex_center().ID == 2
    assert t_5.get_vertex_left().ID == 2
    assert t_5.get_vertex_right().ID == 3
    assert t_5.get_vertex_center().ID == 7
    assert universe.vertex_pool.get(index=2).get_triangle_left().ID == 2
    assert universe.vertex_pool.get(index=2).get_triangle_right().ID == 5
    assert universe.vertex_pool.get(index=3).get_triangle_left().ID == 5
    assert universe.vertex_pool.get(index=3).get_triangle_right().ID == 6

    # Check if bags are updated correctly
    assert universe.is_four_vertex(t_4.get_vertex_left())
    assert universe.four_vertices_bag.contains(t_4.get_vertex_left().ID)

    for id in [3, 9, 17]:
        assert not universe.triangle_flip_bag.contains(id)
