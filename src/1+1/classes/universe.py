# universe.py
#
# Author: Seda den Boer
# Date: 02-01-2024
# 
# Description: The Universe class represents the current state of the triangulation
# and stores properties of the geometry in a convenient matter. It also
# provides member functions that carry out changes on the geometry.


from __future__ import annotations
from typing import cast
from vertex import Vertex
from triangle import Triangle
from pool import Pool
from bag import Bag
import pickle
import sys
sys.setrecursionlimit(10**8)


class Universe:
    """
    The Universe class represents the current state of the triangulation
    and stores properties of the geometry in a convenient matter. It also
    provides member functions that carry out changes on the geometry.

    Args (Attributes):
        total_time (int): Total number of time slices.
        initial_slice_size (int): Initial size of the time slices.

    Attributes:
        vertex_pool (Pool): Pool of vertices.
        triangle_pool (Pool): Pool of triangles.
        triangle_add_bag (Bag): Bag of triangles that can be added.
        four_vertices_bag (Bag): Bag of vertices that have degree 4.
        triangle_flip_bag (Bag): Bag of triangles that can be flipped.
        n_vertices (int): Total number of vertices in the triangulation.
        n_triangles (int): Total number of triangles in the triangulation.
        slice_sizes (dict): Dictionary with the size of each time slice.
        triangle_up_count (int): Number of triangles with an upwards orientation.
        triangle_down_count (int): Number of triangles with a downwards orientation.
    """
    VERTEX_CAPACITY = 1000000
    TRIANGLE_CAPACITY = 2 * VERTEX_CAPACITY
    MINIMAL_TIME = 3
    MINIMAL_SLICE_SIZE = 3
    MINIMAL_VERTEX_CAPACITY = 9

    def __init__(self, total_time: int, initial_slice_size: int):
        if total_time < self.MINIMAL_TIME:
            raise ValueError("Total time must be greater than 3.")
        if initial_slice_size < self.MINIMAL_SLICE_SIZE:
            raise ValueError("Initial slice size must be greater than 3.")
        if self.VERTEX_CAPACITY < self.MINIMAL_VERTEX_CAPACITY:
            raise ValueError("Vertex capacity must be greater than 9.")
        
        self.total_time = total_time
        self.initial_slice_size = initial_slice_size

        # Create pools for vertices and triangles
        self.vertex_pool = Pool(capacity=self.VERTEX_CAPACITY)
        self.triangle_pool = Pool(capacity=self.TRIANGLE_CAPACITY)

        # Create bags
        self.triangle_add_bag = Bag(pool_capacity=self.TRIANGLE_CAPACITY)
        self.four_vertices_bag = Bag(pool_capacity=self.VERTEX_CAPACITY)
        self.triangle_flip_bag = Bag(pool_capacity=self.TRIANGLE_CAPACITY)

        # Total size of the triangulation   
        self.n_vertices = total_time * initial_slice_size
        self.n_triangles = (total_time - 1) * (initial_slice_size - 1) * 2

        # Create a dictionary to store the size of each time slice
        self.slice_sizes = {t: initial_slice_size for t in range(total_time)}

        # Number of triangles with a certain orientation
        self.triangle_up_count = 0
        self.triangle_down_count = 0

        # Initialise the triangulation, store the vertices and triangles
        self.initialise_triangulation()

    def initialise_triangulation(self):
        """
        Initialise the grid with vertices and triangles. Creates
        vertices and triangles and set their connectivity.
        """
        total_time = self.total_time
        width = self.initial_slice_size

        # Create initial vertices
        initial_vertices = []

        for i in range(total_time * width):
            time = i // width
            vertex = Vertex(time)
            vertex_id = self.vertex_pool.occupy(vertex)
            initial_vertices.append(vertex)

        # Create initial triangles, 2 in a square
        initial_triangles = []

        for i in range(total_time):  
            for j in range(width): 
                # Left triangle
                tl = Triangle()
                tl_id = self.triangle_pool.occupy(tl)
                tl.set_vertices(
                    vl=initial_vertices[i * width + j],
                    vr=initial_vertices[i * width + (j + 1) % width],
                    vc=initial_vertices[((i + 1) % total_time) * width + j]
                )

                # Right triangle
                tr = Triangle()
                tr_id = self.triangle_pool.occupy(tr)
                tr.set_vertices(
                    vl=initial_vertices[((i + 1) % total_time) * width + j],
                    vr=initial_vertices[((i + 1) % total_time) * width + (j + 1) % width],
                    vc=initial_vertices[i * width + (j + 1) % width]
                )

                # Add triangles to add and flip bag
                self.triangle_add_bag.add(tl_id)
                self.triangle_add_bag.add(tr_id)
                self.triangle_flip_bag.add(tl_id)
                self.triangle_flip_bag.add(tr_id)

                # Add the triangles to the list of triangles
                initial_triangles.append(tl)
                initial_triangles.append(tr)

                # Update count of up and down oriented triangles
                self.triangle_up_count += 1
                self.triangle_down_count += 1
        
        # Set triangle connectivity
        for i in range(total_time):
            for j in range(width):
                row = 2 * i * width
                column = 2 * j

                tl = initial_triangles[row + column]
                tr = initial_triangles[row + column + 1]

                # Connect left triangle
                tl.set_triangles(
                    tl=initial_triangles[row + (column - 1 + 2 * width) % (2 * width)],
                    tr=initial_triangles[row + column + 1],
                    tc=initial_triangles[(row + column - 2 * width + 1 + 2 * total_time * width) % (2 * total_time * width)]
                )

                # Connect right triangle
                tr.set_triangles(
                    tl=initial_triangles[row + column],
                    tr=initial_triangles[row + (column + 2) % (2 * width)],
                    tc=initial_triangles[(row + column + 2 * width) % (2 * total_time * width)]
                )
        
    def insert_vertex(self, triangle_id: int) -> tuple[Vertex, Triangle, Triangle]:
        """
        Insert a vertex into the triangulation.

        Args:
            triangle (Triangle): Triangle to insert vertex into.

        Returns:
            tuple[Vertex, Triangle, Triangle]: The inserted vertex and the two triangles
            that are created with it.
        """
        # Get the triangle from the triangle pool
        triangle: Triangle = cast(Triangle, self.triangle_pool.get(triangle_id))
        
        # Get the center and right vertex of the triangle
        tc: Triangle = triangle.get_triangle_center()
        vr: Vertex = triangle.get_vertex_right()

        time = triangle.time

        # Create a new vertex and add it to the vertex pool and the bag with 4-vertices
        new_vertex = Vertex(time)
        new_vertex_id = self.vertex_pool.occupy(new_vertex)
        self.four_vertices_bag.add(new_vertex_id)

        # Update the size of this triangulation layer
        self.slice_sizes[triangle.time] += 1

        # Update the original triangles' vertices
        triangle.set_vertex_right(new_vertex)
        tc.set_vertex_right(new_vertex)

        # Update the new vertex's past and future neighbours, and vice versa
        if triangle.is_upwards():
            new_vertex.add_future_neighbour(triangle.get_vertex_center())
            new_vertex.add_past_neighbour(tc.get_vertex_center())
        else:
            new_vertex.add_past_neighbour(triangle.get_vertex_center())
            new_vertex.add_future_neighbour(tc.get_vertex_center())

        # Create two new triangles
        triangle1 = Triangle()
        triangle2 = Triangle()

        # Generate two new triangle IDs and add them to the triangle pool
        triangle_id_1 = self.triangle_pool.occupy(triangle1)
        triangle_id_2 = self.triangle_pool.occupy(triangle2)

        # Add the ids to the bag with all triangles
        self.triangle_add_bag.add(triangle_id_1)
        self.triangle_add_bag.add(triangle_id_2)

        # Set vertices for new triangles
        triangle1.set_vertices(vl=new_vertex, vr=vr, vc=triangle.get_vertex_center())
        triangle2.set_vertices(vl=new_vertex, vr=vr, vc=tc.get_vertex_center())

        # Set triangles for new triangles
        triangle1.set_triangles(tl=triangle, tr=triangle.get_triangle_right(), tc=triangle2)
        triangle2.set_triangles(tl=tc, tr=tc.get_triangle_right(), tc=triangle1)

        # Update flip bag
        if triangle1.type != triangle1.get_triangle_right().type:
            self.triangle_flip_bag.remove(triangle_id)
            self.triangle_flip_bag.add(triangle1.ID)
        if triangle2.type != triangle2.get_triangle_right().type:
            self.triangle_flip_bag.remove(tc.ID)
            self.triangle_flip_bag.add(triangle2.ID)
            
        # Update count of up and down oriented triangles
        self.triangle_up_count += 1
        self.triangle_down_count += 1

    def remove_vertex(self, vertex_id: int) -> tuple[Vertex, Triangle, Triangle]:
        """
        Remove a vertex from the triangulation.

        Args:
            vertex (Vertex): Vertex to remove.

        Returns:
            tuple[Vertex, Triangle, Triangle]: The removed vertex and the two triangles
            that are removed with it.
        """
        # Get the vertex from the vertex pool
        vertex: Vertex = cast(Vertex, self.vertex_pool.get(vertex_id))

        # Get the left and right triangles of the vertex
        tl: Triangle = vertex.get_triangle_left()
        tr: Triangle = vertex.get_triangle_right()

        # Get the center triangles of the left and right triangles of the vertex
        tlc: Triangle = tl.get_triangle_center()
        trc: Triangle = tr.get_triangle_center()

        # Get the right triangles of tr and trc
        trn: Triangle = tr.get_triangle_right()
        trcn: Triangle = trc.get_triangle_right()

        # Update the right triangles of the left- and left center- triangles and vice versa
        tl.set_triangle_right(trn)
        tlc.set_triangle_right(trcn)

        # Update the right vertices of the left- and left center- triangles
        vr: Vertex = tr.get_vertex_right()
        tl.set_vertex_right(vr)
        tlc.set_vertex_right(vr)

        # Update the size of this triangulation layer
        self.slice_sizes[vertex.time] -= 1

        # Update the future and past neighbours of the center vertices
        if tl.is_upwards():
            vertex.delete_future_neighbour(tl.get_vertex_center())
            vertex.delete_past_neighbour(tlc.get_vertex_center())
        else:
            vertex.delete_past_neighbour(tl.get_vertex_center())
            vertex.delete_future_neighbour(tlc.get_vertex_center())

        # Remove deleted triangles from the add bag
        self.triangle_add_bag.remove(tr.ID)
        self.triangle_add_bag.remove(trc.ID)

        # Update the bag with triangles that can be flipped
        if self.triangle_flip_bag.contains(tr.ID):
            self.triangle_flip_bag.remove(tr.ID)
            self.triangle_flip_bag.add(tl.ID)
        
        if self.triangle_flip_bag.contains(trc.ID):
            self.triangle_flip_bag.remove(trc.ID)
            self.triangle_flip_bag.add(tlc.ID)

        # Clear references to avoid memory leaks
        vertex.clear_references()
        tr.clear_references()
        trc.clear_references()

        # Remove deleted triangles from the triangle pool
        self.triangle_pool.free(tr.ID)
        self.triangle_pool.free(trc.ID)

        # Remove vertices that previously had degree 4 from the four vertices bag
        self.four_vertices_bag.remove(vertex.ID)

        # Free the deleted vertex in the vertex pool
        self.vertex_pool.free(vertex.ID)

        # Update count of up and down oriented triangles
        self.triangle_up_count -= 1
        self.triangle_down_count -= 1

    def flip_edge(self, triangle_id : int) -> tuple[Triangle, Triangle]:
        """
        Flip an edge in the triangulation.

        Args:
            triangle (Triangle): Triangle to flip edge in.

        Returns:
            tuple[Triangle, Triangle]: The two triangles that are flipped (t, tr).
        """
        # Get the triangle from the triangle pool
        triangle: Triangle = cast(Triangle, self.triangle_pool.get(triangle_id))

        # Get relevant neighbours
        tr: Triangle = triangle.get_triangle_right()
        tc: Triangle = triangle.get_triangle_center()
        trc: Triangle = tr.get_triangle_center()

        # Swap triangle orientation
        if triangle.is_upwards():
            triangle.get_vertex_left().set_triangle_right(tr)
            triangle.get_vertex_right().set_triangle_left(tr)
        else:
            tr.get_vertex_left().set_triangle_right(triangle)
            tr.get_vertex_right().set_triangle_left(triangle)
            
        # Update center triangles
        triangle.set_triangle_center(trc)
        tr.set_triangle_center(tc)

        # Get relevant vertices
        vl: Vertex = triangle.get_vertex_left()
        vr: Vertex = triangle.get_vertex_right()
        vc: Vertex = triangle.get_vertex_center()
        vrr: Vertex = tr.get_vertex_right()

        # Update vertices of triangles
        triangle.set_vertices(vl=vc, vr=vrr, vc=vl)
        tr.set_vertices(vl=vl, vr=vr, vc=vrr)

        # Update the future and past neighbours of the vertices
        if triangle.is_upwards():
            triangle.get_vertex_left().delete_future_neighbour(tr.get_vertex_right())
        else:
            triangle.get_vertex_left().delete_past_neighbour(tr.get_vertex_right())

        # Remove vertices that previously had degree 4 from the four vertices bag and add new ones
        if self.four_vertices_bag.contains(vl.ID):
            self.four_vertices_bag.remove(vl.ID)
        if self.is_four_vertex(vr):
            self.four_vertices_bag.add(vr.ID)
        if self.is_four_vertex(vc):
            self.four_vertices_bag.add(vc.ID)
        if self.four_vertices_bag.contains(vrr.ID):
            self.four_vertices_bag.remove(vrr.ID)

        # Remove triangles that have the same type as their right neighbour from the flip bag
        if self.triangle_flip_bag.contains(triangle.get_triangle_left().ID) and (
            triangle.type == triangle.get_triangle_left().type):
            self.triangle_flip_bag.remove(triangle.get_triangle_left().ID)
        if self.triangle_flip_bag.contains(tr.ID) and tr.type == tr.get_triangle_right().type:
            self.triangle_flip_bag.remove(tr.ID)

        # Add triangles that have a opposite type than their right neighbour to the flip bag
        if not self.triangle_flip_bag.contains(triangle.get_triangle_left().ID) and (
            triangle.type != triangle.get_triangle_left().type):
            self.triangle_flip_bag.add(triangle.get_triangle_left().ID)
        if not self.triangle_flip_bag.contains(tr.ID) and tr.type != tr.get_triangle_right().type:
            self.triangle_flip_bag.add(tr.ID)

    def is_four_vertex(self, vertex: Vertex) -> bool:
        """
        Checks if a vertex is of degree 4.

        Args:
            vertex (Vertex): Vertex to be checked.

        Returns:
            bool: True if vertex is of degree 4, otherwise False.
        """
        return (vertex.get_triangle_left().get_triangle_right() 
                == vertex.get_triangle_right()) and (
                vertex.get_triangle_left().get_triangle_center().get_triangle_right() 
                == vertex.get_triangle_right().get_triangle_center())

    def get_total_size(self) -> int:
        """
        Get the total size of the triangulation.

        Returns:
            int: Total size of the triangulation.
        """
        return self.vertex_pool.get_number_occupied()
    
    def sort_vertices_periodic(self, vertices):
        # Check if the list is empty or has only one element
        if len(vertices) <= 1:
            return vertices

        # Start with any vertex
        start_vertex = vertices[0]
        current_vertex = start_vertex
        sorted_vertices = [start_vertex]

        # Traverse the neighbors until reaching the starting vertex
        while current_vertex.get_neighbour_right() != start_vertex:
            next_vertex = current_vertex.get_neighbour_right()
            sorted_vertices.append(next_vertex)
            current_vertex = next_vertex

        return sorted_vertices

    def get_triangulation_state(self):
        """
        Get the current state of the triangulation.
        """
        vertex_sheet = {}

        for id in self.vertex_pool.used_indices:
            vertex = self.vertex_pool.get(id)

            # Add vertex to dictionary in key with time as index
            if vertex.time in vertex_sheet:
                vertex_sheet[vertex.time].append(vertex)
            else:
                vertex_sheet[vertex.time] = [vertex]
        
        # Sort the vertices in each time slice in the dictionary based on their right neighbour
        # loop through the dictionary and sort the vertices in each time slice
        for key in vertex_sheet:
            vertex_sheet[key] = self.sort_vertices_periodic(vertex_sheet[key])

        return vertex_sheet

    def print_state(self):
        """
        Print the current state of the triangulation.
        """
        state = self.get_triangulation_state()
        vertex_sheet = list(state.values())

        # Add the edges from the space, future, and past neighbors in the vertex objects
        for y, row in enumerate(vertex_sheet):
            for x, col in enumerate(row):
                space_neighbours = [neighbour.ID for neighbour in col.get_space_neighbours()]
                future_neighbours = [neighbour.ID for neighbour in col.get_future_neighbours()]
                past_neighbours = [neighbour.ID for neighbour in col.get_past_neighbours()]
                print(
                    f"VERTEX {col.ID}: "
                    f"Space neighbours: {space_neighbours}, "
                    f"Future neighbours: {future_neighbours}, "
                    f"Past neighbours: {past_neighbours}"
                )
            print()

    def check_validity(self):
        """
        Check the validity of the triangulation.

        Checks among others:
        - If each vertex has a left and right neighbour.
        - If the time and connectivity of future neighbours is valid.
        - If the time and connectivity of past neighbours is valid.
        - If each triangle has three neighbouring triangles.
        - If each triangle has  base vertices from the same time slice.
        - If each triangle has the correct orientation and time.
        """
        state = self.get_triangulation_state()
        vertex_sheet = list(state.values())

        all_triangles_indices = self.triangle_pool.used_indices
        all_triangles = [self.triangle_pool.elements[i] for i in all_triangles_indices]

        # Check for validity of connections
        for y, row in enumerate(vertex_sheet):
            assert len(row) >= 3, f"Row {y} has {len(row)} vertices"
            for x, col in enumerate(row):
                # Assert that each vertex has a left and right neighbour
                space_neighbours = col.get_space_neighbours()
                assert len(space_neighbours) == 2, f"Vertex {col.ID} has {len(space_neighbours)} space neighbours"

                # Check validity of time and connectivity of future neighbours
                for future_neighbour in col.get_future_neighbours():
                    assert col in future_neighbour.get_past_neighbours(), (
                        f"Future neighbour {future_neighbour.ID} does not have {col.ID} as past neighbour"
                    )
                    if col.time != self.total_time - 1:
                        assert future_neighbour.time == col.time + 1, (
                            f"Future neighbour {future_neighbour.ID} time {future_neighbour.time} is not {col.time + 1}"
                        )
                    else:
                        assert future_neighbour.time == 0, (
                            f"Future neighbour {future_neighbour.ID} time {future_neighbour.time} is not 0"
                        )

                # Check validity of time and connectivity of past neighbours
                for past_neighbour in col.get_past_neighbours():
                    assert col in past_neighbour.get_future_neighbours(), (
                        f"Past neighbour {past_neighbour.ID} does not have {col.ID} as future neighbour"
                    )
                    if col.time != 0:
                        assert past_neighbour.time == col.time - 1, (
                            f"Past neighbour {past_neighbour.ID} time {past_neighbour.time} is not {col.time - 1}"
                        )
                    else:
                        assert past_neighbour.time == self.total_time - 1, (
                            f"Past neighbour {past_neighbour.ID} time {past_neighbour.time} is not {self.total_time - 1}"
                        )

        # Check for validity of triangles
        for triangle in all_triangles:
            assert len(triangle.get_triangles()) == 3, (
                f"Triangle {triangle.ID} has {len(triangle.get_triangles())} neighbouring triangles"
            )
            assert triangle.get_vertex_left().time == triangle.get_vertex_right().time, (
                f"Triangle {triangle.ID} has vertices from different time slices"
            )

            # Check for validity of triangle orientation and time
            if triangle.is_upwards():
                if triangle.time != self.total_time - 1:
                    assert triangle.get_vertex_center().time == triangle.get_vertex_left().time + 1, (
                        f"Triangle {triangle.ID} is upwards but has left vertex with higher time than center vertex."
                    )
                else:
                    assert triangle.get_vertex_center().time == 0, (
                        f"Triangle {triangle.ID} is upwards but has left vertex with higher time than center vertex."
                    )
            else:
                if triangle.time != 0:
                    assert triangle.get_vertex_center().time == triangle.get_vertex_left().time - 1, (
                        f"Triangle {triangle.ID} is downwards but has left vertex with lower time than center vertex."
                    )
                else:
                    assert triangle.get_vertex_center().time == self.total_time - 1, (
                        f"Triangle {triangle.ID} is downwards but has left vertex with lower time than center vertex."
                    )

    def save_to_file(self, filename):
        """
        Save the state of the Universe to a file using pickle.

        Args:
            filename (str): The name of the file to save the state to.
        """
        with open(filename + '.pkl', 'wb') as file:
            pickle.dump(self.__dict__, file, protocol=pickle.HIGHEST_PROTOCOL)

    def load_from_file(self, filename):
        """
        Load the state of the Universe from a file using pickle.

        Args:
            filename (str): The name of the file to load the state from.
        """
        with open(filename, 'rb') as file:
            state = pickle.load(file)
        self.__dict__.update(state)