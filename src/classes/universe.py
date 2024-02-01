# universe.py
#
# Author: Seda den Boer
# Date: 02-01-2024
# 
# Description:

from __future__ import annotations
from typing import cast
from vertex import Vertex
from triangle import Triangle
from pool import Pool
from bag import Bag


class Universe:
    """
    The Universe class represents the current state of the triangulation
    and stores properties of the geometry in a convenient matter. It also
    provides member functions that carry out changes on the geometry.
    """

    def __init__(self, total_time: int, initial_slice_size: int, VERTEX_CAPACITY: int = 100000):
        self.total_time = total_time
        self.initial_slice_size = initial_slice_size

        VERTEX_CAPACITY = VERTEX_CAPACITY
        TRIANGLE_CAPACITY = 2 * VERTEX_CAPACITY

        # Create pools for vertices and triangles
        self.vertex_pool = Pool(capacity=VERTEX_CAPACITY)
        self.triangle_pool = Pool(capacity=TRIANGLE_CAPACITY)

        # Create bags
        self.triangle_add_bag = Bag(pool_capacity=TRIANGLE_CAPACITY)
        self.four_vertices_bag = Bag(pool_capacity=VERTEX_CAPACITY)
        self.triangle_flip_bag = Bag(pool_capacity=TRIANGLE_CAPACITY)

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

        # # Print vertices in grid
        # for i in range(total_time):
        #     for j in range(width):
        #         print(initial_vertices[i * width + j].ID, end=" ")
        #     print()
        
        # # Print the triangles and their left, right, center vertices
        # for t in initial_triangles:
        #     print(f"triangle: {t.ID}")
        #     print(f"vertex left, right, center: {t.vl_.ID, t.vr_.ID, t.vc_.ID}")
        #     print(f"triangle left, right, center: {t.tl_.ID, t.tr_.ID, t.tc_.ID}")
        #     print()
        #     print()
        
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
            # print(f"ADD removed {triangle_id} from flip bag and added {triangle1.ID}")
        if triangle2.type != triangle2.get_triangle_right().type:
            self.triangle_flip_bag.remove(tc.ID)
            self.triangle_flip_bag.add(triangle2.ID)
            # print(f"ADD removed {tc.ID} from flip bag and added {triangle2.ID}\n")
    
        # Update count of up and down oriented triangles
        self.triangle_up_count += 1
        self.triangle_down_count += 1
            
        return new_vertex, triangle1, triangle2
    
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

        # Remove deleted triangles from the add bag
        self.triangle_add_bag.remove(tr.ID)
        self.triangle_add_bag.remove(trc.ID)

        # Update the bag with triangles that can be flipped
        if self.triangle_flip_bag.contains(tr.ID):
            self.triangle_flip_bag.remove(tr.ID)
            self.triangle_flip_bag.add(tl.ID)
            # print(f"DEL removed {tr.ID} from flip bag and added {tl.ID}")
        
        if self.triangle_flip_bag.contains(trc.ID):
            self.triangle_flip_bag.remove(trc.ID)
            self.triangle_flip_bag.add(tlc.ID)
            # print(f"DEL removed {trc.ID} from flip bag and added {tlc.ID}\n")

        # Remove deleted triangles from the triangle pool
        self.triangle_pool.free(tr.ID)
        self.triangle_pool.free(trc.ID)

        # Remove vertices that previously had degree 4 from the four vertices bag
        self.four_vertices_bag.remove(vertex.ID)

        # Free the deleted vertex ib the vertex pool
        self.vertex_pool.free(vertex.ID)

        # Update count of up and down oriented triangles
        self.triangle_up_count -= 1
        self.triangle_down_count -= 1

        return vertex, tr, trc

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

        return triangle, tr
    
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

    def get_triangulation_state(self):
        """
        Get the current state of the triangulation.
        """
        pass


# Example usage:
if __name__ == "__main__":
    universe = Universe(total_time=4, initial_slice_size=4)
    
    print("triangle pool:", universe.triangle_pool.used_indices)
    print("vertex pool:", universe.vertex_pool.used_indices)
    print("triangle_add_bag:", universe.triangle_add_bag)
    print("triangle_flip_bag:", universe.triangle_flip_bag)
    print("four vertices bag:", universe.four_vertices_bag)
    print()
    
    inserted_vertex, new_t1, new_t2 = universe.insert_vertex(triangle_id=10)

    print("triangle pool:", universe.triangle_pool.used_indices)
    print("vertex pool:", universe.vertex_pool.used_indices)
    print("triangle_add_bag:", universe.triangle_add_bag)
    print("triangle_flip_bag:", universe.triangle_flip_bag)
    print("four vertices bag:", universe.four_vertices_bag)
    print()
    
    removed_vertex, removed_t1, removed_t2 = universe.remove_vertex(vertex_id=16)
    
    print("triangle pool:", universe.triangle_pool.used_indices)
    print("vertex pool:", universe.vertex_pool.used_indices)
    print("triangle_add_bag:", universe.triangle_add_bag)
    print("triangle_flip_bag:", universe.triangle_flip_bag)
    print("four vertices bag:", universe.four_vertices_bag)
    print()
    # print(universe.vertex_pool.used_indices)
    # print(universe.triangle_pool.used_indices)
    # print()
    # print("triangle_add_bag:", universe.triangle_add_bag)
    # print("triangle_flip_bag:", universe.triangle_flip_bag)
    # print("four vertices bag:", universe.four_vertices_bag)
    # print()
    # print(universe.triangle_pool.used_indices)
    

