# universe.py
#
# Author: Seda den Boer
# Date: 02-01-2024
# 
# Description:

from vertex import Vertex
from triangle import Triangle
from pool import Pool
from bag import Bag
import numpy as np


class Universe:
    """
    The Universe class represents the current state of the triangulation
    and stores properties of the geometry in a convenient matter. It also
    provides member functions that carry out changes on the geometry.
    """

    def __init__(self, total_time: int, initial_slice_size: int) -> None:
        self.total_time = total_time
        self.initial_slice_size = initial_slice_size

        self.VERTEX_CAPACITY = 50
        self.TRIANGLE_CAPACITY = 2 * self.VERTEX_CAPACITY

        # Create pools for vertices and triangles
        self.vertex_pool = Pool(capacity=self.VERTEX_CAPACITY)
        self.triangle_pool = Pool(capacity=self.TRIANGLE_CAPACITY)

        print("Done creating pools.")

        # Create bags
        self.triangle_add_bag = Bag(pool_capacity=self.TRIANGLE_CAPACITY)
        self.four_vertices_bag = Bag(pool_capacity=self.VERTEX_CAPACITY)
        self.triangle_flip_bag = Bag(pool_capacity=self.TRIANGLE_CAPACITY)

        print("Done creating bags.")

        # Total size of the triangulation   
        self.n_vertices = total_time * initial_slice_size
        self.n_triangles = (total_time - 1) * (initial_slice_size - 1) * 2

        # Create a dictionary to store the size of each time slice
        self.slice_sizes = {t: initial_slice_size for t in range(total_time)}

        # Initialise the triangulation, store the vertices and triangles
        self.triangulation = self.initialise_triangulation()

    def initialise_triangulation(self) -> None:
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
        
        # Create initial triangles
        initial_triangles = []

        for i in range(total_time - 1):  # Only 2 time slices
            triangles = []
            for j in range(width - 1):  # Only 2 vertices per slice
                # Create two triangles for each square
                # Left triangle
                tl = Triangle()
                tl_id = self.triangle_pool.occupy(tl)
                tl.set_vertices(
                    vl=initial_vertices[i * width + j],
                    vr=initial_vertices[i * width + (j + 1)],
                    vc=initial_vertices[(i + 1) * width + j]
                )

                # Right triangle
                tr = Triangle()
                tr_id = self.triangle_pool.occupy(tr)
                tr.set_vertices(
                    vl=initial_vertices[(i + 1) * width + j],
                    vr=initial_vertices[(i + 1) * width + (j + 1)],
                    vc=initial_vertices[i * width + (j + 1)]
                )

                # # Add triangles to flip bag
                # self.triangle_flip_bag.add(tl_id)
                # self.triangle_flip_bag.add(tr_id)

                # Add the triangles to the list of triangles
                triangles.append(tl)
                triangles.append(tr)

            # Add the list of triangles to the list of time slices
            initial_triangles.append(triangles)

        # Create triangle connectivity
        for i in range(len(initial_triangles)):
            for j in range(len(initial_triangles[i])):
                t = initial_triangles[i][j]
                
                # Set the triangle to the left and right
                if t.is_upwards():
                    t.set_triangle_right(initial_triangles[i][j + 1])
                    # If j is not the first triangle in the row
                    if j != 0:
                        t.set_triangle_left(initial_triangles[i][j - 1])
                    
                    # Add the triangle to the flip bag
                    self.triangle_flip_bag.add(t.ID)
                else:
                    t.set_triangle_left(initial_triangles[i][j - 1])
                    # If j is not the last triangle in the row
                    if j != len(initial_triangles[i]) - 1:
                        t.set_triangle_right(initial_triangles[i][j + 1])

                        # Add the triangle to the flip bag 
                        # (triangles with no neighbour to the right should be skipped)
                        self.triangle_flip_bag.add(t.ID)
                    
                # Set the triangle to the center
                if i != 0 and t.is_upwards():
                    t.set_triangle_center(initial_triangles[i - 1][j + 1])
                    self.triangle_add_bag.add(t.ID)
                elif i != total_time - 2 and t.is_downwards():
                    t.set_triangle_center(initial_triangles[i + 1][j - 1])
                    self.triangle_add_bag.add(t.ID)
        
        triangle_array = np.array(initial_triangles)
        
        # Add the left and right triangles to vertex objects
        for triangle in triangle_array.flatten():
            if triangle.is_upwards():
                triangle.vl_.set_triangle_right(triangle)
                triangle.vr_.set_triangle_left(triangle)
        
        # Print vertices in grid
        for i in range(total_time):
            for j in range(width):
                print(initial_vertices[i * width + j].ID, end=" ")
            print()
        
       
    def insert_vertex(self, triangle_id: int) -> tuple["Vertex", "Triangle", "Triangle"]:
        """
        Insert a vertex into the triangulation.

        Args:
            triangle (Triangle): Triangle to insert vertex into.
        """
        # Get the triangle from the triangle pool
        triangle: Triangle = self.triangle_pool.get(triangle_id)

        # Get the center and right vertex of the triangle
        tc: Triangle = triangle.get_triangle_center()
        vr: Vertex = triangle.get_vertex_right()

        # Create a new vertex and add it to the vertex pool
        new_vertex = Vertex(triangle.time)
        new_vertex_id = self.vertex_pool.occupy(new_vertex)

        # Add vertex to bag with all vertices and the four vertices bag
        self.four_vertices_bag.add(new_vertex_id)

        # Update the size of this triangulation layer
        self.slice_sizes[triangle.time] += 1

        # Create triangles
        triangle1 = Triangle()
        triangle2 = Triangle()

        # Generate two new triangle IDs
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

        # Set the triangles for original neighbouring triangles
        triangle1.get_triangle_right().set_triangle_left(triangle1)
        triangle2.get_triangle_right().set_triangle_left(triangle2)

        # Update the original triangles' vertices
        triangle.set_vertex_right(new_vertex)
        triangle.set_triangle_right(triangle1)
        tc.set_vertex_right(new_vertex)
        tc.set_triangle_right(triangle2)

        # Update the vertex triangles
        if triangle.is_upwards():
            new_vertex.set_triangle_left(triangle)
            new_vertex.set_triangle_right(triangle1)
        else:
            new_vertex.set_triangle_left(tc)
            new_vertex.set_triangle_right(triangle2)

        # Update flip bag
        if triangle1.get_triangle_right() and triangle1.type != triangle1.get_triangle_right().type:
            self.triangle_flip_bag.add(triangle1.ID)
        
        if triangle2.get_triangle_right() and triangle2.type != triangle2.get_triangle_right().type:
            self.triangle_flip_bag.add(triangle2.ID)

        if not triangle.get_triangle_left() and triangle.type == triangle.get_triangle_right().type:
            self.triangle_flip_bag.remove(triangle.ID)

        if not tc.get_triangle_left() and tc.type == tc.get_triangle_right().type:
            self.triangle_flip_bag.remove(tc.ID)
        
        return new_vertex, triangle1, triangle2
    
    def remove_vertex(self, vertex_id: int) -> tuple["Vertex", "Triangle", "Triangle"]:
        """
        Remove a vertex from the triangulation.

        Args:
            vertex (Vertex): Vertex to remove.
        """
        # Get the vertex from the vertex pool
        vertex: Vertex = self.vertex_pool.get(vertex_id)

        # Get the left and right triangles of the vertex
        tl: Vertex = vertex.get_triangle_left()
        tr: Vertex = vertex.get_triangle_right()

        # Get the center triangles of the left and right triangles of the vertex
        tlc: Triangle = tl.get_triangle_center()
        trc: Triangle = tr.get_triangle_center()

        # Get the right triangle of the right triangle of the vertex
        trn: Triangle = tr.get_triangle_right()

        # Get the right triangle of the center right triangle of the vertex
        trcn: Triangle = trc.get_triangle_right()

        # Update the right triangles of the left- and left center- triangles
        # and vice versa
        tl.set_triangle_right(trn)
        tlc.set_triangle_right(trcn)
        trn.set_triangle_left(tl)
        trcn.set_triangle_left(tlc)

        # Update the right vertices of the left- and left center- triangles
        vr: Vertex = tr.get_vertex_right()
        tl.set_vertex_right(vr)
        tlc.set_vertex_right(vr)

        # Update the left triangle in the right vertex
        vr.set_triangle_left(tl)

        # Update the size of this triangulation layer
        self.slice_sizes[vertex.time] -= 1

        # Update the triangle pool and bag with deleted triangles
        self.triangle_pool.free(tr.ID)
        self.triangle_pool.free(trc.ID)
        self.triangle_add_bag.remove(tr.ID)
        self.triangle_add_bag.remove(trc.ID)

        # Remove the vertex from the vertex pool and bag with vertices of degree 4
        self.vertex_pool.free(vertex.ID)
        self.four_vertices_bag.remove(vertex.ID)

        # Update the bag with triangles that can be flipped
        if self.triangle_flip_bag.contains(tr.ID):
            self.triangle_flip_bag.remove(tr.ID)
            self.triangle_flip_bag.add(tl.ID)
        
        if self.triangle_flip_bag.contains(trc.ID):
            self.triangle_flip_bag.remove(trc.ID)
            self.triangle_flip_bag.add(tlc.ID)

        return vertex, tr, trc

    def flip_edge(self, triangle_id : int) -> tuple["Triangle", "Triangle"]:
        """
        Flip an edge in the triangulation.

        Args:
            triangle (Triangle): Triangle to flip edge in.
        """
        # Get the triangle from the triangle pool
        triangle: Triangle = self.triangle_pool.get(triangle_id)

        # Get neighbouring triangles and vertices
        tr: Triangle = triangle.get_triangle_right()
        tl: Triangle = triangle.get_triangle_left()
        tc: Triangle = triangle.get_triangle_center()

        vl: Vertex = triangle.get_vertex_left()
        vr: Vertex = triangle.get_vertex_right()
        vc: Vertex = triangle.get_vertex_center()

        if (triangle.is_upwards() and vr == tr.get_vertex_center()) or (
            triangle.is_downwards() and vc == tr.get_vertex_left()):
            # Get right center triangle and right vertex of right triangle
            trc: Triangle = tr.get_triangle_center()
            vrr: Vertex = tr.get_vertex_right()

            # Update left and right triangles of vertices of the original triangle
            if triangle.is_upwards():
                vl.set_triangle_right(tr)
                vr.set_triangle_left(tr)
            else:
                vrr.set_triangle_left(triangle)
                vc.set_triangle_right(triangle)
                
            # Update center triangles
            triangle.set_triangle_center(trc)
            trc.set_triangle_center(triangle)
            tr.set_triangle_center(tc)
            tc.set_triangle_center(tr)

            # Update vertices of triangles
            triangle.set_vertices(vl=vc, vr=vrr, vc=vl)
            tr.set_vertices(vl=vl, vr=vr, vc=vrr)

            # Remove vertices that previously had degree 4 from the four vertices bag
            if self.four_vertices_bag.contains(vl.ID):
                self.four_vertices_bag.remove(vl.ID)
            if self.four_vertices_bag.contains(vrr.ID):
                self.four_vertices_bag.remove(vrr.ID)

            # Add vertices that now have degree 4 to the four vertices bag
            if self.is_four_vertex(vr):
                self.four_vertices_bag.add(vr.ID)
            if self.is_four_vertex(vc.ID):
                self.four_vertices_bag.add(vc.ID)

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

    def is_four_vertex(self, vertex: "Vertex") -> bool:
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
    print(universe.triangle_flip_bag.elements)
    print(universe.triangle_pool.used_indices)
    
    # print("triangle_add_bag:", universe.triangle_add_bag)
    # print("triangle_flip_bag:", universe.triangle_flip_bag)
    # print("four vertices bag:", universe.four_vertices_bag)
    # print()
    
    inserted_vertex, new_t1, new_t2 = universe.insert_vertex(8)
    inserted_vertex, new_t1, new_t2 = universe.insert_vertex(14)
    removed_vertex, removed_t1, removed_t2 = universe.remove_vertex(17)

    # print(universe.vertex_pool.used_indices)
    # print(universe.triangle_pool.used_indices)
    # print()
    # print("triangle_add_bag:", universe.triangle_add_bag)
    # print("triangle_flip_bag:", universe.triangle_flip_bag)
    # print("four vertices bag:", universe.four_vertices_bag)
    # print()
    print(universe.triangle_pool.used_indices)
    

