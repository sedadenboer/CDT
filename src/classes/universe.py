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

        # Create pools for vertices and triangles
        self.vertex_pool = Pool(capacity=Vertex.capacity)
        self.triangle_pool = Pool(capacity=Triangle.capacity)

        # Create a bag for all triangles
        self.triangle_bag = Bag(pool_capacity=Triangle.capacity)

        # Create bags for vertices with degree 4 and triangles that can be flipped
        self.four_vertices_bag = Bag(pool_capacity=Vertex.capacity)
        self.triangle_flip_bag = Bag(pool_capacity=Triangle.capacity)

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
            vertex_id = self.vertex_pool.occupy()
            time = i // width
            vertex = Vertex(vertex_id, time)
            initial_vertices.append(vertex)
        
        # Create initial triangles
        initial_triangles = []

        for i in range(total_time - 1):  # Only 2 time slices
            triangles = []
            for j in range(width - 1):  # Only 2 vertices per slice
                # Create two triangles for each square
                # Left triangle
                tl_id = self.triangle_pool.occupy()
                tl = Triangle(tl_id)
                tl.set_vertices(
                    vl=initial_vertices[i * width + j],
                    vr=initial_vertices[i * width + (j + 1)],
                    vc=initial_vertices[(i + 1) * width + j]
                )

                # Right triangle
                tr_id = self.triangle_pool.occupy()
                tr = Triangle(tr_id)
                tr.set_vertices(
                    vl=initial_vertices[(i + 1) * width + j],
                    vr=initial_vertices[(i + 1) * width + (j + 1)],
                    vc=initial_vertices[i * width + (j + 1)]
                )

                # Add the triangles to relevant bags
                self.triangle_bag.add(tl_id)
                self.triangle_bag.add(tr_id)
                self.triangle_flip_bag.add(tl_id)
                self.triangle_flip_bag.add(tr_id)

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
                    if j != 0:
                        t.set_triangle_left(initial_triangles[i][j - 1])
                else:
                    t.set_triangle_left(initial_triangles[i][j - 1])
                    # If j is not the last triangle in the row
                    if j != len(initial_triangles[i]) - 1:
                        t.set_triangle_right(initial_triangles[i][j + 1])
                
                # Set the triangle to the center
                if i != 0 and t.is_upwards():
                    t.set_triangle_center(initial_triangles[i - 1][j + 1])
                elif i != total_time - 2 and t.is_downwards():
                    t.set_triangle_center(initial_triangles[i + 1][j - 1])
        
        triangle_array = np.array(initial_triangles)
        
        # Add the left and right triangles to vertex objects
        for triangle in triangle_array.flatten():
            if triangle.is_upwards():
                triangle.vl.set_triangle_right(triangle)
                triangle.vr.set_triangle_left(triangle)
        
        # # Print initial vertices
        # for i in range(total_time):
        #     for j in range(width):
        #         print(initial_vertices[i * width + j].ID, end=" ")
        #     print()
        # print()
        
        # # Print the connectivity of the triangles
        # for triangle in triangle_array.flatten():
        #     print("triangle:", triangle.ID, triangle.type)
        #     print("l, r, c:", triangle.tl.ID if triangle.tl else None, triangle.tr.ID if triangle.tr else None, triangle.tc.ID if triangle.tc else None)
        #     print("vertices:", triangle.vl.ID if triangle.vl else None, triangle.vr.ID if triangle.vr else None, triangle.vc.ID if triangle.vc else None)
        #     print()
                
        # # Print the left and right triangles of the vertices
        # for vertex in initial_vertices:
        #     print("vertex:", vertex.ID)
        #     print("l, r:", vertex.tl.ID if vertex.tl else None, vertex.tr.ID if vertex.tr else None)
        #     print()
       
        
    def insert_vertex(self, triangle: "Triangle") -> None:
        """
        Insert a vertex into the triangulation.

        Args:
            triangle (Triangle): Triangle to insert vertex into.
        """
        # Get the center and right vertex of the triangle
        tc = triangle.get_triangle_center()
        vr = triangle.get_vertex_right()
        
        # generate a new vertex ID and create a new vertex
        new_vertex_id = self.vertex_pool.occupy()
        new_vertex = Vertex(new_vertex_id, triangle.time)

        # Add vertex to bag with all vertices and the four vertices bag
        self.four_vertices_bag.add(new_vertex_id)

        # Update the size of this triangulation layer
        self.slice_sizes[triangle.time] += 1

        # Generate two new triangle IDs
        triangle_id_1 = self.triangle_pool.occupy()
        triangle_id_2 = self.triangle_pool.occupy()

        # Add the ids to the bag with all triangles
        self.triangle_bag.add(triangle_id_1)
        self.triangle_bag.add(triangle_id_2)

        # Create triangles
        triangle1 = Triangle(triangle_id_1)
        triangle2 = Triangle(triangle_id_2)

        # Set vertices for new triangles
        triangle1.set_vertices(vl=new_vertex, vr=vr, vc=triangle.get_vertex_center())
        triangle2.set_vertices(vl=new_vertex, vr=vr, vc=tc.get_vertex_center())

        # Set triangles for new triangles
        triangle1.set_triangles(tl=triangle, tr=triangle.get_triangle_right(), tc=triangle2)
        triangle2.set_triangles(tl=tc, tr=tc.get_triangle_right(), tc=triangle1)

        # Update the original triangles
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
        if triangle.type != triangle1.get_triangle_right().type:
            self.triangle_flip_bag.remove(triangle.ID)
            self.triangle_flip_bag.add(triangle1.ID)
        
        if triangle2.type != triangle2.get_triangle_right().type:
            self.triangle_flip_bag.remove(tc.ID)
            self.triangle_flip_bag.add(triangle2.ID) 
    
    def remove_vertex(self, vertex: "Vertex") -> None:
        """
        Remove a vertex from the triangulation.

        Args:
            vertex (Vertex): Vertex to remove.
        """
        # Get the left and right triangles of the vertex
        tl = vertex.get_triangle_left()
        tr = vertex.get_triangle_right()

        # Get the center triangles of the left and right triangles of the vertex
        tlc = tl.get_triangle_center()
        trc = tr.get_triangle_center()

        # Get the right triangle of the right triangle of the vertex
        trn = tr.get_triangle_right()

        # Get the right triangle of the center right triangle of the vertex
        trcn = trc.get_triangle_right()

        # Update the right triangles of the left- and left center- triangles
        tl.set_triangle_right(trn)
        tlc.set_triangle_right(trcn)

        # Update the right vertices of the left- and left center- triangles
        vr = tr.get_vertex_right()
        tl.set_vertex_right(vr)
        tlc.set_vertex_right(vr)

        # Update the left triangle in the right vertex
        vr.set_triangle_left(tl)

        # Update the size of this triangulation layer
        self.slice_sizes[vertex.time] -= 1

        # Update the triangle pool and bag with deleted triangles
        self.triangle_pool.free(tr.ID)
        self.triangle_pool.free(trc.ID)
        self.triangle_bag.remove(tr.ID)
        self.triangle_bag.remove(trc.ID)

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

    def flip_edge(self):
        pass

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


# Example usage:
universe = Universe(total_time=3, initial_slice_size=3)

