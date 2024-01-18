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

        self.vertex_pool = Pool(capacity=Vertex.capacity)
        self.triangle_pool = Pool(capacity=Triangle.capacity)

        self.vertex_bag = Bag(pool_capacity=Vertex.capacity)
        self.triangle_bag = Bag(pool_capacity=Triangle.capacity)
        self.four_vertices_bag = Bag(pool_capacity=Vertex.capacity)

        self.size = total_time * initial_slice_size
        self.slice_sizes = {t: initial_slice_size for t in range(total_time)}

        self.triangulation = self.initialise_triangulation()

    def initialise_triangulation(self) -> list["Vertex"]:
        """
        Initialise the grid with vertices and triangles.

        Returns:
            List["Vertex"]: List of time slices with vertices.
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
            self.vertex_bag.add(vertex_id)
        
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
                    initial_vertices[i * width + j],
                    initial_vertices[i * width + (j + 1)],
                    initial_vertices[(i + 1) * width + j]
                )

                # Right triangle
                tr_id = self.triangle_pool.occupy()
                tr = Triangle(tr_id)
                tr.set_vertices(
                    initial_vertices[(i + 1) * width + j],
                    initial_vertices[(i + 1) * width + (j + 1)],
                    initial_vertices[i * width + (j + 1)]
                )

                # Add the triangles to the bag with all triangles
                self.triangle_bag.add(tl_id)
                self.triangle_bag.add(tr_id)

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
                    t.set_triangle_left(initial_triangles[i][j + 1])
                    if j != 0:
                        t.set_triangle_right(initial_triangles[i][j - 1])
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

        # Add the triangles to vertex objects
        for triangle in triangle_array.flatten():
            if triangle.vl:
                triangle.vl.add_triangle(triangle)
            if triangle.vr:
                triangle.vr.add_triangle(triangle)
            if triangle.vc:
                triangle.vc.add_triangle(triangle)

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
                
        # # Print the triangles of the vertices
        # for vertex in initial_vertices:
        #     print("vertex:", vertex.ID)
        #     for triangle in vertex.triangles:
        #         print(triangle.ID, end=" ")
        #     print()
        

    def insert_vertex(self, triangle: "Triangle") -> None:
        """
        Insert a vertex into the triangulation.

        Args:
            vertex (Vertex): Vertex to insertriangle.
        """
        tc = triangle.get_triangle_center()
        vr = triangle.get_vertex_right()
        
        # Create vertex
        new_vertex_id = self.vertex_pool.occupy()
        vertex = Vertex(new_vertex_id, triangle.time)
        self.slice_sizes[triangle.time] += 1

        # Add vertex to bag with all vertices and the four vertices bag
        self.vertex_bag.add(new_vertex_id)
        self.four_vertices_bag.add(new_vertex_id)

        # Generate two new triangle IDs
        triangle_id_1 = self.triangle_pool.occupy()
        triangle_id_2 = self.triangle_pool.occupy()

        # Add the ids to the bag with all triangles
        self.triangle_bag.add(triangle_id_1)
        self.triangle_bag.add(triangle_id_2)

        # Create triangles
        triangle1 = Triangle(triangle_id_1)
        triangle2 = Triangle(triangle_id_2)

        # Set vertices for original triangles
        triangle.set_vertex_right(vertex)


    def get_random_vertex(self, move_type: str) -> Vertex:
        pass
        
    def get_total_volume(self) -> int:
        """
        Returns the total volume of the triangulation.

        Returns:
            int: Total volume.
        """
        return self.size
    
    def update_links(self, old_vertex_id: int):
        """
        Update the links of the vertices and triangles after adding a vertex.
        """
        pass

# Example usage:
universe = Universe(total_time=3, initial_slice_size=3)

