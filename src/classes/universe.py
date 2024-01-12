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


class Universe(object):
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
        self.triangulation = self.initialize_triangulation()

    def initialize_triangulation(self) -> list["Vertex"]:
        """
        Initialize the grid with vertices and triangles.

        Returns:
            List[[Vertex]]: List of time slices with vertices.
        """
        time = self.total_time
        width = self.initial_slice_size

        # ---- Initialize vertices and their connections ----
        # Create empty list to store time slices with spatial nodes
        vertex_grid = []

        # Loop over all time slices
        for t in range(time):
            # Create empty list to store vertices
            vertices = []

            # Loop over all spatial points
            for _ in range(width):
                # Occupy a space in the vertex pool and get the id
                id = self.vertex_pool.occupy()

                # Add the id to the bag with all candidates
                self.vertex_bag.add(id)

                # Create vertex
                vertex = Vertex(id, t)

                # Add vertex to list of vertices of current time slice
                vertices.append(vertex)

            # Add list of vertices to list of time slices
            vertex_grid.append(vertices)

        # Create spatial edges
        for t in range(time):
            for i in range(width - 1):
                # Get vertex
                vertex = vertex_grid[t][i]
                vertex.right_neighbour = vertex_grid[t][i + 1]
            
            vertex_grid[t][-1].right_neighbour = vertex_grid[t][0]

        # Create time-like edges
        for t in range(time - 1):
            for i in range(width - 1):
                # Get vertices
                current_vertex = vertex_grid[t][i]
                future_vertex1 = vertex_grid[t + 1][i]
                future_vertex2 = vertex_grid[t + 1][i + 1]

                # Add future neighbour
                current_vertex.future_vertices.append(future_vertex1)
                current_vertex.future_vertices.append(future_vertex2)

                # Add past neighbour
                future_vertex1.past_vertices.append(current_vertex)
                future_vertex2.past_vertices.append(current_vertex)

            # Handle the last vertex differently for periodic boundary conditions
            current_vertex = vertex_grid[t][width - 1]
            future_vertex1 = vertex_grid[t + 1][width - 1]
            future_vertex2 = vertex_grid[t + 1][0]

            # Add future neighbour
            current_vertex.future_vertices.append(future_vertex1)
            current_vertex.future_vertices.append(future_vertex2)

            # Add past neighbour
            future_vertex1.past_vertices.append(current_vertex)
            future_vertex2.past_vertices.append(current_vertex)

        # # print vertex_grid by id
        # for time_slice in vertex_grid:
        #     print([vertex.ID for vertex in time_slice])

        # Print all past and future neighbours
        print("id, right_neighbour, future_vertices, past_vertices")
        for t in vertex_grid:
            for v in t:
                print(v.ID, v.right_neighbour.ID, [vertex.ID for vertex in v.future_vertices], [vertex.ID for vertex in v.past_vertices])

        # print(self.vertex_pool.elements)
        # print(self.all_candidates.elements)
        # print(self.all_candidates.indices)
        # self.all_candidates.log()

        # ---- Initialize triangles ----
        
        # Create empty list to store triangles
        triangle_layers = []

        # Loop over all time slices
        for t in range(time - 1):
            triangles = []
            # Loop over all spatial points
            for i in range(width - 1):
                # Occupy two spaces in the triangle pool and get the ids
                id1 = self.triangle_pool.occupy()
                id2 = self.triangle_pool.occupy()

                # Add the ids to the bag with all candidates
                self.triangle_bag.add(id1)
                self.triangle_bag.add(id2)

                # Create triangles
                triangle1 = Triangle(id1)
                triangle2 = Triangle(id2)

                # Set vertices for triangles
                triangle1.set_vertices(
                    vl=vertex_grid[t + 1][i],
                    vr=vertex_grid[t + 1][i + 1],
                    vc=vertex_grid[t][i]
                    )
                
                triangle2.set_vertices(
                    vl=vertex_grid[t][i],
                    vr=vertex_grid[t][i + 1],
                    vc=vertex_grid[t + 1][i + 1]
                    )

                # Add the triangles to the list of triangles of the vertices
                triangle1.vl.add_triangle(triangle1)
                triangle1.vr.add_triangle(triangle1)
                triangle1.vc.add_triangle(triangle1)

                triangle2.vl.add_triangle(triangle2)
                triangle2.vr.add_triangle(triangle2)
                triangle2.vc.add_triangle(triangle2)

                # Add triangles to list of triangles
                triangles.append(triangle1)
                triangles.append(triangle2)

            # Handle the last vertex differently for periodic boundary conditions
            # Occupy two spaces in the triangle pool and get the ids
            id1 = self.triangle_pool.occupy()
            id2 = self.triangle_pool.occupy()

            # Add the ids to the bag with all candidates
            self.triangle_bag.add(id1)
            self.triangle_bag.add(id2)

            # Create triangles
            triangle1 = Triangle(id1)
            triangle2 = Triangle(id2)

            # Set vertices for triangles
            triangle1.set_vertices(
                vl=vertex_grid[t + 1][width - 1],
                vr=vertex_grid[t + 1][0],
                vc=vertex_grid[t][width - 1]
                )
            
            triangle2.set_vertices(
                vl=vertex_grid[t][width - 1],
                vr=vertex_grid[t][0],
                vc=vertex_grid[t + 1][0]
                )

            # Add triangles to list of triangles
            triangles.append(triangle1)
            triangles.append(triangle2)

            # Add list of triangles to list of time slices
            triangle_layers.append(triangles)

        # Set the right, left and center triangles for each triangle
        for i, layer in enumerate(triangle_layers):
            for j, triangle in enumerate(layer):

                # if not on the right edge
                if j != len(layer) - 1:
                    triangle.set_triangle_left(layer[j - 1])
                    triangle.set_triangle_right(layer[j + 1])
                else:
                    triangle.set_triangle_left(layer[j - 1])
                    triangle.set_triangle_right(layer[0])

                if i == 0:
                    if triangle.is_downwards():
                        triangle.set_triangle_center(triangle_layers[i + 1][j + 1])
                elif i == len(triangle_layers) - 1:
                    if triangle.is_upwards():
                        triangle.set_triangle_center(triangle_layers[i - 1][j - 1])
                else:
                    if triangle.is_downwards():
                        triangle.set_triangle_center(triangle_layers[i + 1][j + 1])
                    else:
                        triangle.set_triangle_center(triangle_layers[i - 1][j - 1])

        # # print triangles list with Id
        # print("\nTriangles:")
        # for layer in triangle_layers:
        #     print([triangle.ID for triangle in layer])
        
        # for layer in triangle_layers:
        #     print([triangle.type for triangle in layer])

        # # Print all triangles and their neighbour triangles
        # # if none then print none
        # print("\nid, tl, tr, tc")
        # for t in triangle_layers:
        #     for triangle in t:
        #         print(triangle.ID, triangle.tl.ID if triangle.tl else None, triangle.tr.ID if triangle.tr else None, triangle.tc.ID if triangle.tc else None)

        # # Print all triangles
        # print("\nid, vl, vr, vc")
        # for t in triangle_layers:
        #     for triangle in t:
        #         print(triangle.ID, triangle.vl.ID, triangle.vr.ID, triangle.vc.ID)

        # print()
        # print(self.triangle_pool.elements)
        # print(self.triangle_bag.elements)
        # print(self.triangle_bag.indices)

    def delete_check(self, vertex):
        """
        Check if a vertex can be deleted.
        """
        pass
    
    def get_random_vertex(self, move_type: str) -> Vertex:
        
        if move_type == "add":
            # Pick a random vertex from the bag of all candidates
            return self.all_candidates.pick() # Returns its pool index
        elif move_type == "delete":
            # Pick a random vertex from the bag of deletable candidates
            return self.deletable_candidates.pick()
        else:
            raise ValueError("Invalid move type.")

    def get_total_volume(self) -> int:
        """
        Returns the total volume of the triangulation.

        Returns:
            int: Total volume.
        """
        return self.size
    
# Example usage:
universe = Universe(total_time=3, initial_slice_size=3)

# pool = Pool(capacity=5, silent=True)
# obj1 = pool.occupy()
# obj2 = pool.occupy()
# obj3 = pool.occupy()

# pool.free(obj1)

# obj4 = pool.occupy()

