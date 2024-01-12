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
        self.triangulation = self.initialise_triangulation()

    def initialise_triangulation(self) -> list["Vertex"]:
        """
        initialise the grid with vertices and triangles.

        Returns:
            List["Vertex"]: List of time slices with vertices.
        """
        time = self.total_time
        width = self.initial_slice_size

        # ---- Initialise triangles ----
        # Create grid of blanco vertices
        vertex_grid = [
            [Vertex(self.vertex_pool.occupy(), t) for _ in range(width)] for t in range(time)
        ]

        # Create spatial edges
        for t in range(time):
            for i in range(width - 1):
                # Set right neighbour
                vertex_grid[t][i].right_neighbour = vertex_grid[t][i + 1]

            # Set right neighbour for last vertex
            vertex_grid[t][-1].right_neighbour = vertex_grid[t][0]

        # Create time-like edges
        for t in range(time - 1):
            for i in range(width - 1):
                # Set future vertices
                current_vertex = vertex_grid[t][i]
                future_vertex1 = vertex_grid[t + 1][i]
                future_vertex2 = vertex_grid[t + 1][i + 1]
                
                # Add future vertices to current vertex
                current_vertex.future_vertices.extend([future_vertex1, future_vertex2])
                future_vertex1.past_vertices.append(current_vertex)
                future_vertex2.past_vertices.append(current_vertex)

            # Set future vertices for last vertex
            current_vertex = vertex_grid[t][width - 1]
            future_vertex1 = vertex_grid[t + 1][width - 1]
            future_vertex2 = vertex_grid[t + 1][0]
            
            # Add future vertices to current vertex
            current_vertex.future_vertices.extend([future_vertex1, future_vertex2])
            future_vertex1.past_vertices.append(current_vertex)
            future_vertex2.past_vertices.append(current_vertex)

        # ---- CHECKING ----
        # print vertex_grid by id
        for time_slice in vertex_grid:
            print([vertex.ID for vertex in time_slice])

        # Print all past and future neighbours
        print("id, right_neighbour, future_vertices, past_vertices")
        for t in vertex_grid:
            for v in t:
                print(v.ID, v.right_neighbour.ID, [vertex.ID for vertex in v.future_vertices], [vertex.ID for vertex in v.past_vertices])

        print(self.vertex_pool.elements)
        print(self.vertex_bag.elements)
        print(self.vertex_bag.indices)
        self.vertex_bag.log()

        # ---- Initialise triangles ----
        
        # Create empty list to store triangle layers
        triangle_layers = []

        # Create triangles for each time slice
        for t in range(time - 1):
            # Create empty list to store triangles
            triangles = []
            for i in range(width):
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
                    vl=vertex_grid[(t + 1) % time][i],
                    vr=vertex_grid[(t + 1) % time][(i + 1) % width],
                    vc=vertex_grid[t][i]
                )

                triangle2.set_vertices(
                    vl=vertex_grid[t][i],
                    vr=vertex_grid[t][(i + 1) % width],
                    vc=vertex_grid[(t + 1) % time][(i + 1) % width]
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

            # Add list of triangles to list of time slices
            triangle_layers.append(triangles)

        # Set the right, left and center triangles for each triangle
        for i, layer in enumerate(triangle_layers):
            layer_length = len(layer)
            for j, triangle in enumerate(layer):
                # Get the triangles to the left and right of the current triangle
                triangle_left = layer[(j - 1) % layer_length]
                triangle_right = layer[(j + 1) % layer_length]

                # If the triangle is on the first or last layer, the center triangle is None
                if i == 0:
                    triangle_center = triangle_layers[(i + 1) % time][(j + 1) % layer_length] if triangle.is_downwards() else None
                # If the triangle is on the second or second to last layer, the center triangle is None
                elif i == time - 2:
                    triangle_center = triangle_layers[(i - 1) % time][(j - 1) % layer_length] if triangle.is_upwards() else None
                # If the triangle is on any other layer, the center triangle is the triangle in the next layer
                else:
                    triangle_center = triangle_layers[(i + 1) % time][(j + 1) % layer_length] if triangle.is_downwards() else triangle_layers[(i - 1) % time][(j - 1) % layer_length]

                triangle.set_triangle_left(triangle_left)
                triangle.set_triangle_right(triangle_right)
                triangle.set_triangle_center(triangle_center)

        # ---- CHECKING ----
        # Print triangles list with Id
        print("\nTriangles:")
        for layer in triangle_layers:
            print([triangle.ID for triangle in layer])
        
        for layer in triangle_layers:
            print([triangle.type for triangle in layer])

        # Print all triangles and their neighbour triangles
        # if none then print none
        print("\nid, tl, tr, tc")
        for t in triangle_layers:
            for triangle in t:
                print(triangle.ID,";", triangle.tl.ID if triangle.tl else None, triangle.tr.ID if triangle.tr else None, triangle.tc.ID if triangle.tc else None)

        # Print all triangles
        print("\nid, vl, vr, vc")
        for t in triangle_layers:
            for triangle in t:
                print(triangle.ID, ";", triangle.vl.ID, triangle.vr.ID, triangle.vc.ID)

        print()
        print(self.triangle_pool.elements)
        print(self.triangle_bag.elements)
        print(self.triangle_bag.indices)

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
universe = Universe(total_time=5, initial_slice_size=3)

