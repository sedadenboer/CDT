# vertex.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:

class Vertex:
    """
    Represents a node in the universe, corresponding to a
    vertex in the triangulation.

    Attributes:
    ID: Unique identifier for the node.
    time: Time slice to which the node belongs.
    future_vertices: List of nodes in the future (next time slice) connected to this node.
    past_vertices: List of nodes in the past (previous time slice) connected to this node.
    right_neighbour: The spatially adjacent node to the right (connected in space).
    triangles: List of triangles that contain this node.
    """
    capacity = 100
    
    def __init__(self, ID: int, time: int) -> None:
        self.ID = ID
        self.time = time
        self.future_vertices = []
        self.past_vertices = []
        self.right_neighbour = None
        self.triangles = []
    
    def get_number_of_future_neighbours(self) -> int:
        """
        Returns the number of future neighbours of this node.

        Returns:
            int: Number of future neighbours.
        """
        return len(self.future_vertices)
    
    def get_number_of_past_neighbours(self) -> int:
        """
        Returns the number of past neighbours of this node.

        Returns:
            int: Number of past neighbours.
        """
        return len(self.past_vertices)
    
    def get_total_number_of_neighbours(self) -> int:
        """
        Returns the total number of neighbours of this node.

        Returns:
            int: Total number of neighbours.
        """
        return self.get_number_of_future_neighbours() + self.get_number_of_past_neighbours()

    def add_triangle(self, triangle):
        """
        Add triangle to list of triangles that contain this node.
        """
        self.triangles.append(triangle)

    def __str__(self) -> str:
        return f"{self.ID}"