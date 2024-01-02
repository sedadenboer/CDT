# universe.py
#
# Author: Seda den Boer
# Date: 02-01-2024
# 
# Description:


class SpatialNode(object):
    """
    Represents a spatial node in the universe, corresponding to a
    vertex in the triangulation.

    Attributes:
    ID: Unique identifier for the node.
    time_slice: Time slice to which the node belongs.
    future_nodes: List of nodes in the future (next time slice) connected to this node.
    past_nodes: List of nodes in the past (previous time slice) connected to this node.
    right_neighbour: The spatially adjacent node to the right (connected in space).
    boundary: Indicates whether the node is on the boundary (0 or 1).
    """
    def __init__(self, ID, time_slice) -> None:
        self.ID = ID
        self.time_slice = time_slice
        self.future_nodes = []
        self.past_nodes = []
        self.right_neighbour = None 
        self.boundary = 0  
    
    def __str__(self) -> str:
        return f"{self.ID}"


class Universe(object):
    """
    The Universe class represents the current state of the triangulation
    and stores properties of the geometry in a convenient matter. It also
    provides member functions that carry out changes on the geometry.
    """
    def __init__(self, time_slices, spatial_points) -> None:
        self.grid = []
        self.initialize_grid(time_slices, spatial_points)

    def initialize_grid(self, time_slices, spatial_points):
        """
        """
        # ---- Intialize static grid ----

        # Create empty list to store time slices with spatial nodes
        time_slice_list = []

        for time_slice in range(time_slices):
            # Create empty list to store spatial nodes
            spatial_point_list = []

            for point in range(spatial_points):
                # Create spatial node object and append to list
                spatial_point = SpatialNode(point, time_slice)
                spatial_point_list.append(spatial_point)

            # Append list of spatial nodes to list of time slices
            time_slice_list.append(spatial_point_list)

        # Set the universe's grid to the list of time slices
        self.grid = time_slice_list

        # ---- Create spatial edges ----

        
        # ---- Create temporal edges ----


    def __str__(self) -> str:
        return "\n".join([" ".join([str(node) for node in time_slice]) for time_slice in self.grid])
        
