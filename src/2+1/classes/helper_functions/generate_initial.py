# generate_initial.py
#
# Author: Seda den Boer
# Date: 17-02-2024
#
# Description: this script generates an initial
# configuration for the 2+1D CDT model of S1xS2
# topology, where S1 is a circle and S2 is a sphere.

import numpy as np
import argparse
from typing import Tuple
import os

VERTICES_SPHERE = 5


def get_blocks_ur(v1: int, v2: int, v3: int, v4: int, v5: int, v6: int
                  ) -> Tuple[np.ndarray[int], np.ndarray[int], np.ndarray[int]]:
    """
    Create three NumPy arrays based on the given vertices for upper-right tetrahedra.

    Args:
        v1, v2, v3, v4, v5, v6 (int): Vertices of the simplices.

    Returns:
        Tuple[np.ndarray[int], np.ndarray[int], np.ndarray[int]]: The simplices of the upper right blocks.
    """
    return np.array([v1, v2, v3, v4]), np.array([v2, v3, v4, v6]), np.array([v2, v4, v5, v6])

def get_blocks_ul(v1: int, v2: int, v3: int, v4: int, v5: int, v6: int
                  ) -> Tuple[np.ndarray[int], np.ndarray[int], np.ndarray[int]]:
    """
    Create three NumPy arrays based on the given vertices for upper-left tetrahedra.

    Args:
        v1, v2, v3, v4, v5, v6 (int): Vertices of the simplices.

    Returns:
        Tuple[np.ndarray[int], np.ndarray[int], np.ndarray[int]]: The simplices of the upper left blocks.
    """
    return np.array([v1, v2, v3, v4]), np.array([v2, v3, v4, v5]), np.array([v3, v4, v5, v6])

def generate_sphere(t) -> np.ndarray[int]:
    """
    Generate the initial universe for the 2+1D CDT model with time t.
    The universe has S1xS2 topology, where S1 is a circle and S2 is a sphere.

    Args:
        t (int): Number of time slices in the initial universe.

    Returns:
        np.ndarray[int]: The simplices (IDs) of the initial universe.
    """
    start = 0
    vertices = [] 
    simplices = []

    # Generate the vertices of the sphere
    for i in range(t):
        v = np.arange(start, start + VERTICES_SPHERE, 1)
        start = max(v) + 1
        vertices.append(v)

    vertices.append(np.arange(0, VERTICES_SPHERE, 1))
    vertices = np.asarray(vertices)

    # Generate the simplices of the sphere
    for i in range(t):
        current_vertices = vertices[i]
        next_vertices = vertices[i + 1]

        v0, v1, v2, v3, v4 = current_vertices
        v0_next, v1_next, v2_next, v3_next, v4_next = next_vertices
        
        # Generate the simplices for the current time slice
        s = get_blocks_ur(v0, v1, v2, v0_next, v1_next, v2_next)
        simplices.extend(s)
        s = get_blocks_ur(v0, v3, v1, v0_next, v3_next, v1_next)
        simplices.extend(s)
        s = get_blocks_ur(v0, v2, v3, v0_next, v2_next, v3_next)
        simplices.extend(s)
        
        # Generate the simplices for the next time slice
        s = get_blocks_ul(v4, v2, v1, v4_next, v2_next, v1_next)
        simplices.extend(s)
        s = get_blocks_ul(v4, v1, v3, v4_next, v1_next, v3_next)
        simplices.extend(s)
        s = get_blocks_ul(v4, v3, v2, v4_next, v3_next, v2_next)
        simplices.extend(s)
        
    return np.asarray(simplices)

def find_pairs(simplices: np.ndarray[int]):
    """
    Find the adjacent simplices in the universe.

    Args:
        simplices (np.ndarray[int]): The simplices of the universe, where each simplex 
                                     is represented by an array of vertex IDs.

    Returns:
        np.ndarray[int]: An array where each element contains the indices of adjacent 
                         simplices for the corresponding simplex.
    """
    pairs = []

    # Iterate over each simplex to find its adjacent simplices
    for i in range(len(simplices)):
        triangle_pairs = []
        simplex = simplices[i]

        # Extract vertices of the current simplex
        p0, p1, p2, p3 = simplex

        # Generate the four possible triangle faces of the tetrahedron (simplex)
        triangle_pairs.append([p0, p1, p2])
        triangle_pairs.append([p0, p1, p3])
        triangle_pairs.append([p0, p2, p3])
        triangle_pairs.append([p1, p2, p3])

        adjacent_triangles = []

        # For each triangle face, find simplices sharing that face
        for triangle in triangle_pairs:
            tp0, tp1, tp2 = triangle

            # Find indices of simplices containing each vertex of the triangle
            indices0 = np.argwhere(simplices.T == tp0).T[1]
            indices1 = np.argwhere(simplices.T == tp1).T[1]
            indices2 = np.argwhere(simplices.T == tp2).T[1]

            # Find overlapping simplices that share the first two vertices
            overlap01 = list(set(indices0).intersection(indices1))
            # Find overlapping simplices that share the last two vertices
            overlap12 = list(set(indices1).intersection(indices2))

            # Find simplices that share all three vertices of the triangle
            adjacent_triangle = list(set(overlap01).intersection(overlap12))

            # Collect all simplices adjacent to the current triangle
            adjacent_triangles.append(adjacent_triangle)

        # Find unique indices of adjacent simplices and filter out the current simplex
        unique_pairs = np.unique(adjacent_triangles)
        filtered_pairs = [item for item in unique_pairs if item != i]
        filtered_pairs = np.asarray(filtered_pairs)
        pairs.append(filtered_pairs)

    return np.asarray(pairs)

def vertex_array_sphere(t):
    """
    Generate the vertex array for the sphere.

    Args:
        t (int): Number of time slices in the initial universe.

    Returns:
        np.ndarray[int]: The vertex array for the sphere.
    """
    # Generate an array of vertices for the sphere with t time slices
    vertices_sphere = np.arange(0, t * VERTICES_SPHERE, 1)
    
    # Create an array to store the vertex values 
    vertices = np.zeros(t * VERTICES_SPHERE, dtype=int)
    
    # Assign the time slice index to each vertex
    vertices[vertices_sphere] = np.repeat(np.arange(t), VERTICES_SPHERE)

    return vertices

def prepare_data(vertices: np.ndarray[int], pairs: np.ndarray[int], simplices: np.ndarray[int]) -> np.ndarray[int]:
    """
    Prepare the data for the initial universe, by appending the number of vertices and simplices
    in the correct order.

    Args:
        vertices (np.ndarray[int]): The vertices of the universe.
        pairs (np.ndarray[int]): The adjacent simplices in the universe.
        simplices (np.ndarray[int]): The simplices of the universe.

    Returns:
        np.ndarray[int]: The prepared data for the initial universe.
    """
    data = []

    # Add the number of vertices
    data.append(len(vertices))

    # Add the vertices
    for item in vertices:
        data.append(item)

    # Add the number of vertices and simplices
    data.append(len(vertices))
    data.append(len(pairs))

    # Add the simplices of the sphere and the universe
    for i in range(len(pairs)):
        for item in simplices[i]:
            data.append(item)
        
        for item in pairs[i]:
            data.append(item)

    # Add the number of simplices again
    data.append(len(pairs))

    return np.asarray(data)


if __name__ == "__main__":
    # Argument parser for t
    parser = argparse.ArgumentParser(description='Generate an initial universe for the 2+1D CDt model with time t.')
    parser.add_argument('--t', type=int, help='Number of time slices in the initial universe.')
    args = parser.parse_args()
    if args.t is None:
        raise ValueError("Please provide the number of time slices in the initial universe with --t [int].")
    elif args.t < 3:
        raise ValueError("The number of time slices must be greater than 3.")

    # Generate the initial universe
    simplices = generate_sphere(args.t)
    pairs = find_pairs(simplices)
    vertices = vertex_array_sphere(args.t)
    data = prepare_data(vertices, pairs, simplices)

    # Make a directory to store the initial universes
    PATH = '../initial_universes'
    if not os.path.exists(PATH):
        os.makedirs(PATH)
    filename = f'{PATH}/initial_t{args.t}.CDT'

    # Write the initial universe to a file
    with open(filename, 'w') as f:
        f.write("0\n")
        for item in data:
            f.write("{}\n".format(int(item)))