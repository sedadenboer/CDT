# universe.py
#
# Author: Seda den Boer
# Date: 17-02-2024	
# 
# Description: The Universe class that represents the current state of the triangulation.

from __future__ import annotations
from typing import cast
from vertex import Vertex
from triangle import Triangle
from tetra import Tetrahedron
from pool import Pool
from bag import Bag
import pickle
import resource
import sys
import random

max_rec = 0x100000

# May segfault without this line. 0x100 is a guess at the size of each stack frame.
resource.setrlimit(resource.RLIMIT_STACK, [0x100 * max_rec, resource.RLIM_INFINITY])
sys.setrecursionlimit(max_rec)


class Universe:
    """
    The Universe class represents the current state of the triangulation
    and stores properties of the geometry in a convenient matter. It also
    provides member functions that carry out changes on the geometry.
    """

    def __init__(self, geometry_infilename: str, VERTEX_CAPACITY: int = 1000000):
        VERTEX_CAPACITY = VERTEX_CAPACITY
        TRIANGLE_CAPACITY = 2 * VERTEX_CAPACITY
        TETRAHEDRON_CAPACITY = 2 * VERTEX_CAPACITY

        # Initialize the pools
        self.vertex_pool = Pool(VERTEX_CAPACITY)
        self.triangle_pool = Pool(TRIANGLE_CAPACITY)
        self.tetrahedron_pool = Pool(TETRAHEDRON_CAPACITY)

        # Initialize the bags
        self.tetras_all = Bag(TETRAHEDRON_CAPACITY)
        self.tetras_31 = Bag(TETRAHEDRON_CAPACITY)
        self.vertices_all = Bag(VERTEX_CAPACITY)
        self.vertices_six = Bag(VERTEX_CAPACITY)

        self.initialise(geometry_infilename)

    def initialise(self, geometry_infilename: str) -> bool:
        """
        Initializes the Universe with a given geometry.

        Args:
            geometry_filename (str): Name of the file with the geometry.
            fID_ (int): Identifier of the file.

        Returns:
            bool: True if the initialization was successful, otherwise False.
        """
        # Read the geometry file
        with open(geometry_infilename, 'r') as infile:
            assert infile

            # First line is a switch indicating whether tetrahedron data is ordered by convention
            ordered = bool(infile.readline().strip())

            # Read the number of vertices
            n0 = int(infile.readline())
            print(f"n0: {n0}")
            max_time = 0
            vs = []

            # Read each vertex
            for _ in range(n0):
                line = infile.readline().strip()
                vertex = Vertex(time=int(line))
                vertex.ID = self.vertex_pool.occupy(vertex)
                vs.append(vertex)
                if vertex.time > max_time:
                    max_time = vertex.time

            # Check consistency of the number of vertices
            line = int(infile.readline())
            if line != n0:
                print(f"Error: n0 = {n0}, line = {line}")
                return False

            # Calculate the number of time slices based on the maximum time
            self.n_slices = max_time + 1
            self.slab_sizes = [0] * self.n_slices
            self.slice_sizes = [0] * self.n_slices

            # Read the number of tetrahedra
            n3 = int(infile.readline()) 
            print(f"n3: {n3}")

            # Make n3 tetrahedra
            for _ in range(n3):
                tetra = Tetrahedron()
                tetra.ID = self.tetrahedron_pool.occupy(tetra)

            # Add vertices and neighbouring tetrahedra to the tetrahedra
            for i in range(n3):
                tetra = self.tetrahedron_pool.get(i)

                # Read the vertices and neighboring tetrahedra of the tetrahedron
                tetra_vs = []
                tetra_ts = []
                
                # Read the vertices of the tetrahedron
                for _ in range(4):
                    line = infile.readline().strip()
                    tetra_vs.append(self.vertex_pool.get(int(line)))

                tetra.set_vertices(tetra_vs[0], tetra_vs[1], tetra_vs[2], tetra_vs[3])
          
                # Read the neighboring tetrahedra of the tetrahedron
                for _ in range(4):
                    line = infile.readline().strip()
                    tetra_ts.append(self.tetrahedron_pool.get(int(line)))
                
                # If the tetrahedron is a 3-1 tetrahedron, update the size of the slices
                if tetra.is_31():
                    # Set tetrahedron neighbors for each vertex
                    for j in range(3):
                        vertex = tetra_vs[j]
                        vertex.set_tetra(tetra)

                tetra.set_tetras(tetra_ts[0], tetra_ts[1], tetra_ts[2], tetra_ts[3])

                # Update bags
                self.tetras_all.add(tetra.ID)
                if tetra.is_31():
                    self.tetras_31.add(tetra.ID)

                # Update the sizes
                self.slab_sizes[tetra.get_vertices()[1].time] += 1
                if tetra.is_31():
                    self.slice_sizes[tetra.get_vertices()[0].time] += 1

            # Check consistency of the number of tetrahedra
            line = int(infile.readline())
            if line != n3:
                print(f"Error: n3 = {n3}, line = {line}")
                return False
            
            print(f"Read {geometry_infilename}.")

            # If the tetrahedra are not ordered by convention, order them
            if not ordered:
                for tetra in self.tetrahedron_pool.get_objects():
                    tnbr = tetra.get_tetras()
                    t012, t013, t023, t123 = -1, -1, -1, -1
                    for tn in tnbr:
                        if not tn.has_vertex(tetra.get_vertices()[3]):
                            t012 = tn
                            continue
                        if not tn.has_vertex(tetra.get_vertices()[2]):
                            t013 = tn
                            continue
                        if not tn.has_vertex(tetra.get_vertices()[1]):
                            t023 = tn
                            continue
                        if not tn.has_vertex(tetra.get_vertices()[0]):
                            t123 = tn
                            continue
                    
                    assert t012 >= 0 and t013 >= 0 and t023 >= 0 and t123 >= 0
                    tetra.set_tetras(t012, t013, t023, t123)
            
            # Calculate the number of tetrahedra that contain each vertex and the spatial coordination number
            for vertex in vs:
                cnum = scnum = 0

                for tetra in self.tetrahedron_pool.get_objects():
                    if tetra.has_vertex(vertex):
                        cnum += 1
                    if not tetra.is_31():
                        continue
                    if tetra.get_vertices()[0] == vertex or tetra.get_vertices()[1] == vertex or tetra.get_vertices()[2] == vertex:
                        scnum += 1
                
                vertex.cnum = cnum
                vertex.scnum = scnum
        
        return True


    def insert_vertex(self, triangle_id: int) -> tuple[Vertex, Triangle, Triangle]:
        """
        Insert a vertex into the triangulation.

        Args:
            triangle (Triangle): Triangle to insert vertex into.

        Returns:
            tuple[Vertex, Triangle, Triangle]: The inserted vertex and the two triangles
            that are created with it.
        """
        pass
    
    def remove_vertex(self, vertex_id: int) -> tuple[Vertex, Triangle, Triangle]:
        """
        Remove a vertex from the triangulation.

        Args:
            vertex (Vertex): Vertex to remove.

        Returns:
            tuple[Vertex, Triangle, Triangle]: The removed vertex and the two triangles
            that are removed with it.
        """
        pass

    def flip_edge(self, triangle_id : int) -> tuple[Triangle, Triangle]:
        """
        Flip an edge in the triangulation.

        Args:
            triangle (Triangle): Triangle to flip edge in.

        Returns:
            tuple[Triangle, Triangle]: The two triangles that are flipped (t, tr).
        """
        pass

    def get_total_size(self) -> int:
        """
        Get the total size of the triangulation.

        Returns:
            int: Total size of the triangulation.
        """
        return self.vertex_pool.get_number_occupied()
    
    def sort_vertices_periodic(self, vertices):
        pass

    def get_triangulation_state(self):
        """
        Get the current state of the triangulation.
        """
        pass

    def print_state(self):
        """
        Print the current state of the triangulation.
        """
        pass

    def check_validity(self):
        """
        Check the validity of the triangulation.
        """
        # Check that each tetrahedron has 4 vertices and 4 neighbouring tetrahedra         
        for i in self.tetrahedron_pool.get_objects():
            assert len(i.get_vertices()) == 4
            assert len(i.get_tetras()) == 4

    def save_to_file(self, filename):
        """
        Save the state of the Universe to a file using pickle.

        Args:
            filename (str): The name of the file to save the state to.
        """
        with open(filename, 'wb') as file:
            pickle.dump(self.__dict__, file)


if __name__ == "__main__":
    u = Universe(geometry_infilename="initial_universes/output.txt")