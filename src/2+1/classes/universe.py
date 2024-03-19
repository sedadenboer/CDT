# universe.py
#
# Author: Seda den Boer
# Date: 17-02-2024	
# 
# Description: The Universe class that represents the current state of the triangulation.

from __future__ import annotations
from typing import Any, List
from vertex import Vertex
from halfedge import HalfEdge
from triangle import Triangle
from tetra import Tetrahedron
from pool import Pool
from bag import Bag
import resource
import sys

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
    class Capacity:
        VERTEX = 1000000
        HALFEDGE = 2 * VERTEX
        TRIANGLE = 2 * VERTEX
        TETRAHEDRON = 2 * VERTEX
    
    class Constants:
        CNUM_ADD = 6
        SCNUM_ADD = 3
        N_VERTICES_TETRA = N_TETRA_NEIGHBOURS = 4
        N_VERTICES_TRIANGLE = N_HALFEDGES_TETRA = 3
        APEX_INDEX_31 = OPPOSITE_TETRA_INDEX_31 = 3
        APEX_INDEX_13 = OPPOSITE_TETRA_INDEX_13 = 0
        
    def __init__(self, geometry_infilename: str, strictness: int = 3):
        """
        Args:
            geometry_infilename (str): Name of the file with the geometry.
            strictness (int): The strictness of the manifold conditions.
        """
        self.strictness = strictness

        # Initialize the pools
        self.vertex_pool = Pool(Universe.Capacity.VERTEX)
        self.halfedge_pool = Pool(Universe.Capacity.HALFEDGE)
        self.triangle_pool = Pool(Universe.Capacity.TRIANGLE)
        self.tetrahedron_pool = Pool(Universe.Capacity.TETRAHEDRON)

        # Initialize the bags
        self.tetras_31 = Bag(Universe.Capacity.TETRAHEDRON)

        self.vertex_neighbours = []
        self.triangle_neighbours = []

        self.initialise(geometry_infilename)

    def initialise(self, geometry_infilename: str) -> bool:
        """
        Initializes the Universe with a given geometry.

        Args:
            geometry_filename (str): Name of the file with the geometry.

        Returns:
            bool: True if the initialization was successful, otherwise False.
        """
        # Read the geometry file
        with open(geometry_infilename, 'r') as infile:
            assert infile

            # First line is a switch indicating whether tetrahedron data is ordered by convention
            ordered = int(infile.readline().strip())

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

            # Add vertices and neighboring tetrahedra to the tetrahedra
            for tetra in self.tetrahedron_pool.get_objects():
                # Read the vertices and neighboring tetrahedra of the tetrahedron
                tetra_vs = [self.vertex_pool.get(int(infile.readline().strip())) for _ in range(self.Constants.N_VERTICES_TETRA)]
                tetra_ts = [self.tetrahedron_pool.get(int(infile.readline().strip())) for _ in range(self.Constants.N_TETRA_NEIGHBOURS)]
                
                # Set the vertices and neighboring tetrahedra of the tetrahedron
                tetra.set_vertices(*tetra_vs)
                tetra.set_tetras(*tetra_ts)

                # Update the size of the slabs
                self.slab_sizes[tetra.get_vertices()[1].time] += 1

                # If the tetrahedron is a (3,1)-tetrahedron, update the size of the slices
                if tetra.is_31():
                    # Set base (3,1)-tetrahedron for each vertex
                    for i, vertex in enumerate(tetra_vs):
                        if i != self.Constants.APEX_INDEX_31:
                            vertex.set_tetra(tetra)

                    # Add to the bag and update the slice size
                    self.tetras_31.add(tetra.ID)
                    self.slice_sizes[tetra.get_vertices()[0].time] += 1
    
            # Check consistency of the number of tetrahedra
            line = int(infile.readline())
            if line != n3:
                print(f"Error: n3 = {n3}, line = {line}")
                return False
            
            print(f"Read {geometry_infilename}.")

        # Make sure to order the tetrahedra neighbours by convention
        if ordered == 0:
            for tetra in self.tetrahedron_pool.get_objects():
                tnbr = tetra.get_tetras()
                vs = tetra.get_vertices()

                t012, t013, t023, t123 = -1, -1, -1, -1
                
                # Find tetra opposite to vertices (ordered)
                for tn in tnbr:
                    if not tn.has_vertex(vs[3]):
                        t012 = tn
                        continue
                    if not tn.has_vertex(vs[2]):
                        t013 = tn
                        continue
                    if not tn.has_vertex(vs[1]):
                        t023 = tn
                        continue
                    if not tn.has_vertex(vs[0]):
                        t123 = tn
                        continue
    
                assert t012.ID >= 0 and t013.ID >= 0 and t023.ID >= 0 and t123.ID >= 0

                tetra.set_tetras(t123, t023, t013, t012)
        
        # Calculate the number of tetrahedra that contain each vertex and the spatial coordination number
        for vertex in vs:
            cnum = 0
            scnum = 0

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

    def add(self, tetra31_id: int) -> bool:
        """
        (2,6)-move in the triangulation. Adds a vertex, two (3,1)-tetrahedra,
        and two (1,3)-tetrahedra.

        Args:
            tetra31_id (int): ID of the (3,1)-tetrahedron to perform add move in.

        Returns:
            bool: True if the move was added successfully, otherwise False.
        """
        # Get the tetrahedron object, the time, and its opposite neighbor
        t = self.tetrahedron_pool.get(tetra31_id)
        if not t.is_31():
            return False
        time = t.get_vertices()[0].time
        tv = t.get_tetras()[3]
        if not tv.is_13():
            return False

        # print(f"Add changes t: {t.ID} and tv: {tv.ID}.")

        # Create a new vertex and update its cnum and scnum
        new_vertex = Vertex(time=time)
        self.vertex_pool.occupy(new_vertex)
        new_vertex.cnum = Universe.Constants.CNUM_ADD
        new_vertex.scnum = Universe.Constants.SCNUM_ADD

        # Get relevant vertices of the two tetrahedra, vt=vtop, vb=vbottom
        v0 = t.get_vertices()[0]
        v1 = t.get_vertices()[1]
        v2 = t.get_vertices()[2]
        vt = t.get_vertices()[3]
        vb = tv.get_vertices()[0]
       
        # Create the new tetrahedra for top and bottom
        tn01 = Tetrahedron()
        tn12 = Tetrahedron()
        tn20 = Tetrahedron()
        tvn01 = Tetrahedron()
        tvn12 = Tetrahedron()
        tvn20 = Tetrahedron()
        
        # Add to pool and bag
        self.tetrahedron_pool.occupy(tn01)
        self.tetrahedron_pool.occupy(tn12)
        self.tetrahedron_pool.occupy(tn20)
        self.tetrahedron_pool.occupy(tvn01)
        self.tetrahedron_pool.occupy(tvn12)
        self.tetrahedron_pool.occupy(tvn20)
        self.tetras_31.add(tn01.ID)
        self.tetras_31.add(tn12.ID)
        self.tetras_31.add(tn20.ID)
    
        # Get neighbouring tetrahedra of the two opposing tetrahedra
        to0 = t.get_tetra_opposite(v0)
        to1 = t.get_tetra_opposite(v1)
        to2 = t.get_tetra_opposite(v2)
        tvo0 = tv.get_tetra_opposite(v0)
        tvo1 = tv.get_tetra_opposite(v1)
        tvo2 = tv.get_tetra_opposite(v2)

        # Set the vertices of the new tetrahedra
        tn01.set_vertices(v0, v1, new_vertex, vt)
        tn12.set_vertices(v1, v2, new_vertex, vt)
        tn20.set_vertices(v2, v0, new_vertex, vt)
        tvn01.set_vertices(vb, v0, v1, new_vertex)
        tvn12.set_vertices(vb, v1, v2, new_vertex)
        tvn20.set_vertices(vb, v2, v0, new_vertex)

        # Set neighbouring tetrahedra
        tn01.set_tetras(tn12, tn20, to2, tvn01)
        tn12.set_tetras(tn20, tn01, to0, tvn12)
        tn20.set_tetras(tn01, tn12, to1, tvn20)
        tvn01.set_tetras(tn01, tvn12, tvn20, tvo2)
        tvn12.set_tetras(tn12, tvn20, tvn01, tvo0)
        tvn20.set_tetras(tn20, tvn01, tvn12, tvo1)
        
        # Update the opposite tetra of neighbouring tetras
        to0.exchange_tetra_opposite(t.get_vertex_opposite(v0), tn12)
        to1.exchange_tetra_opposite(t.get_vertex_opposite(v1), tn20)
        to2.exchange_tetra_opposite(t.get_vertex_opposite(v2), tn01)
        tvo0.exchange_tetra_opposite(tv.get_vertex_opposite(v0), tvn12)
        tvo1.exchange_tetra_opposite(tv.get_vertex_opposite(v1), tvn20)
        tvo2.exchange_tetra_opposite(tv.get_vertex_opposite(v2), tvn01)

        # Update the sizes
        self.slab_sizes[time] += 2
        self.slab_sizes[(time - 1 + self.n_slices) % self.n_slices] += 2
        self.slice_sizes[time] += 2

        # Remove original tetrahedra from pool and bag
        self.tetrahedron_pool.free(t.ID)
        self.tetras_31.remove(t.ID)
        self.tetrahedron_pool.free(tv.ID)

        # Update the vertices with the tetrahedra they are part of, and their cnum and scnum
        new_vertex.set_tetra(tn01)
        v0.set_tetra(tn01)
        v1.set_tetra(tn12)
        v2.set_tetra(tn20)

        v0.scnum += 1
        v1.scnum += 1
        v2.scnum += 1
        v0.cnum += 2
        v1.cnum += 2
        v2.cnum += 2
        vt.cnum += 2
        vb.cnum += 2

        # print(f"Added vertex {new_vertex.ID}.")
        # print(f"Added tetras: tn01.ID={tn01.ID}, tn12.ID={tn12.ID}, tn20.ID={tn20.ID}, tvn01.ID={tvn01.ID}, tvn12.ID={tvn12.ID}, tvn20.ID={tvn20.ID}.")
        # print(f"types of new tetras: tn01: {tn01.type}, tn12: {tn12.type}, tn20: {tn20.type}, tvn01: {tvn01.type}, tvn12: {tvn12.type}, tvn20: {tvn20.type}.")
        return True

    def delete(self, vertex_id: int) -> bool:
        """
        (6,2)-move in the triangulation. Deletes a vertex, two (1,3)-tetrahedra,
        and two (3,1)-tetrahedra.

        Args:
            vertex (int): ID of the vertex to perform delete move on.

        Returns:
            bool: True if the move was deleted successfully, otherwise False.
        """
         # Get the vertex object, its time and the two opposing tetrahedra it is part of
        vertex = self.vertex_pool.get(vertex_id)
        if vertex.cnum != self.Constants.CNUM_ADD:
            return False
        time = vertex.time
        t01 = vertex.get_tetra()
        tv01 = t01.get_tetras()[3]

        if not t01.is_31():
            return False
        if not tv01.is_13():
            return False

        # Get the vertex index in the tetrahedron
        vpos = t01.get_vertices().index(vertex)
        assert vpos >= 0

        # Get the base vertices of the two tetrahedra
        v0 = t01.get_vertices()[(vpos + 1) % 3]
        v1 = t01.get_vertices()[(vpos + 2) % 3]
        v2 = t01.get_vertex_opposite(v0)

        # Get all tetrahedra that need to be updated
        t12 = t01.get_tetra_opposite(v0)
        t20 = t01.get_tetra_opposite(v1)
        tv12 = tv01.get_tetra_opposite(v0)
        tv20 = tv01.get_tetra_opposite(v1)
        
        # print(f"Deleting vertex {vertex.ID}.")

        if not t01.is_31():
            return False
        if not t12.is_31():
            return False
        if not t20.is_31():
            return False
        if not tv01.is_13():
            return False
        if not tv12.is_13():
            return False
        if not tv20.is_13():
            return False

        # Get the neighbours of tetras that will be adjusted
        to01 = t01.get_tetra_opposite(vertex)
        to12 = t12.get_tetra_opposite(vertex)
        to20 = t20.get_tetra_opposite(vertex)
        tvo01 = tv01.get_tetra_opposite(vertex)
        tvo12 = tv12.get_tetra_opposite(vertex)
        tvo20 = tv20.get_tetra_opposite(vertex)

        # Disallow tadpole insertions
        if self.strictness == 1:
            if v0.scnum < 3:
                return False
            if v1.scnum < 3:
                return False
            if v2.scnum < 3:
                return False
        # Disallow self-energy insertions
        elif self.strictness >= 2:
            if v0.scnum < 4:
                return False
            if v1.scnum < 4:
                return False
            if v2.scnum < 4:
                return False
       
        # Create new tetrahedra
        tn = Tetrahedron()
        tvn = Tetrahedron()
        self.tetrahedron_pool.occupy(tn)
        self.tetrahedron_pool.occupy(tvn)
        self.tetras_31.add(tn.ID)

        # Get apex of oppposing tetrahedra
        vt = t01.get_vertices()[3]
        vb = tv01.get_vertices()[0]

        # Set the vertices and neighbouring tetrahedra of the new tetrahedra
        tn.set_vertices(v0, v1, v2, vt)
        tvn.set_vertices(vb, v0, v1, v2)
        tn.set_tetras(to12, to20, to01, tvn)
        tvn.set_tetras(tn, tvo12, tvo20, tvo01)

        # Update the vertices with the tetrahedra they are part of, and their cnum and scnum
        v0.set_tetra(tn)
        v1.set_tetra(tn)
        v2.set_tetra(tn)

        v0.scnum -= 1
        v1.scnum -= 1
        v2.scnum -= 1
        v0.cnum -= 2
        v1.cnum -= 2
        v2.cnum -= 2
        vt.cnum -= 2
        vb.cnum -= 2

        # Exchange given tetrahedron for the given vertex with opposite tetrahedron
        to01.exchange_tetra_opposite(t01.get_vertex_opposite(vertex), tn)
        to12.exchange_tetra_opposite(t12.get_vertex_opposite(vertex), tn)
        to20.exchange_tetra_opposite(t20.get_vertex_opposite(vertex), tn)
        tvo01.exchange_tetra_opposite(tv01.get_vertex_opposite(vertex), tvn)
        tvo12.exchange_tetra_opposite(tv12.get_vertex_opposite(vertex), tvn)
        tvo20.exchange_tetra_opposite(tv20.get_vertex_opposite(vertex), tvn)
            
        # Update pool and bag
        self.vertex_pool.free(vertex.ID)
        self.tetrahedron_pool.free(t01.ID)
        self.tetrahedron_pool.free(t12.ID)
        self.tetrahedron_pool.free(t20.ID)
        self.tetras_31.remove(t01.ID)
        self.tetras_31.remove(t12.ID)
        self.tetras_31.remove(t20.ID)
        self.tetrahedron_pool.free(tv01.ID)
        self.tetrahedron_pool.free(tv12.ID)
        self.tetrahedron_pool.free(tv20.ID)

        # Update the sizes
        self.slab_sizes[time] -= 2
        self.slab_sizes[(time - 1 + self.n_slices) % self.n_slices] -= 2
        self.slice_sizes[time] -= 2

        # print(f"Deleted vertex {vertex.ID}, and tetras t01: {t01.ID}, t12: {t12.ID}, t20: {t20.ID}, tv01: {tv01.ID}, tv12: {tv12.ID}, tv20: {tv20.ID}.")

        return True

    def flip(self, tetra012_id: int, tetra230_id: int) -> bool:
        """
        (4,4)-move in the triangulation. Flips shared spatial edge base triangles 
        of two tetrahedra. This move checks has to adhere to the simplicial manifold
        conditions:
        1. The vertices 0 and 2 both have spatial coordination number (i.e. the
        number of vertices connected to the vertex by spacelike links) of at
        least 4. Vertex 0 and 2 are connected by the spacelike link that 
        has to be flipped.
        2. There is no link connecting vertex 1 to vertex 3. (Note that such a link
        is always spacelike.) Vertex 1 and 3 are the vertices are the vertices
        that are supposed to be connected by the flipped spacelink.

        Args:
            tetra012 (int): ID of the first tetrahedron.
            tetra230 (int): ID of the second tetrahedron.

        Returns:
            bool: True if the move was flipped successfully, otherwise False.
        """
        t012 = self.tetrahedron_pool.get(tetra012_id)
        t230 = self.tetrahedron_pool.get(tetra230_id)
        if not t012.is_31():
            return False
        if not t230.is_31():
            return False
        
        # print(f"Flipping t012: {t012.ID} and t230: {t230.ID}.")

        # Get the vertices of the base triangles that are going to be linked
        v1 = t012.get_vertex_opposite_tetra(t230)
        v3 = t230.get_vertex_opposite_tetra(t012)

        # Get the remaining base vertices
        v1pos = t012.get_vertices().index(v1)
        v2 = t012.get_vertices()[(v1pos + 1) % 3]
        v0 = t012.get_vertices()[(v1pos + 2) % 3]

        # Manifold conditions 1: v0 and v2 should have spatial coordination number of at least 4
        if self.strictness >= 1:
            if v1 == v3:
                return False
        if self.strictness >= 2:
            if v0.scnum  < 4:
                return False
            if v2.scnum < 4:
                return False
        if self.strictness >= 3:
            if v1.check_vertex_neighbour(v3):
                return False
            
        # # Manifold conditions 2: there can be no link between v1 and v3
        # for tetra in self.tetrahedron_pool.get_objects():
        #     if v0 in tetra.get_vertices() and v3 in tetra.get_vertices():
        #         return False

        # Get the opposite (1,3)-tetrahedra
        tv012 = t012.get_tetras()[3]
        tv230 = t230.get_tetras()[3]

        if not tv012.is_13():
            return False
        if not tv230.is_13():
            return False
        if not tv012.check_neighbours_tetra(tv230):
            return False

        # Manifold conditions 2: there can be no link between v1 and v3
        if tv012.has_vertex(v3):
            return False
        if tv230.has_vertex(v1):
            return False
        if tv012.has_vertex(v3):
            return False
        if tv230.has_vertex(v1):
            return False
        
        # Get apex of the opposing tetrahedra
        vt = t012.get_vertices()[3]
        vb = tv012.get_vertices()[0]

        # Get opposite neighbouring tetrahedra
        ta01 = t012.get_tetra_opposite(v2)
        ta12 = t012.get_tetra_opposite(v0)
        ta23 = t230.get_tetra_opposite(v0)
        ta30 = t230.get_tetra_opposite(v2)
        tva01 = tv012.get_tetra_opposite(v2)
        tva12 = tv012.get_tetra_opposite(v0)
        tva23 = tv230.get_tetra_opposite(v0)
        tva30 = tv230.get_tetra_opposite(v2)
        
        # Make sure the move is valid
        if ta01 == t230:
            return False
        if ta23 == t012:
            return False
        if tva01 == tv230:
            return False
        if tva23 == tv012:
            return False
        
        # Opposite vertices in opposite tetrahedra
        t012vo2 = t012.get_vertex_opposite(v2)
        t230vo0 = t230.get_vertex_opposite(v0)
        tv012vo2 = tv012.get_vertex_opposite(v2)
        tv230vo0 = tv230.get_vertex_opposite(v0)

        # Assign new names to tetras to be updated
        tn013 = t230
        tn123 = t012
        tvn013 = tv230
        tvn123 = tv012

        # Set the vertices and neighbouring tetrahedra of the new tetrahedra
        tn013.set_vertices(v0, v1, v3, vt)
        tn123.set_vertices(v1, v2, v3, vt)
        tvn013.set_vertices(vb, v0, v1, v3)
        tvn123.set_vertices(vb, v1, v2, v3)
        tn013.set_tetras(tn123, ta30, ta01, tvn013)
        tn123.set_tetras(ta23, tn013, ta12, tvn123)
        tvn013.set_tetras(tn013, tvn123, tva30, tva01)
        tvn123.set_tetras(tn123, tva23, tvn013, tva12)

        # Exchange given tetra for the given vertex with opposite tetra
        ta01.exchange_tetra_opposite(t012vo2, tn013)
        ta23.exchange_tetra_opposite(t230vo0, tn123)
        tva01.exchange_tetra_opposite(tv012vo2, tvn013)
        tva23.exchange_tetra_opposite(tv230vo0, tvn123)

        v0.scnum -= 1
        v1.scnum += 1
        v2.scnum -= 1
        v3.scnum += 1

        if self.strictness >= 2:
            assert v0.scnum >= 3
            assert v2.scnum >= 3

        v0.cnum -= 2
        v1.cnum += 2
        v2.cnum -= 2
        v3.cnum += 2

        # Update vertices' tetrahedra
        v0.set_tetra(tn013)
        v2.set_tetra(tn123)

        # print(f"Flipped, updated tn013: {tn013.ID}, tn123: {tn123.ID}, tvn013: {tvn013.ID}, tvn123: {tvn123.ID}.")

        return True

    def shift_u(self, tetra31_id: int, tetra22_id: int) -> bool:
        """
        (2,3)-move in the triangulation. The shift move selects a (3,1)-simplex
        from the triangulation and replaces it with a neighboring (2,2)-simplex
        if conditions are met. This leads to a rearrangement of the local configuration
        without adding vertices.

        Args:
            tetra31 (int): ID of the (3,1)-tetrahedron.
            tetra22 (int): ID of the (2,2)-tetrahedron.

        Returns:
            bool: True if the move was shifted successfully, otherwise False.
        """
        # Get the 31 and 22 tetra neighbours
        t31 = self.tetrahedron_pool.get(tetra31_id)
        t22 = self.tetrahedron_pool.get(tetra22_id)
        if not t31.is_31():
            return False
        if not t22.is_22():
            return False
        if not t31.check_neighbours_tetra(t22):
            return False

        # print(f"Shifting t31: {t31.ID} and t22: {t22.ID}.")

        # Get the vertices that will be linked
        v0 = t31.get_vertex_opposite_tetra(t22)
        v1 = t22.get_vertex_opposite_tetra(t31)
        
        # The remaining vertices
        v3 = t31.get_vertices()[3]
        v0pos = t31.get_vertices().index(v0)
        v2 = t31.get_vertices()[(v0pos + 1) % 3]
        v4 = t31.get_vertices()[(v0pos + 2) % 3]

        # Get neighbouring tetrahedra that need to be updated after the move
        ta023 = t31.get_tetra_opposite(v4)
        ta034 = t31.get_tetra_opposite(v2)
        ta123 = t22.get_tetra_opposite(v4)
        ta124 = t22.get_tetra_opposite(v3)
        ta134 = t22.get_tetra_opposite(v2)

        # Check if the move is valid
        if ta023.has_vertex(v1):
            return False
        if ta123.has_vertex(v0):
            return False
        if ta034.has_vertex(v1):
            return False
        if ta134.has_vertex(v0):
            return False
        if v0.check_vertex_neighbour(v1):
            return False

        # Create new tetrahedra
        tn31 = Tetrahedron()
        tn22l = Tetrahedron()
        tn22r = Tetrahedron()

        # Add to pool and bag
        self.tetrahedron_pool.occupy(tn31)
        self.tetras_31.add(tn31.ID)
        self.tetrahedron_pool.occupy(tn22l)
        self.tetrahedron_pool.occupy(tn22r)

        # Update connectivity
        tn31.set_vertices(v0, v2, v4, v1)
        tn22l.set_vertices(v0, v2, v1, v3)
        tn22r.set_vertices(v0, v4, v1, v3)
        tn31.set_tetras(ta124, tn22r, tn22l, t31.get_tetras()[3])
        tn22l.set_tetras(ta123, tn22r, ta023, tn31)
        tn22r.set_tetras(ta134, tn22l, ta034, tn31)
        
        time = tn31.get_vertices()[0].time
        self.slab_sizes[time] += 1

        t31.get_tetras()[3].exchange_tetra_opposite(t31.get_tetras()[3].get_vertices()[0], tn31)
        ta023.exchange_tetra_opposite(t31.get_vertex_opposite(v4), tn22l)
        ta034.exchange_tetra_opposite(t31.get_vertex_opposite(v2), tn22r)
        ta123.exchange_tetra_opposite(t22.get_vertex_opposite(v4), tn22l)
        ta124.exchange_tetra_opposite(t22.get_vertex_opposite(v3), tn31)
        ta134.exchange_tetra_opposite(t22.get_vertex_opposite(v2), tn22r)

        v0.cnum += 2
        v1.cnum += 2

        # Remove the original tetrahedra from the neighbours
        self.tetrahedron_pool.free(t31.ID)
        self.tetras_31.remove(t31.ID)
        self.tetrahedron_pool.free(t22.ID)

        # Make sure the vertices are updated with the new 31 tetra
        tn31.get_vertices()[0].set_tetra(tn31)
        tn31.get_vertices()[1].set_tetra(tn31)
        tn31.get_vertices()[2].set_tetra(tn31)

        # print(f"Shifted, created tn31: {tn31.ID}, tn22l: {tn22l.ID}, tn22r: {tn22r.ID}.")

        return True
    
    def ishift_u(self, tetra31_id: int, tetra22l_id: int, tetra22r_id: int) -> bool:
        """
        (3,2)-move in the triangulation. The ishift move removes a (2,2)-simplex
        from the triangulation by replacing its shared timelike edge with a dual
        timelike triangle. This results in the rearrangement of the three-simplices
        in an inverse manner to the shift move.

        Args:
            tetra31 (int): ID of the (3,1)-tetrahedron.
            tetra22 (int): ID of the (2,2)-tetrahedron.
            tetra22r (int): ID of the (2,2)-tetrahedron.

        Returns:
            bool: True if the move was ishifted successfully, otherwise False.
        """
        t31 = self.tetrahedron_pool.get(tetra31_id)
        t22l = self.tetrahedron_pool.get(tetra22l_id)
        t22r = self.tetrahedron_pool.get(tetra22r_id)
        if not t31.is_31():
            return False
        if not t22l.is_22():
            return False
        if not t22r.is_22():
            return False
        if not t31.check_neighbours_tetra(t22l):
            return False
        if not t31.check_neighbours_tetra(t22r):
            return False
        
        # print(f"iShifting t31: {t31.ID}, t22l: {t22l.ID} and t22r: {t22r.ID}.")
        
        # Get the vertices of the interior triangle
        v1 = t31.get_vertices()[3]
        v3 = t22l.get_vertex_opposite_tetra(t31)
        v4 = t31.get_vertex_opposite_tetra(t22l)

        # The remaining vertices
        v4pos = t31.get_vertices().index(v4)
        v0 = t31.get_vertices()[(v4pos + 1) % 3]
        v2 = t31.get_vertices()[(v4pos + 2) % 3]

        # Get neighbouring tetrahedra that need to be updated after the move
        ta023 = t22l.get_tetra_opposite(v1)
        ta034 = t22r.get_tetra_opposite(v1)
        ta123 = t22l.get_tetra_opposite(v0)
        ta124 = t31.get_tetra_opposite(v0)
        ta134 = t22r.get_tetra_opposite(v0)

        # Make sure the move is valid
        if ta023.has_vertex(v4):
            return False
        if ta123.has_vertex(v4):
            return False
        if ta034.has_vertex(v2):
            return False
        if ta124.has_vertex(v3):
            return False
        if ta134.has_vertex(v2):
            return False

        # Create new tetrahedra
        tn31 = Tetrahedron()
        tn22 = Tetrahedron()
        self.tetrahedron_pool.occupy(tn31)
        self.tetras_31.add(tn31.ID)
        self.tetrahedron_pool.occupy(tn22)

        # Update connectivity
        tn31.set_vertices(v0, v2, v4, v3)
        tn22.set_vertices(v2, v4, v1, v3)
        tn31.set_tetras(tn22, ta034, ta023, t31.get_tetras()[3])
        tn22.set_tetras(ta134, ta123, tn31, ta124)

        t31.get_tetras()[3].exchange_tetra_opposite(t31.get_tetras()[3].get_vertices()[0], tn31)
        ta023.exchange_tetra_opposite(t22l.get_vertex_opposite(v1), tn31)
        ta034.exchange_tetra_opposite(t22r.get_vertex_opposite(v1), tn31)
        ta123.exchange_tetra_opposite(t22l.get_vertex_opposite(v0), tn22)
        ta124.exchange_tetra_opposite(t31.get_vertex_opposite(v0), tn22)
        ta134.exchange_tetra_opposite(t22r.get_vertex_opposite(v0), tn22)

        v0.cnum -= 2
        v1.cnum -= 2

        # Free the old tetrahedra
        self.tetrahedron_pool.free(t31.ID)
        self.tetras_31.remove(t31.ID)
        self.tetrahedron_pool.free(t22l.ID)
        self.tetrahedron_pool.free(t22r.ID)

        time = tn31.get_vertices()[0].time
        self.slab_sizes[time] -= 1

        # Make sure the base vertices are updated with the new 31 tetra
        tn31.get_vertices()[0].set_tetra(tn31)
        tn31.get_vertices()[1].set_tetra(tn31)
        tn31.get_vertices()[2].set_tetra(tn31)

        # print(f"iShifted, created tn31: {tn31.ID}, tn22: {tn22.ID}.")

        return True

    def shift_d(self, tetra13_id: int, tetra22_id: int) -> bool:
        """
        (2,3)-move in the triangulation. The shift move selects a (1,3)-simplex
        from the triangulation and replaces it with a neighboring (2,2)-simplex
        if conditions are met. This leads to a rearrangement of the local configuration
        without adding vertices.

        Args:
            tetra13 (int): ID of the (1,3)-tetrahedron.
            tetra22 (int): ID of the (2,2)-tetrahedron.

        Returns:
            bool: True if the move was shifted successfully, otherwise False.
        """
        # Get the tetrahedra
        t13 = self.tetrahedron_pool.get(tetra13_id)
        t22 = self.tetrahedron_pool.get(tetra22_id)
        t31 = t13.get_tetras()[0]
        if not t13.is_13():
            return False
        if not t22.is_22():
            return False
        if not t13.check_neighbours_tetra(t22):
            return False

        # print(f"Shifting t13: {t13.ID} and t22: {t22.ID}.")

        # Get the vertices that will be linked
        v0 = t13.get_vertex_opposite_tetra(t22)
        v1 = t22.get_vertex_opposite_tetra(t13)

        # The remaining vertices
        v3 = t13.get_vertices()[0] # Top
        v0pos = t31.get_vertices().index(v0)
        v2 = t31.get_vertices()[(v0pos + 1) % 3]
        v4 = t31.get_vertices()[(v0pos + 2) % 3]

        # Get the neighbouring tetrahedra
        ta023 = t13.get_tetra_opposite(v4)
        ta034 = t13.get_tetra_opposite(v2)
        ta123 = t22.get_tetra_opposite(v4)
        ta124 = t22.get_tetra_opposite(v3)
        ta134 = t22.get_tetra_opposite(v2)

        # Make sure the move is valid	
        if ta023.has_vertex(v1):
            return False
        if ta123.has_vertex(v0):
            return False
        if ta034.has_vertex(v1):
            return False
        if ta134.has_vertex(v0):
            return False
        if v0.check_vertex_neighbour(v1):
            return False

        # Create new tetrahedra
        tn13 = Tetrahedron()
        tn22l = Tetrahedron()
        tn22r = Tetrahedron()
        self.tetrahedron_pool.occupy(tn13)
        self.tetrahedron_pool.occupy(tn22l)
        self.tetrahedron_pool.occupy(tn22r)

        # Update connectivity
        tn13.set_vertices(v1, v0, v2, v4)
        tn22l.set_vertices(v1, v3, v0, v2)
        tn22r.set_vertices(v1, v3, v0, v4)
        tn13.set_tetras(t13.get_tetras()[0], ta124, tn22r, tn22l)
        tn22l.set_tetras(ta023, tn13, ta123, tn22r)
        tn22r.set_tetras(ta034, tn13, ta134, tn22l)

        time = t31.get_vertices()[0].time
        self.slab_sizes[time] += 1

        # Update tetrahedra connectivity
        t13.get_tetras()[0].exchange_tetra_opposite(t13.get_tetras()[0].get_vertices()[3], tn13)
        ta023.exchange_tetra_opposite(t13.get_vertex_opposite(v4), tn22l)
        ta034.exchange_tetra_opposite(t13.get_vertex_opposite(v2), tn22r)
        ta123.exchange_tetra_opposite(t22.get_vertex_opposite(v4), tn22l)
        ta124.exchange_tetra_opposite(t22.get_vertex_opposite(v3), tn13)
        ta134.exchange_tetra_opposite(t22.get_vertex_opposite(v2), tn22r)

        v0.cnum += 2
        v1.cnum += 2

        # Free the old tetrahedra
        self.tetrahedron_pool.free(t13.ID)
        self.tetrahedron_pool.free(t22.ID)

        # print(f"Shifted, created tn13: {tn13.ID}, tn22l: {tn22l.ID}, tn22r: {tn22r.ID}.")

        return True

    def ishift_d(self, tetra13_id: int, tetra22l_id: int, tetra22r_id: int) -> bool:
        """
        (3,2)-move in the triangulation. The ishift move removes a (2,2)-simplex
        from the triangulation by replacing its shared timelike edge with a dual
        timelike triangle. This results in the rearrangement of the three-simplices
        in an inverse manner to the shift move.

        Args:
            tetra13 (int): ID of the (1,3)-tetrahedron.
            tetra22l (int): ID of the (2,2)-tetrahedron.
            tetra22r (int): ID of the (2,2)-tetrahedron.

        Returns:
            bool: True if the move was ishifted successfully, otherwise False.
        """
        # Get the tetrahedra
        t13 = self.tetrahedron_pool.get(tetra13_id)
        t22l = self.tetrahedron_pool.get(tetra22l_id)
        t22r = self.tetrahedron_pool.get(tetra22r_id)
        t31 = t13.get_tetras()[0]
        if not t13.is_13():
            return False
        if not t22l.is_22():
            return False
        if not t22r.is_22():
            return False
        if not t13.check_neighbours_tetra(t22l):
            return False
        if not t13.check_neighbours_tetra(t22r):
            return False
        if not t31.is_31():
            return False
        
        # print(f"iShifting t13: {t13.ID} and t22l: {t22l.ID} and t22r: {t22r.ID}.")

        # Get the vertices of the inner triangle
        v1 = t13.get_vertices()[0]
        v3 = t22l.get_vertex_opposite_tetra(t13)
        v4 = t13.get_vertex_opposite_tetra(t22l)

        # Get the remaining vertices
        v4pos = t31.get_vertices().index(v4)
        v0 = t31.get_vertices()[(v4pos + 1) % 3]
        v2 = t31.get_vertices()[(v4pos + 2) % 3]

        # Get the neighbouring tetrahedra
        ta023 = t22l.get_tetra_opposite(v1)
        ta034 = t22r.get_tetra_opposite(v1)
        ta123 = t22l.get_tetra_opposite(v0)
        ta124 = t13.get_tetra_opposite(v0)
        ta134 = t22r.get_tetra_opposite(v0)

        # Make sure the move is valid
        if ta023.has_vertex(v4):
            return False
        if ta123.has_vertex(v4):
            return False
        if ta034.has_vertex(v2):
            return False
        if ta124.has_vertex(v3):
            return False
        if ta134.has_vertex(v2):
            return False

        # Create new tetrahedra
        tn13 = Tetrahedron()
        tn22 = Tetrahedron()
        self.tetrahedron_pool.occupy(tn13)
        self.tetrahedron_pool.occupy(tn22)

        # Update connectivity
        tn13.set_vertices(v3, v0, v2, v4)
        tn22.set_vertices(v1, v3, v2, v4)
        tn13.set_tetras(t13.get_tetras()[0], tn22, ta034, ta023)
        tn22.set_tetras(tn13, ta124, ta134, ta123)
    
        # Update tetrahedra connectivity
        t13.get_tetras()[0].exchange_tetra_opposite(t13.get_tetras()[0].get_vertices()[3], tn13)
        ta023.exchange_tetra_opposite(t22l.get_vertex_opposite(v1), tn13)
        ta034.exchange_tetra_opposite(t22r.get_vertex_opposite(v1), tn13)
        ta123.exchange_tetra_opposite(t22l.get_vertex_opposite(v0), tn22)
        ta124.exchange_tetra_opposite(t13.get_vertex_opposite(v0), tn22)
        ta134.exchange_tetra_opposite(t22r.get_vertex_opposite(v0), tn22)

        v0.cnum -= 2
        v1.cnum -= 2

        # Free the old tetrahedra
        self.tetrahedron_pool.free(t13.ID)
        self.tetrahedron_pool.free(t22l.ID)
        self.tetrahedron_pool.free(t22r.ID)

        time = tn13.get_vertices()[3].time
        self.slab_sizes[time] -= 1

        # print(f"iShifted, created tn13: {tn13.ID}, tn22: {tn22.ID}.")
   
        return True

    def update_vertices(self):
        """
        Update the vertices of the Universe.
        """
        # Determine the maximum vertex ID
        max_vertex = 0
        for v in self.vertex_pool.get_objects():
            if v.ID > max_vertex:
                max_vertex = v.ID

        # Clear the vertex neighbours
        self.vertex_neighbours.clear()
        vertex_neighbours = [[] for _ in range(max_vertex + 1)]

        # Update the vertex neighbours
        for v in self.vertex_pool.get_objects():
            nbr = []
            current_tetra = v.get_tetra()
            current_tetra_list = [current_tetra]
            next_tetra_list = []
            done_tetra_list = []

            # Iterate until all tetras are checked
            while current_tetra_list:
                # Get the current tetra to check and its neighbours
                for current_tetra_check in current_tetra_list:
                    for tetra_neighbour in current_tetra_check.get_tetras():
                        # IF the neighbour does not have the current vertex, continue
                        if not tetra_neighbour.has_vertex(v):
                            continue
                    
                    # Mark the current tetrahedron as done
                    if tetra_neighbour not in done_tetra_list:
                        done_tetra_list.append(tetra_neighbour)
                        next_tetra_list.append(tetra_neighbour)

                # Prepare current and next tetra lists for next iteration
                current_tetra_list = next_tetra_list
                next_tetra_list = []

            # Get the vertices of the done tetras
            for done_tetra in done_tetra_list:
                for vertex_done_tetra in done_tetra.get_vertices():
                    # Add the vertex to the neighbours list if it is not already there
                    if vertex_done_tetra not in nbr and vertex_done_tetra != v:
                        nbr.append(vertex_done_tetra)

            # Add the neighbours to the vertex neighbours list
            vertex_neighbours[v.ID] = nbr
        
        self.vertex_neighbours = vertex_neighbours

    def update_triangles(self):
        """
        Update the triangles of the Universe.
        """
        # Free all triangles in the pool
        self.triangle_pool.free_all()
        assert self.triangle_pool.get_number_occupied() == 0
        
        # Get tetras31 
        tetras31_objs = [self.tetrahedron_pool.get(i) for i in self.tetras_31.used_indices]

        for t in tetras31_objs:
            # Create a new triangle
            triangle = Triangle()
            self.triangle_pool.occupy(triangle)
            
            # Triangle is defined by base vertices of the tetrahedron
            triangle.set_vertices(t.get_vertices()[0], t.get_vertices()[1], t.get_vertices()[2])
            
            # Set the triangle of the tetrahedron
            for he in t.get_half_edges():
                he.set_triangle(triangle)

    def update_halfedges(self):
        """
        Update the halfedges of the Universe.
        """
        self.halfedge_pool.free_all()
        assert self.halfedge_pool.get_number_occupied() == 0

        # Get tetras31
        tetras31_objs = [self.tetrahedron_pool.get(i) for i in self.tetras_31.used_indices]

        # Create halfedges for each tetrahedron
        for t in tetras31_objs:
            halfedge_triples = []

            # Three halfedges for each tetrahedron
            for i in range(self.Constants.N_HALFEDGES_TETRA):
                he = HalfEdge()
                self.halfedge_pool.occupy(he)
                he.set_vertices(t.get_vertices()[i], t.get_vertices()[(i + 1) % 3])
                he.set_tetra(t)
                halfedge_triples.append(he)

            t.set_half_edges(halfedge_triples[0], halfedge_triples[1], halfedge_triples[2])

            # Connect the halfedges in a circular list
            for i in range(self.Constants.N_HALFEDGES_TETRA):
                halfedge_triples[i].set_next(halfedge_triples[(i + 1) % 3])
                halfedge_triples[i].set_previous(halfedge_triples[(i + 2) % 3])
        
        # # Set the adjacent opposite halfedges
        # for t in tetras31_objs:
        #     opposite = t.get_tetras()[0]
        #     for he in t.get_half_edges():
        #         for he_opposite in opposite.get_half_edges():
        #             if he.get_vertices() == he_opposite.get_vertices()[::-1] or he.get_vertices() == he_opposite.get_vertices():
        #                 he.set_adjacent(he_opposite)
                            
    def get_total_size(self) -> int:
        """
        Get the total size of the triangulation.

        Returns:
            int: Total size of the triangulation.
        """
        return self.vertex_pool.get_number_occupied()

    def update_geometry(self):
        """
        Update the geometry of the Universe.
        """
        self.update_vertices()
        self.update_halfedges()
        self.update_triangles()

    def has_duplicate_lists(self, list_of_lists: List[List[Any]]) -> bool:
        """
        Check if a list of lists has duplicates, where the order
        of the values in the inner lists don't matter.

        Args:
            list_of_lists (List[List[Any]]): list with lists to be checked

        Returns:
            bool: True if duplicate found, otherwise False
        """
        set_of_sets = set()

        for i, lst in enumerate(list_of_lists):
            set_lst = frozenset(lst)
            # Create a copy of the set without the current set
            temp_set_of_sets = set_of_sets.copy()
            # Remove the current set if present
            temp_set_of_sets.discard(set_lst)

            # Check if the current set has a duplicate
            if set_lst in temp_set_of_sets:
                print(f"i: {i}, duplicate list: {lst}")
                duplicate = self.tetrahedron_pool.get(i)
                duplicate.log()
                return True
            
            set_of_sets.add(set_lst)

        return False

    def check_validity(self):
        """
        Check the validity of the triangulation.
        """
        print("====================================================")
        print(f"Checking validity of the triangulation...")
        print("----------------------------------------------------")
        print(f"Number of tetrahedra: {self.tetrahedron_pool.get_number_occupied()}")

        # Check that each tetrahedron has 4 vertices and 4 neighbouring tetrahedra         
        for t in self.tetrahedron_pool.get_objects():
            assert len(t.get_vertices()) == self.Constants.N_VERTICES_TETRA, \
                f"Error: Tetrahedron {t.ID} has a wrong number of vertices."
            
            assert len(t.get_tetras()) == self.Constants.N_TETRA_NEIGHBOURS, \
                f"Error: Tetrahedron {t.ID} has a wrong number of neighbours. Vertices: {[v.ID for v in t.get_vertices()]}, Neighbours: {[n.ID for n in t.get_tetras()]}."

        for t in self.tetrahedron_pool.get_objects():
            # Make sure that each vertex is contained in the vertex pool
            for i in range(self.Constants.N_VERTICES_TETRA):
                assert self.vertex_pool.contains(t.get_vertices()[i].ID), \
                    f"Error: Tetrahedron {t.ID} has a missing vertex."

                # Make sure that vertices are unique
                for j in range(i + 1, self.Constants.N_VERTICES_TETRA):
                    assert t.get_vertices()[i] != t.get_vertices()[j], \
                        f"Error: Tetrahedron {t.ID} has duplicate vertices."

            # Check if all neighbours still exist
            for i in range(self.Constants.N_TETRA_NEIGHBOURS):             
                neighbour = t.get_tetras()[i] 

                # Make sure that each tetrahedron is contained in the tetrahedron pool
                if not self.tetrahedron_pool.contains(neighbour.ID):
                    print(f"Error: Tetrahedron {t.ID} has a missing neighbour.")
                    t.log()
                    neighbour.log()
                    exit()
 
                # Make sure that tetrahedron neighbours have eachother as neighbours
                assert neighbour.check_neighbours_tetra(t), f"Error: Tetrahedron neighbours {t.ID, neighbour.ID} are not mutual."
                # Make sure that tetrahedron neighbours are unique
                assert neighbour != t, f"Error: Tetrahedron {t.ID} has itself as neighbour."
                # ID of tetrahedron neighbours should be non-negative
                assert neighbour.ID >= 0, f"Error: Tetrahedron {t.ID} has negative neighbour ID."

                # Check for shared vertices between neighbours
                shared_vertices = sum(1 for vertex in neighbour.get_vertices() if t.has_vertex(vertex))
      
                # Between neighbours there should be at least 3 shared vertices
                assert shared_vertices >= 3, \
                    f"Error: Tetrahedron {t.ID} with type {t.type}has a wrong number of shared vertices with neighbour {t.get_tetras()[i].ID} with type {t.get_tetras()[i].type}.\n shared vertices = {shared_vertices}."

                # Check if the tetrahedra neighbours are of the correct type
                if t.is_31():
                    # If the tetrahedron is a (3,1)-simplex, then neighbour 3 should be a (1,3)-simplex
                    if i == self.Constants.APEX_INDEX_31:
                        assert neighbour.is_13(), \
                        f"Error: Tetrahedron {t.ID} has a wrong neighbour type, neigbour: {neighbour.ID}, type: {neighbour.type}."
                    else:
                        assert neighbour.is_22() or t.get_tetras()[i].is_31(), \
                        f"Error: Tetrahedron {t.ID} has a wrong neighbour type, neigbour: {neighbour.ID}, type: {neighbour.type}."
                if t.is_13():
                    # If the tetrahedron is a (1,3)-simplex, then neighbour 0 should be a (3,1)-simplex
                    if i == self.Constants.APEX_INDEX_13:
                        assert neighbour.is_31(), \
                        f"Error: Tetrahedron {t.ID} has a wrong neighbour type, neigbour: {neighbour.ID}, type: {neighbour.type}."
                    else:
                        assert neighbour.is_22() or t.get_tetras()[i].is_13(), \
                        f"Error: Tetrahedron {t.ID} has a wrong neighbour type, neigbour: {neighbour.ID}, type: {neighbour.type}."
            
            # Check if all opposite tetrahedra still exist and are correct
            for i in range(self.Constants.N_TETRA_NEIGHBOURS):
                assert t.get_tetra_opposite(t.get_vertices()[i]) == t.get_tetras()[i], \
                    f"Error: Vertex {neighbour.ID} in tetra {t.ID} has the wrong opposite tetra {t.get_tetras()[i].ID}. {t.get_tetra_opposite(neighbour.ID)} != {t.get_tetras()[i].ID}"
        
        # Check if each tetrahedron has a unique set of vertices
        tetra_vertices = [[v.ID for v in t.get_vertices()] for t in self.tetrahedron_pool.get_objects()]
        assert not self.has_duplicate_lists(tetra_vertices), "Error: Tetrahedra have duplicate vertices."
        
        # Check if vertex-tetrahedron connections are correct
        for v in self.vertex_pool.get_objects():
            assert self.tetrahedron_pool.contains(v.get_tetra().ID), \
                  f"Error: Vertex {v.ID} has a missing tetrahedron. {v.get_tetra().ID} not in tetrahedron pool."

            # Tadpole restriction
            if self.strictness == 1:
                assert v.scnum >= 2
            # Self-energy restriction
            if self.strictness == 2:
                assert v.scnum >= 3
        
        print("Valid! :)")
        print("====================================================")

    def export_geometry(self, geometry_outfilename: str = "output") -> bool:
        """
        Export the geometry of the Universe.
        """
        self.update_geometry()

        vertex_map = {}
        int_v_map = [None] * self.vertex_pool.get_number_occupied()

        for i, v in enumerate(self.vertex_pool.get_objects()):
            vertex_map[v] = i
            int_v_map[i] = v

        tetra_map = {}
        int_t_map = [None] * self.tetrahedron_pool.get_number_occupied()

        for i, t in enumerate(self.tetrahedron_pool.get_objects()):
            tetra_map[t] = i
            int_t_map[i] = t

        out = "1\n"
        out += f"{self.vertex_pool.get_number_occupied()}\n"

        for v in int_v_map:
            out += f"{v.time}\n"

        out += f"{self.vertex_pool.get_number_occupied()}\n"
        out += f"{self.tetrahedron_pool.get_number_occupied()}\n"

        for t in int_t_map:
            for v in t.get_vertices():
                out += f"{vertex_map[v]}\n"
            for tn in t.get_tetras():
                out += f"{tetra_map[tn]}\n"

        out += f"{self.tetrahedron_pool.get_number_occupied()}"

        with open(f"saved_universes/{geometry_outfilename}.txt", "w") as file:
            file.write(out + "\n")

        print(f"Geometry exported to {geometry_outfilename}.")

        return True


if __name__ == "__main__":
    u = Universe(geometry_infilename="initial_universes/sample-g0-T3.cdt")
    u.strictness = 3
    # u.update_geometry()
    # u.check_validity()

    # NOW HARDCODED CHECK BUT SEED CAN RESULT IN DIFFERENT OUTPUT
    # Check delete and add move
    random31 = u.tetras_31.pick()
    size_n3 = u.tetrahedron_pool.size
    size_n31 = u.tetras_31.size
    u.add(random31)
    print(f"Size n3: {size_n3}, Size n31: {size_n31}\nSize n3 after add: {u.tetrahedron_pool.size}, Size n31 after add: {u.tetras_31.size}")
    print(u.delete(15))
    print(f"Size n3 after delete: {u.tetrahedron_pool.size}, Size n31 after delete: {u.tetras_31.size}")

    # # Set of neighbours
    # t31 = u.tetrahedron_pool.get(u.tetras_31.pick())
    # other_t31 = t31.get_tetras()[1]
    # t22 = t31.get_tetras()[0]

    # # For the flip move
    # t012 = u.tetrahedron_pool.get(u.tetras_31.pick())
    # tetra230 = t012.get_tetras()[1]

    # # Check if the tetrahedron is actually flippable (opposite tetras should also be neighbours)
    # while not tetra230.is_31() or not t31.get_tetras()[3].check_neighbours_tetra(tetra230.get_tetras()[3]):
    #     tetra230 = u.tetrahedron_pool.get(u.tetras_31.pick())
    #     tetra230 = t012.get_tetras()[1]
        
    # # Check flip, shift_u and ishift_u move
    # u.flip(t012.ID, tetra230.ID)
    # u.shift_u(t31.ID, t22.ID)
    # u.ishift_u(360, 361, 362)

    # # Other set of neighbours
    # t31_2 = u.tetrahedron_pool.get(u.tetras_31.pick())
    # t13 = t31_2.get_tetras()[3]
    # t22_2 = t13.get_tetras()[3]

    # # Check shift_d and ishift_d move
    # u.shift_d(t13.ID, t22_2.ID)
    # u.ishift_d(360, 361, 362)

    # u.check_validity()
    # u.export_geometry("output_yeah")