# universe_mt.py
#
# Author: Seda den Boer
# Date: 15-05-2024

# Description: Defines the multiple-try version
# of the original Universe class.

from __future__ import annotations
from typing import Any, List, Dict
from vertex import Vertex
from tetra import Tetrahedron
from pool import Pool
from bag import Bag
import numpy as np
import os


class Universe:
    """
    The Universe class represents the current state of the triangulation
    and stores properties of the geometry in a convenient matter. It also
    provides member functions that carry out changes on the geometry.

    Args:
        geometry_infilename (str): Name of the file with the geometry.
        strictness (int): The strictness of the manifold conditions.
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
        self.strictness = strictness

        # Initialize the pools
        self.vertex_pool = Pool(Universe.Capacity.VERTEX)
        self.tetrahedron_pool = Pool(Universe.Capacity.TETRAHEDRON)

        # Initialize the bags
        self.tetras_31 = Bag(Universe.Capacity.TETRAHEDRON)
        self.tetras_22 = Bag(Universe.Capacity.TETRAHEDRON)

        self.vertex_neighbours = {}

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
            self.slab_sizes = np.zeros(self.n_slices)
            self.slice_sizes = np.zeros(self.n_slices)

            # Read the number of tetrahedra
            n3 = int(infile.readline()) 
            print(f"n3: {n3}")

            # Make n3 tetrahedra
            for _ in range(n3):
                tetra = Tetrahedron()
                self.tetrahedron_pool.occupy(tetra)

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
                elif tetra.is_22():
                    self.tetras_22.add(tetra.ID)

                # Clear memory
                del tetra_vs
                del tetra_ts
                
            # Check consistency of the number of tetrahedra
            line = int(infile.readline())
            if line != n3:
                print(f"Error: n3 = {n3}, line = {line}")
                return False

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
        for vertex in self.vertex_pool.get_objects():
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
        
        print(f"Read and loaded {geometry_infilename}.")

        return True

    def add(self, tetra31_id: int, perform: bool) -> bool:
        """
        (2,6)-move in the triangulation. Adds a vertex, two (3,1)-tetrahedra,
        and two (1,3)-tetrahedra.

        Args:
            tetra31_id (int): ID of the (3,1)-tetrahedron to perform add move in.
            perform (bool): Whether to perform the move or not.

        Returns:
            bool: True if the move was added successfully, otherwise False.
        """
        # Get the tetrahedron object, the time, and its opposite neighbor
        t = self.tetrahedron_pool.get(tetra31_id)
        time = t.get_vertices()[0].time
        tv = t.get_tetras()[3]
        # Double check if the tetrahedron is of the right type
        # if not t.is_31():
        #     return False
        # if not tv.is_13():
        #     return False

        if perform:
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

            # Clear references to avoid circular references and memory leaks
            t.clear_references()
            tv.clear_references()

            # Remove original tetrahedra from pool and bag
            self.tetrahedron_pool.free(t.ID)
            self.tetras_31.remove(t.ID)
            self.tetrahedron_pool.free(tv.ID)
            del t
            del tv

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

        return True

    def delete(self, vertex_id: int, perform: bool) -> bool:
        """
        (6,2)-move in the triangulation. Deletes a vertex, two (1,3)-tetrahedra,
        and two (3,1)-tetrahedra.

        Args:
            vertex (int): ID of the vertex to perform delete move on.
            perform (bool): Whether to perform the move or not.

        Returns:
            bool: True if the move was deleted successfully, otherwise False.
        """
         # Get the vertex object, its time and the two opposing tetrahedra it is part of
        vertex = self.vertex_pool.get(vertex_id)
        # if vertex.cnum != self.Constants.CNUM_ADD:
        #     return False
        # if vertex.scnum != self.Constants.SCNUM_ADD:
        #     return False
        time = vertex.time
        t01 = vertex.get_tetra()
        tv01 = t01.get_tetras()[3]

        # # Double check the tetrahedra are of the right type
        # if not t01.is_31():
        #     return False
        # if not tv01.is_13():
        #     return False

        # Get the vertex index in the tetrahedron
        vpos = np.where(t01.get_vertices() == vertex)[0][0]
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

        # if not t01.is_31():
        #     return False
        # if not t12.is_31():
        #     return False
        # if not t20.is_31():
        #     return False
        # if not tv01.is_13():
        #     return False
        # if not tv12.is_13():
        #     return False
        # if not tv20.is_13():
        #     return False
  
        # Get the neighbours of tetras that will be adjusted
        to01 = t01.get_tetra_opposite(vertex)
        to12 = t12.get_tetra_opposite(vertex)
        to20 = t20.get_tetra_opposite(vertex)
        tvo01 = tv01.get_tetra_opposite(vertex)
        tvo12 = tv12.get_tetra_opposite(vertex)
        tvo20 = tv20.get_tetra_opposite(vertex)

        # # Disallow tadpole insertions
        # if self.strictness == 1:
        #     if v0.scnum < 3:
        #         return False
        #     if v1.scnum < 3:
        #         return False
        #     if v2.scnum < 3:
        #         return False
        # # Disallow self-energy insertions
        # elif self.strictness >= 2:
        #     if v0.scnum < 4:
        #         return False
        #     if v1.scnum < 4:
        #         return False
        #     if v2.scnum < 4:
        #         return False
       
        if perform:
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
            
            # Clear references to avoid circular references and memory leaks
            vertex.clear_references()
            t01.clear_references()
            t12.clear_references()
            t20.clear_references()
            tv01.clear_references()
            tv12.clear_references()
            tv20.clear_references()

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
            del vertex
            del t01
            del t12
            del t20
            del tv01
            del tv12
            del tv20

            # Update the sizes
            self.slab_sizes[time] -= 2
            self.slab_sizes[(time - 1 + self.n_slices) % self.n_slices] -= 2
            self.slice_sizes[time] -= 2

        return True

    def flip(self, tetra012_id: int, tetra230_id: int, perform: bool) -> bool:
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
            perform (bool): Whether to perform the move or not.

        Returns:
            bool: True if the move was flipped successfully, otherwise False.
        """
        t012 = self.tetrahedron_pool.get(tetra012_id)
        t230 = self.tetrahedron_pool.get(tetra230_id)

        # Get the opposite (1,3)-tetrahedra
        tv012 = t012.get_tetras()[3]
        tv230 = t230.get_tetras()[3]
        
        # if not t012.is_31():
        #     return False
        # if not t230.is_31():
        #     return False
        # if not tv012.is_13():
        #     return False
        # if not tv230.is_13():
        #     return False
        # if not tv012.check_neighbours_tetra(tv230):
        #     return False
        
        # Get apex of the opposing tetrahedra
        vt = t012.get_vertices()[3]
        vb = tv012.get_vertices()[0]

        # Get the vertices of the base triangles that are going to be linked
        v1 = t012.get_vertex_opposite_tetra(t230)
        v3 = t230.get_vertex_opposite_tetra(t012)
        
        # Get the remaining base vertices
        v1pos = np.where(t012.get_vertices() == v1)[0][0]
        v2 = t012.get_vertices()[(v1pos + 1) % 3]
        v0 = t012.get_vertices()[(v1pos + 2) % 3]

        # # Manifold conditions
        # if self.strictness >= 1:
        #     if v1 == v3:
        #         return False
        # if self.strictness >= 2:
        #     if v0.scnum < 4:
        #         return False
        #     if v2.scnum < 4:
        #         return False
        # if self.strictness >= 3:
        #     if v1.check_vertex_neighbour(v3):
        #         return False

        # Get opposite neighbouring tetrahedra
        ta01 = t012.get_tetra_opposite(v2)
        ta12 = t012.get_tetra_opposite(v0)
        ta23 = t230.get_tetra_opposite(v0)
        ta30 = t230.get_tetra_opposite(v2)
        tva01 = tv012.get_tetra_opposite(v2)
        tva12 = tv012.get_tetra_opposite(v0)
        tva23 = tv230.get_tetra_opposite(v0)
        tva30 = tv230.get_tetra_opposite(v2)

        # # Make sure the move is valid
        # if ta01 == t230:
        #     return False
        # if ta23 == t012:
        #     return False
        # if tva01 == tv230:  
        #     return False
        # if tva23 == tv012:
        #     return False

        if perform:
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

        return True

    def shift_u(self, tetra31_id: int, tetra22_id: int, perform: bool) -> bool:
        """
        (2,3)-move in the triangulation. The shift move selects a (3,1)-simplex
        from the triangulation and replaces it with a neighboring (2,2)-simplex
        if conditions are met. This leads to a rearrangement of the local configuration
        without adding vertices.

        Args:
            tetra31 (int): ID of the (3,1)-tetrahedron.
            tetra22 (int): ID of the (2,2)-tetrahedron.
            perform (bool): Whether to perform the move or not.

        Returns:
            bool: True if the move was shifted successfully, otherwise False.
        """
        # Get the 31 and 22 tetra neighbours
        t31 = self.tetrahedron_pool.get(tetra31_id)
        t22 = self.tetrahedron_pool.get(tetra22_id)

        # Get the vertices that will be linked
        v0 = t31.get_vertex_opposite_tetra(t22)
        v1 = t22.get_vertex_opposite_tetra(t31)
        
        # The remaining vertices
        v3 = t31.get_vertices()[3]
        v0pos = np.where(t31.get_vertices() == v0)[0][0]
        v2 = t31.get_vertices()[(v0pos + 1) % 3]
        v4 = t31.get_vertices()[(v0pos + 2) % 3]

        # Get neighbouring tetrahedra that need to be updated after the move
        ta023 = t31.get_tetra_opposite(v4)
        ta034 = t31.get_tetra_opposite(v2)
        ta123 = t22.get_tetra_opposite(v4)
        ta124 = t22.get_tetra_opposite(v3)
        ta134 = t22.get_tetra_opposite(v2)

        # # Check if the move is valid
        # if ta023.has_vertex(v1):
        #     return False
        # if ta123.has_vertex(v0):
        #     return False
        # if ta034.has_vertex(v1):
        #     return False
        # if ta134.has_vertex(v0):
        #     return False
        # if v0.check_vertex_neighbour(v1):
        #     return False

        if perform:
            # Create new tetrahedra
            tn31 = Tetrahedron()
            tn22l = Tetrahedron()
            tn22r = Tetrahedron()

            # Add to pool and bag
            self.tetrahedron_pool.occupy(tn31)
            self.tetrahedron_pool.occupy(tn22l)
            self.tetrahedron_pool.occupy(tn22r)
            self.tetras_31.add(tn31.ID)
            self.tetras_22.add(tn22l.ID)
            self.tetras_22.add(tn22r.ID)

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

            # Clear references to avoid circular references and memory leaks
            t31.clear_references()
            t22.clear_references()

            # Remove the original tetrahedra from the neighbours
            self.tetrahedron_pool.free(t31.ID)
            self.tetrahedron_pool.free(t22.ID)
            self.tetras_31.remove(t31.ID)
            self.tetras_22.remove(t22.ID)
            del t31
            del t22

            # Make sure the vertices are updated with the new 31 tetra
            tn31.get_vertices()[0].set_tetra(tn31)
            tn31.get_vertices()[1].set_tetra(tn31)
            tn31.get_vertices()[2].set_tetra(tn31)

        return True
    
    def ishift_u(self, tetra31_id: int, tetra22l_id: int, tetra22r_id: int, perform: bool) -> bool:
        """
        (3,2)-move in the triangulation. The ishift move removes a (2,2)-simplex
        from the triangulation by replacing its shared timelike edge with a dual
        timelike triangle. This results in the rearrangement of the three-simplices
        in an inverse manner to the shift move.

        Args:
            tetra31 (int): ID of the (3,1)-tetrahedron.
            tetra22 (int): ID of the (2,2)-tetrahedron.
            tetra22r (int): ID of the (2,2)-tetrahedron.
            perform (bool): Whether to perform the move or not.

        Returns:
            bool: True if the move was ishifted successfully, otherwise False.
        """
        t31 = self.tetrahedron_pool.get(tetra31_id)
        t22l = self.tetrahedron_pool.get(tetra22l_id)
        t22r = self.tetrahedron_pool.get(tetra22r_id)
        
        # Get the vertices of the interior triangle
        v1 = t31.get_vertices()[3]
        v3 = t22l.get_vertex_opposite_tetra(t31)
        v4 = t31.get_vertex_opposite_tetra(t22l)

        # The remaining vertices
        v4pos = np.where(t31.get_vertices() == v4)[0][0]
        v0 = t31.get_vertices()[(v4pos + 1) % 3]
        v2 = t31.get_vertices()[(v4pos + 2) % 3]

        # Get neighbouring tetrahedra that need to be updated after the move
        ta023 = t22l.get_tetra_opposite(v1)
        ta034 = t22r.get_tetra_opposite(v1)
        ta123 = t22l.get_tetra_opposite(v0)
        ta124 = t31.get_tetra_opposite(v0)
        ta134 = t22r.get_tetra_opposite(v0)

        # # Make sure the move is valid
        # if ta023.has_vertex(v4):
        #     return False
        # if ta123.has_vertex(v4):
        #     return False
        # if ta034.has_vertex(v2):
        #     return False
        # if ta124.has_vertex(v3):
        #     return False
        # if ta134.has_vertex(v2):
        #     return False

        if perform:
            # Create new tetrahedra
            tn31 = Tetrahedron()
            tn22 = Tetrahedron()
            self.tetrahedron_pool.occupy(tn31)
            self.tetrahedron_pool.occupy(tn22)
            self.tetras_31.add(tn31.ID)
            self.tetras_22.add(tn22.ID)

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

            # # Clear references to avoid circular references and memory leaks
            t31.clear_references()
            t22l.clear_references()
            t22r.clear_references()

            # Free the old tetrahedra
            self.tetrahedron_pool.free(t31.ID)
            self.tetrahedron_pool.free(t22l.ID)
            self.tetrahedron_pool.free(t22r.ID)
            self.tetras_31.remove(t31.ID)
            self.tetras_22.remove(t22l.ID)
            self.tetras_22.remove(t22r.ID)
            del t31
            del t22l
            del t22r

            time = tn31.get_vertices()[0].time
            self.slab_sizes[time] -= 1

            # Make sure the base vertices are updated with the new 31 tetra
            tn31.get_vertices()[0].set_tetra(tn31)
            tn31.get_vertices()[1].set_tetra(tn31)
            tn31.get_vertices()[2].set_tetra(tn31)

        return True

    def shift_d(self, tetra13_id: int, tetra22_id: int, perform: bool) -> bool:
        """
        (2,3)-move in the triangulation. The shift move selects a (1,3)-simplex
        from the triangulation and replaces it with a neighboring (2,2)-simplex
        if conditions are met. This leads to a rearrangement of the local configuration
        without adding vertices.

        Args:
            tetra13 (int): ID of the (1,3)-tetrahedron.
            tetra22 (int): ID of the (2,2)-tetrahedron.
            perform (bool): Whether to perform the move or not.

        Returns:
            bool: True if the move was shifted successfully, otherwise False.
        """
        # Get the tetrahedra
        t13 = self.tetrahedron_pool.get(tetra13_id)
        t22 = self.tetrahedron_pool.get(tetra22_id)
        t31 = t13.get_tetras()[0]

        # Get the vertices that will be linked
        v0 = t13.get_vertex_opposite_tetra(t22)
        v1 = t22.get_vertex_opposite_tetra(t13)

        # The remaining vertices
        v3 = t13.get_vertices()[0] # Top
        v0pos = np.where(t31.get_vertices() == v0)[0][0]
        v2 = t31.get_vertices()[(v0pos + 1) % 3]
        v4 = t31.get_vertices()[(v0pos + 2) % 3]

        # Get the neighbouring tetrahedra
        ta023 = t13.get_tetra_opposite(v4)
        ta034 = t13.get_tetra_opposite(v2)
        ta123 = t22.get_tetra_opposite(v4)
        ta124 = t22.get_tetra_opposite(v3)
        ta134 = t22.get_tetra_opposite(v2)

        # # Make sure the move is valid	
        # if ta023.has_vertex(v1):
        #     return False
        # if ta123.has_vertex(v0):
        #     return False
        # if ta034.has_vertex(v1):
        #     return False
        # if ta134.has_vertex(v0):
        #     return False
        # if v0.check_vertex_neighbour(v1):
        #     return False
        
        if perform:
            # Create new tetrahedra
            tn13 = Tetrahedron()
            tn22l = Tetrahedron()
            tn22r = Tetrahedron()
            self.tetrahedron_pool.occupy(tn13)
            self.tetrahedron_pool.occupy(tn22l)
            self.tetrahedron_pool.occupy(tn22r)
            self.tetras_22.add(tn22l.ID)
            self.tetras_22.add(tn22r.ID)

            # Update connectivity
            tn13.set_vertices(v1, v0, v2, v4)
            tn22l.set_vertices(v1, v3, v0, v2)
            tn22r.set_vertices(v1, v3, v0, v4)
            tn13.set_tetras(t13.get_tetras()[0], ta124, tn22r, tn22l)
            tn22l.set_tetras(ta023, tn13, ta123, tn22r)
            tn22r.set_tetras(ta034, tn13, ta134, tn22l)

            time = t31.get_vertices()[0].time
            self.slab_sizes[time] += 1

            t13.get_tetras()[0].exchange_tetra_opposite(t13.get_tetras()[0].get_vertices()[3], tn13)
            ta023.exchange_tetra_opposite(t13.get_vertex_opposite(v4), tn22l)
            ta034.exchange_tetra_opposite(t13.get_vertex_opposite(v2), tn22r)
            ta123.exchange_tetra_opposite(t22.get_vertex_opposite(v4), tn22l)
            ta124.exchange_tetra_opposite(t22.get_vertex_opposite(v3), tn13)
            ta134.exchange_tetra_opposite(t22.get_vertex_opposite(v2), tn22r)

            v0.cnum += 2
            v1.cnum += 2

            # Clear references to avoid circular references and memory leaks
            t13.clear_references()
            t22.clear_references()

            # Free the old tetrahedra
            self.tetrahedron_pool.free(t13.ID)
            self.tetrahedron_pool.free(t22.ID)
            self.tetras_22.remove(t22.ID)
            del t13
            del t22

        return True

    def ishift_d(self, tetra13_id: int, tetra22l_id: int, tetra22r_id: int, perform: bool) -> bool:
        """
        (3,2)-move in the triangulation. The ishift move removes a (2,2)-simplex
        from the triangulation by replacing its shared timelike edge with a dual
        timelike triangle. This results in the rearrangement of the three-simplices
        in an inverse manner to the shift move.

        Args:
            tetra13 (int): ID of the (1,3)-tetrahedron.
            tetra22l (int): ID of the (2,2)-tetrahedron.
            tetra22r (int): ID of the (2,2)-tetrahedron.
            perform (bool): Whether to perform the move or not.

        Returns:
            bool: True if the move was ishifted successfully, otherwise False.
        """
        # Get the tetrahedra
        t13 = self.tetrahedron_pool.get(tetra13_id)
        t22l = self.tetrahedron_pool.get(tetra22l_id)
        t22r = self.tetrahedron_pool.get(tetra22r_id)
        t31 = t13.get_tetras()[0]

        # Get the vertices of the inner triangle
        v1 = t13.get_vertices()[0]
        v3 = t22l.get_vertex_opposite_tetra(t13)
        v4 = t13.get_vertex_opposite_tetra(t22l)

        # Get the remaining vertices
        v4pos = np.where(t31.get_vertices() == v4)[0][0]
        v0 = t31.get_vertices()[(v4pos + 1) % 3]
        v2 = t31.get_vertices()[(v4pos + 2) % 3]

        # Get the neighbouring tetrahedra
        ta023 = t22l.get_tetra_opposite(v1)
        ta034 = t22r.get_tetra_opposite(v1)
        ta123 = t22l.get_tetra_opposite(v0)
        ta124 = t13.get_tetra_opposite(v0)
        ta134 = t22r.get_tetra_opposite(v0)

        # # Make sure the move is valid
        # if ta023.has_vertex(v4):
        #     return False
        # if ta123.has_vertex(v4):
        #     return False
        # if ta034.has_vertex(v2):
        #     return False
        # if ta124.has_vertex(v3):
        #     return False
        # if ta134.has_vertex(v2):
        #     return False

        if perform:
            # Create new tetrahedra
            tn13 = Tetrahedron()
            tn22 = Tetrahedron()
            self.tetrahedron_pool.occupy(tn13)
            self.tetrahedron_pool.occupy(tn22)
            self.tetras_22.add(tn22.ID)

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

            # Clear references to avoid circular references and memory leaks
            t13.clear_references()
            t22l.clear_references()
            t22r.clear_references()
            
            # Free the old tetrahedra
            self.tetrahedron_pool.free(t13.ID)
            self.tetrahedron_pool.free(t22l.ID)
            self.tetrahedron_pool.free(t22r.ID)
            self.tetras_22.remove(t22l.ID)
            self.tetras_22.remove(t22r.ID)
            del t13
            del t22l
            del t22r

            time = tn13.get_vertices()[3].time
            self.slab_sizes[time] -= 1

        return True
    
    def update_vertices(self):
        """
        Update the vertices of the Universe.
        """
        # Clear the vertex neighbours
        self.vertex_neighbours.clear()

        # First make sets for each vertex 
        vertex_links = {}

        for vertex in self.vertex_pool.get_objects():
            vertex_links[vertex.ID] = set()

        # Fill the sets
        for tetra in self.tetrahedron_pool.get_objects():
            for vertex in tetra.get_vertices():
              
                # Save the spacelinks of the vertex
                for neighbour in tetra.get_vertices():
                        # Timelinks
                        vertex_links[vertex.ID].add(neighbour.ID)

        # Convert sets to lists
        self.vertex_neighbours = {vertex: list(neighbours) for vertex, neighbours in vertex_links.items()}
        del vertex_links

    def get_total_size(self) -> int:
        """
        Get the total size of the triangulation.

        Returns:
            int: Total size of the triangulation.
        """
        return self.vertex_pool.get_number_occupied()

    def get_curvature_profile(self) -> Dict[int, List[int]]:
        """
        Get the curvature profile of the Universe.

        Returns:
            Dict[int, List[int]]: Curvature profile of the Universe.
        """
        # Make a dict for each time slice with the vertices
        scnums_per_slice = {t: [] for t in range(self.n_slices)}
    
        for v in self.vertex_pool.get_objects():
            scnums_per_slice[v.time].append(v.scnum)
        
        return scnums_per_slice

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
                f"Error: Tetrahedron {t.ID} has a wrong number of vertices.\n{t.log()}\n{self.log()}"
            
            assert len(t.get_tetras()) == self.Constants.N_TETRA_NEIGHBOURS, \
                f"Error: Tetrahedron {t.ID} has a wrong number of neighbours. Vertices: {[v.ID for v in t.get_vertices()]}, Neighbours: {[n.ID for n in t.get_tetras()]}\n{self.log()}."

        for t in self.tetrahedron_pool.get_objects():
            # Make sure that each vertex is contained in the vertex pool
            for i in range(self.Constants.N_VERTICES_TETRA):
                assert self.vertex_pool.contains(t.get_vertices()[i].ID), \
                    f"Error: Tetrahedron {t.ID} has a missing vertex.\n{t.log()}\n{self.log()}"

                # Make sure that vertices are unique
                for j in range(i + 1, self.Constants.N_VERTICES_TETRA):
                    assert t.get_vertices()[i] != t.get_vertices()[j], \
                        f"Error: Tetrahedron {t.ID} has duplicate vertices.\n{t.log()}\n{self.log()}"

            # Check if all neighbours still exist
            for i in range(self.Constants.N_TETRA_NEIGHBOURS):             
                neighbour = t.get_tetras()[i] 

                # Make sure that each tetrahedron is contained in the tetrahedron pool
                if not self.tetrahedron_pool.contains(neighbour.ID):
                    print(f"Error: Tetrahedron {t.ID} has a missing neighbour.\n")
                    t.log()
                    neighbour.log()
                    print()
                    self.log()
                    exit()
 
                # Make sure that tetrahedron neighbours have eachother as neighbours
                assert neighbour.check_neighbours_tetra(t), f"Error: Tetrahedron neighbours {t.ID, neighbour.ID} are not mutual.\n{t.log()}\n{neighbour.log()}\n{self.log()}"
                # Make sure that tetrahedron neighbours are unique
                assert neighbour != t, f"Error: Tetrahedron {t.ID} has itself as neighbour.\n{t.log()}\n{self.log()}"
                # ID of tetrahedron neighbours should be non-negative
                assert neighbour.ID >= 0, f"Error: Tetrahedron {t.ID} has negative neighbour ID.\n{t.log()}\n{self.log()}"

                # Check for shared vertices between neighbours
                shared_vertices = sum(1 for vertex in neighbour.get_vertices() if t.has_vertex(vertex))
      
                # Between neighbours there should be at least 3 shared vertices
                assert shared_vertices == 3, \
                    f"Error: Tetrahedron {t.ID} with type {t.type} has a wrong number of shared vertices with neighbour {t.get_tetras()[i].ID} with type {t.get_tetras()[i].type}.\n shared vertices = {shared_vertices}.\n{t.log()}\n{self.log()}"

                # Check if the tetrahedra neighbours are of the correct type
                if t.is_31():
                    # If the tetrahedron is a (3,1)-simplex, then neighbour 3 should be a (1,3)-simplex
                    if i == self.Constants.APEX_INDEX_31:
                        assert neighbour.is_13(), \
                        f"Error: Tetrahedron {t.ID} has a wrong neighbour type, neigbour: {neighbour.ID}, type: {neighbour.type}.\n{t.log()}\n{neighbour.log()}\n{self.log()}"
                    else:
                        assert neighbour.is_22() or t.get_tetras()[i].is_31(), \
                        f"Error: Tetrahedron {t.ID} has a wrong neighbour type, neigbour: {neighbour.ID}, type: {neighbour.type}.\n{t.log()}\n{neighbour.log()}\n{self.log()}"
                if t.is_13():
                    # If the tetrahedron is a (1,3)-simplex, then neighbour 0 should be a (3,1)-simplex
                    if i == self.Constants.APEX_INDEX_13:
                        assert neighbour.is_31(), \
                        f"Error: Tetrahedron {t.ID} has a wrong neighbour type, neigbour: {neighbour.ID}, type: {neighbour.type}.\n{t.log()}\n{neighbour.log()}\n{self.log()}"
                    else:
                        assert neighbour.is_22() or t.get_tetras()[i].is_13(), \
                        f"Error: Tetrahedron {t.ID} has a wrong neighbour type, neigbour: {neighbour.ID}, type: {neighbour.type}.\n{t.log()}\n{neighbour.log()}\n{self.log()}"
            
            # Check if all opposite tetrahedra still exist and are correct
            for i in range(self.Constants.N_TETRA_NEIGHBOURS):
                assert t.get_tetra_opposite(t.get_vertices()[i]) == t.get_tetras()[i], \
                    f"Error: Vertex {neighbour.ID} in tetra {t.ID} has the wrong opposite tetra {t.get_tetras()[i].ID}. {t.get_tetra_opposite(neighbour.ID)} != {t.get_tetras()[i].ID}. \n{t.log()}\n{self.log()}"
        
        # Check if each tetrahedron has a unique set of vertices
        tetra_vertices = [[v.ID for v in t.get_vertices()] for t in self.tetrahedron_pool.get_objects()]
        assert not self.has_duplicate_lists(tetra_vertices), f"Error: Tetrahedra have duplicate vertices.\n{self.log()}"
        
        # Check if vertex-tetrahedron connections are correct
        for v in self.vertex_pool.get_objects():
            assert self.tetrahedron_pool.contains(v.get_tetra().ID), \
                  f"Error: Vertex {v.ID} has a missing tetrahedron. {v.get_tetra().ID} not in tetrahedron pool.\n{v.get_tetra().log()}\n{self.log()}"

            # Tadpole restriction
            if self.strictness == 1:
                assert v.scnum >= 2
            # Self-energy restriction
            if self.strictness == 2:
                assert v.scnum >= 3
        
        # Check correctness vertex_neighbours
        for v, neighbours in self.vertex_neighbours.items():
            for n in neighbours:
                # Make sure the neighbour is in the vertex pool
                assert self.vertex_pool.contains(n), f"Error: Vertex {v} has a missing neighbour. {n} not in vertex pool.\n{self.log()}"
                # Make sure the neighbour also has the vertex as neighbour
                assert self.vertex_neighbours[n].count(v) == 1, f"Error: Vertex {v} is not a neighbour of vertex {n}.\n{self.log()}"
        
        print("Valid! :)")
        print("====================================================")

    def export_geometry(self, geometry_outfilename: str = "output", k0: float = -1) -> bool:
        """
        Export the geometry of the Universe.

        Args:
            geometry_outfilename (str): Name of the output file.
            as_pickle (bool): Save as pickle or text file.
            k0 (float): k0 value for the Universe.
        """
        if k0 >= 0:
            pathname = f"saved_universes/k0={k0}/{geometry_outfilename}"
        else:
            pathname = f"saved_universes/{geometry_outfilename}"

        if not os.path.exists(os.path.dirname(pathname)):
            os.makedirs(os.path.dirname(pathname))

        # Save as text
        vertex_map = {}
        int_v_map = [None] * len(self.vertex_pool.get_objects())

        i = 0
        for vertex in self.vertex_pool.get_objects():
            vertex_map[vertex] = i
            int_v_map[i] = vertex
            i += 1

        tetra_map = {}
        int_t_map = [None] * self.tetrahedron_pool.get_number_occupied()

        i = 0
        for tetra in self.tetrahedron_pool.get_objects():
            tetra_map[tetra] = i
            int_t_map[i] = tetra
            i += 1

        out = "1\n"  # indicating well-orderedness

        out += str(self.vertex_pool.get_number_occupied()) + "\n"

        for j in range(len(int_v_map)):
            out += str(int_v_map[j].time) + "\n"

        out += str(self.vertex_pool.get_number_occupied()) + "\n"

        out += str(self.tetrahedron_pool.get_number_occupied()) + "\n"

        for j in range(len(int_t_map)):
            for v in int_t_map[j].get_vertices():
                out += str(vertex_map[v]) + "\n"
            for t in int_t_map[j].get_tetras():
                out += str(tetra_map[t]) + "\n"

        out += str(self.tetrahedron_pool.get_number_occupied())

        with open(pathname + ".txt", "w") as file:
            file.write(out + "\n")

        print(f"Geometry exported to {pathname}.txt.")

        return True

    def log(self):
        """
        Log the Universe.
        """
        self.update_vertices()
        print("====================================================")
        print(f"Universe with {self.get_total_size()} vertices and {self.tetrahedron_pool.get_number_occupied()} tetrahedra.")
        print("----------------------------------------------------")
        print(f"Vertices:")
        for vertex in self.vertex_pool.get_objects():
            print(vertex.log())
            print(f"Vertex neighbours: {[nbr for nbr in self.vertex_neighbours[vertex.ID]]}")
        print("----------------------------------------------------")
        print(f"Tetrahedra:")
        for tetra in self.tetrahedron_pool.get_objects():
            tetra.log()
        print("====================================================")


if __name__ == "__main__":
    # u = Universe(geometry_infilename="initial_universes/sample-g0-T3.cdt")
    u = Universe(geometry_infilename="initial_universes/output_g=0_T=200.txt")
    # u.update_vertices()
    # u.check_validity()
    # u.log()