# universe.py
#
# Author: Seda den Boer
# Date: 17-02-2024	
# 
# Description: The Universe class that represents the current state of the triangulation.

from __future__ import annotations
from typing import cast
from vertex import Vertex
from halfedge import HalfEdge
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

    class Capacity:
        VERTEX = 1000000
        HALFEDGE = 2 * VERTEX
        TRIANGLE = 2 * VERTEX
        TETRAHEDRON = 2 * VERTEX
    
    class Constants:
        CNUM_ADD = 6
        SCNUM_ADD = 3

    def __init__(self, geometry_infilename: str, strictness: int = 0, volfix_switch: int = 0):
        self.strictness = strictness
        self.volfix_switch = volfix_switch

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

                # Update bag
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

    def add_move(self, tetra31_id: int) -> bool:
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
        time = t.get_vertices()[0].time
        assert t.is_31()
        tv = t.get_tetras()[3]
        assert tv.is_13()

        # Create a new vertex and update its cnum and scnum
        vertex = Vertex(time=time)
        vertex.ID = self.vertex_pool.occupy(vertex)
        vertex.cnum = Universe.Constants.CNUM_ADD
        vertex.scnum = Universe.Constants.SCNUM_ADD

        # Get relevant vertices of the two tetrahedra
        v0, v1, v2, vt = t.get_vertices()
        vb = tv.get_vertices()[0]

        # Create the new tetrahedra
        tn01 = Tetrahedron()
        tn12 = Tetrahedron()
        tn20 = Tetrahedron()
        tvn01 = Tetrahedron()
        tvn12 = Tetrahedron()
        tvn20 = Tetrahedron()

        # Add to pool and bag
        for t in [tn01, tn12, tn20, tvn01, tvn12, tvn20]:
            self.tetrahedron_pool.occupy(t)
        for t in [tn01, tn12, tn20]:
            self.tetras_31.add(t.ID)

        # Get neighbouring tetrahedra of the two oppsoing tetrahedra
        to0 = t.get_tetra_opposite(v0)
        to1 = t.get_tetra_opposite(v1)
        to2 = t.get_tetra_opposite(v2)
        tvo0 = tv.get_tetra_opposite(v0)
        tvo1 = tv.get_tetra_opposite(v1)
        tvo2 = tv.get_tetra_opposite(v2)

        # Set the vertices of the new tetrahedra
        tn01.set_vertices(v0, v1, vertex, vt)
        tn12.set_vertices(v1, v2, vertex, vt)
        tn20.set_vertices(v2, v0, vertex, vt)
        tvn01.set_vertices(vb, v0, v1, vertex)
        tvn12.set_vertices(vb, v1, v2, vertex)
        tvn20.set_vertices(vb, v2, v0, vertex)

        # Set neighbouring tetrahedra
        tn01.set_tetras(tn12, tn20, to2, tvn01)
        tn12.set_tetras(tn20, tn01, to0, tvn12)
        tn20.set_tetras(tn01, tn12, to1, tvn20)
        
        # Exchange given tetrahedron for the given vertex with opposite tetrahedron
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
        self.tetrahedron_pool.free(tv.ID)
        self.tetras_31.remove(t.ID)

        # Update the vertices with the tetrahedra they are part of, and their cnum and scnum
        vertex.set_tetra(tn01)
        v0.set_tetra(tn01)
        v1.set_tetra(tn12)
        v2.set_tetra(tn20)
        v0.scnum += 1; v1.scnum += 1; v2.scnum += 1
        v0.cnum += 2; v1.cnum += 2; v2.cnum += 2; vt.cnum += 2; vb.cnum += 2

        return True

    def delete_move(self, vertex_id: int) -> bool:
        """
        (6,2)-move in the triangulation. Deletes a vertex, two (1,3)-tetrahedra,
        and two (3,1)-tetrahedra.

        Args:
            vertex (int): ID of the vertex to perform delete move in.

        Returns:
            bool: True if the move was deleted successfully, otherwise False.
        """
        # Get the vertex object, its time and the two opposing tetrahedra it is part of
        vertex = self.vertex_pool.get(vertex_id)
        assert vertex.scnum == Universe.Constants.SCNUM_ADD
        time = vertex.time
        t01 = vertex.get_tetra()
        tv01 = t01.get_tetras()[3]

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
        assert t01.is_31() and t12.is_31() and t20.is_31() and tv01.is_13() and tv12.is_13() and tv20.is_13()
        to01 = t01.get_tetra_opposite(vertex)
        to12 = t12.get_tetra_opposite(vertex)
        to20 = t20.get_tetra_opposite(vertex)
        tvo01 = tv01.get_tetra_opposite(vertex)
        tvo12 = tv12.get_tetra_opposite(vertex)
        tvo20 = tv20.get_tetra_opposite(vertex)

        # Disallow tadpole insertions
        if self.strictness == 1 and (v0.scnum < 3 or v1.scnum < 3 or v2.scnum < 3):
            return False
        # Disallow self-energy insertions
        elif self.strictness >= 2 and (v0.scnum < 4 or v1.scnum < 4 or v2.scnum < 4):
            return False
        
        # Create new tetrahedra
        tn = Tetrahedron()
        tvn = Tetrahedron()
        self.tetrahedron_pool.occupy(tn)
        self.tetrahedron_pool.occupy(tvn)
        self.tetras_31.add(tn.ID)

        # Get apex and base vertex
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
        v0.scnum -= 1; v1.scnum -= 1; v2.scnum -= 1
        v0.cnum -= 2; v1.cnum -= 2; v2.cnum -= 2; vt.cnum -= 2; vb.cnum -= 2

        # Exchange given tetrahedron for the given vertex with opposite tetrahedron
        to01.exchange_tetra_opposite(t01.get_vertex_opposite(vertex), tn)
        to12.exchange_tetra_opposite(t12.get_vertex_opposite(vertex), tn)
        to20.exchange_tetra_opposite(t20.get_vertex_opposite(vertex), tn)
        tvo01.exchange_tetra_opposite(tv01.get_vertex_opposite(vertex), tvn)
        tvo12.exchange_tetra_opposite(tv12.get_vertex_opposite(vertex), tvn)
        tvo20.exchange_tetra_opposite(tv20.get_vertex_opposite(vertex), tvn)

        # Update pool and bag
        for t in [t01, t12, t20, tv01, tv12, tv20]:
            self.tetrahedron_pool.free(t.ID)
        for t in [t01, t12, t20]:
            self.tetras_31.remove(t.ID)
        self.vertex_pool.free(vertex.ID)

        # Update the sizes
        self.slab_sizes[time] -= 2
        self.slab_sizes[(time - 1 + self.n_slices) % self.n_slices] -= 2
        self.slice_sizes[time] -= 2

        return True

    def flip_move(self, tetra012_id: int, tetra230_id: int) -> bool:
        """
        (4,4)-move in the triangulation. Flips shared spatial edge base triangles 
        of two tetrahedra.

        Args:
            tetra012 (int): ID of the first tetrahedron.
            tetra230 (int): ID of the second tetrahedron.

        Returns:
            bool: True if the move was flipped successfully, otherwise False.
        """
        t012 = self.tetrahedron_pool.get(tetra012_id)
        t230 = self.tetrahedron_pool.get(tetra230_id)

        # Get the vertices of the base triangles
        v1 = t012.get_vertex_opposite_tetra(t230)
        v3 = t230.get_vertex_opposite_tetra(t012)
        for i in range(3):
            if t012.get_vertices()[i] == v1:
                v2 = t012.get_vertices()[(i + 1) % 3]
                v0 = t012.get_vertices()[(i + 2) % 3]
                break
        
        # Get the opposite tetrahedra
        tv012 = t012.get_tetras()[3]
        tv230 = t230.get_tetras()[3]

        if self.strictness >= 1 and v1 == v3:
            return False
        # Disallow self-energt insertions in dual graph
        if self.strictness >= 2 and (v0.scnum == 3 or v2.scnum == 3):
            return False
        if self.strictness >= 3 and v1.check_vertex_neighbour(v3):
            return False

        # Get apex and base vertex
        vt = t012.get_vertices()[3]
        vb = tv012.get_vertices()[0]

        # Get opposite neighbouring tetrahedra for initial connected vertices
        ta01 = t012.get_tetra_opposite(v2)
        ta12 = t012.get_tetra_opposite(v0)
        ta23 = t230.get_tetra_opposite(v0)
        ta30 = t230.get_tetra_opposite(v2)
        tva01 = tv012.get_tetra_opposite(v2)
        tva12 = tv012.get_tetra_opposite(v0)
        tva23 = tv230.get_tetra_opposite(v0)
        tva30 = tv230.get_tetra_opposite(v2)
        # Make sure the move is valid
        if ta01 == t230 or ta23 == t012 or tva01 == t230 or tva23 == t012:
            return False
        t012vo2 = t012.get_vertex_opposite(v2)
        t230vo0 = t230.get_vertex_opposite(v0)
        tv012vo2 = tv012.get_vertex_opposite(v2)
        tv230vo0 = tv230.get_vertex_opposite(v0)

        # Create new tetrahedra
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

        # Exchange given tetrahedron for the given vertex with opposite tetrahedron
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

    def shift_u_move(self, tetra31_id: int, tetra22_id: int) -> bool:
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
        t31 = self.tetrahedron_pool.get(tetra31_id)
        t22 = self.tetrahedron_pool.get(tetra22_id)

        v0 = t31.get_vertex_opposite_tetra(t22)
        v1 = t22.get_vertex_opposite_tetra(t31)
        
        v0pos = t31.get_vertices().index(v0)
        assert v0pos >= 0

        v2 = t31.get_vertices()[(v0pos + 1) % 3]
        v4 = t31.get_vertices()[(v0pos + 2) % 3]
        v3 = t31.get_vertices()[3]

        ta023 = t31.get_tetra_opposite(v4)
        ta034 = t31.get_tetra_opposite(v2)
        ta123 = t22.get_tetra_opposite(v4)
        ta124 = t22.get_tetra_opposite(v3)
        ta134 = t22.get_tetra_opposite(v2)

        if ta023.has_vertex(v1) or ta123.has_vertex(v0) or ta034.has_vertex(v1) or ta134.has_vertex(v0):
            return False
        if v0.check_vertex_neighbour(v1):
            return False

        tn31 = Tetrahedron()
        tn22l = Tetrahedron()
        tn22r = Tetrahedron()
        for t in [tn31, tn22l, tn22r]:
            self.tetrahedron_pool.occupy(t)
        self.tetras_31.add(tn31.ID)

        tn31.set_vertices(v0, v2, v4, v1)
        tn22l.set_vertices(v0, v2, v1, v3)
        tn22r.set_vertices(v0, v4, v1, v3)

        tn31.set_tetras(ta124, tn22r, tn22l, t31.tnbr[3])
        tn22l.set_tetras(ta123, tn22r, ta023, tn31)
        tn22r.set_tetras(ta134, tn22l, ta034, tn31)

        time = tn31.get_vertices()[0].time
        self.slab_sizes[time] += 1

        t31.get_tetras()[3].exchange_tetra_opposite(
            t31.get_tetras()[3].get_vertices()[0], tn31
        )

        ta023.exchange_tetra_opposite(t31.get_vertex_opposite(v4), tn22l)
        ta034.exchange_tetra_opposite(t31.get_vertex_opposite(v2), tn22r)
        ta123.exchange_tetra_opposite(t22.get_vertex_opposite(v4), tn22l)
        ta124.exchange_tetra_opposite(t22.get_vertex_opposite(v3), tn31)
        ta134.exchange_tetra_opposite(t22.get_vertex_opposite(v2), tn22r)

        v0.cnum += 2
        v1.cnum += 2

        self.tetrahedron_pool.free(t31.ID)
        self.tetrahedron_pool.free(t22.ID)
        self.tetras_31.remove(t31.ID)

        for v in tn31.get_vertices():
            v.set_tetra(tn31)

        return True
    
    def ishift_u_move(self, tetra31_id: int, tetra22l_id: int, tetra22r_id: int) -> bool:
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

        v1 = t31.get_vertices()[3]
        v3 = t22l.get_vertex_opposite_tetra(t31)
        v4 = t31.get_vertex_opposite_tetra(t22l)

        v4pos = t31.get_vertices().index(v4)
        assert v4pos >= 0

        v0 = t31.get_vertices()[(v4pos + 1) % 3]
        v2 = t31.get_vertices()[(v4pos + 2) % 3]

        ta023 = t22l.get_tetra_opposite(v1)
        ta034 = t22r.get_tetra_opposite(v1)
        ta123 = t22l.get_tetra_opposite(v0)
        ta124 = t31.get_tetra_opposite(v0)
        ta134 = t22r.get_tetra_opposite(v0)

        if ta023.has_vertex(v4) or ta123.has_vertex(v4) or ta034.has_vertex(v2) or ta124.has_vertex(v3) or ta134.has_vertex(v2):
            return False

        tn31 = Tetrahedron()
        tn22 = Tetrahedron()
        self.tetrahedron_pool.occupy(tn31)
        self.tetrahedron_pool.occupy(tn22)
        self.tetras_31.add(tn31.ID)

        tn31.set_vertices(v0, v2, v4, v3)
        tn22.set_vertices(v2, v4, v1, v3)
        tn31.set_tetras(tn22, ta034, ta023, t31.get_tetras()[3])
        tn22.set_tetras(ta134, ta123, tn31, ta124)

        t31.get_tetras()[3].exchange_tetra_opposite(
            t31.get_tetras()[3].get_vertices()[0], tn31
        )
        ta023.exchange_tetra_opposite(t22l.get_vertex_opposite(v1), tn31)
        ta034.exchange_tetra_opposite(t22r.get_vertex_opposite(v1), tn31)
        ta123.exchange_tetra_opposite(t22l.get_vertex_opposite(v0), tn22)
        ta124.exchange_tetra_opposite(t31.get_vertex_opposite(v0), tn22)
        ta134.exchange_tetra_opposite(t22r.get_vertex_opposite(v0), tn22)

        v0.cnum -= 2
        v1.cnum -= 2

        self.tetrahedron_pool.free(t31.ID)
        self.tetrahedron_pool.free(t22l.ID)
        self.tetrahedron_pool.free(t22r.ID)
        self.tetras_31.remove(t31.ID)

        time = tn31.get_vertices()[0].time
        self.slab_sizes[time] -= 1

        for v in tn31.get_vertices():
            v.set_tetra(tn31)

        return True

    def shift_d_move(self, tetra13_id: int, tetra22_id: int) -> bool:
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
        t13 = self.tetrahedron_pool.get(tetra13_id)
        t22 = self.tetrahedron_pool.get(tetra22_id)
        # Only because spatial ordering of (1,3)-tetrahedra is currently not guaranteed
        t31 = t13.get_tetras()[0]

        v0 = t13.get_vertex_opposite_tetra(t22)
        v1 = t22.get_vertex_opposite_tetra(t13)

        v0pos = t13.get_vertices().index(v0)
        assert v0pos >= 0

        v2 = t31.vs[(v0pos + 1) % 3]
        v4 = t31.vs[(v0pos + 2) % 3]
        v3 = t13.vs[0]

        ta023 = t13.get_tetra_opposite(v4)
        ta034 = t13.get_tetra_opposite(v2)
        ta123 = t22.get_tetra_opposite(v4)
        ta124 = t22.get_tetra_opposite(v3)
        ta134 = t22.get_tetra_opposite(v2)

        # Make sure the move is valid	
        if ta023.has_vertex(v1) or ta123.has_vertex(v0) or ta034.has_vertex(v1) or ta134.has_vertex(v0):
            return False
        if v0.check_vertex_neighbour(v1):
            return False

        tn13 = Tetrahedron()
        tn22l = Tetrahedron()
        tn22r = Tetrahedron()
        for t in [tn13, tn22l, tn22r]:
            self.tetrahedron_pool.occupy(t)

        tn13.set_vertices(v1, v0, v2, v4)
        tn22l.set_vertices(v1, v3, v0, v2)
        tn22r.set_vertices(v1, v3, v0, v4)
        tn13.set_tetras(t31.get_tetras()[0], ta124, tn22r, tn22l)
        tn22l.set_tetras(ta023, tn13, ta123, tn22r)
        tn22r.set_tetras(ta034, tn13, ta134, tn22l)

        time = t31.get_vertices()[0].time
        self.slab_sizes[time] += 1

        t13.get_tetras()[3].exchange_tetra_opposite(
            t13.get_tetras()[0].get_vertices()[3], tn13
        )

        ta023.exchange_tetra_opposite(t13.get_vertex_opposite(v4), tn22l)
        ta034.exchange_tetra_opposite(t13.get_vertex_opposite(v2), tn22r)
        ta123.exchange_tetra_opposite(t22.get_vertex_opposite(v4), tn22l)
        ta124.exchange_tetra_opposite(t22.get_vertex_opposite(v3), tn13)
        ta134.exchange_tetra_opposite(t22.get_vertex_opposite(v2), tn22r)

        v0.cnum += 2
        v1.cnum += 2

        self.tetrahedron_pool.remove(t13)
        self.tetrahedron_pool.remove(t22)

        return True

    def ishift_d_move(self, tetra13_id: int, tetra22_id: int, tetra22r_id: int) -> bool:
        """
        (3,2)-move in the triangulation. The ishift move removes a (2,2)-simplex
        from the triangulation by replacing its shared timelike edge with a dual
        timelike triangle. This results in the rearrangement of the three-simplices
        in an inverse manner to the shift move.

        Args:
            tetra13 (int): ID of the (1,3)-tetrahedron.
            tetra22 (int): ID of the (2,2)-tetrahedron.
            tetra22r (int): ID of the (2,2)-tetrahedron.

        Returns:
            bool: True if the move was ishifted successfully, otherwise False.
        """
        t13 = self.tetrahedron_pool.get(tetra13_id)
        t22l = self.tetrahedron_pool.get(tetra22_id)
        t22r = self.tetrahedron_pool.get(tetra22r_id)
        t31 = t13.get_tetras()[0]
        
        v1 = t13.get_vertices()[0]
        v3 = t22l.get_vertex_opposite_tetra(t13)
        v4 = t13.get_vertex_opposite_tetra(t22l)

        v4pos = t13.get_vertices().index(v4)
        assert v4pos >= 0

        v0 = t31.get_vertices()[(v4pos + 1) % 3]
        v2 = t31.get_vertices()[(v4pos + 2) % 3]

        ta023 = t22l.get_tetra_opposite(v1)
        ta034 = t22r.get_tetra_opposite(v1)
        ta123 = t22l.get_tetra_opposite(v0)
        ta124 = t13.get_tetra_opposite(v0)
        ta134 = t22r.get_tetra_opposite(v0)

        if ta023.has_vertex(v4) or ta123.has_vertex(v4) or ta034.has_vertex(v2) or ta124.has_vertex(v3) or ta134.has_vertex(v2):
            return False

        tn13 = Tetrahedron()
        tn22 = Tetrahedron()
        self.tetrahedron_pool.occupy(tn13)
        self.tetrahedron_pool.occupy(tn22)

        tn13.set_vertices(v3, v0, v2, v4)
        tn22.set_vertices(v1, v3, v2, v4)
        tn13.set_tetras(t13.get_tetras()[0], tn22, ta034, ta023)
        tn22.set_tetras(tn13, ta124, ta134, ta123)

        t13.get_tetras()[0].exchange_tetra_opposite(
            t13.get_tetras()[0].get_vertices()[3], tn13
        )
        ta023.exchange_tetra_opposite(t22l.get_vertex_opposite(v1), tn13)
        ta034.exchange_tetra_opposite(t22r.get_vertex_opposite(v1), tn13)
        ta123.exchange_tetra_opposite(t22l.get_vertex_opposite(v0), tn22)
        ta124.exchange_tetra_opposite(t13.get_vertex_opposite(v0), tn22)
        ta134.exchange_tetra_opposite(t22r.get_vertex_opposite(v0), tn22)

        v0.cnum -= 2
        v1.cnum -= 2

        self.tetrahedron_pool.remove(t13)
        self.tetrahedron_pool.remove(t22l)
        self.tetrahedron_pool.remove(t22r)

        time = tn13.get_vertices()[3].time
        self.slab_sizes[time] -= 1

        return True

    def update_vertices(self):
        """
        Update the vertices of the Universe.
        """
        self.vertex_neighbours.clear()
        self.vertex_neighbours = [None] * self.vertex_pool.get_number_occupied()

        for v in self.vertex_pool.get_objects():
            nbr = []
            t = v
            current = [t]
            next = []
            done = []

            while current:
                for tc in current:
                    for tcn in tc.get_tetras():
                        if not tcn.has_vertex(v):
                            continue
                        if tcn not in done:
                            next.append(tcn)
                            done.append(tcn)
                
                current = next
                next = []
            
            for td in done:
                for vd in td.get_vertices():
                    if vd not in nbr and vd != v:
                        nbr.append(vd)
            
            self.vertex_neighbours[v.ID] = nbr

    def update_triangles(self):
        """
        Update the triangles of the Universe.
        """
        # Free all triangles in the pool
        self.triangle_pool.free_all()
        
        # Get tetras31 
        tetras31_objs = [self.tetrahedron_pool.get(i) for i in self.tetras_31.used_indices]

        for t in tetras31_objs:
            # Create a new triangle
            triangle = Triangle()
            self.triangle_pool.occupy(triangle)
            
            triangle.set_vertices(t.get_vertices()[0], t.get_vertices()[1], t.get_vertices()[2])
            triangle.set_half_edges(t.get_half_edges()[0], t.get_half_edges()[1], t.get_half_edges()[2])

            # Set the triangle of the tetrahedron
            for he in t.get_half_edges():
                he.set_triangle(triangle)
        
        self.triangle_neighbours = [None] * self.triangle_pool.get_number_occupied()

        for triangle in self.triangle_pool.get_objects():
            triangle.set_triangle_neighbours(triangle.get_half_edges()[0].get_adjacent().get_triangle(),
                                             triangle.get_half_edges()[1].get_adjacent().get_triangle(),
                                             triangle.get_half_edges()[2].get_adjacent().get_triangle())

            self.triangle_neighbours[triangle.ID] = triangle.get_triangle_neighbours()

    def update_halfedges(self):
        """
        Update the halfedges of the Universe.
        """
        self.halfedge_pool.free_all()
        assert self.halfedge_pool.get_number_occupied() == 0

        # Get tetras31
        tetras31_objs = [self.tetrahedron_pool.get(i) for i in self.tetras_31.used_indices]

        for t in tetras31_objs:
            these = []
            for i in range(3):
                he = HalfEdge()
                self.halfedge_pool.occupy(he)
                he.set_vertices(t.get_vertices()[i], t.get_vertices()[(i + 1) % 3])
                he.set_tetra(t)
                these.append(he)
            
            t.set_half_edges(these[0], these[1], these[2])

            for i in range(3):
                these[i].set_next(these[(i + 1) % 3])
                these[i].set_prev(these[(i + 2) % 3])
        
        for t in tetras31_objs:
            for i in range(3):
                v = t.get_vertices()[i]
                vt = t.get_vertices()[3]
                tc = t.get_tetra_opposite(v)
                tn = Tetrahedron()
                v = vt
                while tc.is_22():
                    tn = tc.get_tetra_opposite(v)

                    vo = v
                    v = tc.get_vertices()[2] if tc.get_vertices()[2] == v else tc.get_vertices()[3]

                    if tn.is_22():
                        if vo == tc.get_vertices()[2]:
                            assert v == tc.get_vertices()[3]
                        if vo == tc.get_vertices()[3]:
                            assert v == tc.get_vertices()[2]
                    
                    tc = tn
        
                assert tc.is_31()
                
                hthis = t.get_half_edges()[(i + 1) % 3]
                hthat = tc.get_half_edge_to(t.get_vertices()[(i + 1) % 3])
                hthis.set_adj(hthat)
                hthat.set_adj(hthis)

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

    def check_validity(self):
        """
        Check the validity of the triangulation.
        """
        print("====================================================")
        print(f"Checking validity of the triangulation...")
        print("====================================================")
        print(f"Number of tetrahedra: {self.tetrahedron_pool.get_number_occupied()}")

        # Check that each tetrahedron has 4 vertices and 4 neighbouring tetrahedra         
        for i in self.tetrahedron_pool.get_objects():
            assert len(i.get_vertices()) == 4

        count = 0
        for t in self.tetrahedron_pool.get_objects():
            count += 1

            for i in range(4):
                assert self.vertex_pool.contains(t.get_vertices()[i].ID)

                for j in range(i + 1, 4):
                    assert t.get_vertices()[i] != t.get_vertices()[j]
            
            # Check if all neighbours still exist
            for i in range(4):
                if not self.tetrahedron_pool.contains(t.get_tetras()[i].ID):
                    print(f"Error: Tetrahedron {t.ID} has a missing neighbour.")
                    t.log()
                    t.get_tetras()[i].log()
                
                assert self.tetrahedron_pool.contains(t.get_tetras()[i].ID)
                assert t.get_tetras()[i].check_neighbours_tetra(t)
                assert t.get_tetras()[i] != t
                assert t.get_tetras()[i] >= 0

                sv = 0
                for j in range(4):
                    vv = t.get_tetras()[i].get_vertices()[j]
                    if t.has_vertex(vv):
                        sv += 1
                assert sv >= 3

                if t.is_31():
                    if i < 3:
                        assert t.get_tetras()[i].is_31() or t.get_tetras()[i].is_22()
                    else:
                        assert t.get_tetras()[i].is_13()
                elif t.is_13():
                    if i == 0:
                        assert t.get_tetras()[i].is_31()
                    else:
                        assert t.get_tetras()[i].is_13() or t.get_tetras()[i].is_22()
            
            for i in range(4):
                assert t.get_tetra_opposite(t.get_vertices()[i]) == t.get_tetras()[i]
                assert t.get_tetras()[i].get_tetra_opposite(t.get_vertex_opposite(t.get_vertices()[i])) == t
                
        for v in self.vertex_pool.get_objects():
            assert self.tetrahedron_pool.contains(v.get_tetra().ID)

            # Tadpole restriction
            if self.strictness == 1:
                assert v.scnum >= 2
            # Self-energy restriction
            if self.strictness == 2:
                assert v.scnum >= 3
        
        print("====================================================")

    def export_geometry(self, geometry_outfilename: str = "output.txt") -> bool:
        """
        Export the geometry of the Universe.
        """
        self.update_geometry()

        vertex_map = {v.ID: i for i, v in enumerate(self.vertex_pool.get_objects())}
        int_v_map = [v.ID for v in self.vertex_pool.get_objects()]

        tetra_map = {t.ID: i for i, t in enumerate(self.tetrahedron_pool.get_objects())}
        int_t_map = [t.ID for t in self.tetrahedron_pool.get_objects()]
        
        # Indicating well-orderedness
        out = "1\n"
        out += f"{self.vertex_pool.get_number_occupied()}\n"
        out += '\n'.join(str(v.time) for v in int_v_map) + '\n'
        out += f"{self.vertex_pool.get_number_occupied()}\n"
        out += f"{self.tetrahedron_pool.get_number_occupied()}\n"
        out += '\n'.join(str(vertex_map[v.ID]) for t in int_t_map for v in t.get_vertices()) + '\n'
        out += '\n'.join(str(tetra_map[t.ID]) for t in int_t_map for t in t.get_tetras())
        out += f"{self.tetrahedron_pool.get_number_occupied()}"

        with open(f"saved_universes/{geometry_outfilename}", "w") as file:
            file.write(out + "\n")
        
        print(f"Geometry exported to {geometry_outfilename}.")
        
        return True

if __name__ == "__main__":
    u = Universe(geometry_infilename="initial_universes/output.txt")