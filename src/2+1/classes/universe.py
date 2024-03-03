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
                    # Set tetrahedron neighbors for each vertex
                    for vertex in tetra_vs:
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
        
        # # TEST
        # types_31 = {}
        # types_13 = {}
        # types_22 = {}
        # # Check which vertices dont have a time assigned
        # for tetra in self.tetrahedron_pool.get_objects():
        #     if tetra.type == '31':
        #         assert [t.type for t in tetra.get_tetras()] == ['13', '31', '31', '22']
        #         if str([v.time for v in tetra.get_vertices()]) in types_31:
        #             types_31[str([v.time for v in tetra.get_vertices()])] += 1
        #         else:
        #             types_31[str([v.time for v in tetra.get_vertices()])] = 1
        #     elif tetra.type == '13':
        #         if str([v.time for v in tetra.get_vertices()]) in types_13:
        #             types_13[str([v.time for v in tetra.get_vertices()])] += 1
        #         else:
        #             types_13[str([v.time for v in tetra.get_vertices()])] = 1
        #         assert [t.type for t in tetra.get_tetras()] == ['22', '22', '13', '31']
        #     elif tetra.type == '22':
        #         if str([v.time for v in tetra.get_vertices()]) in types_22:
        #             types_22[str([v.time for v in tetra.get_vertices()])] += 1
        #         else:
        #             types_22[str([v.time for v in tetra.get_vertices()])] = 1
        #         assert [t.type for t in tetra.get_tetras()] == ['31', '22', '13', '13']

        #     print(f"tetra: {tetra.ID}, type: {tetra.type}")
        #     print(f"vertices: {[v.ID for v in tetra.get_vertices()]}")
        #     print(f"vertex times: {[v.time for v in tetra.get_vertices()]}")
        #     print(f"tetras: {[t.ID for t in tetra.get_tetras()]}")
        #     print(f"tetra times: {[t.time for t in tetra.get_tetras()]}")
        #     print(f"tetra types: {[t.type for t in tetra.get_tetras()]}")
        #     print()

        # print(f"types_31: {types_31}\n")
        # print(f"types_13: {types_13}\n")
        # print(f"types_22: {types_22}\n")
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

    def delete(self, vertex_id: int) -> bool:
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

    def ishift_d(self, tetra13_id: int, tetra22_id: int, tetra22r_id: int) -> bool:
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
        # Free all vertices in the pool
        self.vertex_pool.free_all()
        self.vertex_neighbours.clear()
        self.vertex_neighbours = [None] * self.vertex_pool.get_number_occupied()

        for vertex in self.vertex_pool.get_objects():
            neighbours = []
            current_tetra = vertex.get_tetra()
            current_tetra_list = [current_tetra]
            next_tetra_list = []
            done_tetra_list = []

            # Use a breadth-first search to find all neighbouring vertices
            while current_tetra_list:
                for current_tetra_check in current_tetra_list:
                    for tetra_neighbour in current_tetra_check.get_tetras():
                        # Check if the neighbour tetrahedron contains the vertex
                        if not tetra_neighbour.has_vertex(vertex):
                            continue

                        # Check if the neighbour tetra has been visited
                        if tetra_neighbour not in done_tetra_list:
                            # Mark the neighbour tetra as visited and add it to the next level
                            done_tetra_list.append(tetra_neighbour)
                            next_tetra_list.append(tetra_neighbour)

                # Move to the next level of tetrahedra
                current_tetra_list = next_tetra_list
                next_tetra_list = []

            # Add all vertices from the visited tetrahedra to the neighbours list
            for done_tetra in done_tetra_list:
                for vertex_done_tetra in done_tetra.get_vertices():
                    # Add neighbouring vertices that are not the current vertex
                    if vertex_done_tetra not in neighbours and vertex_done_tetra != vertex:
                        neighbours.append(vertex_done_tetra)

            # Store the list of neighbours for the current vertex
            self.vertex_neighbours[vertex.ID] = neighbours

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
                halfedge_triples[i].set_prev(halfedge_triples[(i + 2) % 3])
        
        # Connect the halfedges of the tetrahedra
        for tetra in tetras31_objs:
            he_id_str = []
            for he in tetra.get_half_edges():
                if he:
                    he_id_str.append(str(he.ID))
                else:
                    he_id_str.append('None')
            print(f"tetra: {tetra.ID}, he: {he_id_str}")
            for i in range(self.Constants.N_VERTICES_TRIANGLE):
                # Get the current vertex and the opposite vertex of the face
                vertex = tetra.get_vertices()[i]
                vertex_t = tetra.get_vertices()[3]

                # Get the opposite tetrahedron along the edge
                tetra_opposite = tetra.get_tetra_opposite(vertex)
                # Update the current vertex to the opposite vertex
                vertex = vertex_t
                
                # Iterate until the adjacent tetrahedron is not a (2,2)-simplex
                while tetra_opposite.is_22():
                    # Get the new opposite tetrahedron
                    new_tetra = tetra_opposite.get_tetra_opposite(vertex)

                    # Update the current vertex to the next vertex along the edge
                    previous_vertex = vertex
                    vertex = tetra_opposite.get_vertices()[2] if tetra_opposite.get_vertices()[2] == vertex \
                        else tetra_opposite.get_vertices()[3]

                    # Check if the next tetrahedron is also a (2,2)-simplex
                    if new_tetra.is_22():
                        # Assert conditions to ensure consistent edge connections
                        if previous_vertex == tetra_opposite.get_vertices()[2]:
                            assert vertex == tetra_opposite.get_vertices()[3], \
                                f"previous_vertex: {previous_vertex.ID}, vertex: {vertex.ID}, \
                                    tetra_opposite: \{tetra_opposite.ID}, new_tetra: {new_tetra.ID}"
                        if previous_vertex == tetra_opposite.get_vertices()[3]:
                            assert vertex == tetra_opposite.get_vertices()[2], \
                                f"previous_vertex: {previous_vertex.ID}, vertex: {vertex.ID}, \
                                    tetra_opposite: {tetra_opposite.ID}, new_tetra: {new_tetra.ID}"
                    
                    # Update the current tetrahedron to the next tetrahedron
                    tetra_opposite = new_tetra

                # assert tetra_opposite.is_13(), f"tetra_opposite: {tetra_opposite.ID}, type: {tetra_opposite.type}"
                
                # Link the current half edge and the half edge in the adjacent tetrahedron
                current_half_edge = tetra.get_half_edges()[(i + 1) % 3]
                adjacent_half_edge = tetra_opposite.get_half_edge_to(tetra.get_vertices()[(i + 1) % 3])
                current_half_edge.set_adjacent(adjacent_half_edge)
                adjacent_half_edge.set_adjacent(current_half_edge)

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
        # self.update_triangles()

    def check_validity(self):
        """
        Check the validity of the triangulation.
        """
        print("====================================================")
        print(f"Checking validity of the triangulation...")
        print("====================================================")
        print(f"Number of tetrahedra: {self.tetrahedron_pool.get_number_occupied()}")

        # Check that each tetrahedron has 4 vertices and 4 neighbouring tetrahedra         
        for t in self.tetrahedron_pool.get_objects():
            assert len(t.get_vertices()) == self.Constants.N_VERTICES_TETRA

        for t in self.tetrahedron_pool.get_objects():
            # Make sure that each vertex is contained in the vertex pool
            for i in range(self.Constants.N_VERTICES_TETRA):
                assert self.vertex_pool.contains(t.get_vertices()[i].ID)

                # Make sure that vertices are unique
                for j in range(i + 1, self.Constants.N_VERTICES_TETRA):
                    assert t.get_vertices()[i] != t.get_vertices()[j]
            
            # Check if all neighbours still exist
            for i in range(self.Constants.N_TETRA_NEIGHBOURS):
                # Make sure that each tetrahedron is contained in the tetrahedron pool
                if not self.tetrahedron_pool.contains(t.get_tetras()[i].ID):
                    print(f"Error: Tetrahedron {t.ID} has a missing neighbour.")
                    t.log()
                    t.get_tetras()[i].log()
                
                assert self.tetrahedron_pool.contains(t.get_tetras()[i].ID)
                # Make sure that tetrahedron neighbours have eachother as neighbours
                assert t.get_tetras()[i].check_neighbours_tetra(t)
                # Make sure that tetrahedron neighbours are unique
                assert t.get_tetras()[i] != t
                # ID of tetrahedron neighbours should be non-negative
                assert t.get_tetras()[i].ID >= 0

                # Check for shared vertices between neighbours
                shared_vertices = 0
                for j in range(self.Constants.N_VERTICES_TETRA):
                    shared_vertex = t.get_tetras()[i].get_vertices()[j]

                    # Keep track of number of shared vertices between neighbours
                    if t.has_vertex(shared_vertex):
                        shared_vertices += 1

                # Between neighbours there should be 3 shared vertices
                assert shared_vertices == 3

                # Check if the tetrahedra neighbours are of the correct type
                if t.is_31():
                    # If the tetrahedron is a (3,1)-simplex, then neighbour 0 should be a (1,3)-simplex
                    if i == 0:
                        assert t.get_tetras()[i].is_13(), f"Error: Tetrahedron {t.ID} has a wrong neighbour type."
                    else:
                        assert t.get_tetras()[i].is_22() or t.get_tetras()[i].is_31(), f"Error: Tetrahedron {t.ID} has a wrong neighbour type."
                if t.is_13():
                    # If the tetrahedron is a (1,3)-simplex, then neighbour 3 should be a (3,1)-simplex
                    if i == 3:
                        assert t.get_tetras()[i].is_31(), f"Error: Tetrahedron {t.ID} has a wrong neighbour type."
                    else:
                        assert t.get_tetras()[i].is_22() or t.get_tetras()[i].is_13(), f"Error: Tetrahedron {t.ID} has a wrong neighbour type."
            
            # Check if all opposite tetrahedra still exist and are correct
            for i in range(self.Constants.N_TETRA_NEIGHBOURS):
                assert t.get_tetra_opposite(t.get_vertices()[i]) == t.get_tetras()[i], \
                    f"Error: Vertex {t.get_vertices()[i].ID} in tetra {t.ID} has the wrong opposite tetra {t.get_tetras()[i].ID}.\
                        {t.get_tetra_opposite(t.get_vertices()[i].ID)} != {t.get_tetras()[i].ID}"
                
        for v in self.vertex_pool.get_objects():
            assert self.tetrahedron_pool.contains(v.get_tetra().ID)

            # Tadpole restriction
            if self.strictness == 1:
                assert v.scnum >= 2
            # Self-energy restriction
            if self.strictness == 2:
                assert v.scnum >= 3
        
        print("====================================================")

    def export_geometry(self, geometry_outfilename: str = "output") -> bool:
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

        with open(f"saved_universes/{geometry_outfilename}.txt", "w") as file:
            file.write(out + "\n")
        
        print(f"Geometry exported to {geometry_outfilename}.")
        
        return True

if __name__ == "__main__":
    u = Universe(geometry_infilename="initial_universes/output.txt")
    # u.update_geometry()
    u.check_validity()