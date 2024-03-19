# observable.py
#
# Author: Seda den Boer
# Date: 29-02-2024
# 
# Description: The Observable class that forms the base class for all observables.


import os
import random
from typing import List
from universe import Universe 


class Observable:
    data_dir: str = ""
    rng: random.Random = random.Random(0)
    done_list: List[bool] = []

    def __init__(self, identifier: str):
        self.name: str = ""
        self.identifier: str = identifier
        self.extension: str = ".dat"
        self.output: str = ""

    def reset(self):
        """
        Reset the observable to its initial state.
        """
        pass

    def measure(self):
        """
        Measure the observable by processing the data and writing the output to a file.	
        """
        self.process()
        self.write()

    def clear(self):
        """
        Clear the output file and reinitialize the observable.
        """
        filename = os.path.join(self.data_dir, f"{self.name}-{self.identifier}{self.extension}")
        with open(filename, "w") as file:
            file.close()
        
        # Reinitialize the observable
        self.reset()

    def write(self):
        """
        Write the output to the observable file.
        """
        filename = os.path.join(self.data_dir, f"{self.name}-{self.identifier}{self.extension}")
        with open(filename, "a") as file:
            file.write(f"{self.output}\n")
            file.close()

    def sphere(origin: int, radius: int) -> List[int]:
        """
        Find vertices within a sphere of a given radius around a given origin.

        Args:
            origin (int): Origin Vertex label.
            radius (int): Radius of the sphere.

        Returns:
            _type_: _description_
        """
        this_depth = []
        next_depth = []

        vertex_list = []
        flipped_vertices = []
        Observable.done_list[origin] = True
        this_depth.append(origin)
        flipped_vertices.append(origin)

        for current_depth in range(radius):
            for v in this_depth:
                for neighbor in Universe.vertex_neighbours[v]:
                    if not Observable.done_list[neighbor]:
                        flipped_vertices.append(neighbor)
                        next_depth.append(neighbor)
                        Observable.done_list[neighbor] = True
                        if current_depth == radius - 1:
                            vertex_list.append(neighbor)
            this_depth = next_depth
            next_depth = []

        for v in flipped_vertices:
            Observable.done_list[v] = False

        return vertex_list

    def sphere2d(origin, radius):
        """
        Find triangles within a sphere of a given radius of the origin in the 2D dual lattice.

        Args:
            origin (_type_): _description_
            radius (_type_): _description_

        Returns:
            _type_: _description_
        """
        this_depth = []
        next_depth = []

        vertex_list = []
        flipped_vertices = []
        Observable.done_list[origin] = True
        this_depth.append(origin)
        flipped_vertices.append(origin)

        for current_depth in range(radius):
            for v in this_depth:
                for neighbor in Universe.vertex_neighbors[v]:
                    if neighbor.time != origin.time:
                        continue
                    if not Observable.done_list[neighbor]:
                        flipped_vertices.append(neighbor)
                        next_depth.append(neighbor)
                        Observable.done_list[neighbor] = True
                        if current_depth == radius - 1:
                            vertex_list.append(neighbor)
            this_depth = next_depth
            next_depth = []

        for v in flipped_vertices:
            Observable.done_list[v] = False

        return vertex_list

    def sphere_dual(origin, radius):
        done = []
        this_depth = []
        next_depth = []

        done.append(origin)
        this_depth.append(origin)

        tetra_list = []

        for current_depth in range(radius):
            for t in this_depth:
                for neighbor in t.tnbr:
                    if neighbor not in done:
                        next_depth.append(neighbor)
                        done.append(neighbor)
                        if current_depth == radius - 1:
                            tetra_list.append(neighbor)
            this_depth = next_depth
            next_depth = []

        return tetra_list

    def sphere2d_dual(origin, radius):
        done = []
        this_depth = []
        next_depth = []

        done.append(origin)
        this_depth.append(origin)

        triangle_list = []

        for current_depth in range(radius):
            for t in this_depth:
                for neighbor in t.trnbr:
                    if neighbor not in done:
                        next_depth.append(neighbor)
                        done.append(neighbor)
                        if current_depth == radius - 1:
                            triangle_list.append(neighbor)
            this_depth = next_depth
            next_depth = []

        return triangle_list
