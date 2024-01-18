# vertex.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:


class Vertex:
    capacity = 100

    def __init__(self, ID, time):
        self.ID = ID
        self.time = time
        self.triangles = []

    def add_triangle(self, triangle):
        self.triangles.append(triangle)
    
    def get_number_of_triangles(self):
        return len(self.triangles)
    
    def get_linked_triangles(self):
        # Get the triangles that intersect with their base
        linked_triangles = []
        for triangle in self.triangles:
            if triangle.vl == self:
                linked_triangles.append(triangle)
