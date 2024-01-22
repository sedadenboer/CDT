# pool.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:

from vertex import Vertex


class Pool(object):
    """
    The Pool class represents a pool of indices that can be used to create
    and destroy objects. It is used to keep track of the objects that are
    currently in use and the objects that are available for use.

    Attributes:
        capacity (int): Maximum number of objects in pool.
        elements (list): List of objects.
        used_indices (set): Set of indices of used objects.
        first (int): Index of first free object.
    """
    def __init__(self, capacity: int) -> None:
        self.capacity = capacity
        self.elements = [None for _ in range(self.capacity)]
        self.x = [0] * self.capacity
        self.p = [i + 1 for i in range(self.capacity)]
        self.p[-1] = None
        self.used_indices = set()
        self.first = 0
    
    def occupy(self, obj: object) -> int:
        """
        Occupies a space in the pool and returns the index of the object.

        Returns:
            int: Index of object.
        """
        # Check if pool is full
        if self.first is None:
            raise Exception("Pool is full.")

        # Get index of first free object
        index = self.first

        # Set index of object
        obj.ID = index

        # Mark object as used
        self.elements[index] = obj
        self.x[index] = 1

        # Add index to set of used indices
        self.used_indices.add(index)

        # Update first free object
        self.first = self.p[index]

        return index

    def free(self, index: int) -> None:
        """
        Frees a space in the pool.

        Args:
            index (int): Index of object to free.

        Raises:
            Exception: If index is not valid.
        """
        # Check if index is valid
        if index not in self.used_indices:
            raise Exception("Index is not valid.")

        # Mark object as unused
        self.elements[index] = None
        self.x[index] = 0

        # Remove index from set of used indices
        self.used_indices.remove(index)

        # Update first free object
        self.p[index] = self.first
        self.first = index

    def get(self, index: int) -> object:
        """
        Gets an object from the pool.

        Args:
            index (int): Index of object to get.

        Returns:
            object: Object.
        """
        return self.elements[index]
    
    def log(self) -> None:
        """
        Prints the state of the pool, including the IDs of the elements.
        """
        print(f"elements: {self.elements}")
        print(f"x: {self.x}")
        print(f"p: {self.p}")
        print(f"\nused indices: {self.used_indices}\nfirst: {self.first}\n")

if __name__ == "__main__":
    vertex = Vertex(0)
    vertex2 = Vertex(0)
    vertex3 = Vertex(0)
    vertex4 = Vertex(0)
    pool = Pool(10)
    index = pool.occupy(vertex)
    index2 = pool.occupy(vertex2)
    pool.free(index)
    index3 = pool.occupy(vertex3)
    index4 = pool.occupy(vertex4)
    pool.free(index2)
    pool.occupy(vertex2)
    

   