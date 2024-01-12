# pool.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description:

class Pool:
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
    def __init__(self, capacity: int, silent: bool = True) -> None:
        self.silent = silent
        self.capacity = capacity
        self.elements = [{'index': i, 'used': False, 'next': i + 1} for i in range(self.capacity)]
        self.elements[-1]['next'] = None
        self.used_indices = set()
        self.first = 0

    def occupy(self) -> int:
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

        # Mark object as used
        self.elements[index]['used'] = True

        # Add index to set of used indices
        self.used_indices.add(index)

        # Update first free object
        self.first = self.elements[index]['next']

        if not self.silent:
            print(self.elements)
            print("added:", index)
            print("used_indices:", self.used_indices)
            print("first:", self.first)
            print()
        
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
        self.elements[index]['used'] = False

        # Remove index from set of used indices
        self.used_indices.remove(index)

        # Update first free object
        self.elements[index]['next'] = self.first
        self.first = index

        if not self.silent:
            print(self.elements)
            print("destroyed:", index)
            print("used_indices:", self.used_indices)
            print("first:", self.first)
            print()