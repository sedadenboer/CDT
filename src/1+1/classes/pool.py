# pool.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description: A pool of indices that can be used to create and destroy objects.

from __future__ import annotations
from typing import TYPE_CHECKING, List, Optional, Set, Union
if TYPE_CHECKING:
    from vertex import Vertex
    from triangle import Triangle


class Pool:
    """
    The Pool class represents a pool of indices that can be used to create
    and destroy objects. It is used to keep track of the objects that are
    currently in use and the objects that are available for use.

    Attributes:
        capacity (int): Maximum number of objects in pool.
        elements (List[Union[Triangle, Vertex, None]]): List of objects.
        used_indices (Set[int]): Set of indices of used objects.
        first (int): Index of first free object.
    """
    def __init__(self, capacity: int):
        self.capacity: int = capacity
        self.elements: List[Union[Triangle, Vertex, None]] = [None for _ in range(self.capacity)]
        self.x: List[int] = [0] * self.capacity
        self.p: List[Optional[int]] = [i + 1 for i in range(self.capacity)]
        self.p[-1] = None
        self.used_indices: Set[int] = set()
        self.first: Optional[int] = 0
    
    def occupy(self, obj: Union[Triangle, Vertex]) -> int:
        """
        Occupies a space in the pool and returns the index of the object.

        Returns:
            int: Index of object.
        """
        # Check if pool is full
        if self.first is None:
            raise Exception("Pool is full.")

        # Get index of first free object
        index: int = self.first

        # Mark object as used
        self.elements[index] = obj
        self.x[index] = 1

        # Add index to set of used indices
        self.used_indices.add(index)

        # Update first free object
        self.first = self.p[index]

        # Set index of object
        obj.ID = index

        return index

    def free(self, index: int):
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

    def get(self, index: int) -> Union[Triangle, Vertex, None]:
        """
        Gets an object from the pool.

        Args:
            index (int): Index of object to get.

        Returns:
            Union[Triangle, Vertex, None]: Object.
        """
        return self.elements[index]
    
    def contains(self, index: int) -> bool:
        """
        Checks if the pool contains an object with the given index.

        Args:
            index (int): Index of object to check.

        Returns:
            bool: True if pool contains object, False otherwise.
        """
        if index < 0 or index >= self.capacity:
            return False
        
        if not self.elements[index]:
            return False
        
        return True
    
    def get_number_occupied(self) -> int:
        """
        Gets the number of occupied spaces in the pool.

        Returns:
            int: Number of occupied spaces.
        """
        return len(self.used_indices)
    
    def log(self):
        """
        Prints the state of the pool, including the IDs of the elements.
        """
        print(f"elements: {self.elements}")
        print(f"x: {self.x}")
        print(f"p: {self.p}")
        print(f"\nused indices: {self.used_indices}\nfirst: {self.first}\n")
