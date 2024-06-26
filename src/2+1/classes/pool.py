# pool.py
#
# Author: Seda den Boer
# Date: 17-02-2024
# 
# Description: Defines a pool of indices that can be used to create and destroy objects.

from __future__ import annotations
import numpy as np
import random
from typing import TYPE_CHECKING, List, Optional, Union
if TYPE_CHECKING:
    from vertex import Vertex
    from tetra import Tetrahedron


class Pool:
    """
    The Pool class represents a pool of indices that can be used to create
    and destroy objects. It is used to keep track of the objects that are
    currently in use and the objects that are available for use.

    Args (Attributes):
        capacity (int): Maximum number of objects in pool.

    Attributes:
        elements (np.ndarray): Array of objects.
        x (np.ndarray): Array of flags indicating if object is in use.
        p (np.ndarray): Array of pointers to next free object.
        first (int): Index of first free object.
        size (int): Number of objects in pool.
    """
    def __init__(self, capacity: int):
        self.capacity: int = capacity
        self.elements: np.ndarray = np.empty(self.capacity, dtype=object)
        self.x: np.ndarray = np.zeros(self.capacity, dtype=int)
        self.p: np.ndarray = np.arange(1, self.capacity + 1, dtype=object)
        self.p[-1] = None
        self.first: Optional[int] = 0
        self.size = 0
    
    def occupy(self, obj: Union[Vertex, Tetrahedron]) -> int:
        """
        Occupies a space in the pool and returns the index of the object.

        Args:
            obj (Union[Vertex, Tetrahedron]): Object to occupy space in pool.

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
        self.size += 1

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
        if self.x[index] == 0:
            raise Exception(f"Object with index {index} is not in the pool.")

        # Mark object as unused
        self.elements[index] = None
        self.x[index] = 0
        self.size -= 1

        # Update first free object
        self.p[index] = self.first
        self.first = index

    def free_all(self):
        """
        Frees all spaces in the pool.
        """
        self.elements = np.empty(self.capacity, dtype=object)
        self.x = np.zeros(self.capacity, dtype=int)
        self.p = np.arange(1, self.capacity + 1, dtype=object)
        self.p[-1] = None
        self.first = 0
        self.size = 0

    def get(self, index: int) -> Union[Vertex, Tetrahedron, None]:
        """
        Gets an object from the pool.

        Args:
            index (int): Index of object to get.

        Returns:
            Union[Vertex, Tetrahedron, None]: Object.
        """
        return self.elements[index]
    
    def get_used_indices(self) -> List[int]:
        """
        Gets the indices of the used objects in the pool.

        Returns:
            List[int]: Set of indices of used objects.
        """
        return np.where(self.x == 1)[0]
    
    def pick(self) -> int:
        """
        Picks a random object from the pool.

        Raises:
            Exception: If pool is empty.

        Returns:
            int: Index of object.
        """
        if self.size == 0:
            raise Exception("Pool is empty.")
        else:
            # Pick random element from pool
            random_element = self.elements[random.randint(0, self.size - 1)]

            # If random element is empty, pick another element
            while not random_element:
                random_element = self.elements[random.randint(0, self.size - 1)]

            return random_element.ID
 
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
        return self.size
    
    def get_objects(self) -> np.darray[Union[Vertex, Tetrahedron]]:
        """
        Gets the objects in the pool without the empty spaces.

        Returns:
            np.darray[Union[Vertex, Tetrahedron]]: Array of objects.
        """
        return self.elements[self.x == 1]
    
    def log(self):
        """
        Prints the state of the pool, including the IDs of the elements.
        """
        print(f"elements: {self.elements}")
        print(f"x: {self.x}")
        print(f"p: {self.p}")
        print(f"first: {self.first}\n")