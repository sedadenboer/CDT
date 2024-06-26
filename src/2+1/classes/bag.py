# bag.py
#
# Author: Seda den Boer
# Date: 17-02-2024
# 
# Description: Defines a bag of objects that can be added to and removed from.

import random
import numpy as np
from typing import Optional, List


class Bag:
    """
    The Bag class represents a bag of objects (pool indices) that can
    be added to and removed from. It is used to be able to pick random
    objects from a set of objects with certain properties.

    Args (Attributes): 
        pool_capacity (int): Maximum number of objects in bag.
        
    Attributes:
        pool_capacity (int): Maximum number of objects in bag.
        elements (np.ndarray): Array of pool indices.
        indices (np.ndarray): Array of indices of pool indices.
        size (int): Number of pool indices in bag.
    """
    EMPTY = -1

    def __init__(self, pool_capacity: int):
        self.pool_capacity = pool_capacity
        self.elements = np.full(pool_capacity, self.EMPTY, dtype=int)
        self.indices = np.full(pool_capacity, self.EMPTY, dtype=int)
        self.size = 0

    def add(self, pool_index: int):
        """
        Add pool index to Bag.

        Args:
            pool_index (int): Object to add to bag.

        Raises:
            ValueError: If bag is full.
        """
        if self.size >= self.pool_capacity:
            raise ValueError("Bag is full.")
        
        self.elements[self.size] = pool_index
        self.indices[pool_index] = self.size
        self.size += 1

    def remove(self, pool_index: int):
        """
        Remove pool index from Bag.

        Args:
            pool_index (int): Object to remove from bag.

        Raises:
            ValueError: If pool index is not in bag.
        """
        if self.indices[pool_index] == self.EMPTY:
            raise ValueError(f"Object with pool_index {pool_index} is not in the bag.")
  
        index = self.indices[pool_index]
        last_pool_index = self.elements[self.size - 1]

        # Replace pool_index at index with last_pool_index
        self.elements[index] = last_pool_index
        self.indices[last_pool_index] = index

        # Clear last_pool_index from Bag
        self.elements[self.size - 1] = self.EMPTY
        self.indices[pool_index] = self.EMPTY
        self.size -= 1

    def clear(self):
        """
        Clear Bag.
        """
        self.elements = np.full(self.pool_capacity, self.EMPTY, dtype=int)
        self.indices = np.full(self.pool_capacity, self.EMPTY, dtype=int)
        self.size = 0
        
    def pick(self) -> Optional[int]:
        """
        Pick a random object from Bag.

        Raises:
            ValueError: If Bag is empty.

        Returns:
            int: Random pool index from Bag.
        """
        if self.size <= 0:
            raise ValueError("Bag is empty.")
        elif self.size > 0:
            # Pick random element from Bag
            random_element_id = self.elements[random.randint(0, self.size - 1)]

            # If random element is empty, pick another random element
            while random_element_id == self.EMPTY:
                random_element_id = self.elements[random.randint(0, self.size - 1)]

            return random_element_id
    
    def contains(self, pool_index: int) -> bool:
        """
        Check if Bag contains pool index.

        Args:
            pool_index (int): Object to check if in bag.

        Returns:
            bool: True if pool index is in bag, False otherwise.
        """
        if pool_index < 0 or pool_index >= self.pool_capacity:
            return False
        
        return self.indices[pool_index] != self.EMPTY

    def get_number_occupied(self) -> int:
        """
        Get number of occupied spaces in Bag.

        Returns:
            int: Number of occupied spaces in Bag.
        """
        return self.size
    
    def get_used_indices(self) -> List[int]:
        """
        Gets the indices of the used objects in the Bag.

        Returns:
            List[int]: List of indices of used objects in Bag.
        """
        return np.where(self.indices != self.EMPTY)[0]
    
    def log(self):
        """
        Print Bag.
        """
        print(f"Bag size: {self.size}")
        print("Elements:")
        for i in range(self.size):
            print(f"size index {i}: {self.elements[i]}")
        print(self.elements)
        print(self.indices)
        print("--")
    
    def __str__(self) -> str:
        return str(self.elements)


if __name__ == "__main__":
    bag = Bag(5)
    bag.log()
    bag.add(0)
    bag.add(1)
    bag.add(2)
    bag.add(3)
    bag.add(4)
    bag.log()
    bag.remove(2)
    bag.remove(0)
    bag.log()
    print(bag.contains(2))
    print(bag.get_number_occupied())
    print(bag.get_used_indices())
    bag.clear()
    bag.log()