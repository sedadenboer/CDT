# bag.py
#
# Author: Seda den Boer
# Date: 03-01-2024
# 
# Description: A bag of indices that can be used to pick random objects
# from a set of objects with certain properties.

import random
from typing import Optional, Set


class Bag:
    """
    The Bag class represents a bag of objects (pool indices) that can
    be added to and removed from. It is used to be able to pick random
    objects from a set of objects with certain properties.

    Attributes:
        pool_capacity (int): Maximum number of objects in bag.
        elements (List[int]): List of pool indices.
        indices (List[int]): List of indices of pool indices.
        size (int): Number of pool indices in bag.
    """
    EMPTY = -1

    def __init__(self, pool_capacity: int):
        self.pool_capacity = pool_capacity
        self.elements = [self.EMPTY] * pool_capacity
        self.indices = [self.EMPTY] * pool_capacity
        self.used_indices: Set[int] = set()
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
        self.used_indices.add(pool_index)
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

        self.used_indices.remove(pool_index)
        self.size -= 1

    def pick(self) -> Optional[int]:
        """
        Pick a random object from Bag.

        Returns:
            int: Random pool index from Bag.
        """
        if self.size == 0:
            raise ValueError("Bag is empty.")
        elif self.size > 0:
            return random.choice(list(self.used_indices))
        
        # Bag is empty
        return None 
    
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
        Returns:
            int: Number of occupied spaces in Bag.
        """
        return self.size
    
    def log(self) -> None:
        """
        Print Bag.
        """
        print("elements")
        for i in range(self.size):
            print(f"size index {i}: {self.elements[i]}")

        print("--")
    
    def __str__(self) -> str:
        return str(self.elements)