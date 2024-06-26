# observable.py
#
# Author: Seda den Boer
# Date: 29-02-2024
# 
# Description: The Observable class that forms the base class for all observables.

import numpy as np
import os
from copy import deepcopy
from typing import Any


class Observable:
    """
    Represents an observable in the simulation.

    Args (Attributes):
        observable (str): The name of the observable.
        thermal_sweeps (int): The number of thermal sweeps.
        main_sweeps (int): The number of main sweeps.
        k0 (float): The k0 value.
        measuring_interval (int): The interval at which to measure the observable.

    Attributes:
        data (list): The data of the observable.
        next (int): The next index to measure.
    """
    def __init__(self, observable: str, thermal_sweeps: int, main_sweeps: int, k0: float, measuring_interval: int):
        self.observable = observable
        self.thermal_sweeps = thermal_sweeps
        self.main_sweeps = main_sweeps
        self.k0 = k0
        self.measuring_interval = measuring_interval
        self.data = [None] * (((thermal_sweeps + main_sweeps) // measuring_interval) + 1)
        self.next = 0

    def measure(self, data_point: Any):
        """
        Measure a data point.

        Args:
            data_point (Any): The data point to measure.
        """
        if isinstance(data_point, (int, float, str)):
            self.data[self.next] = data_point
        else:
            self.data[self.next] = deepcopy(data_point)
        self.next += 1
        assert self.next <= len(self.data)
    
    def get_data(self) -> np.ndarray:
        """
        Retrieve the observable data.

        Returns:
            np.ndarray: The observable data as a numpy array.
        """
        return self.data
    
    def clear_data(self):
        """
        Clear the observable data.
        """
        self.data = [None] * (((self.thermal_sweeps + self.main_sweeps) // self.measuring_interval) + 1)

    def save_data(self, filename: str):
        """
        Save the observable data to a file.

        Args:
            filename (str): The filename to save the data to.
        """
        # Create the directory if it does not exist
        pathname = f'measurements/k0={self.k0}/'
        if not os.path.exists(pathname):
            os.makedirs(pathname)

        np.save(pathname + filename, self.data)
