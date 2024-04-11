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

    Attributes:
        observable (str): The name of the observable.
        thermal_sweeps (int): The number of thermal sweeps.
        main_sweeps (int): The number of main sweeps.
        k_steps (int): The number of k steps.
        data (list): The data of the observable.
        k0 (float): The k0 value.
    """
    def __init__(self, observable: str, thermal_sweeps: int, main_sweeps: int, k_steps: int, k0: float):
        self.observable = observable
        self.thermal_sweeps = thermal_sweeps
        self.main_sweeps = main_sweeps
        self.k_steps = k_steps
        self.k0 = k0
        self.data = []

    def measure(self, data_point: Any):
        """
        Measure a data point.

        Args:
            data_point (Any): The data point to measure.
        """
        self.data.append(deepcopy(data_point))
    
    def get_data(self) -> np.ndarray:
        """
        Retrieve the observable data.

        Returns:
            np.ndarray: The observable data as a numpy array.
        """
        return np.array(self.data)
    
    def clear_data(self):
        """
        Clear the observable data.
        """
        self.data.clear()

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
