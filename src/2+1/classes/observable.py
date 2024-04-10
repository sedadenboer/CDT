# observable.py
#
# Author: Seda den Boer
# Date: 29-02-2024
# 
# Description: The Observable class that forms the base class for all observables.

import numpy as np
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
    """
    def __init__(self, observable: str, thermal_sweeps: int, main_sweeps: int, k_steps: int):
        self.obserable = observable
        self.thermal_sweeps = thermal_sweeps
        self.main_sweeps = main_sweeps
        self.k_steps = k_steps
        self.data = []

    def measure(self, data_point: Any):
        """
        Measure a data point.

        Args:
            data_point (Any): The data point to measure.
        """
        self.data.append(data_point)
    
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
        np.save(filename, self.data)
