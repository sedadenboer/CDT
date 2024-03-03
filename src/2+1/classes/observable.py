# observable.py
#
# Author: Seda den Boer
# Date: 29-02-2024
# 
# Description: The Observable class that forms the base class for all observables.


class Observable:
    """
    The Observable class provides some of the basic functionality
    that one may require when implementing certain observables.
    """
    def __init__(self):
        pass

    def update(self):
        pass

    def measure(self):
        pass

    def set_data_dir(self, data_dir: str):
        """
        Set the data directory for the observable.

        Args:
            data_dir (str): The directory to store the data in.
        """
        self.data_dir = data_dir
        self.data = []