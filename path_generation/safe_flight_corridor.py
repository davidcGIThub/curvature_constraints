import numpy as np
from dataclasses import dataclass

@dataclass
class SFC:
    """Safe Flight Corridor Data Class"""
    dimensions: np.ndarray(dtype=np.float64)
    translation: np.ndarray(dtype=np.float64)
    rotation: np.ndarray(dtype=np.float64)

    def getRotatedBounds(self):
        center = self.rotation.T @ self.translation
        max_bounds = center + self.dimensions/2
        min_bounds = center - self.dimensions/2
        return min_bounds, max_bounds