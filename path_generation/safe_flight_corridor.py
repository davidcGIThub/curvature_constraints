import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt

@dataclass
class SFC:
    """Safe Flight Corridor Data Class"""
    dimensions: np.ndarray(dtype=np.float64, shape=(2,1))
    translation: np.ndarray(dtype=np.float64, shape=(2,1))
    rotation: np.ndarray(dtype=np.float64, shape=(2,2))

    def getRotatedBounds(self):
        center = self.rotation.T @ self.translation
        max_bounds = center + self.dimensions/2
        min_bounds = center - self.dimensions/2
        return min_bounds, max_bounds
    
def plot_2D_sfc(sfc: SFC):
    min_bounds, max_bounds = sfc.getRotatedBounds()
    x_min = min_bounds.item(0)
    x_max = max_bounds.item(0)
    y_min = min_bounds.item(1)
    y_max = max_bounds.item(1)
    points = np.array([[x_min, x_min, x_max, x_max, x_min],
                        [y_min, y_max, y_max, y_min, y_min]])
    points = sfc.rotation @ points
    plt.plot(points[0,:], points[1,:])

def plot_2D_sfcs(sfcs: list):
    if sfcs != None:
        for sfc_index in range(len(sfcs)):
            plot_2D_sfc(sfcs[sfc_index])

