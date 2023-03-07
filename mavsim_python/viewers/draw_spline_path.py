"""
mavsim_python: drawing tools
    - Beard & McLain, PUP, 2012
    - Update history:
        4/15/2019 - BGM
"""

import numpy as np
import pyqtgraph.opengl as gl
from bsplinegenerator.bsplines import BsplineEvaluation


class DrawSplinePath:
    def __init__(self, control_points, order, window):
        self._order = order
        bspline = BsplineEvaluation(control_points, self._order, 0,1)
        points_xyz, time_data = bspline.get_spline_data(100)
        self._R = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        points = self._R @ points_xyz
        self.path_plot_object = gl.GLLinePlotItem(pos=points.T,
                                                  width=1,
                                                  antialias=True,
                                                  mode='line_strip')
        window.addItem(self.path_plot_object)

    def update(self, control_points):
        bspline = BsplineEvaluation(control_points, self._order, 0,1)
        points_xyz, time_data = bspline.get_spline_data(100)
        points = self._R @ points_xyz
        self.path_plot_object.setData(pos=points.T)


