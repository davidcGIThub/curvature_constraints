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
        self._R = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        bspline = BsplineEvaluation(control_points, self._order, 0,1)
        points_xyz, time_data = bspline.get_spline_data(100)
        start_point = points_xyz[:,0][:,None]
        final_point = points_xyz[:,-1][:,None]
        end_points = np.concatenate((start_point,final_point),1)
        end_points = self._R @ end_points
        print("end_points: " , end_points)
        start_vel_vector = bspline.get_derivative_at_time_t(0,1)
        end_time = bspline.get_end_time()
        end_vel_vector = bspline.get_derivative_at_time_t(end_time,1)
        points = self._R @ points_xyz
        self._blue = np.array([[30, 144, 255, 255]])/255.
        self._red = (1., 0., 0., 1)
        self.path_plot_object = gl.GLLinePlotItem(pos=points.T,
                                                  width=1,
                                                  antialias=True,
                                                  mode='line_strip')
        self.endpoint_plot_object = gl.GLScatterPlotItem(pos=end_points.T,
                                                        size=10,
                                                        color = self._red)
        window.addItem(self.path_plot_object)
        window.addItem(self.endpoint_plot_object)

    def update(self, control_points):
        bspline = BsplineEvaluation(control_points, self._order, 0,1)
        points_xyz, time_data = bspline.get_spline_data(100)
        start_point = points_xyz[:,0][:,None]
        final_point = points_xyz[:,-1][:,None]
        end_points = np.concatenate((start_point,final_point),1)
        end_points = self._R @ end_points
        points = self._R @ points_xyz
        self.path_plot_object.setData(pos=points.T)
        self.endpoint_plot_object.setData(pos=end_points.T,
                                          size=10,
                                          color = self._red)


