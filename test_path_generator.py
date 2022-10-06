
import numpy as np
import matplotlib.pyplot as plt
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator import PathGenerator

waypoints = np.array([[1,4],[2,4]])

max_curvature = .1
dimension = np.shape(waypoints)[0]
order = 3
path_gen = PathGenerator(order, dimension)
control_points = path_gen.generate_trajectory(waypoints)
scale_factor = 1
start_time = 0
bspline = BsplineEvaluation(control_points, order, start_time, scale_factor, False)
number_data_points = 1000
spline_data, time_data = bspline.get_spline_data(number_data_points)

bspline.plot_spline(number_data_points)
# plt.plot(spline_data[0,:], spline_data[1,:])
# plt.scatter(waypoints[0,:], waypoints[1,:])
# plt.plot(control_points[0,:], control_points[1,:])
# plt.scatter(control_points[0,:], control_points[1,:])
# plt.show()
# bspline.plot_derivative(number_data_points, 1)
# bspline.plot_spline_vs_time(number_data_points)
# bspline.plot_curvature(number_data_points)
# bspline.plot_derivative_magnitude(number_data_points, 1)
# bspline.plot_derivative_magnitude(number_data_points, 2)