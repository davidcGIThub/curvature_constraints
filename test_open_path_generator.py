
import numpy as np
import matplotlib.pyplot as plt
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator_open import PathGenerator
import time

## try different initial conditions
## move the max velocity around
## try direction constraint


waypoints = np.array([[0,5],[0,0]])
# directions = np.array([0,0])
velocities = np.array([[-1,-1],[0,0]])
velocities = velocities/np.linalg.norm(velocities,2,0)
print("velocities: " , velocities)
max_velocity = 1
max_curvature = 0.5

dimension = np.shape(waypoints)[0]
order = 3
curvature_method = "roots_of_curvature_derivative"
path_gen = PathGenerator(order, dimension, curvature_method)
start_time = time.time()
control_points, scale_factor = path_gen.generate_trajectory(waypoints, velocities, max_velocity, max_curvature)
end_time = time.time()
print("control_points: " ,  control_points)
print("time: " , end_time - start_time)
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
bspline.plot_derivative_vs_time(number_data_points, 1)
# bspline.plot_spline_vs_time(number_data_points)
bspline.plot_curvature(number_data_points)
bspline.plot_derivative_magnitude(number_data_points, 1)
# bspline.plot_derivative_magnitude(number_data_points, 2)