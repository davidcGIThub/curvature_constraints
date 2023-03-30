
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator_lin import PathGenerator
import time

## try different initial conditions
## move the max velocity around
## try direction constraint
## add function to get path length
max_curvature = 1.5
order = 3
waypoints = np.array([[0,6],[0,0]])
dimension = np.shape(waypoints)[0]
# velocities = np.array([[1,1],[0,0]]) # 0
# velocities = np.array([[1,0],[0,1]]) # 1
velocities = np.array([[0,0],[1,-1]]) # 2
# velocities = np.array([[-1,-1],[0,0]]) # 3
# velocities = np.array([[0,0],[1,1]]) # 4
# velocities = np.array([[-1,1],[0,0]]) # 5
velocities = velocities/np.linalg.norm(velocities,2,0) # normalize velocities


colors = np.array(["r", "c", "m", "y"])
fig, ax = plt.subplots(1,2)
max_x = 0
max_y = 0
min_x = 0
min_y = 0
path_gen = PathGenerator(order, dimension)
start_time = time.time()
control_points = path_gen.generate_path(waypoints, velocities, max_curvature)
end_time = time.time()
spline_start_time = 0
bspline = BsplineEvaluation(control_points, order, spline_start_time, 1, False)
number_data_points = 1000000
spline_data, time_data = bspline.get_spline_data(number_data_points)
spline_at_knot_points, knot_points = bspline.get_spline_at_knot_points()
curvature_data, time_data = bspline.get_spline_curvature_data(number_data_points)
velocity_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,1)
    # centripetal_acceleration_data, time_data = bspline.get_centripetal_acceleration_data(number_data_points)
acceleration_magnitude_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,2)
ax[0].plot(spline_data[0,:],spline_data[1,:],label="path",color=colors[0])
ax[1].plot(time_data,velocity_data,label="path",color=colors[1])
evaluation_time = end_time - start_time
path_length = bspline.get_arc_length(number_data_points)
acceleration_integral = np.sum(acceleration_magnitude_data)*time_data[1]
curvature_extrema = np.max(curvature_data)
print("max curvature: " , np.max(max_curvature))
print("curvature extrema: " , np.max(curvature_data))
print("path length: " , path_length)
print("evaluation time: " , evaluation_time)
ax[0].scatter(waypoints[0,:],waypoints[1,:],color='k', s=8)
ax[0].set_xlabel("Path",weight='bold')
ax[1].set_xlabel("Velocity",weight='bold')
plt.show()

