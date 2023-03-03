
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator_open import PathGenerator
import time

## try different initial conditions
## move the max velocity around
## try direction constraint
## add function to get path length

waypoints = np.array([[0,5],[0,0]])
# velocities = np.array([[1,1],[0,0]]) # 0
velocities = np.array([[1,0],[0,1]]) # 1
# velocities = np.array([[0,0],[1,-1]]) # 2
# velocities = np.array([[-1,-1],[0,0]]) # 3
# velocities = np.array([[0,0],[1,1]]) # 4
# velocities = np.array([[-1,1],[0,0]]) # 5

##### For Path optimization!! #####
# minimize length with directions, not velocities
# minimize acceleration with max vel 2 times velocities at the start points

initial_control_points = None
# initial_control_points = np.array([[2,0,-2,-2,2.5,7,7,5,3],[0,0,0,4,4,4,0,0,0]]) # case 3 top
# initial_control_points = np.array([[2,0,-2,-2,2.5,7,7,5,3],[0,0,0,-4,-4,-4,0,0,0]]) # case 3 bottom
# initial_control_points = np.array([[5,-1,-2,3,2,7,6,0],[3,-2,3,4,-3,-3,1,-3]]) # case 3 middle
# initial_control_points = np.array([[0,0,0,2.5,2.5,2.5,5,5,5],[-3,0,3,3,0,-3,-3,0,3]]) # case 4
# initial_control_points = np.array([[2,0,-2,-2,2,2,3,5,7],[0,0,0,-2,-2,-2,0,0,0]]) # case 5

velocities = velocities/np.linalg.norm(velocities,2,0)

# objective_function_type = "minimize_control_point_distance_and_time"
objective_function_type = "minimize_distance_and_time"
# objective_function_type = "minimize_acceleration"
# objective_function_type = "minimize_velocity"
# objective_function_type = "maximize_velocity"

max_velocity = 2
# max_acceleration = 10
max_curvature = 1
# center_sfc = np.array([[2.5],[2]]) # top
center_sfc = np.array([[2.5],[-2]]) # bottom
# center_sfc = np.array([[2.5],[0]])
dimensions_sfc = np.array([[12],[4]])

dimension = np.shape(waypoints)[0]
order = 3
# curvature_method = "roots_of_curvature_derivative"
curvature_method = "roots_numerator_and_denominator"
# curvature_method = "control_point_derivatives"
# curvature_method = "curvature_at_min_velocity"

path_gen = PathGenerator(order, dimension, curvature_method, objective_function_type)
start_time = time.time()
control_points, scale_factor = path_gen.generate_trajectory(waypoints, velocities, \
    max_velocity, max_curvature, center_sfc, dimensions_sfc,initial_control_points)
end_time = time.time()

spline_start_time = 0
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)
spline_at_knot_points, knot_points = bspline.get_spline_at_knot_points()
curvature_data, time_data = bspline.get_spline_curvature_data(number_data_points)
velocity_data, time_data = bspline.get_spline_derivative_data(5,1)
velocity_magnitude_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,1)
acceleration_magnitude_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,2)

# bspline.plot_spline(number_data_points)
# bspline.plot_curvature(number_data_points)
# bspline.plot_derivative_magnitude(number_data_points, 2)
fig, ax = plt.subplots()
ax.plot(control_points[0,:],control_points[1,:],color='tab:orange')
ax.plot(spline_data[0,:],spline_data[1,:],label="path",color="tab:blue")
ax.scatter(spline_at_knot_points[0,:],spline_at_knot_points[1,:],label="interval endpoints",color="tab:blue")
ax.scatter(control_points[0,:],control_points[1,:],label="control points",facecolors='none', edgecolors='tab:orange')
plt.title("Optimized Path")
# ax.add_patch(Rectangle((center_sfc[0]-dimensions_sfc[0]/2, center_sfc[1]-dimensions_sfc[1]/2),\
#      dimensions_sfc[0], dimensions_sfc[1], edgecolor = 'tab:green', fill=False))
# plt.legend()
plt.show()

dynamics_figure, axis = plt.subplots(1,3)
axis[0].set_title("Curvature")
axis[0].plot(time_data, curvature_data, color="tab:green")
axis[0].set_xlabel("time")
axis[0].tick_params(width=1.00)
axis[1].set_title("Velocity Magnitude")
axis[1].plot(time_data, velocity_magnitude_data, color="b")
axis[1].tick_params(width=1.00)
axis[1].set_xlabel("time")
axis[2].set_title("Acceleration Magnitude")
axis[2].plot(time_data, acceleration_magnitude_data, color="tab:red")
axis[2].tick_params(width=1.00)
axis[2].set_xlabel("time")
plt.show()

print("control points: " ,  control_points)
print("method: ", curvature_method)
print("objective: " , objective_function_type)
print("evaluation time: " , np.round(end_time - start_time,3) , "sec")
print("path length: ", np.round(bspline.get_arc_length(),3), "units")
print("path time: ", np.round(bspline.get_end_time(),3), "sec")
print("highest curvature: " , np.round(np.max(curvature_data),2))
print("V0: ", np.round(velocity_data[:,0],2), "units/sec")
print("Vf: ", np.round(velocity_data[:,-1],2), "units/sec")
print("vel mean: " , np.round(np.mean(velocity_magnitude_data),3))
print("vel std: " , np.round(np.std(velocity_magnitude_data),3))
print("acceleration_integral: " , round(np.sum(acceleration_magnitude_data)*time_data[1],2))
