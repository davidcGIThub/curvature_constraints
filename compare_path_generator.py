
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator_compare import PathGenerator
from path_generation.path_generator_2nd_order import PathGenerator2nd
from path_generation.spline_order_converter import SmoothingSpline
import time

## try different initial conditions
## move the max velocity around
## try direction constraint
## add function to get path length
max_curvature = 1.5
order = 3
waypoints = np.array([[0,6],[0,0]])
dimension = np.shape(waypoints)[0]
velocities = np.array([[1,1],[0,0]]) # 0
velocities = np.array([[1,0],[0,1]]) # 1
# velocities = np.array([[0,0],[1,-1]]) # 2
# velocities = np.array([[-1,-1],[0,0]]) # 3
# velocities = np.array([[0,0],[1,1]]) # 4
# velocities = np.array([[-1,1],[0,0]]) # 5
velocities = velocities/np.linalg.norm(velocities,2,0) # normalize velocities

initial_control_points = None
# initial_control_points = np.array([[2,0,-2,-2,2.5,7,7,5,3],[0,0,0,4,4,4,0,0,0]]) # case 3 top
# initial_control_points = np.array([[2,0,-2,-2,2.5,7,7,5,3],[0,0,0,-4,-4,-4,0,0,0]]) # case 3 bottom
# initial_control_points = np.array([[5,-1,-2,3,2,7,6,0],[3,-2,3,4,-3,-3,1,-3]]) # case 3 middle
# initial_control_points = np.array([[0,0,0,2.5,2.5,2.5,5,5,5],[-3,0,3,3,0,-3,-3,0,3]]) # case 4
# initial_control_points = np.array([[2,0,-2,-2,2,2,3,5,7],[0,0,0,-2,-2,-2,0,0,0]]) # case 5

# curvature_method = "roots_of_curvature_derivative"
# curvature_method = "roots_numerator_and_denominator"
# curvature_method = "control_point_derivatives"
curvature_method = "constrain_max_acceleration_and_min_velocity"
# curvature_methods = ["roots_of_curvature_derivative", "roots_numerator_and_denominator", "control_point_derivatives_mdm", "control_point_derivatives_rotate", "constrain_max_acceleration_and_min_velocity"]
curvature_methods = ["geometric","roots_numerator_and_denominator", "control_point_derivatives_mdm", "constrain_max_acceleration_and_min_velocity"]
# curvature_methods = ["roots_numerator_and_denominator", "control_point_derivatives_rotate", "constrain_max_acceleration_and_min_velocity"]
# curvature_methods = ["roots_numerator_and_denominator","geometric","geometric","geometric"]
# colors = np.array(["r", "c"])
num_methods = len(curvature_methods)
colors = np.array(["r", "c", "m", "y", "b"])
fig, ax = plt.subplots(3,num_methods)
max_x = 0
max_y = 0
min_x = 0
min_y = 0

# smoother = SmoothingSpline(new_order, dimension, 1000)
# new_control_points, new_scale_factor = smoother.generate_new_control_points(control_points,scale_factor,order)

for i in range(0,len(curvature_methods)):
    print(i)
    if curvature_methods[i] == "geometric":
        path_gen = PathGenerator2nd(dimension)
        start_time = time.time()
        control_points, scale_factor = path_gen.generate_path(waypoints, velocities, max_curvature)
        end_time = time.time()
        order_ = 2
    else:
        path_gen = PathGenerator(order, dimension, curvature_methods[i])
        start_time = time.time()
        control_points, scale_factor = path_gen.generate_path(waypoints, velocities, max_curvature)
        end_time = time.time()
        order_ = order

    spline_start_time = 0
    bspline = BsplineEvaluation(control_points, order_, spline_start_time, scale_factor, False)
    number_data_points = 1000000
    spline_data, time_data = bspline.get_spline_data(number_data_points)
    spline_at_knot_points, knot_points = bspline.get_spline_at_knot_points()
    curvature_data, time_data = bspline.get_spline_curvature_data(number_data_points)
    velocity_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,1)
    # centripetal_acceleration_data, time_data = bspline.get_centripetal_acceleration_data(number_data_points)
    acceleration_magnitude_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,2)
    ax[0,i].plot(spline_data[0,:],spline_data[1,:],label="path",color=colors[i])
    ax[1,i].plot(time_data,velocity_data,label="path",color=colors[i])
    ratio = 1
    max_x = np.max(np.concatenate((spline_data[0,:].flatten(),np.array([max_x])),0))
    min_x = np.min(np.concatenate((spline_data[0,:].flatten(),np.array([min_x])),0))
    max_y = np.max(np.concatenate((spline_data[1,:].flatten(),np.array([max_y])),0))
    min_y = np.min(np.concatenate((spline_data[1,:].flatten(),np.array([min_y])),0))
    if curvature_methods[i] == "roots_of_curvature_derivative":
        ax[0,i].set_title("Roots of Curvature \n Derivative")
    elif curvature_methods[i] == "roots_numerator_and_denominator":
        ax[0,i].set_title("Analytical Roots of \n Numerator & Denominator")
    elif curvature_methods[i] == "control_point_derivatives_mdm":
        ax[0,i].set_title("Control Point \n Derivatives MDM")
    elif curvature_methods[i] == "constrain_max_acceleration_and_min_velocity":
        ax[0,i].set_title("Max Acceleration \n and Min Velocity")
    elif curvature_methods[i] == "geometric":
        ax[0,i].set_title("Geometric")
    evaluation_time = end_time - start_time
    path_length = bspline.get_arc_length(number_data_points)
    acceleration_integral = np.sum(acceleration_magnitude_data)*time_data[1]
    # centripetal_acceleration_integral = np.sum(centripetal_acceleration_data)*time_data[1]/time_data[-1]
    # curvature_integral = np.sum(curvature_data)*time_data[1]
    curvature_extrema = np.max(curvature_data)
    width = 0.5
    ax[2,0].bar([i], [evaluation_time] , color=colors[i])
    ax[2,1].bar([i], [path_length] , color=colors[i])
    ax[2,2].bar([i], [curvature_extrema] , color=colors[i])
    ax[2,3].bar([i], [acceleration_integral] , color=colors[i])
    head_width = 0.3
    head_length = 0.25
    ax[0,i].arrow(waypoints[0,0], waypoints[1,0], velocities[0,0]/2, velocities[1,0]/2,
                  head_width=head_width, head_length=head_length, color = 'k')
    ax[0,i].arrow(waypoints[0,1], waypoints[1,1], velocities[0,1]/2, velocities[1,1]/2,
                  head_width=head_width, head_length=head_length, color = 'k')
    ax[0,i].scatter(waypoints[0,:],waypoints[1,:],color='k', s=8)


for i in range(0,len(curvature_methods)):
    ax[0,i].set_aspect(abs((max_x-min_x)/(max_y-min_y))*ratio)
    # ax[0,i].set_ylim([-1, 1])
ax[2,0].set_ylabel("Evaluation Time",weight='bold')
ax[2,1].set_ylabel("Path Length",weight='bold')
ax[2,2].set_ylabel("Curvature Extrema",weight='bold')
ax[2,2].plot([-1,4],[max_curvature, max_curvature], color='k')
ax[2,2].text(-0.5,max_curvature+0.15,"max curvature")
ax[0,0].set_ylabel("Optimized Paths",weight='bold')
# ax[1,2].legend()
ax[2,3].set_ylabel("accel_integral",weight='bold')
ax[1,0].set_ylabel("Velocity Mag",weight='bold')

fig.supxlabel("Curvature Constraint Methods")
fig.suptitle("Curvature Contrained Path Optimizations: Case 3", weight='bold')
plt.tight_layout()
plt.subplots_adjust(wspace=None, hspace=None)
# plt.legend()
plt.show()

