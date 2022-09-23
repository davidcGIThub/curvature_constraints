import numpy as np
from bsplinegenerator.bsplines import BsplineEvaluation
from max_curvature_evaluators.helper_files.max_curvature_test_functions import test_performance_of_max_curvature_function, create_control_point_and_max_curvature_list
import random
import matplotlib.pyplot as plt

# isRestricted = True
# method = "maximize_curvature_equation"
# method = "roots_of_curvature_derivative"
# method = "discrete_evaluations"
# method = "curvature_at_min_velocity"
# method = "max_numerator_over_min_denominator"
# method = "control_point_derivatives"
# method = "geometric"
# method = "angle_constrained_control_points"
order = 3
num_iterations = 1000
isRestricted = False
methods = np.array(["discrete_evaluations", "maximize_curvature_equation", "roots_of_curvature_derivative", \
    "curvature_at_min_velocity", "max_numerator_over_min_denominator", "control_point_derivatives"])#, "geometric"])
colors = np.array(["b", "g", "r", "c", "m", "y", "orange",])
# methods = np.array(["discrete_evaluations", "maximize_curvature_equation", \
#     "roots_of_curvature_derivative", "control_point_derivatives"])
# colors = np.array(["b", "g", "r", "y"])
control_point_list, max_curvature_list = create_control_point_and_max_curvature_list(order, num_iterations)
plt.figure(1)
# plt.plot([0,len(methods)],[0,0],color = "gray", linestyle=":") #, label="mode")
plt.plot([0,len(methods)],[0,0],color = "gray") #, label = "mean")
for i in range(len(methods)):
    method = methods[i]
    color = colors[i]
    mode, mean, std, relative_error_array, average_time, std_time = \
        test_performance_of_max_curvature_function(order,method,control_point_list, max_curvature_list)
    print("method: " , method)
    print("mean: " , mean)
    print("std: ", std)
    print("mode ", mode)
    print("high: " , np.max(relative_error_array))
    print("low: " , np.min(relative_error_array))
    print("average_time: ", average_time)
    print("std_time: " , std_time)
    print(" ")
    x_data = np.zeros(len(relative_error_array)) + i
    plt.scatter(x_data, relative_error_array,facecolors='none',edgecolors=color)
    plt.plot([i-0.25,i,i+0.25],np.array([mode,mode,mode]), color = color, linestyle=":")
    plt.plot([i-0.25,i,i+0.25],np.array([mean,mean,mean]), color = color, label = method)

plt.legend()
plt.ylabel("Relative Error")
plt.xlabel("Methods of Finding Curvature Extrema")
plt.title("Curvature Bounding Methods for 5th Order Spline")
plt.yscale("log") 
plt.show()


plt.figure(2)
# plt.plot([0,len(methods)],[0,0],color = "gray", linestyle=":") #, label="mode")
plt.plot([0,len(methods)],[0,0],color = "gray") #, label = "mean")
for i in range(len(methods)):
    method = methods[i]
    color = colors[i]
    mode, mean, std, relative_error_array, average_time, std_time = \
        test_performance_of_max_curvature_function(order,method,control_point_list, max_curvature_list)
    print("method: " , method)
    print("mean: " , mean)
    print("std: ", std)
    print("mode ", mode)
    print("high: " , np.max(relative_error_array))
    print("low: " , np.min(relative_error_array))
    print("average_time: ", average_time)
    print("std_time: " , std_time)
    print(" ")
    x_data = np.zeros(len(relative_error_array)) + i
    plt.scatter(x_data, -relative_error_array,facecolors='none',edgecolors=color)
    plt.plot([i-0.25,i,i+0.25],-np.array([mode,mode,mode]), color = color, linestyle=":")
    plt.plot([i-0.25,i,i+0.25],-np.array([mean,mean,mean]), color = color, label = method)

plt.legend()
plt.ylabel("Relative Error")
plt.xlabel("Methods of Finding Curvature Extrema")
plt.title("Curvature Bounding Methods for 5th Order Spline")
plt.yscale("log") 
plt.show()
