import numpy as np
from bsplinegenerator.bsplines import BsplineEvaluation
from max_curvature_evaluators.helper_files.max_curvature_test_functions import test_performance_of_max_curvature_function, create_control_point_and_max_curvature_list, \
    get_curvature_bound_and_time, create_edge_cases
import matplotlib.pyplot as plt

order = 3
num_iterations = 1000
isRestricted = False
methods = np.array(["discrete_evaluations", "maximize_curvature_equation", "roots_of_curvature_derivative", "roots_of_curvature_numerator_and_denominator", "control_point_derivatives"]) #, "geometric"])
colors = np.array(["b", "g", "r", "c", "m"])
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']

control_point_list, max_curvature_list = create_control_point_and_max_curvature_list(order, num_iterations)
fig = plt.figure()
ax = fig.add_subplot()
ax.set_xticks([1,2,3,4,5])
ax.set_xticklabels(['discrete \n evaluations','maximize curvature \n equation', \
                    'roots of curvature \n derivative',  'roots of curvature \n numerator & denominator',
                    'control point \n derivatives'])

ax.plot([0,len(methods)],[0,0],color = "gray",label="mean")
ax.plot([0,len(methods)],[0,0],color = "gray",label="median", linestyle='dotted') 
for i in range(len(methods)):
    method = methods[i]
    color = colors[i]
    median, mean, std, relative_error_array, average_time, std_time = \
        test_performance_of_max_curvature_function(order,method,control_point_list, max_curvature_list)
    print("method: " , method)
    print("mean: " , mean)
    print("std: ", std)
    print("median ", median)
    print("high: " , np.max(relative_error_array))
    print("low: " , np.min(relative_error_array))
    print("average_time: ", average_time)
    print("std_time: " , std_time)
    print(" ")
    j = i + 1
    x_data = np.zeros(len(relative_error_array)) + j
    ax.scatter(x_data, relative_error_array,facecolors='none',edgecolors=color)
    ax.plot([j-0.25,j,j+0.25],np.array([median,median,median]), color = color, linestyle=":")
    ax.plot([j-0.25,j,j+0.25],np.array([mean,mean,mean]), color = color)

plt.ylabel("Relative Error for Random B-splines", weight='bold')
plt.xlabel("Methods of Evaluating Curvature Extrema",  weight='bold')
plt.title("B-Spline Curvature Extrema Error Calculations", weight='bold')
plt.yscale("symlog") 
plt.yscale('symlog' , linthresh=0.00001)
plt.legend()
plt.show()

# control_point_list, max_curvature_list = create_edge_cases(order)

# for i in range(len(methods)):
# # i = 0
#     method = methods[i]
#     color = colors[i]
#     print("method: " , method)
#     curvature_bound_1, evaluation_time = get_curvature_bound_and_time(method,control_point_list[0],order)
#     curvature_bound_2, evaluation_time = get_curvature_bound_and_time(method,control_point_list[1],order)
#     curvature_bound_3, evaluation_time = get_curvature_bound_and_time(method,control_point_list[2],order)
#     curvature_bound_4, evaluation_time = get_curvature_bound_and_time(method,control_point_list[3],order)
#     curvature_bound_5, evaluation_time = get_curvature_bound_and_time(method,control_point_list[4],order)
#     curvature_bound_6, evaluation_time = get_curvature_bound_and_time(method,control_point_list[5],order)
#     curvature_bound_7, evaluation_time = get_curvature_bound_and_time(method,control_point_list[6],order)

#     fig = plt.figure()
#     ax = fig.add_subplot()
#     # zero curvature, constant velocity
#     ax.set_xticks([1,2,3,4,5,6,7])
#     ax.set_xticklabels(['constant \n velocity','zero \n centripetal \n acceleration','180 degree \n turn','starts at \n zero velocity' ,\
#     "stagnant", "almost \n straight line", "sharp turn \n around"])
#     ax.scatter([1],[curvature_bound_1], facecolors='none', edgecolors=color, marker ="o", label = "constant velocity")
#     ax.scatter([1],[max_curvature_list[0]], color=color, marker="x")
#     # zero curvature, with forward acceleration
#     ax.scatter([2],[curvature_bound_2], facecolors='none',edgecolors=color, marker ="o", label = "zero centripetal acceleration")
#     ax.scatter([2],[max_curvature_list[1]], color=color, marker="x")
#     # infinite curvature middle, hits zero velocity middle
#     # plt.scatter([3],[curvature_bound_3], facecolors='none',edgecolors=color, marker ="o", label = "180 degree turn")
#     # plt.scatter([3],[max_curvature_list[2]], color=color, marker="x")
#     # zero curvature, hits zero velocity at endpoint
#     ax.scatter([4],[curvature_bound_4], facecolors='none',edgecolors=color, marker ="o", label = "starts at zero velocity")
#     ax.scatter([4],[max_curvature_list[3]], color=color, marker="x")
#     # zero velocity whole time
#     ax.scatter([5],[curvature_bound_5], facecolors='none',edgecolors=color, marker ="o", label = "stagnant")
#     ax.scatter([5],[max_curvature_list[4]], color=color, marker="x")
#     # almost straight line
#     ax.scatter([6],[curvature_bound_6], facecolors='none',edgecolors=color, marker ="o", label = "almost straight line")
#     ax.scatter([6],[max_curvature_list[5]], color=color, marker="x")
#     # almost straight line, sharp turn around
#     plt.scatter([7],[curvature_bound_1], facecolors='none',edgecolors=color, marker ="o", label = "sharp turn around")
#     plt.scatter([7],[max_curvature_list[6]], color=color, marker="x")
#     plt.title(method)
#     plt.yscale('symlog', linthresh=0.05)
#     # plt.locator_params(axis='y', numticks=5)
#     plt.show()



