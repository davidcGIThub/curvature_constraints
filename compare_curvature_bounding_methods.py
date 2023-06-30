import numpy as np
from bsplinegenerator.bsplines import BsplineEvaluation
from max_curvature_evaluators.helper_files.max_curvature_test_functions import test_performance_of_max_curvature_function, create_control_point_and_max_curvature_list, \
    get_curvature_bound_and_time, create_edge_cases
import matplotlib.pyplot as plt

order = 2
num_iterations = 100
isRestricted = False
if order == 2:
    methods = np.array(["discrete_evaluations", "maximize_curvature_equation", "roots_of_curvature_derivative", "geometric", "roots_of_curvature_numerator_and_denominator", "control_point_derivatives"]) #, "geometric"])
    colors = np.array(["tab:orange","b", "g", "r", "c", "m"])
    label_list = ['discrete evaluations','newton method' ,
                    'roots of curvature derivative', 'geometric', 'analytical roots',
                    'control points']
    tick_list = [1, 2, 3, 4, 5, 6]
if order == 3:
    methods = np.array(["discrete_evaluations", "maximize_curvature_equation", "roots_of_curvature_derivative", "roots_of_curvature_numerator_and_denominator", "control_point_derivatives"]) #, "geometric"])
    colors = np.array(["tab:orange","b", "g", "c", "m"])
    label_list = ['discrete evaluations','newton method' ,
                    'roots of curvature derivative', 'analytical roots',
                    'control points']
    tick_list = [1, 2, 3, 4, 5]
if order >= 4:
    methods = np.array(["discrete_evaluations", "maximize_curvature_equation", "roots_of_curvature_derivative", "control_point_derivatives"])
    colors  = np.array(["tab:orange" , "b" , "g" , "m"])
    label_list = ['discrete evaluations','newton method' ,
                    'roots of curvature derivative',
                    'control points']
    tick_list = [1, 2, 3, 4]
    
# curvature_methods = ["geometric","roots_numerator_and_denominator", "control_point_derivatives_mdm", "constrain_max_acceleration_and_min_velocity"]
# colors = np.array(["r", "c", "m", "y"])
plt.rcParams['font.family'] = 'DeJavu Serif'
plt.rcParams['font.serif'] = ['Times New Roman']

control_point_list, max_curvature_list = create_control_point_and_max_curvature_list(order, num_iterations)
fig, (ax1, ax2) = plt.subplots(1, 2,gridspec_kw={'width_ratios': [2, 1]})
ax1.set_xticks(tick_list)
# ax.set_xticklabels(label_list)
evaluation_time_array = np.zeros(len(methods))
evaluation_time_std_array = np.zeros(len(methods))
ax1.plot([0,len(methods)],[0,0],color = "gray",label="mean")
ax1.plot([0,len(methods)],[0,0],color = "gray",label="median", linestyle='dotted') 
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
    evaluation_time_array[i] = average_time
    evaluation_time_std_array[i] = std_time
    j = i + 1
    x_data = np.zeros(len(relative_error_array)) + j
    ax1.scatter(x_data, relative_error_array,facecolors='none',edgecolors=color, label=str(i+1) + " " + label_list[i])
    ax1.plot([j-0.25,j,j+0.25],np.array([median,median,median]), color = color, linestyle=":")
    ax1.plot([j-0.25,j,j+0.25],np.array([mean,mean,mean]), color = color)

ax1.set_ylabel("Curvature Error for Random B-splines")
ax1.set_xlabel("Evaluation Methods")
ax1.set_title( str(order) + " Order Curvature Error Calculations")
# ax1.set_yscale("symlog") 
ax1.set_yscale('symlog' , linthresh=0.00001)
ax1.legend()

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




ax2.set_xticks(tick_list)
for i in range(len(methods)):
    ax2.errorbar(i+1,evaluation_time_array[i],evaluation_time_std_array[i], fmt='-o',capsize=7,color=colors[i],label= str(i) + " " + label_list[i])


ax2.set_xlabel("Evaluation Methods")
ax2.set_ylabel("Computationa Time")
ax2.set_title("Time Analysis")
ax2.set_yscale("log") 
ax2.yaxis.set_label_position("right")
# ax2.legend()
plt.show()