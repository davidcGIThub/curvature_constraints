import numpy as np
from scipy.optimize import root_scalar
from helper import calculate_curvature, calculate_curvature_derivative
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.matrix_evaluation import get_M_matrix
import time
import matplotlib.pyplot as plt
"""
Uses the brent_q bracketing method 
"""

order = 3
control_points = np.random.randint(10, size=(2,order+1)) # random
control_points = np.array( [[3, 7, 8, 8], [2, 6, 5, 9]])
# control_points = np.array([[8, 2, 4, 0],
#  [1, 2, 2, 5]])
# control_points = np.array([[0, 1, 1, 4],[8, 2, 5, 1]])
print("control_points: " , control_points)

def find_max_curvature_root_finder(control_points):
    num_control_points = count_number_of_control_points(control_points)
    order = num_control_points - 1
    num_intervals = order+1
    M = get_M_matrix(0, order, [], False)
    def curvature_derivative(t):
        return calculate_curvature_derivative(t,M,control_points)
    interval_start = 0
    interval_length = 1/(num_intervals)
    curvature_start = calculate_curvature(0,M,control_points)
    curvature_end = calculate_curvature(1,M,control_points)
    max_curvature = np.max((curvature_start,curvature_end))
    for i in range(num_intervals):
        try:
            solution = root_scalar(curvature_derivative,method='brentq',bracket=[interval_start,interval_start+interval_length])
            root = solution.root
            curvature = calculate_curvature(root,M,control_points)
            max_curvature = np.max((curvature,max_curvature))
        except:
            pass
        interval_start = interval_start+interval_length
    return max_curvature

total_time = 0
iterations = 100
average_accuracy = 0
# for i in range(iterations):
#     control_points = np.random.randint(10, size=(2,order+1)) # random
#     start_time = time.time()
#     max_curvature = find_max_curvature_root_finder(control_points)
#     end_time = time.time()
#     total_time += end_time - start_time
#     bspline = BsplineEvaluation(control_points,order,0,1)
#     curvature_data, time_data = bspline.get_spline_curvature_data(10000)
#     true_max = np.max(curvature_data)
#     if true_max < 1e-10 and max_curvature < 1e-10:
#         accuracy = 1
#     elif true_max < 1e-10:
#         accuracy = 0
#     elif true_max == np.inf and max_curvature == np.inf:
#         accuracy = 1
#     elif true_max == np.inf:
#         accuracy = 0
#     else:
#         accuracy = max_curvature/true_max
#     average_accuracy += accuracy/iterations
#     print("average_accuracy: " , average_accuracy)
# ave_time = total_time/iterations
# print("ave_accuracy: " , average_accuracy)
# print("ave_time: " , ave_time)

max_curvature = find_max_curvature_root_finder(control_points)
# print("max_curvature: ", max_curvature)
# bspline = BsplineEvaluation(control_points,order,0,1)
# curvature_data, time_data = bspline.get_spline_curvature_data(1000)
# print("true_max_curvature: " , np.max(curvature_data))
# # bspline.plot_spline(1000)
# bspline.plot_curvature(1000)

bspline = BsplineEvaluation(control_points,order,0,1)
curvature_data, time_data = bspline.get_spline_curvature_data(10000)
root_left_arg = np.argmax(curvature_data[0:5000])
root_right_arg = np.argmax(curvature_data[5000:9999]) + 5000
roots = np.array([0,0])
roots_time = np.array([time_data[root_left_arg],time_data[root_right_arg]])
plt.plot(time_data[0:-1],curvature_data[1:]-curvature_data[0:-1], label = "curvature derivative")
plt.scatter(roots_time, roots, label="roots")
plt.legend()
plt.xlabel("time")
plt.ylabel("curvature derivative")
plt.title("Derivative of the Curvature of a 3rd Order Spline")
plt.show()