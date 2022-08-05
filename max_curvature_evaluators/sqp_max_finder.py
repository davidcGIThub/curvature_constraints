import numpy as np
import time
from scipy.optimize import minimize, Bounds
from helper import calculate_curvature, calculate_curvature_derivative
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.matrix_evaluation import get_M_matrix
import matplotlib.pyplot as plt

"""
Uses the slsqp method 
"""

order = 5
control_points = np.random.randint(10, size=(2,order+1)) # random
# control_points = np.array([[8, 8, 4, 0],[1, 1, 4, 5]])
# control_points = np.array([[0, 1, 1, 4],[8, 5, 5, 1]])
# control_points = np.array([[7, 1, 7, 1],[9, 7, 9, 4]]) #produces nan
# control_points = np.array( [[3, 7, 8, 8], [2, 6, 5, 9]])
print("control_points: " , control_points)

def find_max_curvature_sqp_method(control_points):
    num_control_points = count_number_of_control_points(control_points)
    order = num_control_points - 1
    M = get_M_matrix(0, order, [], False)
    num_initial_guesses = order+1
    def curvature_objective_function(time):
        t = time[0]
        curvature = calculate_curvature(t,M,control_points)
        return -curvature
    def curvature_derivative(time):
        t = time[0]
        curvature_derivative = calculate_curvature_derivative(t,M,control_points)
        return -curvature_derivative
    time_bounds = Bounds(lb=0.0, ub =1.0)
    t = 0
    max_curvature = 0
    for i in range(num_initial_guesses):
        t = t + i/(num_initial_guesses-1)
        result = minimize(curvature_objective_function, x0=np.array([t]), \
            method='SLSQP', bounds=time_bounds, jac=curvature_derivative)
        curvature = -result.fun
        max_curvature = np.max([curvature,max_curvature])
    return max_curvature

start_time = time.time()
max_curvature = find_max_curvature_sqp_method(control_points)
end_time = time.time()
total_time = 0
iterations = 100
average_accuracy = 0
variance = 0
for i in range(iterations):
    control_points = np.random.randint(10, size=(2,order+1)) # random
    start_time = time.time()
    max_curvature = find_max_curvature_sqp_method(control_points)
    end_time = time.time()
    total_time += end_time - start_time
    bspline = BsplineEvaluation(control_points,order,0,1)
    curvature_data, time_data = bspline.get_spline_curvature_data(10000)
    true_max = np.max(curvature_data)
    if true_max < 1e-10 and max_curvature < 1e-10:
        accuracy = 1
    elif true_max < 1e-10:
        accuracy = 0
    elif true_max == np.inf and max_curvature == np.inf:
        accuracy = 1
    elif true_max == np.inf:
        accuracy = 0
    else:
        accuracy = max_curvature/true_max
    average_accuracy += accuracy/iterations
    print("accuracy: " , accuracy)
    variance += (1-accuracy)**2
ave_time = total_time/iterations
standard_deviation = np.sqrt(variance/iterations)
print("standard deviation: " , standard_deviation)
print("ave_accuracy: " , average_accuracy)
print("ave_time: " , ave_time)

# print("total time: " , end_time - start_time)
# print("max_curvature: ", max_curvature)
# bspline = BsplineEvaluation(control_points,order,0,1)
# curvature_data, time_data = bspline.get_spline_curvature_data(10000)
# max_left_arg = np.argmax(curvature_data[0:5000])
# max_right_arg = np.argmax(curvature_data[5000:9999]) + 5000
# local_max = np.array([curvature_data[max_left_arg],curvature_data[max_right_arg]])
# local_max_time = np.array([time_data[max_left_arg],time_data[max_right_arg]])
# initial_conditions, intial_cond_time = bspline.get_spline_curvature_data(order+1)
# plt.plot(time_data,curvature_data, label = "curvature")
# plt.scatter(intial_cond_time, initial_conditions, label = "intial conditions")
# plt.scatter(local_max_time, local_max, label = "local maxima")
# plt.plot(time_data,curvature_data*0 + max_curvature , label = "estimated max curvature")
# plt.legend()
# plt.xlabel("time")
# plt.ylabel("curvature")
# plt.title("Curvature of 3rd Order Spline")
# plt.show()
# print("true_max_curvature: " , np.max(curvature_data))
# bspline.plot_spline(1000)
# bspline.plot_spline_vs_time(1000)
# bspline.plot_curvature(1000)
# bspline.plot_derivative_magnitude(1000,1)