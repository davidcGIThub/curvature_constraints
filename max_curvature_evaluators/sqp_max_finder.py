import numpy as np
import time
from scipy.optimize import minimize, Bounds
from helper import calculate_curvature, calculate_curvature_derivative
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.matrix_evaluation import get_M_matrix

"""
Uses the slsqp method 
"""

order = 3
control_points = np.random.randint(10, size=(2,order+1)) # random
control_points = np.array([[8, 8, 4, 0],[1, 1, 4, 5]])
# control_points = np.array([[0, 1, 1, 4],[8, 5, 5, 1]])
# control_points = np.array([[7, 1, 7, 1],[9, 7, 9, 4]]) #produces nan
print("control_points: " , control_points)

def find_max_curvature_sqp_method(control_points):
    num_control_points = count_number_of_control_points(control_points)
    order = num_control_points - 1
    M = get_M_matrix(0, order, [], False)
    num_initial_guesses = 4
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
print("total time: " , end_time - start_time)
print("max_curvature: ", max_curvature)
bspline = BsplineEvaluation(control_points,order,0,1)
curvature_data, time_data = bspline.get_spline_curvature_data(10000)
print("true_max_curvature: " , np.max(curvature_data))
bspline.plot_spline(1000)
bspline.plot_spline_vs_time(1000)
bspline.plot_curvature(1000)
bspline.plot_derivative_magnitude(1000,1)