import numpy as np
from scipy.optimize import root_scalar
from helper import calculate_curvature, calculate_curvature_derivative
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.matrix_evaluation import get_M_matrix

"""
Uses the brent_q bracketing method 
"""

order = 3
control_points = np.random.randint(10, size=(2,order+1)) # random
# control_points = np.array([[8, 2, 4, 0],
#  [1, 2, 2, 5]])
control_points = np.array([[0, 1, 1, 4],[8, 2, 5, 1]])
print("control_points: " , control_points)

def find_max_curvature_root_finder(control_points):
    num_control_points = count_number_of_control_points(control_points)
    order = num_control_points - 1
    num_intervals = 4
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

max_curvature = find_max_curvature_root_finder(control_points)
print("max_curvature: ", max_curvature)
bspline = BsplineEvaluation(control_points,order,0,1)
curvature_data, time_data = bspline.get_spline_curvature_data(1000)
print("true_max_curvature: " , np.max(curvature_data))
# bspline.plot_spline(1000)
bspline.plot_curvature(1000)