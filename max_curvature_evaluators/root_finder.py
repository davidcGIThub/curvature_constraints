import numpy as np
from scipy.optimize import root_scalar
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import calculate_curvature, calculate_curvature_derivative
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.matrix_evaluation import get_M_matrix
import time
import matplotlib.pyplot as plt
"""
Uses the brent_q bracketing method 
"""

def get_m_matrix_root(order):
    M = get_M_matrix(0, order, [], False)
    return M

def find_max_curvature_root_finder(control_points,order,M):
    num_intervals = order+1
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