import numpy as np
import time
from scipy.optimize import minimize, Bounds
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import calculate_curvature, calculate_curvature_derivative
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.matrix_evaluation import get_M_matrix
import matplotlib.pyplot as plt

"""
Uses the slsqp method 
"""

def get_m_matrix_sqp(order):
    M = get_M_matrix(0, order, [], False)
    return M

def find_max_curvature_sqp_method(control_points,order,M):
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