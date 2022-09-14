import numpy as np
import time
from bsplinegenerator.bsplines import BsplineEvaluation 
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import get_matrix
from max_curvature_evaluators.angle_constrained_max_finder import create_random_control_points_greater_than_angles
from max_curvature_evaluators.sqp_max_finder import find_max_curvature_sqp_method, get_m_matrix_sqp
from max_curvature_evaluators.root_finder import find_max_curvature_root_finder, get_m_matrix_root
from max_curvature_evaluators.discrete_evaluations import get_matrices_discrete_evaluations, get_max_curvature_by_checking_n_points
from max_curvature_evaluators.max_at_min_velocity_finder import find_curvature_at_min_velocity_magnitude

def get_speed_and_accuracy_of_max_curvature_function(order,isRestricted,method):
    control_points = create_control_points(isRestricted,order)
    curvature_bound, evaluation_time = get_curvature_bound_and_time(method, control_points,order)
    bspline = BsplineEvaluation(control_points,order,0,1)
    curvature_data, time_data = bspline.get_spline_curvature_data(10000)
    true_max_curvature = np.max(curvature_data)
    if true_max_curvature < 1e-10 and curvature_bound < 1e-10:
        error_rate = 0
    elif true_max_curvature < 1e-10:
        error_rate = 100
    elif true_max_curvature == np.inf and curvature_bound == np.inf:
        error_rate = 0
    elif true_max_curvature == np.inf:
        error_rate = 100
    else:
        error_rate = np.abs(curvature_bound - true_max_curvature)/true_max_curvature*100
    accuracy = 100-error_rate
    return accuracy, evaluation_time

def create_control_points(isRestricted, order):
    if isRestricted:
        control_points = create_random_control_points_greater_than_angles(order+1,order,np.random.rand()*10)
    else:
        control_points = np.random.randint(10, size=(3,order+1)) #random
    return control_points

def test_performance_of_max_curvature_function(order,isRestricted,method,num_iterations):
    accuracy_array = np.zeros(num_iterations)
    evaluation_time_array = np.zeros(num_iterations)
    for i in range(num_iterations):
        accuracy, evaluation_time = get_speed_and_accuracy_of_max_curvature_function(order,isRestricted,method)
        accuracy_array[i] = accuracy
        evaluation_time_array[i] = evaluation_time
    average_accuracy = np.average(accuracy_array)
    average_time = np.average(evaluation_time_array)
    std_accuracy = np.sqrt(np.sum((average_accuracy - accuracy_array)**2)/num_iterations)
    std_time = np.sqrt(np.sum((average_time - evaluation_time_array)**2)/num_iterations)
    return average_accuracy, std_accuracy, average_time, std_time

def get_curvature_bound_and_time(method,control_points,order):
    if method == "maximize_curvature_equation":
        M = get_matrix(order)
        start_time = time.time()
        curvature_bound = find_max_curvature_sqp_method(control_points,order,M)
        evaluation_time = time.time() - start_time
    elif method == "roots_of_curvature_derivative":
        M = get_matrix(order)
        start_time = time.time()
        curvature_bound = find_max_curvature_root_finder(control_points,order,M)
        evaluation_time = time.time() - start_time
    elif method == "discrete_evaluations":
        n_points = 100
        M_vel, M_accel, L_vel, L_accel = get_matrices_discrete_evaluations(n_points,order)
        start_time = time.time()
        curvature_bound = get_max_curvature_by_checking_n_points(control_points, M_vel, M_accel, L_vel, L_accel)
        evaluation_time = time.time() - start_time
    elif method == "curvature_at_min_velocity":
        M = get_matrix(order)
        start_time = time.time()
        curvature_bound = find_curvature_at_min_velocity_magnitude(control_points, order, M)
        evaluation_time = time.time() - start_time
    elif method == "control_point_derivatives":
        pass
    elif method == "max_numerator_over_min_denominator":
        pass
    elif method == "":
        pass
    elif method == "angle_constrained_control_points":
        pass
    return curvature_bound, evaluation_time