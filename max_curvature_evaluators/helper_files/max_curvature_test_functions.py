import numpy as np
import time
from sklearn.neighbors import KernelDensity
from bsplinegenerator.bsplines import BsplineEvaluation 
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import get_matrix
from max_curvature_evaluators.angle_constrained_max_finder import create_random_control_points_greater_than_angles
from max_curvature_evaluators.sqp_max_finder import find_max_curvature_sqp_method, get_m_matrix_sqp
from max_curvature_evaluators.root_finder import find_max_curvature_root_finder, get_m_matrix_root
from max_curvature_evaluators.discrete_evaluations import get_matrices_discrete_evaluations, get_max_curvature_by_checking_n_points
# from max_curvature_evaluators.discrete_evaluations_slow import get_matrices_discrete_evaluations, get_max_curvature_by_checking_n_points
from max_curvature_evaluators.max_at_min_velocity_finder import find_curvature_at_min_velocity_magnitude
from max_curvature_evaluators.control_point_method import get_control_point_curvature_bound
from max_curvature_evaluators.max_numerator_over_min_denominator import find_curvature_using_max_numerator_over_min_denominator
from max_curvature_evaluators.geometric_max_finder import get_max_curvature
from max_curvature_evaluators.angle_constrained_max_finder import get_constant, get_curvature_bound

def get_speed_and_accuracy_of_max_curvature_function(order,method,control_points, true_max_curvature):
    curvature_bound, evaluation_time = get_curvature_bound_and_time(method, control_points,order)
    relative_error = (curvature_bound - true_max_curvature)/(1+abs(true_max_curvature))
    if true_max_curvature == np.inf or curvature_bound == np.inf:
        relative_error = 0
    return relative_error, evaluation_time

def create_control_point_and_max_curvature_list(order, num_iterations, dimension=2):
    control_point_list = []
    max_curvature_list = []
    for i in range(num_iterations):
        control_points, max_curvature = get_control_points_and_max_curvature(order)
        control_point_list.append(control_points)
        max_curvature_list.append(max_curvature)
    return control_point_list, max_curvature_list

def create_edge_cases(order):
    if order == 2:
        control_point_list = create_2nd_order_edge_cases()
    elif order == 3:
        control_point_list = create_3rd_order_edge_cases()
    elif order == 5:
        control_point_list = create_5th_order_edge_cases()
    max_curvature_list = []
    for i in range(len(control_point_list)):
        control_points = control_point_list[i]
        bspline = BsplineEvaluation(control_points,order,0,1)
        curvature_data, time_data = bspline.get_spline_curvature_data(10000)
        max_curvature = np.max(curvature_data)
        max_curvature_list.append(max_curvature)
    return control_point_list, max_curvature_list
    
def create_2nd_order_edge_cases():
    # zero curvature, constant velocity
    points_1 = np.array([[0,1,2],[0,1,2],[0,1,2]])
    # zero curvature, with acceleration
    points_2 = np.array([[0,1,3],[0,1,3],[0,1,3]])
    # infinit curvature middle, hits zero velocity in middle
    points_3 = np.array([[0,2,1],[0,2,1],[0,2,1]])
    # zero curvature, hits zero velocity at endpoint
    points_4 = np.array([[1,1,2],[1,1,2],[1,1,2]])
    # zero velocity whole time, infinite curvature
    points_5 = np.array([[1,1,1],[1,1,1],[1,1,1]])
    # almost straight line
    points_6 = np.array([[0,1.01,2],[0,1,2.03],[0.01,1,2]])
    # almost straight line
    points_7 = np.array([[0.03,1.99,1.01],[0,2.02,1.07],[0.01,2.04,1]])
    return [points_1, points_2, points_3, points_4, points_5, points_6, points_7]

def create_3rd_order_edge_cases():
    # zero curvature, constant velocity
    points_1 = np.array([[0,1,2,3],[0,1,2,3],[0,1,2,3]])
    # zero curvature with some acceleration
    points_2 = np.array([[0,1,3,6],[0,1,3,6],[0,1,3,6]])
    # infinite curvature middle, hits zero velocity middle
    points_3 = np.array([[1,4,3,2],[1,4,3,2],[1,4,3,2]])
    # zero curvature, hits zero velocity at endpoint
    points_4 = np.array([[1,1,1,2],[1,1,1,2],[1,1,1,2]])
    # zero velocity whole time, infinite curvature
    points_5 = np.array([[1,1,1,1],[1,1,1,1],[1,1,1,1]])
    # almost straight line
    points_6 = np.array([[0,1.01,2.01,3.06],[0,1,2.03,3.01],[0.01,1,1.95,2.98]])
    # almost straight line, sharp turn around
    points_7 = np.array([[1,4.02,2.99,2],[1.05,3.94,2.96,2],[1.02,4,3.06,2.02]])
    return [points_1, points_2, points_3, points_4, points_5, points_6, points_7]

def create_5th_order_edge_cases():
    # zero curvature, constant velocity
    points_1 = np.array([[0,1,2,3,4,5],[0,1,2,3,4,5],[0,1,2,3,4,5]])
    # zero curvature, with some acceleration
    points_2 = np.array([[0,1,3,4,7,9],[0,1,3,4,7,9],[0,1,3,4,7,9]])
    # infinite curvature middle, hits zero velocity middle
    points_3 = np.array([[1,3,4,5,3,0],[1,3,4,5,3,0],[1,3,4,5,3,0]])
    # zero curvature, hits zero velocity at endpoint
    points_4 = np.array([[1,1,1,1,1,2],[1,1,1,1,1,2],[1,1,1,1,1,2]])
    # zero velocity whole time, infinite curvature
    points_5 = np.array([[1,1,1,1,1,1],[1,1,1,1,1,1],[1,1,1,1,1,1]])
    # almost straight line
    points_6 = np.array([[0,1.01,2,3.06,4.01,5],[0,1,2.03,3.01,4,5.03],[0.01,1,2,3,4.05,5.02]])
    # almost straight line, infinite curvature middle, hits zero velocity middle
    points_7 = np.array([[1.01,2.99,4.05,4.96,3.01,0],[1,3.01,4,5,3.06,0],[1,3.03,4,5.02,2.97,0.01]])
    return [points_1, points_2, points_3, points_4, points_5, points_6, points_7]

# # zero curvature, constant velocity
# control_points = np.array([[0,1,2,3],[0,1,2,3]]) # correct
# # zero curvature with some acceleration
# control_points = np.array([[0,1,3,6],[0,1,3,6]]) 
# # infinite curvature middle, hits zero velocity middle
# control_points = np.array([[1,4,3,2],[1,4,3,2]]) 
# # zero curvature, hits zero velocity at endpoint
# control_points = np.array([[1,1,1,2],[1,1,1,2]])
# # zero velocity whole time, infinite curvature
# control_points = np.array([[1,1,1,1],[1,1,1,1]])
# # almost straight line
# control_points = np.array([[0,1.01,2.03,3.06],[0,1,2.03,3.01],[0.01,1,2.02,3]])
# # almost straight, infinite curvature middle, hits zero velocity middle
# control_points = np.array([[1.01,3.98,3.02,2.04],[.99,4.03,2.97,2.02]]) 

def get_control_points_and_max_curvature(order,dimension=2):
    valid_control_points = False
    while not valid_control_points:
        control_points = np.random.randint(10, size=(dimension,order+1)) #random
        bspline = BsplineEvaluation(control_points,order,0,1)
        velocity_data, time_data = bspline.get_derivative_magnitude_data(1000,1)
        min_velocity = np.min(velocity_data)
        acceleration_data, time_data = bspline.get_derivative_magnitude_data(1000,1)
        min_acceleration = np.min(acceleration_data)
        curvature_data, time_data = bspline.get_spline_curvature_data(10000)
        max_curvature = np.max(curvature_data)
        cp_max_curvature = get_control_point_curvature_bound(control_points, order)
        relative_error_cp_max = np.abs(cp_max_curvature - max_curvature)/(1+np.abs(max_curvature))
        if min_velocity > 10e-5 and min_acceleration > 10e-5 and relative_error_cp_max < 10000:
            valid_control_points = True
    return control_points, max_curvature

def test_performance_of_max_curvature_function(order,method,control_point_list, max_curvature_list):
    num_iterations = len(control_point_list)
    relative_error_array = np.zeros(num_iterations)
    evaluation_time_array = np.zeros(num_iterations)
    for i in range(num_iterations):
        control_points = control_point_list[i]
        max_curvature = max_curvature_list[i]
        relative_error, evaluation_time = get_speed_and_accuracy_of_max_curvature_function(order,method,control_points, max_curvature)
        relative_error_array[i] = relative_error
        evaluation_time_array[i] = evaluation_time
    average_time = np.average(evaluation_time_array)
    high_error = np.max(relative_error_array)
    low_error = np.min(relative_error_array)
    std_time = np.sqrt(np.sum((average_time - evaluation_time_array)**2)/num_iterations)
    kde = KernelDensity(kernel='gaussian',bandwidth=0.5).fit(relative_error_array.reshape([-1,1]))
    samples = np.linspace(low_error,high_error,10000)
    scores = kde.score_samples(samples.reshape([-1,1]))
    index_highest_score = np.argmax(scores)
    mode = samples[index_highest_score]
    mean = np.mean(relative_error_array)
    std = np.std(relative_error_array)
    return mode, mean, std, relative_error_array, average_time, std_time

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
        n_points = 1000
        M_vel, M_accel, L_vel, L_accel = get_matrices_discrete_evaluations(n_points,order)
        start_time = time.time()
        curvature_bound = get_max_curvature_by_checking_n_points(control_points, M_vel, M_accel, L_vel, L_accel)
        evaluation_time = time.time() - start_time
    elif method == "curvature_at_min_velocity":
        M = get_matrix(order)
        start_time = time.time()
        curvature_bound = find_curvature_at_min_velocity_magnitude(control_points, order, M)
        evaluation_time = time.time() - start_time
    elif method == "max_numerator_over_min_denominator":
        M = get_matrix(order)
        start_time = time.time()
        curvature_bound = find_curvature_using_max_numerator_over_min_denominator(control_points, order, M)
        evaluation_time = time.time() - start_time
    elif method == "control_point_derivatives":
        start_time = time.time()
        curvature_bound = get_control_point_curvature_bound(control_points,order)
        evaluation_time = time.time() - start_time
    elif method == "geometric":
        start_time = time.time()
        curvature_bound = get_max_curvature(control_points)
        evaluation_time = time.time() - start_time
    elif method == "angle_constrained_control_points":
        isBezier = False
        constant = get_constant(order,isBezier)
        start_time = time.time()
        curvature_bound = get_curvature_bound(control_points,constant)
        evaluation_time = time.time() - start_time
    else:
        print( method , " not found")
    return curvature_bound, evaluation_time