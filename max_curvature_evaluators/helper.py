
import numpy as np
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.matrix_evaluation import get_T_derivative_vector, get_M_matrix

def get_matrix(order):
    M = get_M_matrix(0, order, [], False)
    return M

def calculate_velocity_magnitude(t,M,control_points):
    num_control_points = count_number_of_control_points(control_points)
    order = num_control_points - 1
    dT = get_T_derivative_vector(order,t,0,1,1)
    velocity = np.dot(control_points,np.dot(M,dT)).flatten()
    velocity_magnitude = np.linalg.norm(velocity)
    return velocity_magnitude 

def calculate_curvature(t,M,control_points):
    num_control_points = count_number_of_control_points(control_points)
    order = num_control_points - 1
    dT = get_T_derivative_vector(order,t,0,1,1)
    d2T = get_T_derivative_vector(order,t,0,2,1)
    velocity = np.dot(control_points,np.dot(M,dT)).flatten()
    acceleration = np.dot(control_points,np.dot(M,d2T)).flatten()
    numerator = np.linalg.norm(np.cross(velocity,acceleration))
    denominator = np.linalg.norm(velocity)**3
    curvature = numerator/denominator
    return curvature

def calculate_curvature_derivative(t,M,control_points):
    num_control_points = count_number_of_control_points(control_points)
    order = num_control_points - 1
    dT = get_T_derivative_vector(order,t,0,1,1)
    d2T = get_T_derivative_vector(order,t,0,2,1)
    d3T = get_T_derivative_vector(order,t,0,3,1)
    velocity = np.dot(control_points,np.dot(M,dT)).flatten()
    acceleration = np.dot(control_points,np.dot(M,d2T)).flatten()
    jerk = np.dot(control_points,np.dot(M,d3T)).flatten()
    vel_cross_accel = np.cross(velocity,acceleration)
    vel_cross_jerk = np.cross(velocity,jerk)
    vel_dot_accel = np.dot(velocity,acceleration)
    norm_vel_cross_accel = np.linalg.norm(vel_cross_accel)
    norm_velocity = np.linalg.norm(velocity)
    numerator_1 = np.dot(vel_cross_accel,vel_cross_jerk)
    denominator_1 = norm_vel_cross_accel*norm_velocity**3
    numerator_2 = 3*vel_dot_accel*norm_vel_cross_accel
    denominator_2 = norm_velocity**5
    curvature_derivative = numerator_1/denominator_1 - numerator_2/denominator_2
    return curvature_derivative