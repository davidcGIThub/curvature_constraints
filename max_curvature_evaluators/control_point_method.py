import numpy as np
from bsplinegenerator.helper_functions import count_number_of_control_points, get_dimension
from bsplinegenerator.bspline_to_bezier import convert_to_bezier_control_points
from bsplinegenerator.bspline_to_minvo import bezier_to_minvo_control_points
from bsplinegenerator.bspline_to_minvo import convert_to_minvo_control_points
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.bsplines import BsplineEvaluation
import random
import matplotlib.pyplot as plt
from max_curvature_evaluators.helper_files.mdm_algorithm_adapted import MDM
from max_curvature_evaluators.helper_files.cross_term_control_points import \
    get_cross_term_norm_value_from_second_order_spline, get_cross_term_bezier_control_points_from_third_order_spline,\
    get_cross_term_bezier_control_points_from_fourth_order_spline, get_cross_term_bezier_control_points_from_fifth_order_spline

def get_control_point_curvature_bound(control_points,order,scale_factor=1,method="mdm",rotation_angle = None):
    # print("control_points: " , control_points)
    dimension = get_dimension(control_points)
    control_point_velocities = (control_points[:,1:] - control_points[:,0:-1])/scale_factor
    control_point_accelerations = (control_point_velocities[:,1:]-control_point_velocities[:,0:-1])/scale_factor
    if count_number_of_control_points(control_point_velocities) > 1:
        bezier_control_point_velocities = convert_to_bezier_control_points(control_point_velocities)
        # minvo_control_point_velocities = convert_to_minvo_control_points(control_point_velocities)
    if count_number_of_control_points(control_point_accelerations) > 1:
        control_point_accelerations = convert_to_bezier_control_points(control_point_accelerations)
    max_cross_term = get_cross_term_norm_bound(control_points,order,dimension)
    max_acceleration = np.max(np.linalg.norm(control_point_accelerations,2,0))
    if method == "mdm":
        mdm = MDM(bezier_control_point_velocities,dimension,500,1)
        # mdm = MDM(minvo_control_point_velocities,dimension,500,1)
        min_velocity = mdm.get_min_distance()
    if method == "rotate":
        min_velocity = get_min_velocity_rotation_method(bezier_control_point_velocities,rotation_angle)
        # min_velocity = get_min_velocity_rotation_method(minvo_control_point_velocities,rotation_angle)
    if min_velocity < 1e-10:
        # print("control_point_velocities: " , control_point_velocities)
        # print("sign of vel: " , np.sign(control_point_velocities))
        if np.any(control_point_velocities>0) and np.any(control_point_velocities<0):
            max_curvature = np.inf
            # print("case 1")
        else:
            max_curvature = 0
            # print("case 2")
    else:
        max_curvature_cross_term_method = max_cross_term/min_velocity**3
        max_curvature_acceleration_method = max_acceleration/min_velocity**2
        max_curvature = np.min([max_curvature_acceleration_method, \
                                max_curvature_cross_term_method])
    return max_curvature

def get_min_velocity_rotation_method(control_point_velocities,rotation_angle):
    # mean_points = np.mean(control_point_velocities,1)
    # dx = mean_points.item(0)
    # dy = mean_points.item(1)
    # psi = -np.arctan2(dy,dx)
    psi = rotation_angle
    c_psi = np.cos(psi)
    s_psi = np.sin(psi)
    rotation = np.array([[c_psi, -s_psi], 
                         [s_psi,  c_psi]])
    rotated_points = np.dot(rotation, control_point_velocities)
    min_velocity = np.max( [np.min(rotated_points[0,:]), 0])
    return min_velocity


def get_cross_term_norm_bound(control_points, order, dimension):
    if order == 2:
        cross_term_norm_bound = get_cross_term_norm_value_from_second_order_spline(control_points)
    else:
        if order == 3:
            bez_cross_term_cps = get_cross_term_bezier_control_points_from_third_order_spline(control_points, dimension)
        elif order == 4:
            bez_cross_term_cps = get_cross_term_bezier_control_points_from_fourth_order_spline(control_points, dimension)
        elif order == 5:
            bez_cross_term_cps = get_cross_term_bezier_control_points_from_fifth_order_spline(control_points, dimension)
        if dimension == 2:
            cross_term_norm_bound = np.max(np.abs(bez_cross_term_cps))
        else:
            cross_term_norm_bound = np.max(np.linalg.norm(bez_cross_term_cps,2,0))
    return cross_term_norm_bound
