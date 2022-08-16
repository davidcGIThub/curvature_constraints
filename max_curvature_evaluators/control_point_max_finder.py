import numpy as np
from bsplinegenerator.helper_functions import count_number_of_control_points, get_dimension
from bsplinegenerator.bspline_to_bezier import convert_to_bezier_control_points
from bsplinegenerator.helper_functions import count_number_of_control_points
from mdm_algorithm_adapted import MDM

def get_control_point_curvature(control_points,scale_factor):
    dimension = get_dimension(control_points)
    control_point_velocities = (control_points[:,1:] - control_points[:,0:-1])/scale_factor
    control_point_accelerations = (control_point_velocities[:,1:]-control_point_velocities[:,0:-1])/scale_factor
    if count_number_of_control_points(control_point_velocities) > 1:
        control_point_velocities = convert_to_bezier_control_points(control_point_velocities)
    if count_number_of_control_points(control_point_accelerations) > 1:
        control_point_accelerations = convert_to_bezier_control_points(control_point_accelerations)
    max_acceleration = np.max(np.linalg.norm(control_point_accelerations,2,0))
    mdm = MDM(control_point_velocities,dimension,True)
    min_velocity = mdm.get_min_distance()
    if min_velocity == 0:
         max_curvature = 0
    else:
        max_curvature = max_acceleration/min_velocity**2
    return max_curvature
