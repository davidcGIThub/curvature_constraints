from curses.ascii import BS
import numpy as np
import random
from bsplinegenerator.matrix_evaluation import get_M_matrix
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.bsplines import BsplineEvaluation

n_points = 10
scale_factor = 1
dimension = 3
control_points = np.random.randint(10, size=(2, random.randint(4, 6)) ) # random
num_control_points = count_number_of_control_points(control_points)
order = num_control_points - 1
vel_order = order - 1
accel_order = order - 2
L_vel = np.ones((vel_order+1,n_points))
L_accel = np.ones((accel_order+1,n_points))
steps_array = np.linspace(0,1,n_points)
for i in range(vel_order+1):
    L_vel[i,:] = steps_array**(vel_order-i)
for i in range(accel_order+1):
    L_accel[i,:] = steps_array**(accel_order-i)
M_vel = get_M_matrix(0, vel_order, [], False)
M_accel = get_M_matrix(0, accel_order, [], False)

def get_max_curvature_by_checking_n_points(control_points, M_vel, M_accel, L_vel, L_accel, scale_factor):
    control_point_velocities = (control_points[:,1:] - control_points[:,0:-1])/scale_factor
    control_point_accelerations = (control_point_velocities[:,1:]-control_point_velocities[:,0:-1])/scale_factor
    velocity_data = np.dot(np.dot(control_point_velocities,M_vel),L_vel)
    acceleration_data = np.dot(np.dot(control_point_accelerations,M_accel),L_accel)
    velocity_magnitude = np.linalg.norm(velocity_data,2,0)
    curvature = np.abs(np.cross(velocity_data.T, acceleration_data.T).T) / velocity_magnitude**3
    # curvature = np.linalg.norm(np.cross(velocity_data.T, acceleration_data.T).T,2,0) / velocity_magnitude**3
    return np.max(curvature)



