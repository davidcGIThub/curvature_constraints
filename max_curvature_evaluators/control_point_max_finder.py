from email.utils import decode_rfc2231
from matplotlib import scale
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.optimize import minimize, Bounds
from yaml import scan
from helper import calculate_curvature, calculate_curvature_derivative
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.matrix_evaluation import get_M_matrix
from bsplinegenerator.bspline_to_bezier import convert_to_bezier_control_points
from beziercurvegenerator.bezier_curves import BezierCurve
from bsplinegenerator.helper_functions import count_number_of_control_points
from mdm_algorithm_adapted import MDM

order = 3
scale_factor = 1
control_points = np.random.randint(10, size=(2,order+1))*1.0 # random
if order == 2:
    angle = np.pi/2
if order == 3:
    angle = np.pi/2
elif order == 4:
    angle = np.pi/3
elif order == 5:
    angle = np.pi/4

# angle = np.pi/2

def create_random_control_points_greater_than_angles(num_control_points,angle):
    control_points = np.zeros((2,num_control_points))
    len_ave = 3
    for i in range(num_control_points):
        if i == 0:
            control_points[:,i][:,None] = np.array([[0],[0]])
        elif i == 1:
            random_vec = np.random.rand(2,1)
            next_vec = len_ave*random_vec/np.linalg.norm(random_vec)
            control_points[:,i][:,None] = control_points[:,i-1][:,None] + next_vec
        else:
            new_angle = angle*2*(0.5-np.random.rand())
            R = np.array([[np.cos(new_angle), -np.sin(new_angle)],[np.sin(new_angle), np.cos(new_angle)]])
            prev_vec = control_points[:,i-1][:,None] - control_points[:,i-2][:,None]
            unit_prev_vec = prev_vec/np.linalg.norm(prev_vec)
            next_vec = len_ave*np.dot(R,unit_prev_vec)#*np.random.rand()
            control_points[:,i][:,None] = control_points[:,i-1][:,None] + next_vec
    return control_points

# control_points = create_random_control_points_greater_than_angles(order+1,angle)
# control_points = np.array([[1,3,4,8],[1,3,4,8]])

def get_inscribed_points(control_points,order):
    M = get_M_matrix(0, order, [], False)
    times = np.linspace(0,1,order+1)
    L_vectors = np.ones((order+1,order+1))
    for i in range(order):
        L_vectors[i,:] = times**(order-i)
    inscribed_points = np.dot(control_points,np.dot(M,L_vectors))
    return inscribed_points

def get_inscribed_velocities(control_points,scale_factor):
    order = count_number_of_control_points(control_points)-1
    inscribed_points = get_inscribed_points(control_points,order)
    inscribed_velocities = (inscribed_points[:,1:]-inscribed_points[:,0:-1])/(scale_factor/order)
    return inscribed_velocities

def get_control_point_curvatures_1(control_points,scale_factor):
    order = count_number_of_control_points(control_points)-1
    inscribed_points = get_inscribed_points(control_points,order)
    inscribed_velocities = (inscribed_points[:,1:]-inscribed_points[:,0:-1])/(scale_factor/order)
    control_point_velocities = (control_points[:,1:] - control_points[:,0:-1])/scale_factor
    control_point_acceleration = (control_point_velocities[:,1:]-control_point_velocities[:,0:-1])/scale_factor
    control_point_acceleration_norm = np.linalg.norm(control_point_acceleration,2,0)
    inscribed_velocity_norm = np.linalg.norm(inscribed_velocities,2,0)
    curvatures_1 = control_point_acceleration_norm/inscribed_velocity_norm[:-1]**2
    curvatures_2 = control_point_acceleration_norm/inscribed_velocity_norm[1:]**2
    curvatures = np.max((curvatures_1,curvatures_2),0)
    time_cp = np.linspace(0,scale_factor,len(curvatures))
    return curvatures,time_cp

def get_control_point_curvatures_2(control_points,scale_factor):
    dimension = 2
    order = count_number_of_control_points(control_points)-1
    control_point_velocities = (control_points[:,1:] - control_points[:,0:-1])/scale_factor
    control_point_acceleration = (control_point_velocities[:,1:]-control_point_velocities[:,0:-1])/scale_factor
    control_point_acceleration_norm = np.linalg.norm(control_point_acceleration,2,0)
    max_acceleration = np.max(control_point_acceleration_norm)
    b_spline_vel_control_points = convert_to_bezier_control_points(control_point_velocities)
    mdm = MDM(np.transpose(b_spline_vel_control_points),dimension,True)
    result = mdm.solve()
    min_velocity = np.linalg.norm(result)
    max_curvature = max_acceleration/min_velocity**2
    return max_curvature, max_acceleration, min_velocity

num_data_points = 1000
bspline = BsplineEvaluation(control_points,order,0,scale_factor)
curvature , time_data = bspline.get_spline_curvature_data(1000)
curv_cp_1, time_cp_curv = get_control_point_curvatures_1(control_points,scale_factor)


total_time = 0
iterations = 100
average_accuracy = 0
variance = 0
for i in range(iterations):
    # control_points = np.random.randint(10, size=(2,order+1))*1.0 # random
    control_points = create_random_control_points_greater_than_angles(order+1,angle)
    print("here")
    start_time = time.time()
    max_curvature, max_acceleration, min_velocity = get_control_point_curvatures_2(control_points,scale_factor)
    end_time = time.time()
    print("here2")
    total_time += end_time - start_time
    bspline = BsplineEvaluation(control_points,order,0,1)
    curvature_data, time_data = bspline.get_spline_curvature_data(10000)
    velocity_data, time_data = bspline.get_derivative_magnitude_data(10000,1)
    acceleration_data, time_data = bspline.get_derivative_magnitude_data(10000,2)
    true_max = np.max(curvature_data)
    if true_max < 1e-10 and max_curvature < 1e-10:
        accuracy = 1
    elif true_max < 1e-10:
        accuracy = 0
    elif true_max == np.inf and max_curvature == np.inf:
        accuracy = 1
    elif true_max == np.inf:
        accuracy = 0
    else:
        accuracy = max_curvature/true_max
    average_accuracy += accuracy/iterations
    variance += (1-accuracy)**2
    plt.plot(time_data,curvature_data,label="curvature")
    plt.plot(time_data, time_data*0+max_curvature, label="curvature_bound")
    plt.show()
    plt.plot(time_data,velocity_data,label="velocity")
    plt.plot(time_data, time_data*0+min_velocity, label="vel_min_bound")
    plt.legend()
    plt.show()
    # plt.plot(time_data,acceleration_data,label="acceleration")
    # plt.plot(time_data, time_data*0+max_acceleration, label="accel_max_bound")
    # plt.legend()
    # plt.show()
    print("accuracy: " , accuracy)

ave_time = total_time/iterations
standard_deviation = np.sqrt(variance/iterations)
print("standard deviation: " , standard_deviation)
print("ave_accuracy: " , average_accuracy)
print("ave_time: " , ave_time)


