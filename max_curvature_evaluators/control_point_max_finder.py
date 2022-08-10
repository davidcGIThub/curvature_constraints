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

order = 4
dimension = 2
control_points = np.random.randint(10, size=(dimension,order+1)) # random
# print("control_points: " , control_points)

def get_control_point_curvature(control_points,scale_factor):
    control_point_velocities = (control_points[:,1:] - control_points[:,0:-1])/scale_factor
    control_point_accelerations = (control_point_velocities[:,1:]-control_point_velocities[:,0:-1])/scale_factor
    if count_number_of_control_points(control_point_velocities) > 1:
        control_point_velocities = convert_to_bezier_control_points(control_point_velocities)
    if count_number_of_control_points(control_point_accelerations) > 1:
        control_point_accelerations = convert_to_bezier_control_points(control_point_accelerations)
    max_acceleration = np.max(np.linalg.norm(control_point_accelerations,2,0))
    mdm = MDM(control_point_velocities,dimension,True)
    min_velocity = np.linalg.norm(mdm.solve())
    max_curvature = max_acceleration/min_velocity**2
    return max_curvature, max_acceleration, min_velocity

num_data_points = 1000
scale_factor = 1
bspline = BsplineEvaluation(control_points,order,0,scale_factor)
curvature , time_data = bspline.get_spline_curvature_data(1000)
velocity, time_data = bspline.get_derivative_magnitude_data(1000,1)
acceleration, time_data = bspline.get_derivative_magnitude_data(1000,2)
max_curvature, max_acceleration, min_velocity = get_control_point_curvature(control_points,scale_factor)

velocity_control_points = bspline.get_control_point_derivatives(1)
velocity_bspline = BsplineEvaluation(velocity_control_points,order-1,0,scale_factor,False)
bezeir_velocity_control_points = velocity_bspline.get_bezier_control_points()
print("bezeir_velocity_control_points: " , bezeir_velocity_control_points)
velocity_bspline.plot_bezier_curves(1000)
# plt.plot(time_data, curvature, label="curvature")
# plt.plot(time_data,time_data*0 + max_curvature, label="max_curvature")
# plt.legend()
# plt.show()

# plt.plot(time_data, velocity, label="velocity")
# plt.plot(time_data,time_data*0 + min_velocity, label="min_velocity")
# plt.legend()
# plt.show()

# plt.plot(time_data, acceleration, label="acceleration")
# plt.plot(time_data,time_data*0 + max_acceleration, label="max_acceleration")
# plt.legend()
# plt.show()

# total_time = 0
# iterations = 100
# average_accuracy = 0
# variance = 0
# for i in range(iterations):
#     # control_points = np.random.randint(10, size=(2,order+1))*1.0 # random
#     control_points = create_random_control_points_greater_than_angles(order+1,angle)
#     control_points = np.array([[1,2,6,7,8,14],[1,2,4,7,8,14]])
#     print("here")
#     start_time = time.time()
#     max_curvature, max_acceleration, min_velocity = get_control_point_curvatures_2(control_points,scale_factor)
#     end_time = time.time()
#     print("here2")
#     total_time += end_time - start_time
#     bspline = BsplineEvaluation(control_points,order,0,1)
#     curvature_data, time_data = bspline.get_spline_curvature_data(10000)
#     velocity_data, time_data = bspline.get_derivative_magnitude_data(10000,1)
#     acceleration_data, time_data = bspline.get_derivative_magnitude_data(10000,2)
#     true_max = np.max(curvature_data)
#     if true_max < 1e-10 and max_curvature < 1e-10:
#         accuracy = 1
#     elif true_max < 1e-10:
#         accuracy = 0
#     elif true_max == np.inf and max_curvature == np.inf:
#         accuracy = 1
#     elif true_max == np.inf:
#         accuracy = 0
#     else:
#         accuracy = max_curvature/true_max
#     average_accuracy += accuracy/iterations
#     variance += (1-accuracy)**2
#     plt.plot(time_data,curvature_data,label="curvature")
#     plt.plot(time_data, time_data*0+max_curvature, label="curvature_bound")
#     plt.show()
#     plt.plot(time_data,velocity_data,label="velocity")
#     plt.plot(time_data, time_data*0+min_velocity, label="vel_min_bound")
#     plt.legend()
#     plt.show()
#     # plt.plot(time_data,acceleration_data,label="acceleration")
#     # plt.plot(time_data, time_data*0+max_acceleration, label="accel_max_bound")
#     # plt.legend()
#     # plt.show()
#     print("accuracy: " , accuracy)

# ave_time = total_time/iterations
# standard_deviation = np.sqrt(variance/iterations)
# print("standard deviation: " , standard_deviation)
# print("ave_accuracy: " , average_accuracy)
# print("ave_time: " , ave_time)


