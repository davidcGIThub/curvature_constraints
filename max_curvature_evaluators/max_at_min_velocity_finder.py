import numpy as np
from cube_root_solver import analytical_solver, cube_root_plotter
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from helper import calculate_velocity_magnitude, calculate_curvature, get_matrix
import matplotlib.pyplot as plt
import time

order = 3
control_points = np.random.randint(10, size=(2,order+1)) # random
# control_points = np.array([[0, 1, 1, 4],[8, 5, 5, 1]])
# control_points = np.array([[9, 9, 2, 8],[2, 2, 9, 4]])
# control_points = np.array([[9, 9, 2, 2],[2, 2, 9, 9]])
# control_points = np.array([[1, 9, 1, 2],[2, 4, 2, 9]])
# control_points = np.array([[0, 2, 2, 0],[1, 3, 2, 8]])
# control_points = np.array([[9, 8, 9, 6],[9, 8, 4, 0]])
# control_points = np.array([[0, 1, 9, 8],[5, 6, 6, 0]])
control_points = np.array( [[3, 7, 8, 8], [2, 6, 5, 9]])

print("control_points: " , control_points)

def find_time_at_min_velocity(control_points):
    num_control_points = count_number_of_control_points(control_points)
    order = num_control_points -1
    M = get_matrix(order)
    P = control_points
    J = np.dot( np.dot(M.T,P.T) , np.dot(P,M))
    if order == 1:
        A = 0
        B = 0
        C = 0
        D = 0
    elif order == 2:
        A = 0
        B = 0
        C = 8*J[0,0]
        D = 4*J[1,0]
    elif order == 3:
        A = 36*J[0,0]
        B = 12*J[0,1] + 24*J[1,0]
        C = 8*J[1,1] + 12*J[2,0]
        D = 4*J[2,1]
    else:
        raise Exception("Function only capable of 1-3rd order splines")
    # cube_root_plotter(A,B,C,D)
    roots = analytical_solver(A,B,C,D)
    t_array = roots + [0, 1]
    t_min_velocity = 0
    min_velocity = np.inf
    for i in range(len(t_array)):
        t = t_array[i]
        if t >= 0 and t <= 1:
            velocity = calculate_velocity_magnitude(t, M, control_points)
            if velocity < min_velocity:
                t_min_velocity = t
                min_velocity = velocity
    return t_min_velocity

def find_curvature_bound(control_points):
    t_min_vel = find_time_at_min_velocity(control_points)
    order = 3
    bspline = BsplineEvaluation(control_points,order,0,1)
    min_vel = bspline.get_derivative_magnitude_at_time(t_min_vel,1)
    if t_min_vel < 0.5:
        t_accel = 0
    else:
        t_accel = 1
    max_accel_mag = bspline.get_derivative_magnitude_at_time(t_accel,2)
    return max_accel_mag/min_vel**2

def find_centr_accel_bound(control_points):
    bspline = BsplineEvaluation(control_points,3,0,1)
    accel_data, time_data = bspline.get_derivative_magnitude_data(10000,2)
    vel_data, time_data = bspline.get_derivative_magnitude_data(10000,1)
    t_accel = time_data[np.argmax(accel_data)]
    t_vel = time_data[np.argmin(vel_data)]
    centr_accel_1 = bspline.get_centripetal_acceleration_at_time_t(t_accel)
    centr_accel_2 = bspline.get_centripetal_acceleration_at_time_t(t_vel)
    return np.max((centr_accel_1,centr_accel_2))

def find_ang_rate_bound(control_points):
    bspline = BsplineEvaluation(control_points,3,0,1)
    accel_data, time_data = bspline.get_derivative_magnitude_data(10000,2)
    vel_data, time_data = bspline.get_derivative_magnitude_data(10000,1)
    t_accel = time_data[np.argmax(accel_data)]
    t_vel = time_data[np.argmin(vel_data)]
    ang_rate_1 = bspline.get_angular_rate_at_time_t(t_accel)
    ang_rate_2 = bspline.get_angular_rate_at_time_t(t_vel)
    return np.max((ang_rate_1,ang_rate_2))

def find_curvature_at_min_velocity_magnitude(control_points):
    num_control_points = count_number_of_control_points(control_points)
    order = num_control_points -1
    M = get_matrix(order)
    P = control_points
    J = np.dot( np.dot(M.T,P.T) , np.dot(P,M))
    if order == 1:
        A = 0
        B = 0
        C = 0
        D = 0
    elif order == 2:
        A = 0
        B = 0
        C = 8*J[0,0]
        D = 4*J[1,0]
    elif order == 3:
        A = 36*J[0,0]
        B = 12*J[0,1] + 24*J[1,0]
        C = 8*J[1,1] + 12*J[2,0]
        D = 4*J[2,1]
    else:
        raise Exception("Function only capable of 1-3rd order splines")
    # cube_root_plotter(A,B,C,D)
    roots = analytical_solver(A,B,C,D)
    times_to_check = roots
    t_min_velocity = 0
    min_velocity = np.inf
    for i in range(len(times_to_check)):
        t = times_to_check[i]
        if t >= 0 and t <= 1:
            velocity = calculate_velocity_magnitude(t, M, control_points)
            if velocity < min_velocity:
                t_min_velocity = t
                min_velocity = velocity
    max_curvature = calculate_curvature(t_min_velocity,M,control_points)
    curv_end = calculate_curvature(1,M,control_points)
    curv_start = calculate_curvature(0,M,control_points)
    max_curvature = np.max((curv_end,curv_start,max_curvature))
    return max_curvature, t_min_velocity

total_time = 0
iterations = 100
average_accuracy = 0
# for i in range(iterations):
#     control_points = np.random.randint(10, size=(2,order+1)) # random
#     start_time = time.time()
#     max_curvature, t_min_velocity = find_curvature_at_min_velocity_magnitude(control_points)
#     end_time = time.time()
#     total_time += end_time - start_time
#     bspline = BsplineEvaluation(control_points,order,0,1)
#     curvature_data, time_data = bspline.get_spline_curvature_data(10000)
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
#     print("average_accuracy: " , average_accuracy)
# ave_time = total_time/iterations
# print("ave_accuracy: " , average_accuracy)
# print("ave_time: " , ave_time)

# curvature_bound = find_curvature_bound(control_points)
# estimated_max_curvature , t_min= find_curvature_at_min_velocity_magnitude(control_points)
# print("estimated_max_curvature: ", estimated_max_curvature)
# # print("t min vel: " , find_time_at_min_velocity(control_points))
bspline = BsplineEvaluation(control_points,order,0,1)
curvature, time_data = bspline.get_spline_curvature_data(10000)
velocity, time_data = bspline.get_derivative_magnitude_data(10000 , 1)
acceleration, time_data = bspline.get_derivative_magnitude_data(10000 , 2)
bspline.plot_derivative_magnitude(10000,2)

plt.plot(time_data[0:-1],velocity[1:] - velocity[0:-1],label="d/dt ||b(t)'||")
plt.legend()
plt.xlabel("time")
plt.ylabel("velocity magnitude derivative")
plt.title("Derivative of the Velocity Magnitude")
plt.show()

plt.plot(time_data[0:-1],acceleration[1:]-acceleration[0:-1],label="d/dt ||b(t)''||")
# plt.legend()
plt.xlabel("time")
plt.ylabel("acceleration magnitude")
plt.title("Accleration Magnitude")
plt.show()

plt.plot(time_data,velocity,label="velocity")
plt.plot(time_data,curvature,label="curvature")
plt.legend()
plt.xlabel("time")
plt.title("Velocity Magnitude and Curvature")
plt.show()
# curvature_data, time_data = bspline.get_spline_curvature_data(1000)
# print("true_max_curvature: " , np.max(curvature_data))
# centripetal_acceleration = bspline.get_centripetal_acceleration_at_time_t(t_min)
# angular_rate = bspline.get_angular_rate_at_time_t(t_min)
# print("angular_rate: " , angular_rate)
# centr_data, time_data = bspline.get_centripetal_acceleration_data(1000)
# ang_rate_data, time_data = bspline.get_angular_rate_data(1000)
# print("ang_rate_max: " , np.max(ang_rate_data) )
# print("ang_rate_bound: " , find_ang_rate_bound(control_points))
# # print("centripetal_acceleration: " , centripetal_acceleration)
# # print("true_centr_max: " , np.max(centr_data))
# # print("centrl accel bound: " , find_centr_accel_bound(control_points))
# # print("curvature_bound: " , curvature_bound)
# bspline.plot_spline(1000)
# bspline.plot_derivative_magnitude(1000,1)
# bspline.plot_curvature(1000)
# bspline.plot_centripetal_acceleration(1000)
# bspline.plot_derivative_magnitude(1000,1)





