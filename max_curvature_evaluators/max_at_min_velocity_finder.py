import numpy as np
from cube_root_solver import analytical_solver, cube_root_plotter
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from helper import calculate_velocity_magnitude, calculate_curvature, get_matrix

order = 3
control_points = np.random.randint(10, size=(2,order+1)) # random
# control_points = np.array([[0, 1, 1, 4],[8, 2, 5, 1]])
# control_points = np.array([[9, 4, 2, 8],[2, 7, 9, 4]])
control_points = np.array([[0, 2, 2, 0],[1, 3, 2, 8]])

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
    min_vel = bspline.get_derivative_at_time_t(t_min_vel,1).flatten()
    if t_min_vel < 0.5:
        t_accel = 0
    else:
        t_accel = 1
    max_accel = bspline.get_derivative_at_time_t(t_accel,2).flatten()
    curvature_bound = np.linalg.norm(np.cross(min_vel,max_accel))/np.linalg.norm(min_vel)**3
    return curvature_bound

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
    t_min_velocity = 0
    min_velocity = np.inf
    for i in range(len(roots)):
        t = roots[i]
        if t >= 0 and t <= 1:
            velocity = calculate_velocity_magnitude(t, M, control_points)
            if velocity < min_velocity:
                t_min_velocity = t
                min_velocity = velocity
    curvature_start = calculate_curvature(0,M,control_points)
    curvature_end = calculate_curvature(1,M,control_points)
    curvature_t_min_veloclity = calculate_curvature(t_min_velocity,M,control_points)
    max_curvature = np.max([curvature_start, curvature_end, curvature_t_min_veloclity])
    return max_curvature

curvature_bound = find_curvature_bound(control_points)
estimated_max_curvature = find_curvature_at_min_velocity_magnitude(control_points)
print("estimated_max_curvature: ", estimated_max_curvature)
print("t min vel: " , find_time_at_min_velocity(control_points))
bspline = BsplineEvaluation(control_points,order,0,1)
curvature_data, time_data = bspline.get_spline_curvature_data(1000)
print("true_max_curvature: " , np.max(curvature_data))
print("curvature_bound: " , curvature_bound)
bspline.plot_spline(1000)
bspline.plot_derivative_magnitude(1000,1)
bspline.plot_derivative_magnitude(1000,2)
bspline.plot_curvature(1000)





