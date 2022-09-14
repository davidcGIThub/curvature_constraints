import numpy as np
from max_curvature_evaluators.helper_files.cube_root_solver import solver
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import calculate_velocity_magnitude, calculate_curvature, get_matrix
import matplotlib.pyplot as plt
import time


def find_curvature_at_min_velocity_magnitude(control_points, order, M):
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
    roots = solver(A,B,C,D)
    times_to_check = roots
    t_min_velocity = 0
    min_velocity = np.inf
    for i in range(len(times_to_check)):
        t = times_to_check[i]
        if t >= 0 and t <= 1:
            velocity = calculate_velocity_magnitude(t, M, control_points,order)
            if velocity < min_velocity:
                t_min_velocity = t
                min_velocity = velocity
    max_curvature = calculate_curvature(t_min_velocity,M,control_points,order)
    curv_end = calculate_curvature(1,M,control_points,order)
    curv_start = calculate_curvature(0,M,control_points,order)
    max_curvature = np.max((curv_end,curv_start,max_curvature))
    return max_curvature






