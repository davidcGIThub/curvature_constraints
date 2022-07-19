import numpy as np
from cube_root_solver import analytical_solver, cube_root_plotter
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from helper import calculate_velocity_magnitude, calculate_curvature, get_matrix

order = 3
control_points = np.random.randint(10, size=(2,order+1)) # random
# control_points = np.array([[0, 1, 1, 4],[8, 2, 5, 1]])
# control_points = np.array([[9, 4, 2, 8],[2, 7, 9, 4]])


print("control_points: " , control_points)

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
    print("A: ", A)
    print("B: ", B)
    print("C: ", C)
    print("D: ", D)
    # cube_root_plotter(A,B,C,D)
    roots = analytical_solver(A,B,C,D)
    print("roots: " , roots)
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

estimated_max_curvature = find_curvature_at_min_velocity_magnitude(control_points)
print("estimated_max_curvature: ", estimated_max_curvature)
bspline = BsplineEvaluation(control_points,order,0,1)
curvature_data, time_data = bspline.get_spline_curvature_data(10000)
print("true_max_curvature: " , np.max(curvature_data))
# bspline.plot_derivative_magnitude(1000,1)
bspline.plot_curvature(1000)





