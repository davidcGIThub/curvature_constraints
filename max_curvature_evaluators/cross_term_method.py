import numpy as np
from bsplinegenerator.helper_functions import count_number_of_control_points, get_dimension
from bsplinegenerator.bspline_to_bezier import convert_to_bezier_control_points
from bsplinegenerator.helper_functions import count_number_of_control_points
from bsplinegenerator.bsplines import BsplineEvaluation
import random
import matplotlib.pyplot as plt
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
    print("control_point_velocities: " , control_point_velocities)
    mdm = MDM(control_point_velocities,dimension,True)
    min_velocity = mdm.get_min_distance()
    print("min_velocity: " , min_velocity)
    print("max_acceleration: " , max_acceleration)
    if min_velocity == 0:
         max_curvature = 0
    else:
        max_curvature = max_acceleration/min_velocity**2
    return max_curvature

def get_control_point_curvature_cross_method(control_points,order,scale_factor):
    dimension = get_dimension(control_points)
    control_point_velocities = (control_points[:,1:] - control_points[:,0:-1])/scale_factor
    control_point_accelerations = (control_point_velocities[:,1:]-control_point_velocities[:,0:-1])/scale_factor
    if count_number_of_control_points(control_point_velocities) > 1:
        control_point_velocities = convert_to_bezier_control_points(control_point_velocities)
    if count_number_of_control_points(control_point_accelerations) > 1:
        control_point_accelerations = convert_to_bezier_control_points(control_point_accelerations)
    max_cross_term = get_cross_term_norm_bound(control_points,order,dimension)
    mdm = MDM(control_point_velocities,dimension,True)
    min_velocity = mdm.get_min_distance()
    if min_velocity == 0:
         max_curvature = 0
    else:
        max_curvature = max_cross_term/min_velocity**3
    return max_curvature

def get_cross_term_norm_bound(control_points, order, dimension):
    if order == 2:
        cross_term_norm_bound = get_cross_term_norm_value_from_second_order_spline(control_points)
    elif order == 3:
        cross_term_cps = get_cross_term_bezier_control_points_from_third_order_spline(control_points, dimension)
        if dimension == 2:
            cross_term_cp_norm = cross_term_cps
        else: #dimension == 3
            cross_term_cp_norm = np.linalg.norm(cross_term_cps,2,0)
        cross_term_norm_bound = np.max(cross_term_cp_norm)
    return cross_term_norm_bound

def get_cross_term_norm_value_from_second_order_spline(control_points):
    P_0 = control_points[:,0]
    P_1 = control_points[:,1]
    P_2 = control_points[:,2]
    cross_term_norm = np.linalg.norm(np.cross(P_0,P_1)+np.cross(P_1,P_2)+np.cross(P_2,P_0))
    return cross_term_norm

def get_cross_term_bezier_control_points_from_third_order_spline(control_points,dimension):
    P_0 = control_points[:,0][:,None]
    P_1 = control_points[:,1][:,None]
    P_2 = control_points[:,2][:,None]
    P_3 = control_points[:,3][:,None]
    p1 = P_0-2*P_1+P_2
    p2 = P_0/2-3*P_1/2+3*P_2/2-P_3/2
    p3 = P_0/2 - P_2/2
    if dimension == 3:
        Y1 = np.array([[0,1,0],[0,0,1],[1,0,0]])
        Y2 = np.array([[0,0,1],[1,0,0],[0,1,0]])
        Y3 = np.array([[0,0,1],[0,0,-1],[0,1,0]])
        Y4 = np.array([[0,1,0],[1,0,0],[1,0,0]])
    elif dimension == 2:
        Y1 = Y4 = np.array([1,0])
        Y2 = Y3 = np.array([0,1])
    c2 =  np.dot(Y1,p1) * np.dot(Y2,p2) - np.dot(Y2,p1) * np.dot(Y1,p2) \
        + np.dot(Y3,p1) * np.dot(Y4,2*p2) - np.dot(Y4,p1) * np.dot(Y3,2*p2)
    c1 =  np.dot(Y1,p3) * np.dot(Y2,2*p2) - np.dot(Y2,p3) * np.dot(Y1,2*p2)
    c0 =  np.dot(Y3,p3) * np.dot(Y4,p1) - np.dot(Y4,p3) * np.dot(Y3,p1)
    a0 = c0
    a1 = (2*a0+c1)/2
    a2 = 2*a1 - a0 + c2 
    if dimension == 2:
        cross_term_control_points = np.abs(np.array([a0,a1,a2])).flatten()
    else:
        cross_term_control_points = np.concatenate((a0,a1,a2),1)
    return cross_term_control_points

def get_cross_term_bezier_control_points_from_fourth_order_spline(control_points,dimension):
    P_0 = control_points[:,0][:,None]
    P_1 = control_points[:,1][:,None]
    P_2 = control_points[:,2][:,None]
    P_3 = control_points[:,3][:,None]
    P_4 = control_points[:,4][:,None]
    p1 = P_0-2*P_1+P_2
    p2 = P_0/2-3*P_1/2+3*P_2/2-P_3/2
    p3 = P_0/2 - P_2/2
    if dimension == 3:
        Y1 = np.array([[0,1,0],[0,0,1],[1,0,0]])
        Y2 = np.array([[0,0,1],[1,0,0],[0,1,0]])
        Y3 = np.array([[0,0,1],[0,0,-1],[0,1,0]])
        Y4 = np.array([[0,1,0],[1,0,0],[1,0,0]])
    elif dimension == 2:
        Y1 = Y4 = np.array([1,0])
        Y2 = Y3 = np.array([0,1])
    c_4 = np.dot(Y1, 24*P_0 - 72*P_1 + 72*P_2 - 24*P_3)*np.dot(Y2,4*P_0 - 16*P_1 + 24*P_2 - 16*P_3 + 4*P_4) \
         - np.dot(Y2,24*P_0 - 72*P_1 + 72*P_2 - 24*P_3)*np.dot(Y1,4*P_0 - 16*P_1 + 24*P_2 - 16*P_3 + 4*P_4) \
         + np.dot(Y3,12*P_0 - 36*P_1 + 36*P_2 - 12*P_3)*np.dot(Y4,12*P_0 - 48*P_1 + 72*P_2 - 48*P_3 + 12*P_4) \
         - np.dot(Y4,12*P_0 - 36*P_1 + 36*P_2 - 12*P_3)*np.dot(Y3,12*P_0 - 48*P_1 + 72*P_2 - 48*P_3 + 12*P_4)
    
    c_3 = np.dot(Y2,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3)*np.dot(Y1,4*P_0 - 16*P_1 + 24*P_2 - 16*P_3 + 4*P_4) \
        - np.dot(Y1,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3)*np.dot(Y2,4*P_0 - 16*P_1 + 24*P_2 - 16*P_3 + 4*P_4) \
        - np.dot(Y3,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3)*np.dot(Y4,12*P_0 - 48*P_1 + 72*P_2 - 48*P_3 + 12*P_4) \
        + np.dot(Y4,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3)*np.dot(Y3,12*P_0 - 48*P_1 + 72*P_2 - 48*P_3 + 12*P_4) \
        + np.dot(Y4,12*P_0 - 36*P_1 + 36*P_2 - 12*P_3)*np.dot(Y3,24*P_0 - 72*P_1 + 72*P_2 - 24*P_3) \
        - np.dot(Y3,12*P_0 - 36*P_1 + 36*P_2 - 12*P_3)*np.dot(Y4,24*P_0 - 72*P_1 + 72*P_2 - 24*P_3)
    
    c_2 = np.dot(Y2,4*P_0 + 12*P_1 - 12*P_2 - 4*P_3)*np.dot(Y1,12*P_0 - 48*P_1 + 72*P_2 - 48*P_3 + 12*P_4) \
        - np.dot(Y1,4*P_0 + 12*P_1 - 12*P_2 - 4*P_3)*np.dot(Y2,12*P_0 - 48*P_1 + 72*P_2 - 48*P_3 + 12*P_4) \
        + np.dot(Y4,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3)*np.dot(Y3,12*P_0 - 36*P_1 + 36*P_2 - 12*P_3) \
        - np.dot(Y3,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3)*np.dot(Y4,12*P_0 - 36*P_1 + 36*P_2 - 12*P_3) \
        - np.dot(Y4,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3)*np.dot(Y3,24*P_0 - 72*P_1 + 72*P_2 - 24*P_3) \
        + np.dot(Y3,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3)*np.dot(Y4,24*P_0 - 72*P_1 + 72*P_2 - 24*P_3)
    
    c_1 = np.dot(Y1, 4*P_0 + 12*P_1 - 12*P_2 - 4*P_3)*np.dot(Y2,24*P_0 - 72*P_1 + 72*P_2 - 24*P_3) \
         - np.dot(Y2,4*P_0 + 12*P_1 - 12*P_2 - 4*P_3)*np.dot(Y1,24*P_0 - 72*P_1 + 72*P_2 - 24*P_3)

    c_0 = - np.dot(Y4,4*P_0 + 12*P_1 - 12*P_2 - 4*P_3)*np.dot(Y3,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3) \
          + np.dot(Y4,12*P_0 - 12*P_1 - 12*P_2 + 12*P_3)*np.dot(Y3,4*P_0 + 12*P_1 - 12*P_2 - 4*P_3)
    # a0 = c0
    # a1 = (2*a0+c1)/2
    # a2 = 2*a1 - a0 + c2 
    if dimension == 2:
        cross_term_control_points = np.abs(np.array([a0,a1,a2])).flatten()
    else:
        cross_term_control_points = np.concatenate((a0,a1,a2),1)
    return cross_term_control_points

dimension = 2
num_data_points = 1000
order = 3
control_points = np.random.randint(10, size=(dimension,order+1)) # random
# control_points = np.array([[0,1,4,5],[0,1,4,5]])
# control_points = np.array([[0, 4, 2, 2],[4, 5, 7, 2]])
# control_points = np.array([[4, 1, 3, 3],[9, 3, 5, 0]])
control_points = np.array([[1, 4, 3, 8],[3, 4, 2, 1]])
print("control_points: " , control_points)
bspline = BsplineEvaluation(control_points,order,0,1)
acceleration, time_data = bspline.get_spline_derivative_data(num_data_points,2)
velocity, time_data = bspline.get_spline_derivative_data(num_data_points,1)
curvature, time_data = bspline.get_spline_curvature_data(num_data_points)
cross_term = np.cross(velocity.T,acceleration.T).T
if dimension == 2:
    cross_term_norm = np.abs(cross_term)
elif dimension == 3:
    cross_term_norm = np.linalg.norm(cross_term,2,0)


curvature_bound_cross = get_control_point_curvature_cross_method(control_points,order,1)
curvature_bound_control = get_control_point_curvature(control_points,1)
plt.plot(time_data,curvature,label="curvature")
plt.plot(time_data,time_data*0+curvature_bound_cross,label="bound cross")
plt.plot(time_data,time_data*0+curvature_bound_control,label="bound control")
plt.legend()
plt.show()

# bspline.plot_spline(num_data_points)
# bspline.plot_derivative_magnitude(num_data_points,1)
# bspline.plot_derivative_magnitude(num_data_points,2)
# cross_term_norm_bound = get_cross_term_norm_bound(control_points, order, dimension)

# plt.plot(time_data,cross_term_norm,label="cross term norm")
# plt.plot(time_data,time_data*0+cross_term_norm_bound,label="bound")
# plt.show()
# cross_term_cps = get_cross_term_bezier_control_points_from_third_order_spline(control_points,dimension)
# print("cross_term_cps: " , cross_term_cps)
# if dimension == 3:
#     plt.figure()
#     ax = plt.axes(projection='3d')
#     ax.set_box_aspect(aspect =(1,1,1))
#     ax.plot(cross_term[0,:], cross_term[1,:],cross_term[2,:],label="cross_term")
#     ax.scatter(cross_term_cps[0,:], cross_term_cps[1,:],cross_term_cps[2,:],label="cross_term_cps")
#     plt.show()
# elif dimension == 2:
#     plt.plot(time_data, cross_term)
#     plt.scatter(np.linspace(0,1,2*order-3), cross_term_cps)
#     plt.show()

