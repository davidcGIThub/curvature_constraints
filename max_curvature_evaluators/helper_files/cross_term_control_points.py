import numpy as np
from bsplinegenerator.bspline_to_minvo import bezier_to_minvo_control_points

def get_bezier_cross_term_norm_bound(control_points, order, dimension):
    if order == 2:
        cross_term_norm_bound = get_cross_term_norm_value_from_second_order_spline(control_points)
    elif order == 3:
        bez_cross_term_cps = get_cross_term_bezier_control_points_from_third_order_spline(control_points, dimension)
    elif order == 4:
        bez_cross_term_cps = get_cross_term_bezier_control_points_from_fourth_order_spline(control_points, dimension)
    elif order == 5:
        bez_cross_term_cps = get_cross_term_bezier_control_points_from_fifth_order_spline(control_points, dimension)
    if dimension == 2:
        bez_max_cross_term_norm = np.max(np.abs(bez_cross_term_cps))
    else:
        bez_max_cross_term_norm = np.max(np.linalg.norm(bez_cross_term_cps,2,0))
    return bez_max_cross_term_norm

def get_minvo_cross_term_norm_bound(control_points, order, dimension):
    if order == 2:
        cross_term_norm_bound = get_cross_term_norm_value_from_second_order_spline(control_points)
    elif order == 3:
        minv_cross_term_cps = get_cross_term_minvo_control_points_from_third_order_spline(control_points, dimension)
    elif order == 4:
        minv_cross_term_cps = get_cross_term_minvo_control_points_from_fourth_order_spline(control_points, dimension)
    elif order == 5:
        minv_cross_term_cps = get_cross_term_minvo_control_points_from_fifth_order_spline(control_points, dimension)
    if dimension == 2:
        cross_term_norm_bound = np.max(np.abs(minv_cross_term_cps))
    else:
        cross_term_norm_bound = np.max(np.linalg.norm(minv_cross_term_cps,2,0))
    return cross_term_norm_bound

def get_cross_term_norm_value_from_second_order_spline(control_points):
    P_0 = control_points[:,0]
    P_1 = control_points[:,1]
    P_2 = control_points[:,2]
    cross_term_norm = np.linalg.norm(np.cross(P_0,P_1)+np.cross(P_1,P_2)+np.cross(P_2,P_0))
    return cross_term_norm

def get_cross_term_minvo_control_points_from_third_order_spline(control_points,dimension):
    bezier_cross_term_control_points = get_cross_term_bezier_control_points_from_third_order_spline(control_points,dimension)
    minvo_cross_term_control_points = bezier_to_minvo_control_points(bezier_cross_term_control_points)
    return minvo_cross_term_control_points

def get_cross_term_minvo_control_points_from_fourth_order_spline(control_points,dimension):
    bezier_cross_term_control_points = get_cross_term_bezier_control_points_from_fourth_order_spline(control_points,dimension)
    minvo_cross_term_control_points = bezier_to_minvo_control_points(bezier_cross_term_control_points)
    return minvo_cross_term_control_points

def get_cross_term_minvo_control_points_from_fifth_order_spline(control_points,dimension):
    bezier_cross_term_control_points = get_cross_term_bezier_control_points_from_fifth_order_spline(control_points,dimension)
    minvo_cross_term_control_points = bezier_to_minvo_control_points(bezier_cross_term_control_points)
    return minvo_cross_term_control_points

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
        cross_term_control_points = np.array([a0,a1,a2]).flatten()
    else:
        cross_term_control_points = np.concatenate((a0,a1,a2),1)
    return cross_term_control_points

def get_cross_term_bezier_control_points_from_fifth_order_spline(control_points,dimension):
    P_0 = control_points[:,0][:,None]
    P_1 = control_points[:,1][:,None]
    P_2 = control_points[:,2][:,None]
    P_3 = control_points[:,3][:,None]
    P_4 = control_points[:,4][:,None]
    P_5 = control_points[:,5][:,None]
    p1 = P_0/4  - P_1/2      + P_3/2      - P_4/4
    p2 = P_0/24 - (5*P_1)/24 + (5*P_2)/12 - (5*P_3)/12 + (5*P_4)/24 - P_5/24
    p3 = P_0/6  - (2*P_1)/3  + P_2        - (2*P_3)/3  + P_4/6
    p4 = P_0/6  + P_1/3      - P_2        + P_3/3      + P_4/6
    p5 = P_0/24 + (5*P_1)/12 - (5*P_3)/12 - P_4/24
    if dimension == 3:
        Y1 = np.array([[0,1,0],[0,0,1],[1,0,0]])
        Y2 = np.array([[0,0,1],[1,0,0],[0,1,0]])
        Y3 = np.array([[0,0,1],[0,0,-1],[0,1,0]])
        Y4 = np.array([[0,1,0],[1,0,0],[1,0,0]])
    elif dimension == 2:
        Y1 = Y4 = np.array([1,0])
        Y2 = Y3 = np.array([0,1])
    c_6 = np.dot(Y1,3*p3)*np.dot(Y2,p2) - np.dot(Y2,3*p3)*np.dot(Y1,p2) + \
          np.dot(Y3,p3)*np.dot(Y4,4*p2) - np.dot(Y4,p3)*np.dot(Y3,4*p2)
    c_5 = np.dot(Y1,p1)*np.dot(Y2,4*p2) - np.dot(Y2,p1)*np.dot(Y1,4*p2) + \
          np.dot(Y4,p3)*np.dot(Y3,3*p3) - np.dot(Y3,p3)*np.dot(Y4,3*p3) + \
          np.dot(Y3,2*p1)*np.dot(Y4,p2) - np.dot(Y4,2*p1)*np.dot(Y3,p2)
    c_4 = np.dot(Y2,p4)*np.dot(Y1,4*p2) - np.dot(Y1,p4)*np.dot(Y2,4*p2) - \
          np.dot(Y3,p4)*np.dot(Y4,p2) + np.dot(Y4,p4)*np.dot(Y3,p2) + \
          np.dot(Y3,p1)*np.dot(Y4,3*p3) - np.dot(Y4,p1)*np.dot(Y3,3*p3) - \
          np.dot(Y3,2*p1)*np.dot(Y4,p3) + np.dot(Y4,2*p1)*np.dot(Y3,p3)
    c_3 = np.dot(Y4,p4)*np.dot(Y3,3*p3) - np.dot(Y4,3*p3)*np.dot(Y3,p4) - \
          np.dot(Y3,p5)*np.dot(Y4,4*p2) + np.dot(Y4,p5)  *np.dot(Y3,4*p2) - \
          np.dot(Y4,2*p1)*np.dot(Y3,p1) + np.dot(Y3,2*p1)*np.dot(Y4,p1) + \
          np.dot(Y4,p3)*np.dot(Y3,p4) - np.dot(Y3,p3)  *np.dot(Y4,p4)
    c_2 = np.dot(Y1,2*p1)*np.dot(Y2,p4) - np.dot(Y2,2*p1)*np.dot(Y1,p4) + \
          np.dot(Y3,p1)*np.dot(Y4,p4) - np.dot(Y4,p1)*np.dot(Y3,p4) + \
          np.dot(Y3,p5)*np.dot(Y4,3*p3) - np.dot(Y4,p5)*np.dot(Y3,3*p3)
    c_1 = np.dot(Y2,2*p1)*np.dot(Y1,p5) - np.dot(Y1,2*p1)*np.dot(Y2,p5)
    c_0 = np.dot(Y3,p5)*np.dot(Y4,p4) - np.dot(Y4,p5)*np.dot(Y3,p4)
    a_0 = c_0
    a_1 = (c_1 + 6*a_0)/6
    a_2 = (c_2 - 15*a_0 + 30*a_1)/15
    a_3 = (c_3 + 60*a_2 + 20*a_0 - 60*a_1)/20
    a_4 = (c_4 - 15*a_0 + 60*a_1 - 90*a_2 + 60*a_3)/15
    a_5 = (c_5 - 30*a_1 + 6*a_0  + 60*a_2 - 60*a_3 + 30*a_4)/6
    a_6 =  c_6 - a_0    + 6*a_1  - 15*a_2 + 20*a_3 - 15*a_4 + 6*a_5
    if dimension == 2:
        cross_term_control_points = np.array([a_0,a_1,a_2,a_3,a_4,a_5,a_6]).flatten()
    else:
        cross_term_control_points = np.concatenate((a_0,a_1,a_2,a_3,a_4,a_5,a_6),1)
    return cross_term_control_points

def get_cross_term_bezier_control_points_from_fourth_order_spline(control_points,dimension):
    P_0 = control_points[:,0][:,None]
    P_1 = control_points[:,1][:,None]
    P_2 = control_points[:,2][:,None]
    P_3 = control_points[:,3][:,None]
    P_4 = control_points[:,4][:,None]
    p1 = (P_0 - 3*P_1 + 3*P_2 - P_3)/2
    p2 = P_0/6 - 2*P_1/3 + P_2 - 2*P_3/3 + P_4/6
    p3 = (P_0 - P_1 - P_2 + P_3)/2
    p4 = (P_0/3 + P_1 - P_2 - P_3/3)/2
    if dimension == 3:
        Y1 = np.array([[0,1,0],[0,0,1],[1,0,0]])
        Y2 = np.array([[0,0,1],[1,0,0],[0,1,0]])
        Y3 = np.array([[0,0,1],[0,0,-1],[0,1,0]])
        Y4 = np.array([[0,1,0],[1,0,0],[1,0,0]])
    elif dimension == 2:
        Y1 = Y4 = np.array([1,0])
        Y2 = Y3 = np.array([0,1])
    c_4 = np.dot(Y1, 2*p1)*np.dot(Y2,p2) - np.dot(Y2,2*p1)*np.dot(Y1,p2) \
         + np.dot(Y3,p1)*np.dot(Y4,3*p2) - np.dot(Y4,p1)*np.dot(Y3,3*p2)

    c_3 = np.dot(Y2,p3)*np.dot(Y1,p2) - np.dot(Y1,p3)*np.dot(Y2,p2) \
        - np.dot(Y3,p3)*np.dot(Y4,3*p2) + np.dot(Y4,p3)*np.dot(Y3,3*p2) \
        + np.dot(Y4,p1)*np.dot(Y3,2*p1) - np.dot(Y3,p1)*np.dot(Y4,2*p1)
    
    c_2 = np.dot(Y2,p4)*np.dot(Y1,3*p2) - np.dot(Y1,p4)*np.dot(Y2,3*p2) \
        + np.dot(Y4,p3)*np.dot(Y3,p1) - np.dot(Y3,p3)*np.dot(Y4,p1) \
        - np.dot(Y4,p3)*np.dot(Y3,2*p1) + np.dot(Y3,p3)*np.dot(Y4,2*p1)
    
    c_1 = np.dot(Y1, p4)*np.dot(Y2,2*p1) - np.dot(Y2,p4)*np.dot(Y1,2*p1)

    c_0 = - np.dot(Y4,p4)*np.dot(Y3,p3) + np.dot(Y4,p3)*np.dot(Y3,p4)
    a0 = c_0
    a1 = (4*a0+c_1)/4
    a2 = (c_2 + 12*a1 - 6*a0)/6 
    a3 = (c_3 + 12*a2 + 4*a0 - 12*a1)/4
    a4 = c_4 + 4*a1 + 4*a3 - 6*a2 - a0
    if dimension == 2:
        cross_term_control_points = np.array([a0,a1,a2,a3,a4]).flatten()
    else:
        cross_term_control_points = np.concatenate((a0,a1,a2,a3,a4),1)
    return cross_term_control_points
