import numpy as np
from bsplinegenerator.helper_functions import count_number_of_control_points
import time
from helper import calculate_curvature, get_matrix
from bsplinegenerator.bsplines import BsplineEvaluation
from beziercurvegenerator.bezier_curves import BezierCurve

"""
Uses the slsqp method 
"""

order = 3
control_points = np.random.randint(10, size=(2,order+1)) # random
control_points = np.array([[0, 3, 6, -1],[0, 2, 7, 9]])

print("control_points: " , control_points)

def get_max_curvature(control_points):
    order = count_number_of_control_points(control_points) - 1
    if order < 2 or order > 3:
        raise Exception("Cannot compute geometric max curvature for " , order , " spline.")
    elif order == 2:
        max_curvature = get_2nd_order_spline_max_curvature(control_points)
    elif order == 3:
        max_curvature = get_3rd_order_spline_max_curvature(control_points)
    return max_curvature

def get_2nd_order_spline_max_curvature(control_points):
    p0 = control_points[:,0]
    p1 = control_points[:,1]
    p2 = control_points[:,2]
    m = (p0+p2)/2
    leg_start = p1-p0
    leg_middle = p1-m
    leg_end = p1-p2
    A = np.linalg.norm(np.cross(leg_start,leg_end))/2
    if np.dot(leg_start,leg_middle) <= 0:
        max_curvature = 2*A/np.linalg.norm(leg_start)**3
    elif np.dot(leg_middle,leg_end) <= 0:
        max_curvature = 2*A/np.linalg.norm(leg_end)**3
    else:
        max_curvature = 2*np.linalg.norm(leg_middle)**3 / A**2
    return max_curvature

def get_3rd_order_spline_max_curvature(control_points):
    # B0 = control_points[:,0]
    # B1 = control_points[:,1]
    # B2 = control_points[:,2]
    # B3 = control_points[:,3]
    # W0 = (B2-B0)/6
    # W1 = (B2-B1)/3
    # W2 = (B3-B1)/6
    P0 = control_points[:,0]
    P1 = control_points[:,1]
    P2 = control_points[:,2]
    P3 = control_points[:,3]
    W0 = P1-P0
    W1 = P2-P1
    W2 = P3-P2
    print("W0: " , W0 , " W1: " , W1 , " W2: " , W2)
    print("condition 1: " , W1[0] - W0[0] == 0)
    print("condition 2: " , np.linalg.norm(W0-2*W1+W2) > 0)
    if (W1[0] == W0[0]):
        a = 3*(W0[1] - 2*W1[1] + W2[1])
        print("condition 3: " , a != 0)
        alpha_1 = W0[1]/W0[0]
        alpha_2 = (W1[1] - W0[1])/(3*W0[0]**2)
        alpha_3 = a / (27*W0[0]**3)
        if a < 0:
            alpha_1 = -alpha_1
            alpha_2 = -alpha_2
            alpha_3 = -alpha_3
        beta_1 = alpha_1 - (alpha_2**2)/alpha_3
        print("a: " , a)
        print("alpha_1: " , alpha_1 ," alpha_2: " , alpha_2, " alpha_3: " , alpha_3 )
        t1 = np.sqrt((-2*beta_1 + np.sqrt(9*beta_1**2 + 5))/5)/(3*W0[0])
        t2 = -np.sqrt((-2*beta_1 + np.sqrt(9*beta_1**2 + 5))/5)/(3*W0[0])
        print("t1: " , t1 , " t2: " , t2)
        bezier_curve = BezierCurve(control_points,0,1)
        curvature_1 = bezier_curve.evaluate_curvature_at_time(0)
        curvature_2 = bezier_curve.evaluate_curvature_at_time(1)
        max_curvature = np.max([curvature_1, curvature_2])
        if t1 > 0 and t1 < 1:
            curvature = bezier_curve.evaluate_curvature_at_time(t1)
            max_curvature = np.max([curvature, max_curvature])
        if t2 > 0 and t2 < 1:
            curvature = bezier_curve.evaluate_curvature_at_time(t2)
            max_curvature = np.max([curvature, max_curvature])
        return max_curvature
    
start_time = time.time()
max_curvature = get_max_curvature(control_points)
end_time = time.time()
print("total time: " , end_time - start_time)
print("max_curvature: ", max_curvature)
# bspline = BsplineEvaluation(control_points,order,0,1)
# curvature_data, time_data = bspline.get_spline_curvature_data(10000)
# print("true_max_curvature: " , np.max(curvature_data))
bezier_curve = BezierCurve(control_points,0,1)
curvature_data, time_data = bezier_curve.get_curvature_data(1000)
print("true_max_curvature: " , np.max(curvature_data))
bezier_curve.plot_curve_data(1000)
bezier_curve.plot_curvature(1000)
