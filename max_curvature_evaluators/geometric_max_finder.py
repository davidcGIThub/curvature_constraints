import numpy as np
from bsplinegenerator.helper_functions import count_number_of_control_points
import time
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import calculate_curvature, get_matrix
from bsplinegenerator.bsplines import BsplineEvaluation
from beziercurvegenerator.bezier_curves import BezierCurve
import matplotlib.pyplot as plt

"""
Uses the slsqp method 
"""

order = 2
control_points = np.random.randint(10, size=(2,order+1)) # random

print("control_points: " , control_points)

def get_max_curvature(control_points):
    order = count_number_of_control_points(control_points) - 1
    if order != 2:
        raise Exception("Cannot compute geometric max curvature for " , order , " spline.")
    elif order == 2:
        F = np.array([[1,1,0],[0,2,0],[0,1,1]])/2
        bezier_control_points = np.dot(F,control_points.T).T
        max_curvature = get_2nd_order_spline_max_curvature(bezier_control_points)
    return max_curvature

def get_2nd_order_spline_max_curvature(bezier_control_points):
    p0 = bezier_control_points[:,0]
    p1 = bezier_control_points[:,1]
    p2 = bezier_control_points[:,2]
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

bspline = BsplineEvaluation(control_points,order,0,1)
bspline_data, time_data_bspline = bspline.get_spline_data(100)
start_time = time.time()
max_curvature = get_max_curvature(control_points)
end_time = time.time()
print("total time: " , end_time - start_time)
print("max_curvature: ", max_curvature)
F = np.array([[1,1,0],[0,2,0],[0,1,1]])/2
bezier_control_points = np.dot(F,control_points.T).T
bezier_curve = BezierCurve(bezier_control_points,0,1)
curve_data, time_data = bezier_curve.get_curve_data(1000)
curvature_data, time_data = bezier_curve.get_curvature_data(1000)
print("true_max_curvature: " , np.max(curvature_data))
plt.plot(curve_data[0,:],curve_data[1,:])
plt.plot(bspline_data[0,:],bspline_data[1,:])
plt.show()
bezier_curve.plot_curvature(1000)
