import numpy as np
from bsplinegenerator.helper_functions import count_number_of_control_points
import time
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import calculate_curvature, get_matrix
from bsplinegenerator.bsplines import BsplineEvaluation
from beziercurvegenerator.bezier_curves import BezierCurve
import matplotlib.pyplot as plt

def get_max_curvature(control_points):
    order = count_number_of_control_points(control_points) - 1
    if order != 2:
        raise Exception("Cannot compute geometric max curvature for " , order , " spline.")
    elif order == 2:
        F = np.array([[1,1,0],[0,2,0],[0,1,1]])/2
        bezier_control_points = np.dot(F,control_points.T).T
        max_curvature = get_2nd_order_spline_max_curvature(bezier_control_points)
    return max_curvature

def get_max_curvature(control_points, order):
    num_cont_pts = count_number_of_control_points(control_points)
    num_intervals = num_cont_pts - 2
    if order != 2:
        raise Exception("Cannot compute geometric max curvature for " , order , " spline.")
    elif order == 2:
        conversion_matrix = np.array([[1,1,0],
                                      [0,2,0],
                                      [0,1,1]])/2  
        max_curvature = 0
        for i in range(num_intervals):
            interval_control_points = control_points[:,i:i+3]
            bezier_control_points = np.transpose(np.dot(conversion_matrix, np.transpose(interval_control_points)))
            curvature_extrema = get_2nd_order_spline_max_curvature(bezier_control_points)
            if curvature_extrema > max_curvature:
                max_curvature = curvature_extrema
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
        # print("method 1")
        norm_start = np.linalg.norm(leg_start)
        if norm_start == 0:
            max_curvature = 0
        else:
            max_curvature = A/norm_start**3
    elif np.dot(leg_middle,leg_end) <= 0:
        norm_end = np.linalg.norm(leg_end)
        if norm_end == 0:
            max_curvature = 0
        else:
            max_curvature = A/norm_end**3
        # print("method 2")
    else:
        if A == 0:
            max_curvature = np.inf
        else:
            max_curvature = np.linalg.norm(leg_middle)**3 / (A*A)
        # print("method 3")
    return max_curvature
