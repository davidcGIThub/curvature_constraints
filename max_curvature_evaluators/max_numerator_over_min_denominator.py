import numpy as np
from max_curvature_evaluators.helper_files.cube_root_solver import solver
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.helper_functions import count_number_of_control_points
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import calculate_velocity_magnitude, calculate_cross_term_magnitude
import matplotlib.pyplot as plt
import time

def find_curvature_using_max_numerator_over_min_denominator(control_points, order, M):
    dimension = np.shape(control_points)[0]
    min_velocity = find_min_velocity_magnitude(control_points, order, M)
    max_acceleration = find_max_acceleration(control_points, order)
    max_cross_term = find_max_cross_term(control_points, order, M, dimension)
    min_numerator = np.min((max_cross_term, max_acceleration))
    if min_velocity == 0:
        if min_numerator == 0:
            max_curvature = 0
        else:
            max_curvature = np.inf
    else:
        curvature1 = max_acceleration/min_velocity**2
        curvature2 = max_cross_term/min_velocity**3
        max_curvature = np.min((curvature1,curvature2))
    return max_curvature

def find_min_velocity_magnitude(control_points, order, M):
    P = control_points
    J = np.dot(np.dot(M.T,P.T) , np.dot(P,M))
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
    times_to_check = np.concatenate((roots,np.array([0,1])))
    min_velocity = np.inf
    for i in range(len(times_to_check)):
        t = times_to_check[i]
        if t >= 0 and t <= 1:
            velocity = calculate_velocity_magnitude(t, M, control_points,order)
            if velocity < min_velocity:
                min_velocity = velocity
    return min_velocity

def find_max_cross_term(control_points, order, M, dimension):
    if order == 2:
        p0 = control_points[:,0]
        p1 = control_points[:,1]
        p2 = control_points[:,2]
        max_cross_term = np.linalg.norm(np.cross(p0,p1) + np.cross(p1,p2) + np.cross(p2,p0))
    elif order == 3:
        A,B,C,D = get_cross_coeficients(dimension,control_points)
        roots = solver(A,B,C,D)
        times_to_check = np.concatenate((roots,np.array([0,1])))
        max_cross_term = 0
        for i in range(len(times_to_check)):
            t = times_to_check[i]
            if t >= 0 and t <= 1:
                cross_term = calculate_cross_term_magnitude(t,M,control_points,order)
                if cross_term > max_cross_term:
                    max_cross_term = cross_term
    return max_cross_term

def find_max_acceleration(control_points,order):
    if order == 2:
        MT = np.array([[1],[-2],[1]])
        max_acceleration = np.linalg.norm(np.dot(control_points,MT))
    elif order == 3:
        alpha = 1
        MT_start = np.array([[1],[-2],[1],[0]])
        MT_end = np.array([[1-1/alpha],[3/alpha-2],[1-3/alpha],[1/alpha]])
        start_acceleration = np.linalg.norm(np.dot(control_points,MT_start))
        end_acceleration = np.linalg.norm(np.dot(control_points,MT_end))
        max_acceleration = np.max((start_acceleration,end_acceleration))
    return max_acceleration

def get_cross_coeficients(dimension,control_points):
    if dimension == 3:
        p0x = control_points[0,0]
        p0y = control_points[1,0]
        p0z = control_points[2,0]
        p1x = control_points[0,1]
        p1y = control_points[1,1]
        p1z = control_points[2,1]
        p2x = control_points[0,2]
        p2y = control_points[1,2]
        p2z = control_points[2,2]
        p3x = control_points[0,3]
        p3y = control_points[1,3]
        p3z = control_points[2,3]
        c_3 = ((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) \
            *((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - \
            (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2) + ((p0z - 2*p1z + p2z)*\
            (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0z - 3*p1z + 3*p2z - p3z))* \
            ((p0x*p1z)/2 - (p1x*p0z)/2 - p0x*p2z + p2x*p0z + (p0x*p3z)/2 + (3*p1x*p2z)/2 - (3*p2x*p1z)/2 - \
            (p3x*p0z)/2 - p1x*p3z + p3x*p1z + (p2x*p3z)/2 - (p3x*p2z)/2) + ((p0z - 2*p1z + p2z)* \
            (p0y - 3*p1y + 3*p2y - p3y) - (p0y - 2*p1y + p2y)*(p0z - 3*p1z + 3*p2z - p3z))*((p0y*p1z)/2 - \
            (p1y*p0z)/2 - p0y*p2z + p2y*p0z + (p0y*p3z)/2 + (3*p1y*p2z)/2 - (3*p2y*p1z)/2 - (p3y*p0z)/2 - \
            p1y*p3z + p3y*p1z + (p2y*p3z)/2 - (p3y*p2z)/2)
        c_2 = - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y))*\
            ((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - \
            (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2) - ((p0z/2 - p2z/2)*\
            (p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0z - 3*p1z + 3*p2z - p3z))*((p0x*p1z)/2 - \
            (p1x*p0z)/2 - p0x*p2z + p2x*p0z + (p0x*p3z)/2 + (3*p1x*p2z)/2 - (3*p2x*p1z)/2 - (p3x*p0z)/2 - \
            p1x*p3z + p3x*p1z + (p2x*p3z)/2 - (p3x*p2z)/2) - ((p0z/2 - p2z/2)*(p0y - 3*p1y + 3*p2y - p3y) - \
            (p0y/2 - p2y/2)*(p0z - 3*p1z + 3*p2z - p3z))*((p0y*p1z)/2 - (p1y*p0z)/2 - p0y*p2z + p2y*p0z + \
            (p0y*p3z)/2 + (3*p1y*p2z)/2 - (3*p2y*p1z)/2 - (p3y*p0z)/2 - p1y*p3z + p3y*p1z + (p2y*p3z)/2 - \
            (p3y*p2z)/2) - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*\
            (p0y - 3*p1y + 3*p2y - p3y))*((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)* \
            (p0y - 3*p1y + 3*p2y - p3y)) - ((p0z/2 - p2z/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*\
            (p0z - 3*p1z + 3*p2z - p3z))*((p0z - 2*p1z + p2z)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*\
            (p0z - 3*p1z + 3*p2z - p3z)) - ((p0z/2 - p2z/2)*(p0y - 3*p1y + 3*p2y - p3y) - (p0y/2 - p2y/2)*\
            (p0z - 3*p1z + 3*p2z - p3z))*((p0z - 2*p1z + p2z)*(p0y - 3*p1y + 3*p2y - p3y) - (p0y - 2*p1y + p2y)\
            *(p0z - 3*p1z + 3*p2z - p3z))
        c_1 = ((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y - 2*p1y + p2y)*\
            (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) + ((p0z/2 - p2z/2)*\
            (p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0z - 2*p1z + p2z))*((p0z - 2*p1z + p2z)* \
            (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0z - 3*p1z + 3*p2z - p3z)) + ((p0z/2 - p2z/2)*\
            (p0y - 2*p1y + p2y) - (p0y/2 - p2y/2)*(p0z - 2*p1z + p2z))*((p0z - 2*p1z + p2z)* \
            (p0y - 3*p1y + 3*p2y - p3y) - (p0y - 2*p1y + p2y)*(p0z - 3*p1z + 3*p2z - p3z)) + \
            ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y))**2 + \
            ((p0z/2 - p2z/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0z - 3*p1z + 3*p2z - p3z))**2 + \
            ((p0z/2 - p2z/2)*(p0y - 3*p1y + 3*p2y - p3y) - (p0y/2 - p2y/2)*(p0z - 3*p1z + 3*p2z - p3z))**2
        c_0 = -((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y/2 - p2y/2)* \
            (p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y)) - ((p0z/2 - p2z/2)* \
            (p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0z - 2*p1z + p2z))*((p0z/2 - p2z/2)* \
            (p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0z - 3*p1z + 3*p2z - p3z)) - \
            ((p0z/2 - p2z/2)*(p0y - 2*p1y + p2y) - (p0y/2 - p2y/2)*(p0z - 2*p1z + p2z))*((p0z/2 - p2z/2)* \
            (p0y - 3*p1y + 3*p2y - p3y) - (p0y/2 - p2y/2)*(p0z - 3*p1z + 3*p2z - p3z))
    elif dimension == 2:
        p0x = control_points[0,0]
        p0y = control_points[1,0]
        p1x = control_points[0,1]
        p1y = control_points[1,1]
        p2x = control_points[0,2]
        p2y = control_points[1,2]
        p3x = control_points[0,3]
        p3y = control_points[1,3]
        c_3 = ((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) \
            *((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - \
            (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2)
        c_2 = - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y))*\
            ((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - \
            (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2) - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*\
            (p0y - 3*p1y + 3*p2y - p3y))*((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)* \
            (p0y - 3*p1y + 3*p2y - p3y))
        c_1 = ((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y - 2*p1y + p2y)*\
            (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) + \
            ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y))**2
        c_0 = -((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y/2 - p2y/2)* \
            (p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y))
    return c_3, c_2, c_1, c_0