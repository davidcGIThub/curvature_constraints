""" 
This module generates a B-spline path from one point to another
with some given direction vector at each point, and with some
curvature constraint
"""

import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint
from bsplinegenerator.matrix_evaluation import get_M_matrix
# from bsplinegenerator

class PathGenerator:
    """ 
    This class generates a B-spline path from one point to another
    with some given direction vector at each point, and with some
    curvature constraint
    """

    def __init__(self, order, dimension):
        # 3rd order spline 5 intervals - 8 control points
        # 4th order spline 4 intervals - 8 control points
        # 5th order spline 3 intervals - 8 control points
        self._order = order
        self._dimension = dimension
        self._num_control_points = 8

    def generate_trajectory(self, waypoints):
        # create initial conditions
        self._dimension = np.shape(waypoints)[0]
        initial_control_points = self.__create_initial_control_points(waypoints)
        optimization_variables = initial_control_points.flatten()
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints)
        objectiveFunction = self.__minimize_distance_objective_function
        minimize_options = {'disp': True}#, 'maxiter': self.maxiter, 'ftol': tol}
        # perform optimizationd
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            method='SLSQP', 
            # method = 'trust-constr',
            constraints=(waypoint_constraint), 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x ,(self._dimension,self._num_control_points))
        return optimized_control_points


    def __minimize_distance_objective_function(self,variables):
        control_points = np.reshape(variables,(self._dimension,self._num_control_points))
        distances = self.__get_distances_between_points(control_points)
        return np.sum(distances)

    def __get_distances_between_points(self,points):
        number_of_points = np.shape(points)[1]
        first_points = points[:,0:number_of_points-1]
        next_points = points[:,1:number_of_points]
        distances = np.sqrt(np.sum(((next_points - first_points)**2),0))
        return distances

    def __create_initial_control_points(self, waypoints):
        start_waypoint = waypoints[:,0]
        end_waypoint = waypoints[:,1]
        control_points = np.linspace(start_waypoint,end_waypoint,self._num_control_points).T
        return control_points

    def __create_waypoint_constraint(self, waypoints):
        M = get_M_matrix(0, self._order, np.array([]), False)
        num_waypoints = 2
        m = num_waypoints
        n = self._num_control_points
        k = self._order
        d = self._dimension
        constraint_matrix = np.zeros((m*d,n*d))
        Gamma_0 = np.zeros((self._order+1,1))
        Gamma_0[self._order,0] = 1
        Gamma_f = np.ones((self._order+1,1))
        M_Gamma_0_T = np.dot(M,Gamma_0).T
        M_Gamma_f_T = np.dot(M,Gamma_f).T
        for i in range(self._dimension):
            constraint_matrix[i*m   ,  i*n        : i*n+k+1] = M_Gamma_0_T
            constraint_matrix[i*m+1 , (i+1)*n-k-1 : (i+1)*n] = M_Gamma_f_T
        constraint = LinearConstraint(constraint_matrix, lb=waypoints.flatten(), ub=waypoints.flatten())
        return constraint