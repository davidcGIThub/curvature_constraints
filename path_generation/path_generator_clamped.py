""" 
This module generates a B-spline path from one point to another
with some given direction vector at each point, and with some
curvature constraint
"""

import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, Bounds
from bsplinegenerator.matrix_evaluation import get_M_matrix, get_T_derivative_vector
from max_curvature_evaluators.root_finder import find_max_curvature_root_finder
import sys
# from bsplinegenerator

class PathGenerator:
    """ 
    This class generates a clamped B-spline path from one point to another
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
        self._num_intermediate_control_points = self._num_control_points - 2
        self._M = get_M_matrix(0, self._order, np.array([]), False)
        self._current_waypoints = np.array([])

    def generate_trajectory(self, waypoints,directions,max_curvature):
        # create initial conditions
        self._current_waypoints = waypoints
        self._dimension = np.shape(waypoints)[0]
        initial_intermediate_control_points = self.__create_initial_intermediate_control_points(waypoints)
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_intermediate_control_points.flatten(),[initial_scale_factor]))
        # define constraints and objective function and constraints
        # direction_constraint = self.__create_2D_directional_constraints(directions)
        objectiveFunction = self.__minimize_distance_and_time_objective_function
        direction_constraint = self.__create_2D_directional_constraints(directions)
        objective_variable_bounds = self.__create_objective_variable_bounds()
        minimize_options = {'disp': True}#, 'maxiter': self.maxiter, 'ftol': tol}
        # perform optimization
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            method='SLSQP', 
            # method = 'trust-constr',
            bounds=objective_variable_bounds,
            constraints=(direction_constraint), 
            options = minimize_options)
        # retrieve data
        optimized_intermediate_control_points = np.reshape(result.x[0:-1] , \
            (self._dimension,self._num_intermediate_control_points))
        optimized_control_points =  np.concatenate((self._current_waypoints[:,0][:,None], \
            optimized_intermediate_control_points,self._current_waypoints[:,1][:,None]),1)
        optimized_scale_factor = result.x[-1]
        return optimized_control_points, optimized_scale_factor

    def __minimize_distance_and_time_objective_function(self,variables):
        intermediate_control_points = np.reshape(variables[0:self._num_intermediate_control_points*self._dimension], \
            (self._dimension,self._num_intermediate_control_points))
        control_points = np.concatenate((self._current_waypoints[:,0][:,None], \
            intermediate_control_points,self._current_waypoints[:,1][:,None]),1)
        # print("control_points: " , control_points)
        distances = self.__get_distances_between_points(control_points)
        scale_factor = variables[-1]
        return np.sum(distances) + scale_factor

    def __create_objective_variable_bounds(self):
        lower_bounds = np.zeros(self._num_intermediate_control_points*self._dimension + 1) - np.inf
        upper_bounds = np.zeros(self._num_intermediate_control_points*self._dimension + 1) + np.inf
        lower_bounds[self._num_intermediate_control_points*self._dimension] = 0.00001
        return Bounds(lb=lower_bounds, ub = upper_bounds)

    def __get_distances_between_points(self,points):
        number_of_points = np.shape(points)[1]
        first_points = points[:,0:number_of_points-1]
        next_points = points[:,1:number_of_points]
        distances = np.sqrt(np.sum(((next_points - first_points)**2),0))
        return distances

    def __create_initial_intermediate_control_points(self, waypoints):
        start_waypoint = waypoints[:,0]
        end_waypoint = waypoints[:,1]
        control_points = np.linspace(start_waypoint,end_waypoint,self._num_control_points).T
        intermediate_control_points = control_points[:,1:-1]
        return intermediate_control_points



    def __create_2D_directional_constraints(self,directions):
        def direction_constraint_function(variables):
            constraints = np.zeros(2)
            intermediate_control_points = np.reshape(variables[0:self._num_intermediate_control_points*self._dimension], \
                (self._dimension,self._num_intermediate_control_points))
            start_vector = intermediate_control_points[:,0] - self._current_waypoints[:,0]
            end_vector = self._current_waypoints[:,1] - intermediate_control_points[:,-1]
            # constraints = np.zeros(4)
            start_angle = np.arctan2(start_vector[1],start_vector[0])
            end_angle = np.arctan2(end_vector[1],end_vector[0])
            print("start_vector: " , start_vector)
            print("end_vector: " , end_vector)
            constraints[0] = start_angle - directions[0]
            constraints[1] = end_angle - directions[1]
            print("constraints: " , constraints)
            return constraints
        lower_bound = 0
        upper_bound = 0
        direction_constraint = NonlinearConstraint(direction_constraint_function , lb= lower_bound, ub=upper_bound)
        return direction_constraint
