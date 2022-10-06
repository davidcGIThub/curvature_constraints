""" 
This module generates a B-spline path from one point to another
with some given direction vector at each point, and with some
curvature constraint
"""

import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint
from bsplinegenerator.matrix_evaluation import get_M_matrix, get_T_derivative_vector
from bsplinegenerator.bspline_to_bezier import get_bspline_to_bezier_conversion_matrix
from max_curvature_evaluators.root_finder import find_max_curvature_root_finder
import sys
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

    def generate_trajectory(self, waypoints,directions,max_curvature):
        # create initial conditions
        self._dimension = np.shape(waypoints)[0]
        initial_control_points = self.__create_initial_control_points(waypoints)
        optimization_variables = initial_control_points.flatten()
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints)
        direction_constraint = self.__create_2D_directional_constraints(directions)
        # curvature_constraint = self.__create_curvature_constraint(max_curvature)
        objectiveFunction = self.__minimize_distance_objective_function
        minimize_options = {'disp': True}#, 'maxiter': self.maxiter, 'ftol': tol}
        # perform optimizationd
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            method='SLSQP', 
            # method = 'trust-constr',
            constraints=(waypoint_constraint,direction_constraint), 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x ,(self._dimension,self._num_control_points))
        return optimized_control_points


    def __minimize_distance_objective_function(self,variables):
        control_points = np.reshape(variables,(self._dimension,self._num_control_points))
        distances = self.__get_distances_between_points(control_points)
        return np.sum(distances)

    # def __minimize_distance_and_curvature_objective_function(self,variables):
    #     control_points = np.reshape(variables,(self._dimension,self._num_control_points))
    #     distances = self.__get_distances_between_points(control_points)
    #     return np.sum(distances)

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

    # def __create_2D_directional_constraints(self,directions):
    #     F = np.transpose(get_bspline_to_bezier_conversion_matrix(self._order))
    #     num_directions = 2
    #     def direction_constraint_function(variables):
    #         constraints = np.zeros(self._dimension*num_directions)
    #         control_points = np.reshape(variables,(self._dimension,self._num_control_points))
    #         segement_1_control_points = control_points[:,0:self._order+1]
    #         segement_2_control_points = control_points[:,self._num_control_points-self._order-1:]
    #         start_bezier_control_points = np.dot(segement_1_control_points,F)
    #         end_bezier_control_points = np.dot(segement_2_control_points,F)
    #         direction_start = directions[:,0]
    #         direction_end = directions[:,1]
    #         vector_start = start_bezier_control_points[:,1] - start_bezier_control_points[:,0]
    #         vector_end = end_bezier_control_points[:,-1] - end_bezier_control_points[:,-2]
    #         norm_vector_start = np.linalg.norm(vector_start)
    #         norm_vector_end = np.linalg.norm(vector_start)
    #         norm_direction_start = np.linalg.norm(direction_start)
    #         norm_direction_end = np.linalg.norm(direction_end)
    #         constraints[0] = direction_start[0]*norm_vector_start - norm_direction_start*vector_start[0]
    #         constraints[1] = direction_start[1]*norm_vector_start - norm_direction_start*vector_start[1]
    #         constraints[2] = direction_end[0]*norm_vector_end - norm_direction_end*vector_end[0]
    #         constraints[3] = direction_end[1]*norm_vector_end - norm_direction_end*vector_end[1]
    #         print("constraints: " , constraints)
    #         return constraints
    #     lower_bound = 0
    #     upper_bound = 0
    #     direction_constraint = NonlinearConstraint(direction_constraint_function , lb= lower_bound, ub=upper_bound)
    #     return direction_constraint

    def __create_2D_directional_constraints(self,directions):
        T_0 = get_T_derivative_vector(self._order,0,0,1,1)
        T_f = get_T_derivative_vector(self._order,1,0,1,1)
        M = get_M_matrix(0, self._order, np.array([]), False)
        def direction_constraint_function(variables):
            constraints = np.zeros(4)
            control_points = np.reshape(variables,(self._dimension,self._num_control_points))
            segement_1_control_points = control_points[:,0:self._order+1]
            segement_2_control_points = control_points[:,self._num_control_points-self._order-1:]
            start_velocity = np.dot(segement_1_control_points,np.dot(M,T_0))
            end_velocity = np.dot(segement_2_control_points,np.dot(M,T_f))
            # start_angle = np.arctan2(start_velocity[1],start_velocity[0])
            # end_angle = np.arctan2(end_velocity[1],end_velocity[0])
            direction_start = directions[:,0]
            direction_end = directions[:,1]
            norm_start_velocity = np.linalg.norm(start_velocity)
            norm_end_velocity = np.linalg.norm(end_velocity)
            norm_direction_start = np.linalg.norm(direction_start)
            norm_direction_end = np.linalg.norm(direction_end)
            constraints[0] = direction_start[0]*norm_start_velocity - norm_direction_start*start_velocity[0]
            constraints[1] = direction_start[1]*norm_start_velocity - norm_direction_start*start_velocity[1]
            constraints[2] = direction_end[0]*norm_end_velocity - norm_direction_end*end_velocity[0]
            constraints[3] = direction_end[1]*norm_end_velocity - norm_direction_end*end_velocity[1]
            print("start_velocity: " , start_velocity)
            print("end_velocity: " , end_velocity)
            print("constraints: " , constraints)
            return constraints
        lower_bound = 0
        upper_bound = 0
        direction_constraint = NonlinearConstraint(direction_constraint_function , lb= lower_bound, ub=upper_bound)
        return direction_constraint


    # def __create_curvature_constraint(self, max_curvature):
    #     M = get_M_matrix(0, self._order, np.array([]), False)
    #     num_intervals = self._num_control_points - self._order
    #     def curvature_constraint_function(variables):
    #         # constraints = np.zeros(num_intervals)
    #         control_points = np.reshape(variables,(self._dimension,self._num_control_points))
    #         largest_curvature = 0
    #         for i in range(num_intervals):
    #             control_points_per_interval = control_points[:,i:i+self._order+1]
    #             largest_curvature = find_max_curvature_root_finder(control_points_per_interval,self._order,M)
    #         if largest_curvature > sys.maxsize:
    #             largest_curvature = sys.maxsize
    #         constraint = largest_curvature - max_curvature
    #         print("constraint: " , constraint)
    #         return constraint
    #     lower_bound = -np.inf
    #     upper_bound = 0
    #     curvature_constraint = NonlinearConstraint(curvature_constraint_function , lb= lower_bound, ub=upper_bound)
    #     return curvature_constraint
        
    # def __get_max_curvature():

    # def __create_3D_directional_constraints_at_waypoints(self,directions):
    #     F = np.transpose(get_bspline_to_bezier_conversion_matrix(self._order))
    #     num_directions = 2
    #     def direction_constraint_function(variables):
    #         constraints = np.zeros(self._dimension*num_directions)
    #         control_points = np.reshape(variables,(self._dimension,self._num_control_points))
    #         segement_1_control_points = control_points[:,0:self._order+1]
    #         segement_2_control_points = control_points[:,self._num_control_points-self._order-1:]
    #         direction_start = directions[:,0]
    #         direction_end = directions[:,1]
    #         vector_start = segement_1_control_points[:,1] - segement_1_control_points[:,0]
    #         vector_end = segement_2_control_points[:,-1] - segement_2_control_points[:,-2]
    #         for i in range(self._dimension):
    #             direction_start[i-1]*vector_start[i] - direction_start[i]*vector_start[i-1]
    #             direction_end[i-1]*vector_end[i] - direction_end[i]*vector_end[i-1]
