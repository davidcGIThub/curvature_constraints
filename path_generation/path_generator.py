"""
This module generates a 3rd order B-spline path between two waypoints,
waypoint directions, curvature constraint, and adjoining 
safe flight corridors.
"""
import os
import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, Bounds
from path_generation.matrix_evaluation import get_M_matrix, get_T_derivative_vector
from PathObjectivesAndConstraints.python_wrappers.objective_functions import ObjectiveFunctions
from PathObjectivesAndConstraints.python_wrappers.curvature_constraints import CurvatureConstraints
from PathObjectivesAndConstraints.python_wrappers.obstacle_constraints import ObstacleConstraints
from PathObjectivesAndConstraints.python_wrappers.waypoint_velocity_constraints import WaypointConstraints
from max_curvature_evaluators.max_numerator_over_min_denominator import find_curvature_using_max_numerator_over_min_denominator


class PathGenerator:
    """
    This class generates a 3rd order B-spline path between two waypoints,
    waypoint directions, curvature constraint, and adjoining 
    safe flight corridors.
    """

# when minimizing distance and time need constriant over max velocity
# when minimizing acceleration, no max velocity constraint needed

    def __init__(self, dimension):
        self._dimension = dimension
        self._order = 3
        self._M = get_M_matrix(self._order)
        self._control_points_per_corridor = 4
        self._num_control_points = 8
        self._objective_func_obj = ObjectiveFunctions(self._dimension)
        self._curvature_const_obj = CurvatureConstraints(self._dimension)
        self._obstacle_const_obj = ObstacleConstraints(self._dimension)
        self._waypoint_vel_const_obj = WaypointConstraints(self._dimension)
        

    def generate_path(self, waypoints, waypoint_directions, max_curvature):
        # create initial conditions
        self._dimension = np.shape(waypoints)[0]
        initial_control_points = self.__create_initial_control_points(waypoints)
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor]))
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints)
        velocity_constraint = self.__create_waypoint_direction_constraint(waypoint_directions)
        curvature_constraint = self.__create_curvature_constraint(max_curvature)
        objectiveFunction = self.__minimize_acceleration_objective_function
        objective_variable_bounds = self.__create_objective_variable_bounds()
        minimize_options = {'disp': True}#, 'maxiter': self.maxiter, 'ftol': tol}
        # perform optimization
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            method='SLSQP', 
            bounds=objective_variable_bounds,
            constraints=(\
            curvature_constraint, 
            velocity_constraint,
            waypoint_constraint), 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x[0:-1] ,(self._dimension,self._num_control_points))
        optimized_scale_factor = result.x[-1]
        return optimized_control_points

# bounds over just the scale factor
    def __create_objective_variable_bounds(self):
        lower_bounds = np.zeros(self._num_control_points*self._dimension + 1) - np.inf
        upper_bounds = np.zeros(self._num_control_points*self._dimension + 1) + np.inf
        lower_bounds[self._num_control_points*self._dimension] = 0.00001
        return Bounds(lb=lower_bounds, ub = upper_bounds)

    # def __minimize_acceleration_objective_function(self, variables):
    #     # for third order splines only
    #     control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
    #         (self._dimension,self._num_control_points))
    #     scale_factor = variables[-1]
    #     return self._objective_func_obj.minimize_acceleration_and_time(control_points, scale_factor)

    def __minimize_acceleration_objective_function(self, variables):
        # for third order splines only
        control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
        scale_factor = variables[-1]
        num_intervals = self._num_control_points - self._order
        sum_of_integrals = 0
        for i in range(num_intervals):
            p0 = control_points[:,i]
            p1 = control_points[:,i+1]
            p2 = control_points[:,i+2]
            p3 = control_points[:,i+3]
            sum_of_integrals += np.sum((p0 - 3*p1 + 3*p2 - p3)**2) 
        python_objective = sum_of_integrals + scale_factor
        c_objective = self._objective_func_obj.minimize_acceleration_and_time(control_points, scale_factor)
        print("python objective: " , python_objective)
        print("c++ objective: " , c_objective)
        return c_objective

    def __create_initial_control_points(self, waypoints):
        start_waypoint = waypoints[:,0]
        end_waypoint = waypoints[:,1]
        control_points = np.linspace(start_waypoint,end_waypoint,self._num_control_points).T
        return control_points

    def __create_waypoint_constraint(self, waypoints):
        num_waypoints = 2
        m = num_waypoints
        n = self._num_control_points
        k = self._order
        d = self._dimension
        constraint_matrix = np.zeros((m*d,n*d))
        Gamma_0 = np.zeros((self._order+1,1))
        Gamma_0[self._order,0] = 1
        Gamma_f = np.ones((self._order+1,1))
        M_Gamma_0_T = np.dot(self._M,Gamma_0).T
        M_Gamma_f_T = np.dot(self._M,Gamma_f).T
        for i in range(self._dimension):
            constraint_matrix[i*m   ,  i*n        : i*n+k+1] = M_Gamma_0_T
            constraint_matrix[i*m+1 , (i+1)*n-k-1 : (i+1)*n] = M_Gamma_f_T
        constraint_matrix = np.concatenate((constraint_matrix,np.zeros((m*d,1))),1)
        constraint = LinearConstraint(constraint_matrix, lb=waypoints.flatten(), ub=waypoints.flatten())
        return constraint

    # def __create_waypoint_direction_constraint(self, direction_vectors):
    #     def velocity_constraint_function(variables):
    #         control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
    #         (self._dimension,self._num_control_points))
    #         scale_factor = variables[-1]
    #         constraints = self._waypoint_vel_const_obj.velocity_at_waypoints_constraints(control_points,
    #             scale_factor, direction_vectors)
    #         return constraints.flatten()
    #     lower_bound = 0
    #     upper_bound = 0
    #     velocity_vector_constraint = NonlinearConstraint(velocity_constraint_function, lb= lower_bound, ub=upper_bound)
    #     return velocity_vector_constraint

    def __create_waypoint_direction_constraint(self, velocities):
        def velocity_constraint_function(variables):
            constraints = np.zeros(self._dimension*2)
            control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
            segement_1_control_points = control_points[:,0:self._order+1]
            segement_2_control_points = control_points[:,self._num_control_points-self._order-1:]
            scale_factor = variables[-1]
            T_0 = get_T_derivative_vector(self._order,0,0,1,scale_factor)
            T_f = get_T_derivative_vector(self._order,scale_factor,0,1,scale_factor)
            start_velocity = np.dot(segement_1_control_points,np.dot(self._M,T_0)).flatten()
            end_velocity = np.dot(segement_2_control_points,np.dot(self._M,T_f)).flatten()
            desired_start_velocity = velocities[:,0]
            desired_end_velocity = velocities[:,1]
            constraints[0:self._dimension] = start_velocity - desired_start_velocity
            constraints[self._dimension:] = end_velocity - desired_end_velocity
            return constraints
        lower_bound = 0
        upper_bound = 0
        velocity_vector_constraint = NonlinearConstraint(velocity_constraint_function, lb= lower_bound, ub=upper_bound)
        return velocity_vector_constraint

    # def __create_curvature_constraint(self, max_curvature):
    #     def curvature_constraint_function(variables):
    #         control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
    #         (self._dimension,self._num_control_points))
    #         return self._curvature_const_obj.get_spline_curvature_constraint(control_points,max_curvature)
    #     lower_bound = -np.inf
    #     upper_bound = 0
    #     curvature_constraint = NonlinearConstraint(curvature_constraint_function , lb = lower_bound, ub = upper_bound)
    #     return curvature_constraint

    def __create_curvature_constraint(self, max_curvature):
        def curvature_constraint_function(variables):
            control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
            max_curvature_of_spline_intervals = self.__get_max_curvature_of_each_spline_interval(control_points)
            largest_curvature = np.max(max_curvature_of_spline_intervals)
            constraint = largest_curvature - max_curvature
            return constraint
        lower_bound = -np.inf
        upper_bound = 0
        curvature_constraint = NonlinearConstraint(curvature_constraint_function , lb = lower_bound, ub = upper_bound)
        return curvature_constraint
        
    def __get_max_curvature_of_each_spline_interval(self,control_points):
        num_intervals = self._num_control_points - self._order
        max_curvatures = np.zeros(num_intervals)
        for i in range(num_intervals):
            control_points_per_interval = control_points[:,i:i+self._order+1]
            max_curvatures[i] = find_curvature_using_max_numerator_over_min_denominator(control_points_per_interval,self._order,self._M)
        return  max_curvatures