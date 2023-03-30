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
from PathObjectivesAndConstraints.python_wrappers.waypoint_constraints import WaypointConstraints
from max_curvature_evaluators.max_numerator_over_min_denominator import find_curvature_using_max_numerator_over_min_denominator


class PathGenerator:
    """
    This class generates a 3rd order B-spline path between two waypoints,
    waypoint directions, curvature constraint, and adjoining 
    safe flight corridors.
    """

# when minimizing distance and time need constraint over max velocity
# when minimizing acceleration, no max velocity constraint needed

    def __init__(self, dimension, num_obstacles = 1):
        self._dimension = dimension
        self._num_obstacles = num_obstacles
        self._order = 3
        self._M = get_M_matrix(self._order)
        self._control_points_per_corridor = 4
        self._num_control_points = 8
        self._objective_func_obj = ObjectiveFunctions(self._dimension)
        self._curvature_const_obj = CurvatureConstraints(self._dimension)
        self._obstacle_const_obj = ObstacleConstraints(self._dimension, self._num_obstacles)
        self._waypoint_vel_const_obj = WaypointConstraints(self._dimension)
        
    def generate_path(self, waypoints, waypoint_directions, max_curvature, obstacles = None):
        if(obstacles == None):
            return self.generate_obstacle_free_path(waypoints, waypoint_directions, max_curvature)
        else:
            return self.generate_obstructed_path(waypoints, waypoint_directions, max_curvature, obstacles)
    
    def generate_obstacle_free_path(self, waypoints, waypoint_directions, max_curvature):
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
        return optimized_control_points
    
    def generate_obstructed_path(self, waypoints, waypoint_directions, max_curvature, obstacles):
        # create initial conditions
        self._dimension = np.shape(waypoints)[0]
        initial_control_points = self.__create_initial_control_points(waypoints)
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor]))
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints)
        velocity_constraint = self.__create_waypoint_direction_constraint(waypoint_directions)
        curvature_constraint = self.__create_curvature_constraint(max_curvature)
        print("obstacles[0]: ", obstacles[0])
        print("obstacles[1]: ", obstacles[1])
        obstacle_constraint = self.__create_obstacle_constraint(obstacles[0], obstacles[1])
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
            waypoint_constraint,
            obstacle_constraint), 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x[0:-1] ,(self._dimension,self._num_control_points))
        return optimized_control_points

# bounds over just the scale factor
    def __create_objective_variable_bounds(self):
        lower_bounds = np.zeros(self._num_control_points*self._dimension + 1) - np.inf
        upper_bounds = np.zeros(self._num_control_points*self._dimension + 1) + np.inf
        lower_bounds[self._num_control_points*self._dimension] = 0.00001
        return Bounds(lb=lower_bounds, ub = upper_bounds)

    def __minimize_acceleration_objective_function(self, variables):
        # for third order splines only
        control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
        scale_factor = variables[-1]
        return self._objective_func_obj.minimize_acceleration_and_time(control_points, scale_factor)

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

    def __create_waypoint_direction_constraint(self, direction_vectors):
        def velocity_constraint_function(variables):
            control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
            scale_factor = variables[-1]
            constraints = self._waypoint_vel_const_obj.velocity_at_waypoints_constraints(control_points,
                scale_factor, direction_vectors)
            return constraints.flatten()
        lower_bound = 0
        upper_bound = 0
        velocity_vector_constraint = NonlinearConstraint(velocity_constraint_function, lb= lower_bound, ub=upper_bound)
        return velocity_vector_constraint

    def __create_curvature_constraint(self, max_curvature):
        num_intervals = self._num_control_points - self._order
        def curvature_constraint_function(variables):
            control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
            return self._curvature_const_obj.get_interval_curvature_constraints(control_points,max_curvature)
        lower_bound = np.zeros(num_intervals) - np.inf
        upper_bound = np.zeros(num_intervals)
        curvature_constraint = NonlinearConstraint(curvature_constraint_function , lb = lower_bound, ub = upper_bound)
        return curvature_constraint
    
    def __create_obstacle_constraint(self, obstacle_centers, obstacle_radii):
        # num_obstacles = len(obstacle_radii)
        num_obstacles = self._num_obstacles
        if num_obstacles == 1:
            def obstacle_constraint_function(variables):
                control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
                (self._dimension,self._num_control_points))
                return self._obstacle_const_obj.getObstacleDistanceToSpline(control_points, \
                    obstacle_radii.item(0), obstacle_centers)
        else:
            def obstacle_constraint_function(variables):
                control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
                (self._dimension,self._num_control_points))
                return self._obstacle_const_obj.getObstacleDistancesToSpline(control_points, \
                    obstacle_radii, obstacle_centers)
        lower_bound = np.zeros(num_obstacles) - np.inf
        upper_bound = np.zeros(num_obstacles)
        obstacle_constraint = NonlinearConstraint(obstacle_constraint_function , lb = lower_bound, ub = upper_bound)
        return obstacle_constraint