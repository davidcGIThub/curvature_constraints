"""
This module generates a 3rd order B-spline path between two waypoints,
waypoint directions, curvature constraint, and adjoining 
safe flight corridors.
"""
import os
import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, Bounds
from dataclasses import dataclass
import ctypes
import pathlib
from matrix_evaluation import get_M_matrix, get_T_derivative_vector


class PathGenerator:
    """
    This class generates a 3rd order B-spline path between two waypoints,
    waypoint directions, curvature constraint, and adjoining 
    safe flight corridors.
    """

# when minimizing distance and time need constriant over max velocity
# when minimizing acceleration, no max velocity constraint needed

    def __init__(self):
        dimension = 2
        self._dimension = dimension
        self._order = 3
        self._M = get_M_matrix(self._order)
        self._control_points_per_corridor = 4
        self._num_control_points = 8
        script_dir = os.path.abspath(os.path.dirname(__file__))
        libname_str = os.path.join(script_dir)
        libname = pathlib.Path(libname_str)
        self._curvature_lib = ctypes.CDLL(libname / "libCurvatureEvaluator.so")
        self._curvature_lib.find_spline_curvature_bound.restype = ctypes.c_double

    def generate_minimum_acceleration_path(self, waypoints, waypoint_directions, max_curvature, sfcs):
        # create initial conditions
        self._dimension = np.shape(waypoints)[0]
        initial_control_points = self.__create_initial_control_points(waypoints,self._num_control_points)
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor]))
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints, self._num_control_points)
        velocity_constraint = self.__create_waypoint_direction_constraint(waypoint_directions, self._num_control_points)
        curvature_constraint = self.__create_curvature_constraint(max_curvature, self._num_control_points)
        objectiveFunction = self.__minimize_acceleration_objective_function
        objective_variable_bounds = self.__create_objective_variable_bounds(self._num_control_points)
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
    def __create_objective_variable_bounds(self, num_control_points):
        lower_bounds = np.zeros(num_control_points*self._dimension + 1) - np.inf
        upper_bounds = np.zeros(num_control_points*self._dimension + 1) + np.inf
        lower_bounds[num_control_points*self._dimension] = 0.00001
        return Bounds(lb=lower_bounds, ub = upper_bounds)

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
        return sum_of_integrals + scale_factor

    def __create_initial_control_points(self, waypoints):
        start_waypoint = waypoints[:,0]
        end_waypoint = waypoints[:,1]
        control_points = np.linspace(start_waypoint,end_waypoint,self._num_control_points).T
        return control_points

    def __create_waypoint_constraint(self, waypoints, num_control_points):
        num_waypoints = 2
        m = num_waypoints
        n = num_control_points
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

    def __create_waypoint_direction_constraint(self, direction_vectors, num_control_points):
        def velocity_constraint_function(variables):
            constraints = np.zeros(self._dimension*2)
            control_points = np.reshape(variables[0:num_control_points*self._dimension], \
            (self._dimension,num_control_points))
            segement_1_control_points = control_points[:,0:self._order+1]
            segement_2_control_points = control_points[:,num_control_points-self._order-1:]
            scale_factor = variables[-1]
            T_0 = get_T_derivative_vector(self._order,0,0,1,scale_factor)
            T_f = get_T_derivative_vector(self._order,scale_factor,0,1,scale_factor)
            start_direction = np.dot(segement_1_control_points,np.dot(self._M,T_0)).flatten()
            end_direction = np.dot(segement_2_control_points,np.dot(self._M,T_f)).flatten()
            desired_start_direction = direction_vectors[:,0]
            desired_end_direction = direction_vectors[:,1]
            constraints[0:self._dimension] = start_direction - desired_start_direction
            constraints[self._dimension:] = end_direction - desired_end_direction
            return constraints
        lower_bound = 0
        upper_bound = 0
        velocity_vector_constraint = NonlinearConstraint(velocity_constraint_function, lb= lower_bound, ub=upper_bound)
        return velocity_vector_constraint

## C++ code method
    def __create_curvature_constraint(self, max_curvature, num_control_points):
        def curvature_constraint_function(variables):
            array_length = num_control_points*self._dimension
            c_array = (ctypes.c_double * array_length)(*variables[0:num_control_points*self._dimension])
            num_cnt_pts = ctypes.c_int(num_control_points)
            largest_curvature = self._curvature_lib.find_spline_curvature_bound(c_array, num_cnt_pts)
            constraint = largest_curvature - max_curvature
            return constraint
        lower_bound = -np.inf
        upper_bound = 0
        curvature_constraint = NonlinearConstraint(curvature_constraint_function , lb = lower_bound, ub = upper_bound)
        return curvature_constraint