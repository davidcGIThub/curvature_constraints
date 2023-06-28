""" 
This module generates a B-spline path from one point to another
with some given direction vector at each point, and with some
curvature constraint
"""
import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, Bounds
from bsplinegenerator.matrix_evaluation import get_M_matrix, get_T_derivative_vector
from bsplinegenerator.bspline_to_bezier import get_bspline_to_bezier_conversion_matrix, convert_to_bezier_control_points
from max_curvature_evaluators.geometric_max_finder import get_max_curvature
from bsplinegenerator.bspline_to_bezier import get_composite_bspline_to_bezier_conversion_matrix
from max_curvature_evaluators.helper_files.cube_root_solver import solver
# from bsplinegenerator

class PathGenerator2nd:
    """ 
    This class generates a B-spline path from one point to another
    with some given direction vector at each point, and with some
    curvature constraint
    """

    def __init__(self, dimension):
        self._order = 2
        self._dimension = dimension
        self._num_control_points = 9
        self._M = get_M_matrix(0, self._order, np.array([]), False)
        self._objective_function_type = "minimize_acceleration"
        self._F_composite = get_composite_bspline_to_bezier_conversion_matrix(self._num_control_points, self._order)
        
    def generate_path(self, waypoints, velocities, max_curvature, initial_control_points = None):
        self._dimension = np.shape(waypoints)[0]
        optimization_variables = self.__create_optimization_variables(waypoints, initial_control_points)
        objectiveFunction = self.__get_objective_function()
        objective_variable_bounds = self.__create_objective_variable_bounds()
        minimize_options = {'disp': True}#, 'maxiter': self.maxiter, 'ftol': tol}
        constraints = self.__get_constraints(waypoints, velocities, max_curvature)
        # perform optimization
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            method='SLSQP', 
            # method = 'trust-constr',
            bounds=objective_variable_bounds,
            constraints=constraints, 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x[0:self._num_control_points*self._dimension] ,(self._dimension,self._num_control_points))
        optimized_scale_factor = result.x[self._num_control_points*self._dimension]
        return optimized_control_points, optimized_scale_factor
    
    def __create_optimization_variables(self, waypoints, initial_control_points=None):
        if initial_control_points is None:
            initial_control_points = self.__create_initial_control_points(waypoints)
        initial_scale_factor = 1
        start_waypoint_scalar = 1
        end_waypoint_scalar = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor],[start_waypoint_scalar, end_waypoint_scalar]))
        return optimization_variables
    
    def __get_constraints(self, waypoints, velocities, max_curvature):
        waypoint_constraint = self.__create_waypoint_constraint(waypoints)
        velocity_constraint = self.__create_waypoint_velocity_constraint(velocities)
        curvature_constraint = self.__create_curvature_constraint(max_curvature)
        constraints = (waypoint_constraint, velocity_constraint, curvature_constraint)
        return constraints

    def __get_objective_function(self):
        # return self.__minimize_jerk_cps_objective_function
        return self.__minimize_acceleration_control_points_objective_function
        # return self.__minimize_velocity_control_points_objective_function
    
    # def __minimize_jerk_cps_objective_function(self, variables):
    #     # for third order splines only
    #     # print(" ")
    #     # print("Minimize Jerk:")
    #     control_points, scale_factor =  self.__get_control_points_and_scale_factor(variables)
    #     # print("control_points: " , control_points)
    #     jerk_cps = control_points[:,3:] - 3*control_points[:,2:-1] + 3*control_points[:,1:-2] - control_points[:,0:-3]
    #     square_jerk_control_points = np.sum(jerk_cps**2,0)
    #     objective = np.sum(square_jerk_control_points)
    #     # print("objective: " , objective)
    #     return objective
    
    # def __minimize_velocity_control_points_objective_function(self, variables):
    #     # for third order splines only
    #     control_points, scale_factor =  self.__get_control_points_and_scale_factor(variables)
    #     # print(" ")
    #     # print("Minimize Vel:")
    #     # print("control_points: " , control_points)
    #     velocity_cps =  control_points[:,0:-1] - control_points[:,1:]
    #     velocity_control_points_squared_sum = np.sum(velocity_cps**2,0)
    #     objective = np.sum(velocity_control_points_squared_sum)
    #     # print("objective: " , objective)
    #     return objective
    
    def __minimize_acceleration_control_points_objective_function(self, variables):
        # for third order splines only
        control_points, scale_factor =  self.__get_control_points_and_scale_factor(variables)
        # print(" ")
        # print("Minimize Accel:")
        # print("control_points: " , control_points)
        acceleration_cps =  control_points[:,2:] - 2*control_points[:,1:-1] + control_points[:,0:-2]
        accel_control_points_squared_sum = np.sum(acceleration_cps**2,0)
        objective = np.sum(accel_control_points_squared_sum)
        # print("objective: " , objective)
        return objective
    
    def __create_objective_variable_bounds(self):
        num_extra_cols = 3
        lower_bounds = np.zeros(self._num_control_points*self._dimension + num_extra_cols) - np.inf
        upper_bounds = np.zeros(self._num_control_points*self._dimension + num_extra_cols) + np.inf
        lower_bounds[self._num_control_points*self._dimension:] = 0.000001
        return Bounds(lb=lower_bounds, ub = upper_bounds)

    def __create_initial_control_points(self, waypoints):
        start_waypoint = waypoints[:,0]
        end_waypoint = waypoints[:,1]
        control_points = np.linspace(start_waypoint,end_waypoint,self._num_control_points).T
        return control_points

    def __create_waypoint_constraint(self, waypoints):
        num_extra_cols = 3
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
        constraint_matrix = np.concatenate((constraint_matrix,np.zeros((m*d,num_extra_cols))),1)
        constraint = LinearConstraint(constraint_matrix, lb=waypoints.flatten(), ub=waypoints.flatten())
        return constraint
    
    def __get_control_points_and_scale_factor(self, variables):
        control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
        scale_factor = variables[self._num_control_points*self._dimension]
        return control_points, scale_factor
    
    def __get_waypoint_scalars(self, variables):
        start_waypoint_scalar = variables[self._num_control_points*self._dimension+1]
        end_waypoint_scalar = variables[self._num_control_points*self._dimension+2]
        return start_waypoint_scalar, end_waypoint_scalar

##### need to adjust for second order spline
    def __create_waypoint_velocity_constraint(self, velocities):
        def velocity_constraint_function(variables):
            constraints = np.zeros(self._dimension*2)
            control_points, scale_factor =  self.__get_control_points_and_scale_factor(variables)
            start_waypoint_scalar, end_waypoint_scalar = self.__get_waypoint_scalars(variables)
            desired_start_velocity = velocities[:,0]
            desired_end_velocity = velocities[:,1]
            start_velocity = end_waypoint_scalar*(control_points[:,1] - control_points[:,0])
            end_velocity = end_waypoint_scalar*(control_points[:,-1] - control_points[:,-2])
            constraints[0:self._dimension] = start_velocity - desired_start_velocity
            constraints[self._dimension:] = end_velocity - desired_end_velocity
            return constraints
        lower_bound = 0
        upper_bound = 0
        velocity_vector_constraint = NonlinearConstraint(velocity_constraint_function, lb= lower_bound, ub=upper_bound)
        return velocity_vector_constraint

    def __create_curvature_constraint(self, max_curvature):
        def curvature_constraint_function(variables):
            control_points, scale_factor =  self.__get_control_points_and_scale_factor(variables)
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
            max_curvatures[i] = get_max_curvature(control_points_per_interval, self._order)
        return  max_curvatures
    