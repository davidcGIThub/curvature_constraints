""" 
This module generates a B-spline path from one point to another
with some given direction vector at each point, and with some
curvature constraint
"""
import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, Bounds
from bsplinegenerator.matrix_evaluation import get_M_matrix, get_T_derivative_vector
from bsplinegenerator.bspline_to_bezier import get_bspline_to_bezier_conversion_matrix
from max_curvature_evaluators.root_finder import find_max_curvature_root_finder
from max_curvature_evaluators.max_numerator_over_min_denominator import find_curvature_using_max_numerator_over_min_denominator
from max_curvature_evaluators.control_point_method import get_control_point_curvature_bound
from bsplinegenerator.bspline_to_bezier import get_composite_bspline_to_bezier_conversion_matrix
from max_curvature_evaluators.helper_files.cube_root_solver import solver
# from bsplinegenerator

class PathGenerator:
    """ 
    This class generates a B-spline path from one point to another
    with some given direction vector at each point, and with some
    curvature constraint
    """

    def __init__(self, order, dimension, curvature_method):
        self._order = order
        self._dimension = dimension
        self._num_control_points = 7
        self._M = get_M_matrix(0, self._order, np.array([]), False)
        self._curvature_method = curvature_method
        if self._curvature_method == "constrain_max_acceleration_and_min_velocity":
            self._objective_function_type = "minimize_distance_and_time"
        else:
            self._objective_function_type = "minimize_acceleration"
        self._F_composite = get_composite_bspline_to_bezier_conversion_matrix(self._num_control_points, self._order)

    def generate_path(self, waypoints, velocities, max_curvature, initial_control_points = None):
        if self._curvature_method == "constrain_max_acceleration_and_min_velocity":
            return self.generate_trajectory_indirect(waypoints, velocities, max_curvature, initial_control_points)
        else:
            return self.generate_trajectory_direct(waypoints, velocities, max_curvature, initial_control_points)
        
    def generate_trajectory_indirect(self, waypoints, velocities, max_curvature, initial_control_points = None):
        # create initial conditions
        min_velocity = 0.5
        max_acceleration = max_curvature*min_velocity**2
        # max_cross_term_mag = max_curvature*min_velocity**3
        self._dimension = np.shape(waypoints)[0]
        if initial_control_points is None:
            initial_control_points = self.__create_initial_control_points(waypoints)
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor]))
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints)
        velocity_constraint = self.__create_waypoint_velocity_constraint(velocities)
        min_velocity_constraint = self.__create_min_velocity_constraint(min_velocity)
        max_acceleration_constraint = self.__create_maximum_acceleration_constraint(max_acceleration)
        # max_cross_term_constraint = self.__create_max_cross_term_constraint(max_cross_term_mag)
        objectiveFunction = self.__get_objective_function()
        objective_variable_bounds = self.__create_objective_variable_bounds()
        minimize_options = {'disp': True}#, 'maxiter': self.maxiter, 'ftol': tol}
        # perform optimization
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            method='SLSQP', 
            # method = 'trust-constr',
            bounds=objective_variable_bounds,
            constraints=(waypoint_constraint, velocity_constraint, \
                max_acceleration_constraint, \
                # max_cross_term_constraint, \
                min_velocity_constraint), 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x[0:-1] ,(self._dimension,self._num_control_points))
        optimized_scale_factor = result.x[-1]
        return optimized_control_points, optimized_scale_factor

    def generate_trajectory_direct(self, waypoints,velocities, max_curvature, initial_control_points = None):
        # create initial conditions
        self._dimension = np.shape(waypoints)[0]
        if initial_control_points is None:
            initial_control_points = self.__create_initial_control_points(waypoints)
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor]))
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints)
        velocity_constraint = self.__create_waypoint_velocity_constraint(velocities)
        curvature_constraint = self.__create_curvature_constraint(max_curvature)
        objectiveFunction = self.__get_objective_function()
        objective_variable_bounds = self.__create_objective_variable_bounds()
        minimize_options = {'disp': True}#, 'maxiter': self.maxiter, 'ftol': tol}
        # perform optimization
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            method='SLSQP', 
            # method = 'trust-constr',
            bounds=objective_variable_bounds,
            constraints=(waypoint_constraint, velocity_constraint, \
                curvature_constraint
                ), 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x[0:-1] ,(self._dimension,self._num_control_points))
        optimized_scale_factor = result.x[-1]
        return optimized_control_points, optimized_scale_factor

    def __get_objective_function(self):
        # return self.__minimize_acceleration_and_distance_objective_function
        # return self.__minimize_acceleration_and_time_objective_function
        # return self.__minimize_acceleration_objective_function
        if self._objective_function_type == "minimize_distance_and_time":
            # return self.__minimize_distance_and_time_objective_function
            # return self.__minimize_acceleration_objective_function
            return self.__minimize_control_point_distance
        elif self._objective_function_type == "minimize_acceleration":
            # return self.__minimize_distance_and_time_objective_function
            # return self.__minimize_acceleration_objective_function
            return self.__minimize_control_point_distance
            # return self.__minimize_distance_objective_function

    def __minimize_distance_and_time_objective_function(self, variables):
        # for third order splines only
        control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
        scale_factor = variables[-1]
        num_intervals = self._num_control_points - self._order
        sum_of_integrals = 0
        sum_of_integ = 0
        for i in range(num_intervals):
            p0 = control_points[:,i]
            p1 = control_points[:,i+1]
            p2 = control_points[:,i+2]
            p3 = control_points[:,i+3]
            a = (p0/2 - (3*p1)/2 + (3*p2)/2 - p3/2)**2
            b = -2*(p0 - 2*p1 + p2)*(p0/2 - (3*p1)/2 + (3*p2)/2 - p3/2)
            c = (p0 - 2*p1 + p2)**2 + 2*(p0/2 - p2/2)*(p0/2 - (3*p1)/2 + (3*p2)/2 - p3/2)
            d = -2*(p0/2 - p2/2)*(p0 - 2*p1 + p2)
            f =  (p0/2 - p2/2)**2
            integral = np.sum(a/5 + b/4 + c/3 + d/2 + f)
            sum_of_integrals += integral 
            p1x = control_points[0,i]
            p2x = control_points[0,i+1]
            p3x = control_points[0,i+2]
            p4x = control_points[0,i+3]
            p1y = control_points[1,i]
            p2y = control_points[1,i+1]
            p3y = control_points[1,i+2]
            p4y = control_points[1,i+3]
            integ = ((9*((p1x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/6 - (p2x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/2 + (p3x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/2 - (p4x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/6 + (p1y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/6 - (p2y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/2 + (p3y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/2 - (p4y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/6)))/5 + (- (6*((p1x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/2 - p2x*(p1x/6 - p2x/2 + p3x/2 - p4x/6) + (p3x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/2 + (p1y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/2 - p2y*(p1y/6 - p2y/2 + p3y/2 - p4y/6) + (p3y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/2)) - (6*((p1x*(p1x/2 - p2x + p3x/2))/6 - (p2x*(p1x/2 - p2x + p3x/2))/2 + (p3x*(p1x/2 - p2x + p3x/2))/2 - (p4x*(p1x/2 - p2x + p3x/2))/6 + (p1y*(p1y/2 - p2y + p3y/2))/6 - (p2y*(p1y/2 - p2y + p3y/2))/2 + (p3y*(p1y/2 - p2y + p3y/2))/2 - (p4y*(p1y/2 - p2y + p3y/2))/6)))/4 + ((4*((p1x*(p1x/2 - p2x + p3x/2))/2 - p2x*(p1x/2 - p2x + p3x/2) + (p3x*(p1x/2 - p2x + p3x/2))/2 + (p1y*(p1y/2 - p2y + p3y/2))/2 - p2y*(p1y/2 - p2y + p3y/2) + (p3y*(p1y/2 - p2y + p3y/2))/2)) + (3*((p1x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/2 - (p3x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/2 + (p1y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/2 - (p3y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/2)) + (3*((p1x*(p1x/2 - p3x/2))/6 - (p2x*(p1x/2 - p3x/2))/2 + (p3x*(p1x/2 - p3x/2))/2 - (p4x*(p1x/2 - p3x/2))/6 + (p1y*(p1y/2 - p3y/2))/6 - (p2y*(p1y/2 - p3y/2))/2 + (p3y*(p1y/2 - p3y/2))/2 - (p4y*(p1y/2 - p3y/2))/6)))/3 + (- (2*((p1x*(p1x/2 - p3x/2))/2 - p2x*(p1x/2 - p3x/2) + (p3x*(p1x/2 - p3x/2))/2 + (p1y*(p1y/2 - p3y/2))/2 - p2y*(p1y/2 - p3y/2) + (p3y*(p1y/2 - p3y/2))/2)) - (2*((p1x*(p1x/2 - p2x + p3x/2))/2 - (p3x*(p1x/2 - p2x + p3x/2))/2 + (p1y*(p1y/2 - p2y + p3y/2))/2 - (p3y*(p1y/2 - p2y + p3y/2))/2)))/2 + ((p1x*(p1x/2 - p3x/2))/2 - (p3x*(p1x/2 - p3x/2))/2 + (p1y*(p1y/2 - p3y/2))/2 - (p3y*(p1y/2 - p3y/2))/2)
            sum_of_integ += integ
        # print("control_points: " , control_points)
        # print("scale_factor: " , scale_factor)
        # print("answer: " , sum_of_integrals + scale_factor )
        # print(" ")
        return sum_of_integrals
    

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
        num_intervals = self._num_control_points - self._order
        sum_of_integrals = 0
        sum_of_integ = 0
        for i in range(num_intervals):
            p0 = control_points[:,i]
            p1 = control_points[:,i+1]
            p2 = control_points[:,i+2]
            p3 = control_points[:,i+3]
            p1x = control_points[0,i]
            p2x = control_points[0,i+1]
            p3x = control_points[0,i+2]
            p4x = control_points[0,i+3]
            p1y = control_points[1,i]
            p2y = control_points[1,i+1]
            p3y = control_points[1,i+2]
            p4y = control_points[1,i+3]
            p1z = 0
            p2z = 0
            p3z = 0
            p4z = 0
            integ = (36*((p1x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/6 - (p2x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/2 + (p3x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/2 - (p4x*(p1x/6 - p2x/2 + p3x/2 - p4x/6))/6 + (p1y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/6 - (p2y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/2 + (p3y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/2 - (p4y*(p1y/6 - p2y/2 + p3y/2 - p4y/6))/6 + (p1z*(p1z/6 - p2z/2 + p3z/2 - p4z/6))/6 - (p2z*(p1z/6 - p2z/2 + p3z/2 - p4z/6))/2 + (p3z*(p1z/6 - p2z/2 + p3z/2 - p4z/6))/2 - (p4z*(p1z/6 - p2z/2 + p3z/2 - p4z/6))/6))
            sum_of_integ += integ
            integral_vector = (-p0 + 3*p1 - 3*p2 + p3)
            sum_of_integrals += np.sum(integral_vector**2)
        # print("sum_of_integrals: " , sum_of_integrals)
        # print("sum_of_integ: ", sum_of_integ)
        return sum_of_integrals + scale_factor
    
    def __minimize_acceleration_and_time_objective_function(self, variables):
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
        # print("control_points: " , control_points)
        # print("scale_factor: " , scale_factor)
        # print("answer: " , sum_of_integrals + scale_factor )
        # print(" ")
        return sum_of_integrals + scale_factor
    
    def __minimize_control_point_distance(self, variables):
        # for third order splines only
        control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
        scale_factor = variables[-1]
        distance_vectors = control_points[:,1:] - control_points[:,0:-1]
        distances_squared = np.sum(distance_vectors**2,0)**2
        return np.sum(distances_squared)
    
    def __minimize_acceleration_and_distance_objective_function(self, variables):
        # for third order splines only
        control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
        scale_factor = variables[-1]
        num_intervals = self._num_control_points - self._order
        sum_of_acceleration_integrals = 0
        sum_of_distance_integrals = 0
        for i in range(num_intervals):
            p0 = control_points[:,i]
            p1 = control_points[:,i+1]
            p2 = control_points[:,i+2]
            p3 = control_points[:,i+3]
            sum_of_acceleration_integrals += np.sum((p0 - 3*p1 + 3*p2 - p3)**2) 
            a = (p0/2 - (3*p1)/2 + (3*p2)/2 - p3/2)**2
            b = -2*(p0 - 2*p1 + p2)*(p0/2 - (3*p1)/2 + (3*p2)/2 - p3/2)
            c = (p0 - 2*p1 + p2)**2 + 2*(p0/2 - p2/2)*(p0/2 - (3*p1)/2 + (3*p2)/2 - p3/2)
            d = -2*(p0/2 - p2/2)*(p0 - 2*p1 + p2)
            f =  (p0/2 - p2/2)**2
            sum_of_distance_integrals += np.sum(a/5 + b/4 + c/3 + d/2 + f)
        return sum_of_acceleration_integrals + sum_of_distance_integrals

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

    def __create_waypoint_velocity_constraint(self, velocities):
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

    def __create_maximum_acceleration_constraint(self, max_acceleration):
        def maximum_acceleration_constraint_function(variables):
            control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
                (self._dimension,self._num_control_points))
            scale_factor = variables[-1]
            velocity_control_points = (control_points[:,1:] - control_points[:,0:-1]) / scale_factor
            acceleration_control_points = (velocity_control_points[:,1:] - velocity_control_points[:,0:-1]) / scale_factor
            norm_acceleration_control_points = np.linalg.norm(acceleration_control_points,2,0)
            acceleration_constraints = norm_acceleration_control_points - max_acceleration
            return acceleration_constraints
        lower_bound = -np.inf
        upper_bound = 0
        max_velocity_constraint = NonlinearConstraint(maximum_acceleration_constraint_function, \
             lb = lower_bound, ub= upper_bound)
        return max_velocity_constraint

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
            if self._curvature_method == "roots_of_curvature_derivative":
                max_curvatures[i] = find_max_curvature_root_finder(control_points_per_interval,self._order,self._M)
            elif self._curvature_method == "roots_numerator_and_denominator":
                max_curvatures[i] = find_curvature_using_max_numerator_over_min_denominator(control_points_per_interval,self._order,self._M)
            elif self._curvature_method == "control_point_derivatives":
                max_curvatures[i] = get_control_point_curvature_bound(control_points_per_interval,self._order)
        return  max_curvatures

    def __create_min_velocity_constraint(self, min_velocity):
        def min_velocity_constraint_function(variables):
            control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
            scale_factor = variables[-1]
            min_velocity_of_spline = self.__get_min_velocity_of_spline(control_points,scale_factor)
            constraint = min_velocity - min_velocity_of_spline
            return constraint
        lower_bound = -np.inf
        upper_bound = 0
        min_velocity_constraint = NonlinearConstraint(min_velocity_constraint_function , lb = lower_bound, ub = upper_bound)
        return min_velocity_constraint

    def __get_min_velocity_of_spline(self, control_points, scale_factor):
        min_velocity = np.inf
        # print("#### new run #####")
        for i in range(self._num_control_points-self._order):
            interval_control_points = control_points[:,i:i+self._order+1]
            velocity = self.__find_min_velocity_magnitude(interval_control_points, self._order, self._M, scale_factor)
            if velocity < min_velocity:
                min_velocity = velocity
        return min_velocity

    def __create_max_cross_term_constraint(self, max_cross_term_mag):
        def max_cross_term_constraint_function(variables):
            control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
            scale_factor = variables[-1]
            max_cross_term_of_spline = self.__get_max_cross_term_of_spline(control_points,scale_factor)
            constraint = max_cross_term_of_spline - max_cross_term_mag
            return constraint
        lower_bound = -np.inf
        upper_bound = 0
        min_velocity_constraint = NonlinearConstraint(max_cross_term_constraint_function, lb = lower_bound, ub = upper_bound)
        return min_velocity_constraint
    
    def __get_max_cross_term_of_spline(self, control_points, scale_factor):
        max_cross_term_mag = np.inf
        # print("#### new run #####")
        for i in range(self._num_control_points-self._order):
            interval_control_points = control_points[:,i:i+self._order+1]
            cross_term_mag = self.__find_min_velocity_magnitude(interval_control_points, self._order, self._M, scale_factor)
            if cross_term_mag > max_cross_term_mag:
                max_cross_term_mag = cross_term_mag
        return max_cross_term_mag

    def __find_min_velocity_magnitude(self, control_points, order, M, scale_factor):
        P = control_points
        J = np.dot(np.dot(M.T,P.T) , np.dot(P,M))
        A = 36*J[0,0]
        B = 12*J[0,1] + 24*J[1,0]
        C = 8*J[1,1] + 12*J[2,0]
        D = 4*J[2,1]
        roots = solver(A,B,C,D)
        times_to_check = np.concatenate((roots,np.array([0,1])))
        min_velocity = np.inf
        for i in range(len(times_to_check)):
            t = times_to_check[i]
            if t >= 0 and t <= 1:
                velocity = self.__calculate_velocity_magnitude(t, M, control_points,order,scale_factor)
                if velocity < min_velocity:
                    min_velocity = velocity
        return min_velocity
    
    def __calculate_velocity_magnitude(self, t,M,control_points,order,scale_factor):
        dT = get_T_derivative_vector(order,t*scale_factor,0,1,scale_factor)
        velocity = np.dot(control_points,np.dot(M,dT)).flatten()
        velocity_magnitude = np.linalg.norm(velocity)
        return velocity_magnitude

    def __find_max_cross_term(self, control_points, order, M):
        A,B,C,D = self.__get_cross_coeficients(control_points)
        roots = solver(A,B,C,D)
        times_to_check = np.concatenate((roots,np.array([0,1])))
        max_cross_term = 0
        for i in range(len(times_to_check)):
            t = times_to_check[i]
            if t >= 0 and t <= 1:
                cross_term = self.__calculate_cross_term_magnitude(t,M,control_points,order)
                if cross_term > max_cross_term:
                    max_cross_term = cross_term
        return max_cross_term
    
    def __get_cross_coeficients(self, control_points):
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
    

    def __calculate_cross_term_magnitude(self,t,M,control_points,order):
        dT = get_T_derivative_vector(order,t,0,1,1)
        d2T = get_T_derivative_vector(order,t,0,2,1)
        velocity = np.dot(control_points,np.dot(M,dT)).flatten()
        acceleration = np.dot(control_points,np.dot(M,d2T)).flatten()
        cross_term_magnitude = np.linalg.norm(np.cross(velocity,acceleration))
        return cross_term_magnitude