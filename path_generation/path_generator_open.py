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
from max_curvature_evaluators.max_at_min_velocity_finder import find_curvature_at_min_velocity_magnitude
from bsplinegenerator.bspline_to_bezier import get_composite_bspline_to_bezier_conversion_matrix
from max_curvature_evaluators.helper_files.cube_root_solver import solver
# from bsplinegenerator

class PathGenerator:
    """ 
    This class generates a B-spline path from one point to another
    with some given direction vector at each point, and with some
    curvature constraint
    """

    def __init__(self, order, dimension, curvature_method, objective_function_type):
        # 3rd order spline 5 intervals - 8 control points
        # 4th order spline 4 intervals - 8 control points
        # 5th order spline 3 intervals - 8 control points
        self._order = order
        self._dimension = dimension
        self._num_control_points = 8
        self._M = get_M_matrix(0, self._order, np.array([]), False)
        self._curvature_method = curvature_method
        self._objective_function_type = objective_function_type
        self._F_composite = get_composite_bspline_to_bezier_conversion_matrix(self._num_control_points, self._order)

    def generate_trajectory(self, waypoints,velocities, max_velocity, \
                max_curvature, center_sfc, dimensions_sfc, initial_control_points = None):
        # create initial conditions
        min_velocity = 0.5
        max_acceleration = max_curvature*min_velocity**2
        self._dimension = np.shape(waypoints)[0]
        if initial_control_points is None:
            initial_control_points = self.__create_initial_control_points(waypoints)
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor]))
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints)
        velocity_constraint = self.__create_waypoint_velocity_constraint(velocities)
        safe_corridor_constraint = self.__create_safe_corridor_constraint(center_sfc, dimensions_sfc)
        direction_constraint = self.__create_waypoint_direction_constraint(velocities)
        # max_velocity_constraint = self.__create_maximum_velocity_constraint(max_velocity)
        min_velocity_constraint = self.__create_min_velocity_constraint(min_velocity)
        max_acceleration_constraint = self.__create_maximum_acceleration_constraint(max_acceleration)
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
                max_acceleration_constraint, \
                min_velocity_constraint, \
                # curvature_constraint
                ), 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x[0:-1] ,(self._dimension,self._num_control_points))
        optimized_scale_factor = result.x[-1]
        return optimized_control_points, optimized_scale_factor

    def __get_objective_function(self):
        if self._objective_function_type == "minimize_control_point_distance_and_time":
            return self.__minimize_control_point_distance_and_time_objective_function
        elif self._objective_function_type == "minimize_distance_and_time":
            return self.__minimize_distance_objective_function
        elif self._objective_function_type == "minimize_acceleration":
            return self.__minimize_acceleration_objective_function
        elif self._objective_function_type == "minimize_velocity":
            return self.__minimize_velocity_objective_function
        elif self._objective_function_type == "maximize_velocity":
            return self.__maximize_velocity_objective_function
        else:
            return self.__minimize_control_point_distance_and_time_objective_function

    def __minimize_control_point_distance_and_time_objective_function(self,variables):
        control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
        distances = self.__get_distances_between_points(control_points)
        scale_factor = variables[-1]
        return np.sum(distances) + scale_factor

    def __create_objective_variable_bounds(self):
        lower_bounds = np.zeros(self._num_control_points*self._dimension + 1) - np.inf
        upper_bounds = np.zeros(self._num_control_points*self._dimension + 1) + np.inf
        lower_bounds[self._num_control_points*self._dimension] = 0.00001
        return Bounds(lb=lower_bounds, ub = upper_bounds)

    def __get_distances_between_points(self,points):
        number_of_points = np.shape(points)[1]
        first_points = points[:,0:number_of_points-1]
        next_points = points[:,1:number_of_points]
        distances = np.sqrt(np.sum(((next_points - first_points)**2),0))
        return distances

    def __minimize_distance_objective_function(self, variables):
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
            a = (p0/2 - (3*p1)/2 + (3*p2)/2 - p3/2)**2
            b = -2*(p0 - 2*p1 + p2)*(p0/2 - (3*p1)/2 + (3*p2)/2 - p3/2)
            c = (p0 - 2*p1 + p2)**2 + 2*(p0/2 - p2/2)*(p0/2 - (3*p1)/2 + (3*p2)/2 - p3/2)
            d = -2*(p0/2 - p2/2)*(p0 - 2*p1 + p2)
            f =  (p0/2 - p2/2)**2
            integral = np.sum(a/5 + b/4 + c/3 + d/2 + f)
            sum_of_integrals += integral 
        return sum_of_integrals + scale_factor

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
        return sum_of_integrals #+ scale_factor

    def __minimize_velocity_objective_function(self, variables):
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
            a = (p0 - 3*p1 + 3*p2 - p3)**2
            b = (-2*(p0 - 2*p1 + p2)*(p0 - 3*p1 + 3*p2 - p3))
            c = (p0 - 2*p1 + p2)**2
            sum_of_integrals += np.sum(a/3 + b/2 + c)
        return sum_of_integrals + 1/scale_factor

    def __maximize_velocity_objective_function(self, variables):
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
            a = (p0 - 3*p1 + 3*p2 - p3)**2
            b = (-2*(p0 - 2*p1 + p2)*(p0 - 3*p1 + 3*p2 - p3))
            c = (p0 - 2*p1 + p2)**2
            sum_of_integrals += np.sum(a/3 + b/2 + c)
        return -(sum_of_integrals + 1/scale_factor)

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

    def __create_maximum_velocity_constraint(self, max_velocity):
        def maximum_velocity_constraint_function(variables):
            control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
                (self._dimension,self._num_control_points))
            scale_factor = variables[-1]
            velocity_control_points = (control_points[:,1:] - control_points[:,0:-1]) / scale_factor
            norm_velocity_control_points = np.linalg.norm(velocity_control_points,2,0)
            velocity_constraints = norm_velocity_control_points - max_velocity
            return velocity_constraints
        lower_bound = -np.inf
        upper_bound = 0
        max_velocity_constraint = NonlinearConstraint(maximum_velocity_constraint_function, \
             lb = lower_bound, ub= upper_bound)
        return max_velocity_constraint       

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
            elif self._curvature_method == "curvature_at_min_velocity":
                max_curvatures[i] = find_curvature_at_min_velocity_magnitude(control_points_per_interval, self._order, self._M)
        return  max_curvatures

    def __create_waypoint_direction_constraint(self, direction):
        def direction_constraint_function(variables):
            constraints = np.zeros(2)
            control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
            (self._dimension,self._num_control_points))
            segement_1_control_points = control_points[:,0:self._order+1]
            segement_2_control_points = control_points[:,self._num_control_points-self._order-1:]
            scale_factor = variables[-1]
            T_0 = get_T_derivative_vector(self._order,0,0,1,scale_factor)
            T_f = get_T_derivative_vector(self._order,scale_factor,0,1,scale_factor)
            start_velocity = np.dot(segement_1_control_points,np.dot(self._M,T_0)).flatten()
            end_velocity = np.dot(segement_2_control_points,np.dot(self._M,T_f)).flatten()
            start_angle = np.arctan2(start_velocity[1],start_velocity[0])
            end_angle = np.arctan2(end_velocity[1],end_velocity[0])
            desired_start_angle = np.arctan2(direction[1,0],direction[0,0])
            desired_end_angle = np.arctan2(direction[1,1],direction[0,1])
            constraints[0] = start_angle - desired_start_angle
            constraints[1] = end_angle - desired_end_angle
            # print("constraints: " , constraints)
            return constraints
        lower_bound = 0
        upper_bound = 0
        direction_constraint = NonlinearConstraint(direction_constraint_function, lb= lower_bound, ub=upper_bound)
        return direction_constraint

    def __create_safe_corridor_constraint(self, center, dimensions):
        def safe_corridor_constraint_function(variables):
            bspline_control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
                (self._dimension,self._num_control_points))
            bezier_control_points = np.transpose(np.dot(self._F_composite, np.transpose(bspline_control_points)))
            control_points = bezier_control_points
            constraints = (np.abs(control_points - center) - dimensions/2).flatten()
            return constraints
        lower_bound = -np.inf
        upper_bound = 0
        safe_corridor_constraint = NonlinearConstraint(safe_corridor_constraint_function, lb = lower_bound, ub = upper_bound)
        return safe_corridor_constraint

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
            # print("velocity: " , velocity)
            if velocity < min_velocity:
                min_velocity = velocity
        return min_velocity

    def __calculate_velocity_magnitude(self, t,M,control_points,order,scale_factor):
        dT = get_T_derivative_vector(order,t*scale_factor,0,1,scale_factor)
        velocity = np.dot(control_points,np.dot(M,dT)).flatten()
        velocity_magnitude = np.linalg.norm(velocity)
        return velocity_magnitude

    def __find_min_velocity_magnitude(self, control_points, order, M, scale_factor):
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
                velocity = self.__calculate_velocity_magnitude(t, M, control_points,order,scale_factor)
                if velocity < min_velocity:
                    min_velocity = velocity
        return min_velocity