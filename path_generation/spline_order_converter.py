"""
This module generates a 3rd order B-spline path between two waypoints,
waypoint directions, curvature constraint, and adjoining 
safe flight corridors.
"""
import os
import numpy as np
from scipy.optimize import minimize, NonlinearConstraint

class SmoothingSpline:
    """
    This class generates a new spline from a previous one
    """

    def __init__(self, order, dimension, resolution):
        self._dimension = dimension
        self._resolution = resolution # points spline
        self._order = order
        
    def generate_new_control_points(self, old_control_points, old_scale_factor, old_order):
        num_cont_pts = self.__get_num_control_points(old_control_points, old_order)
        initial_control_points = self.create_initial_control_points(old_control_points, old_order, num_cont_pts)
        scale_factor = self.__get_new_scale_factor(initial_control_points, old_scale_factor, old_control_points, old_order)
        objective_function = self.__get_objective_function(num_cont_pts, old_order, old_control_points, old_scale_factor,scale_factor)
        point_constraint = self.__get_point_constraints(num_cont_pts, old_order, old_control_points,old_scale_factor, scale_factor)
        result = minimize(
            objective_function,
            x0=initial_control_points.flatten(),
            constraints=(point_constraint),
            method='SLSQP')
        optimized_control_points = np.reshape(result.x[0:num_cont_pts*self._dimension] ,(self._dimension,num_cont_pts))
        return optimized_control_points, scale_factor
    
    def __get_objective_function(self, num_cont_pts, old_order, old_control_points, old_scale_factor, scale_factor):
        old_points = self.__matrix_bspline_evaluation_for_dataset(old_order, old_control_points, self._resolution)
        def smoother(variables):
            control_points = self.__get_objective_variables(variables, num_cont_pts)
            points = self.__matrix_bspline_evaluation_for_dataset(self._order, control_points, self._resolution)
            points_difference = (old_points - points)**2
            return np.sum(points_difference)
        return smoother
    
    def __get_point_constraints(self, num_cont_pts, old_order, old_control_points,old_scale_factor, scale_factor):
        # constraining the initial and final location, velocity, and acceleration
        num_pts = 2
        old_points = self.__matrix_bspline_evaluation_for_dataset(old_order, old_control_points, num_pts)
        old_velocity_points = self.__matrix_bspline_derivative_evaluation_for_dataset(old_order, 1, old_scale_factor, old_control_points, num_pts)
        old_acceleration_points = self.__matrix_bspline_derivative_evaluation_for_dataset(old_order, 2, old_scale_factor, old_control_points, num_pts)
        def point_constraint_function(variables):
            control_points = self.__get_objective_variables(variables, num_cont_pts)
            points = self.__matrix_bspline_evaluation_for_dataset(self._order, control_points, num_pts)
            points_difference = (old_points - points)
            velocity_points = self.__matrix_bspline_derivative_evaluation_for_dataset(self._order, 1, scale_factor, control_points, num_pts)
            velocity_difference = (old_velocity_points - velocity_points)
            acceleration_points = self.__matrix_bspline_derivative_evaluation_for_dataset(self._order, 2, scale_factor, control_points, num_pts)
            acceleration_difference = (old_acceleration_points - acceleration_points)
            return np.concatenate((points_difference.flatten(), velocity_difference.flatten(), acceleration_difference.flatten()))
        lower_bound = 0
        upper_bound = 0
        point_constraint = NonlinearConstraint(point_constraint_function, lower_bound, upper_bound)
        return point_constraint

    def __get_objective_variables(self, variables, num_cont_pts):
        control_points = np.reshape(variables[0:num_cont_pts*self._dimension], \
                    (self._dimension,num_cont_pts))
        return control_points

    def __get_num_control_points(self, old_control_points, old_order):
        old_num_intervals = self.__count_number_of_control_points(old_control_points) - old_order
        new_num_intervals = int(old_num_intervals*2.5)
        new_num_control_points = new_num_intervals + self._order
        return new_num_control_points
    
    def __get_new_scale_factor(self, control_points, old_scale_factor, old_control_points, old_order):
        old_num_intervals = self.__count_number_of_control_points(old_control_points) - old_order
        end_time = old_num_intervals*old_scale_factor
        num_intervals = self.__count_number_of_control_points(control_points) - self._order
        scale_factor = end_time/ num_intervals
        return scale_factor

    def create_initial_control_points(self, old_pts, old_order, num_cont_pts):
        old_num_pts = self.__count_number_of_control_points(old_pts)
        distances = np.linalg.norm(old_pts[:,1:] - old_pts[:,0:-1],2,0)
        old_num_segments = old_num_pts - 1
        num_segments = num_cont_pts - 1
        for i in range(old_num_segments-1):
            distances[i+1] = distances[i+1] + distances[i]
        distance_between_cont_pts = distances[old_num_segments-1] / (num_segments)
        old_segment_num = 0
        current_distance = 0
        prev_point_location = old_pts[:,0]
        step_distance = 0
        new_cont_pts = np.zeros((self._dimension, num_cont_pts))
        for i in range(num_segments):
            old_interval_start_point = old_pts[:,old_segment_num]
            old_interval_end_point = old_pts[:,old_segment_num+1]
            vector_to_point = old_interval_end_point - old_interval_start_point
            unit_vector_to_point = vector_to_point / (np.linalg.norm(vector_to_point))
            new_cont_pts[:,i] = prev_point_location + unit_vector_to_point*step_distance
            prev_point_location = new_cont_pts[:,i]
            step_distance = distance_between_cont_pts
            current_distance = current_distance + step_distance
            if distances[old_segment_num] < current_distance:
                temp = np.copy(distances)
                temp = temp - current_distance
                temp[temp < 0] = np.inf
                old_segment_num = np.argmin(temp)
                step_distance = current_distance - distances[old_segment_num-1]
                prev_point_location = old_pts[:,old_segment_num]
        new_cont_pts[:,-1] = old_pts[:,-1]
        return new_cont_pts
    
    def __matrix_bspline_evaluation_for_dataset(self, order, control_points, num_points):
        """
        This function evaluates the B spline for a given time data-set
        """
        #initialize variables
        dimension = self.__get_dimension(control_points)
        number_of_control_points = self.__count_number_of_control_points(control_points)
        num_intervals = number_of_control_points - order
        #create steps matrix
        time_data = np.linspace(0,num_intervals,num_points)
        # Find M matrix
        M = self.__get_M_matrix(order)
        #Evaluate spline data
        spline_data = np.zeros((dimension,num_points))
        marker = 0
        for i in range(num_intervals):
            P = control_points[:,i:i+order+1]
            if i == num_intervals - 1:
                steps_array = time_data[(time_data >= i) & (time_data <= i+1)] - i
            else:
                steps_array = time_data[(time_data >= i) & (time_data < i+1)] - i
            num_point_interval = len(steps_array)
            L = np.ones((order+1,num_point_interval))
            for i in range(order+1):
                L[i,:] = steps_array**(order-i)
            spline_data_over_interval = np.dot(np.dot(P,M),L)
            spline_data[:,marker:marker+num_point_interval] = spline_data_over_interval
            marker = marker + num_point_interval
        return spline_data

    def __matrix_bspline_derivative_evaluation_for_dataset(self, order, derivative_order, scale_factor, control_points, num_points):
        """
        This function evaluates the B spline for a given time data-set
        """
        # Initialize variables
        dimension = self.__get_dimension(control_points)
        number_of_control_points = self.__count_number_of_control_points(control_points)
        num_intervals = number_of_control_points - order
        #create steps matrix
        time_data = np.linspace(0,num_intervals,num_points)
        # Find M matrix
        M = self.__get_M_matrix(order)
        K = self.__create_k_matrix(order,derivative_order,scale_factor)
        # Evaluate Spline data
        marker = 0
        spline_derivative_data = np.zeros((dimension,num_points))
        for i in range(num_intervals):
            P = control_points[:,i:i+order+1]
            # Find M matrix if clamped
            if i == num_intervals - 1:
                steps_array = time_data[(time_data >= i) & (time_data <= i+1)] - i
            else:
                steps_array = time_data[(time_data >= i) & (time_data < i+1)] - i
            num_point_interval = len(steps_array)
            L_r = np.zeros((order+1,num_point_interval))
            for i in range(order-derivative_order+1):
                L_r[i,:] = steps_array**(order-derivative_order-i)
            spline_derivative_data_over_interval = np.dot(np.dot(P,M),np.dot(K,L_r))
            spline_derivative_data[:,marker:marker+num_point_interval] = spline_derivative_data_over_interval
            marker = marker + num_point_interval
        return spline_derivative_data

    def __create_k_matrix(self,order,derivative_order,scale_factor):
        K = np.zeros((order+1,order+1))
        for i in range(order-derivative_order+1):
            K[i,i] = np.math.factorial(order-i)/np.math.factorial(order-derivative_order-i)
        K = K/scale_factor**(derivative_order)
        return K

    def __get_M_matrix(self, order):
        if order > 5:
            print("Error: Cannot compute higher than 5th order matrix evaluation")
            return None
        if order == 0:
            return 1
        if order == 1:
            M = self.__get_1_order_matrix()
        if order == 2:
            M = self.__get_2_order_matrix()
        elif order == 3:
            M = self.__get_3_order_matrix()
        elif order == 4:
            M = self.__get_4_order_matrix()
        elif order == 5:
            M = self.__get_5_order_matrix()
        else:
            raise Exception("Cannot return M matrix for spline of order " , order)
        return M

    def __get_T_derivative_vector(self, order,t,tj,rth_derivative,scale_factor):
        T = np.zeros((order+1,1))
        t_tj = t-tj
        for i in range(order-rth_derivative+1):
            T[i,0] = (t_tj**(order-rth_derivative-i))/(scale_factor**(order-i)) * np.math.factorial(order-i)/np.math.factorial(order-i-rth_derivative)
        return T

    def __get_T_vector(self, order,t,tj,scale_factor):
        T = np.ones((order+1,1))
        t_tj = t-tj
        for i in range(order+1):
            if i > order:
                T[i,0] = 0
            else:
                T[i,0] = (t_tj/scale_factor)**(order-i)
        return T

    def __get_1_order_matrix(self):
        M = np.array([[-1,1],
                        [1,0]])
        return M

    def __get_2_order_matrix(self):
        M = .5*np.array([[1,-2,1],
                            [-2,2,1],
                            [1,0,0]])
        return M

    def __get_3_order_matrix(self):
        M = np.array([[-2 ,  6 , -6 , 2],
                        [ 6 , -12 ,  0 , 8],
                        [-6 ,  6 ,  6 , 2],
                        [ 2 ,  0 ,  0 , 0]])/12
        return M

    def __get_4_order_matrix(self):
        M = np.array([[ 1 , -4  ,  6 , -4  , 1],
                        [-4 ,  12 , -6 , -12 , 11],
                        [ 6 , -12 , -6 ,  12 , 11],
                        [-4 ,  4  ,  6 ,  4  , 1],
                        [ 1 ,  0  ,  0 ,  0  , 0]])/24
        return M

    def __get_5_order_matrix(self):
        M = np.array([[-1  ,  5  , -10 ,  10 , -5  , 1],
                        [ 5  , -20 ,  20 ,  20 , -50 , 26],
                        [-10 ,  30 ,  0  , -60 ,  0  , 66],
                        [ 10 , -20 , -20 ,  20 ,  50 , 26],
                        [-5  ,  5  ,  10 ,  10 ,  5  , 1 ],
                        [ 1  ,  0  ,  0  ,  0  ,  0  , 0]])/120
        return M


    def __get_dimension(self, control_points):
        if control_points.ndim == 1:
            dimension = 1
        else:
            dimension = len(control_points)
        return dimension

    def __count_number_of_control_points(self, control_points):
        if control_points.ndim == 1:
            number_of_control_points = len(control_points)
        else:
            number_of_control_points = len(control_points[0])
        return number_of_control_points