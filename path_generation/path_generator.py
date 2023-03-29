"""
This module generates a 3rd order B-spline path between two waypoints,
waypoint directions, curvature constraint, and adjoining 
safe flight corridors.
"""
import os
import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, Bounds
from path_generation.matrix_evaluation import get_M_matrix
from PathObjectivesAndConstraints.python_wrappers.objective_functions import ObjectiveFunctions
from PathObjectivesAndConstraints.python_wrappers.curvature_constraints import CurvatureConstraints
from PathObjectivesAndConstraints.python_wrappers.obstacle_constraints import ObstacleConstraints
from PathObjectivesAndConstraints.python_wrappers.waypoint_velocity_constraints import WaypointConstraints
from path_generation.safe_flight_corridor import SFC_2D
from bsplinegenerator.bspline_to_minvo import get_composite_bspline_to_minvo_conversion_matrix



class PathGenerator:
    """
    This class generates a 3rd order B-spline path between two waypoints,
    waypoint directions, curvature constraint, and adjoining 
    safe flight corridors.
    """

# when minimizing distance and time need constraint over max velocity
# when minimizing acceleration, no max velocity constraint needed

    def __init__(self, dimension: int):
        self._dimension = dimension
        self._order = 3
        self._M = get_M_matrix(self._order)
        self._control_points_per_corridor = 4
        # 0 sfc = 5 intervals - 8
        # 1 sfc = 5 intervals - 8
        # 2 sfc = 6 intervals - 9 cps
        # 3 sfc = 8 intervals - 11 cps
        # 4 sfc = 10 intervals - 13 cps
        # 5 sfc = 12 intervals - 15 cps
        self._objective_func_obj = ObjectiveFunctions(self._dimension)
        self._curvature_const_obj = CurvatureConstraints(self._dimension)
        self._waypoint_vel_const_obj = WaypointConstraints(self._dimension)
        
    def generate_path(self, waypoints: np.ndarray, waypoint_directions: np.ndarray, max_curvature: np.float64, sfcs: list = None):
        if(sfcs == None):
            return self.generate_unbounded_path(waypoints, waypoint_directions, max_curvature)
        else:
            return self.generate_bounded_path(waypoints, waypoint_directions, max_curvature, sfcs)
    
    def generate_unbounded_path(self, waypoints, waypoint_directions, max_curvature):
        # create initial conditions
        num_cont_pts = self.get_num_control_points(0)
        self._dimension = np.shape(waypoints)[0]
        initial_control_points = self.__create_initial_control_points(waypoints, num_cont_pts)
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor]))
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints, num_cont_pts)
        velocity_constraint = self.__create_waypoint_direction_constraint(waypoint_directions, num_cont_pts)
        curvature_constraint = self.__create_curvature_constraint(max_curvature, num_cont_pts)
        objectiveFunction = self.__minimize_acceleration_objective_function
        objective_variable_bounds = self.__create_objective_variable_bounds(num_cont_pts)
        minimize_options = {'disp': True}#, 'maxiter': self.maxiter, 'ftol': tol}
        # perform optimization
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            args=(num_cont_pts,),
            method='SLSQP', 
            bounds=objective_variable_bounds,
            constraints=(\
            curvature_constraint, 
            velocity_constraint,
            waypoint_constraint), 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x[0:-1] ,(self._dimension,num_cont_pts))
        return optimized_control_points
    
    def generate_bounded_path(self, waypoints, waypoint_directions, max_curvature, sfcs):
        # create initial conditions
        num_cont_pts = self.get_num_control_points(len(sfcs))
        self._dimension = np.shape(waypoints)[0]
        initial_control_points = self.__create_initial_control_points(waypoints, num_cont_pts)
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor]))
        # define constraints and objective function and constraints
        waypoint_constraint = self.__create_waypoint_constraint(waypoints, num_cont_pts)
        velocity_constraint = self.__create_waypoint_direction_constraint(waypoint_directions, num_cont_pts)
        curvature_constraint = self.__create_curvature_constraint(max_curvature, num_cont_pts)
        sfc_constraint = self.__create_safe_flight_corridor_constraint(sfcs, num_cont_pts)
        objectiveFunction = self.__minimize_acceleration_objective_function
        objective_variable_bounds = self.__create_objective_variable_bounds(num_cont_pts)
        minimize_options = {'disp': True} #, 'maxiter': self.maxiter, 'ftol': tol}
        # perform optimization
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            args=(num_cont_pts,),
            method='SLSQP', 
            bounds=objective_variable_bounds,
            constraints=(\
            curvature_constraint, 
            velocity_constraint,
            waypoint_constraint,
            sfc_constraint), 
            options = minimize_options)
        # retrieve data
        optimized_control_points = np.reshape(result.x[0:-1] ,(self._dimension,num_cont_pts))
        return optimized_control_points

# bounds over just the scale factor
    def __create_objective_variable_bounds(self, num_cont_pts):
        lower_bounds = np.zeros(num_cont_pts*self._dimension + 1) - np.inf
        upper_bounds = np.zeros(num_cont_pts*self._dimension + 1) + np.inf
        lower_bounds[num_cont_pts*self._dimension] = 0.00001
        return Bounds(lb=lower_bounds, ub = upper_bounds)

    def __minimize_acceleration_objective_function(self, variables, num_cont_pts):
        # for third order splines only
        control_points = np.reshape(variables[0:num_cont_pts*self._dimension], \
            (self._dimension,num_cont_pts))
        scale_factor = variables[-1]
        return self._objective_func_obj.minimize_acceleration_and_time(control_points, scale_factor)

    def __create_initial_control_points(self, waypoints, num_cont_pts):
        start_waypoint = waypoints[:,0]
        end_waypoint = waypoints[:,1]
        control_points = np.linspace(start_waypoint,end_waypoint,num_cont_pts).T
        return control_points

    def __create_waypoint_constraint(self, waypoints, num_cont_pts):
        num_waypoints = 2
        m = num_waypoints
        n = num_cont_pts
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

    def __create_waypoint_direction_constraint(self, direction_vectors, num_cont_pts):
        def velocity_constraint_function(variables):
            control_points = np.reshape(variables[0:num_cont_pts*self._dimension], \
            (self._dimension,num_cont_pts))
            scale_factor = variables[-1]
            constraints = self._waypoint_vel_const_obj.velocity_at_waypoints_constraints(control_points,
                scale_factor, direction_vectors)
            return constraints.flatten()
        lower_bound = 0
        upper_bound = 0
        velocity_vector_constraint = NonlinearConstraint(velocity_constraint_function, lb= lower_bound, ub=upper_bound)
        return velocity_vector_constraint

    def __create_curvature_constraint(self, max_curvature, num_cont_pts):
        num_intervals = num_cont_pts - self._order
        def curvature_constraint_function(variables):
            control_points = np.reshape(variables[0:num_cont_pts*self._dimension], \
            (self._dimension,num_cont_pts))
            return self._curvature_const_obj.get_interval_curvature_constraints(control_points,max_curvature)
        lower_bound = np.zeros(num_intervals) - np.inf
        upper_bound = np.zeros(num_intervals)
        curvature_constraint = NonlinearConstraint(curvature_constraint_function , lb = lower_bound, ub = upper_bound)
        return curvature_constraint

    def __create_safe_flight_corridor_constraint(self, sfcs, num_cont_pts):
        # create the rotation matrix.
        num_corridors = len(sfcs)
        num_minvo_cont_pts = (num_cont_pts - self._order)*(self._order+1)
        intervals_per_corridor = self.get_intervals_per_corridor(num_corridors)
        M_rot = self.get_composite_sfc_rotation_matrix(intervals_per_corridor, sfcs, num_minvo_cont_pts)
        # create the bspline to minvo conversion matrix 
        M_minvo = get_composite_bspline_to_minvo_conversion_matrix(\
            num_cont_pts, self._order)
        zero_block = np.zeros((num_minvo_cont_pts,num_cont_pts))
        zero_col = np.zeros((num_minvo_cont_pts, 1))
        if self._dimension == 2:
            M_minvo = np.block([[M_minvo, zero_block, zero_col],
                                        [zero_block, M_minvo, zero_col]])
        conversion_matrix = M_rot @ M_minvo
        #create bounds
        lower_bounds = np.zeros((self._dimension, num_minvo_cont_pts))
        upper_bounds = np.zeros((self._dimension, num_minvo_cont_pts))
        index = 0
        for corridor_index in range(num_corridors):
            num_intervals = intervals_per_corridor[corridor_index]
            lower_bound, upper_bound = sfcs[corridor_index].getRotatedBounds()
            num_points = num_intervals*(self._order+1)
            lower_bounds[:,index:index+num_points] = lower_bound
            upper_bounds[:,index:index+num_points] = upper_bound
            index = index+num_points
        safe_corridor_constraints = LinearConstraint(conversion_matrix, lb=lower_bounds.flatten(), ub=upper_bounds.flatten())
        return safe_corridor_constraints
    
    def get_composite_sfc_rotation_matrix(self, intervals_per_corridor, sfcs, num_minvo_cont_pts):
        num_corridors = len(intervals_per_corridor)
        M_len = num_minvo_cont_pts*self._dimension
        M_rot = np.zeros((M_len, M_len))
        num_cont_pts_per_interval = self._order + 1
        interval_count = 0
        dim_step = num_minvo_cont_pts
        for corridor_index in range(num_corridors):
            rotation = sfcs[corridor_index].rotation.T
            num_intervals = intervals_per_corridor[corridor_index]
            for interval_index in range(num_intervals):
                for cont_pt_index in range(num_cont_pts_per_interval):
                    index = interval_count*num_cont_pts_per_interval+cont_pt_index
                    M_rot[index, index] = rotation[0,0]
                    M_rot[index, index + dim_step] = rotation[0,1]
                    M_rot[index + dim_step, index] = rotation[1,0]
                    M_rot[index + dim_step, index + dim_step] = rotation[1,1]
                    if self._dimension == 3:

                        M_rot[2*dim_step + index, index] = rotation[2,0]
                        M_rot[2*dim_step + index, index + dim_step] = rotation[2,1]
                        M_rot[2*dim_step + index, index + 2*dim_step] = rotation[2,2]
                        M_rot[dim_step + index, index + 2*dim_step] = rotation[1,2]
                        M_rot[index, index + 2*dim_step] = rotation[0,2]
                interval_count += 1
        return M_rot

    
    def get_intervals_per_corridor(self, num_corridors):
        if num_corridors == 1:
            # 5 intervals
            return (5)
        elif num_corridors == 2:
            # 6 intervals
            return (3,3)
        elif num_corridors == 3:
            # 8 intervals
            return (3,2,3)
        elif num_corridors == 4:
            # 10 intervals
            return (3,2,2,3)
        elif num_corridors == 5:
            # 12 intervals
            return (3,2,2,2,3)
        
    def get_num_control_points(self, num_corridors):
        if num_corridors < 2:
            num_cont_pts = 8
        else:
            num_cont_pts = np.sum(self.get_intervals_per_corridor(num_corridors)) + self._order
        return num_cont_pts
        