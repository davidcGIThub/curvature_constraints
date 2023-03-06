import numpy as np
from bsplinegenerator.bsplines import BsplineEvaluation
from path_following.matrix_evaluation import matrix_bspline_evaluation_for_dataset, get_M_matrix, \
    get_T_derivative_vector

class BsplinePathFollower:

    def __init__(self, order, path_gain = 1, distance_gain = 1):
        self._order = order
        self._path_gain = path_gain
        self._distance_gain = distance_gain
        self._initial_num_points_to_check_per_interval = 100
        self._num_points_to_check_per_interval = 10000

    def get_commands(self, control_points, position, desired_airspeed):
        closest_point, path_vector = self.get_closest_point_and_direction_vector(control_points, position)
        direction_desired = self.get_desired_direction_vector(closest_point, position, path_vector, desired_airspeed)
        course_angle_command = np.arctan2(direction_desired.item(1), direction_desired.item(0))
        climb_rate_command = desired_airspeed * direction_desired.item(2)
        airspeed_command = desired_airspeed
        return [course_angle_command, climb_rate_command, airspeed_command]

    def get_desired_direction_vector(self, closest_point, position, path_vector, desired_airspeed):
        distance_vector = closest_point - position
        desired_direction_vector = distance_vector*self._distance_gain + path_vector*desired_airspeed*self._path_gain
        desired_direction_vector = desired_direction_vector/ np.linalg.norm(desired_direction_vector)
        return desired_direction_vector

    def get_closest_point_and_direction_vector(self, control_points, position):
        control_point_set = self.__get_closest_control_point_set(control_points, position)
        closest_point, t = self.__get_closest_point_and_t(control_point_set, position)
        direction_vector = self.__get_direction_vector(t, control_point_set)
        return closest_point, direction_vector

    def __get_closest_control_point_set(self, control_points, position):
        num_control_points = self.__count_number_of_control_points(control_points)
        number_of_knot_points = num_control_points + self._order + 1
        knot_points = np.arange(number_of_knot_points) - self._order
        dataset = matrix_bspline_evaluation_for_dataset(control_points, \
                knot_points, self._initial_num_points_to_check_per_interval)
        distances = np.linalg.norm(position - dataset,2,0)
        end_knot_point = knot_points[num_control_points]
        intial_knot_point = int(np.argmin(distances)/len(distances)*end_knot_point)
        control_point_set = control_points[:,intial_knot_point:intial_knot_point+self._order+1]
        return control_point_set
    
    def __get_closest_point_and_t(self, control_points, position):
        num_control_points = self._order + 1
        number_of_knot_points = num_control_points + self._order + 1
        knot_points = np.arange(number_of_knot_points) - self._order
        dataset = matrix_bspline_evaluation_for_dataset(control_points, \
                    knot_points, self._num_points_to_check_per_interval)
        distances = np.linalg.norm(position - dataset,2,0)
        closest_point_index = np.argmin(distances)
        closest_point = dataset[:,closest_point_index][:,None]
        t = closest_point_index/len(distances)
        return closest_point, t
    
    def __get_direction_vector(self, t, control_point_set):
        M = get_M_matrix(0, self._order, [], False)
        T = get_T_derivative_vector(self._order,t,0,1,1)
        derivative_at_time_t = np.dot(control_point_set, np.dot(M,T))
        direction_vector = derivative_at_time_t/np.linalg.norm(derivative_at_time_t,2)
        return direction_vector

    def __count_number_of_control_points(self, control_points):
        if control_points.ndim == 1:
            number_of_control_points = len(control_points)
        else:
            number_of_control_points = len(control_points[0])
        return number_of_control_points

    



# add code for when outside the transition region???

#need to trim parths of path when get past certian parts of it

