import ctypes 
import pathlib 
import os 
import numpy as np

script_dir = os.path.abspath(os.path.dirname(__file__))
libname_str = os.path.join(script_dir)
libname = pathlib.Path(libname_str)
lib = ctypes.CDLL(libname / "libPathObjectivesAndConstraints.so")

class WaypointConstraints(object):

    def __init__(self, dimension):
        ND_POINTER_DOUBLE = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags="C")
        ND_POINTER_C_DOUBLE = np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(dimension,2))
        self._dimension = dimension
        if dimension == 2:
            lib.WaypointConstraints_2.argtypes = [ctypes.c_void_p]
            lib.WaypointConstraints_2.restype = ctypes.c_void_p
            lib.velocity_at_waypoints_constraints_2.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ctypes.c_double, ND_POINTER_DOUBLE]
            lib.velocity_at_waypoints_constraints_2.restype = ND_POINTER_C_DOUBLE
            self.obj = lib.WaypointConstraints_2(0)
        else: # value == 3
            lib.WaypointConstraints_3.argtypes = [ctypes.c_void_p]
            lib.WaypointConstraints_3.restype = ctypes.c_void_p
            lib.velocity_at_waypoints_constraints_3.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ctypes.c_double, ND_POINTER_DOUBLE]
            lib.velocity_at_waypoints_constraints_3.restype = ND_POINTER_C_DOUBLE
            self.obj = lib.WaypointConstraints_3(0)

    def velocity_at_waypoints_constraints(self, cont_pts, scale_factor, desired_velocities):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float64')
        des_vel_array = desired_velocities.flatten().astype('float64')
        if self._dimension == 2:
            objective = lib.velocity_at_waypoints_constraints_2(self.obj, cont_pts_array, num_cont_pts, scale_factor, des_vel_array)
        else: # value = 3
            objective = lib.velocity_at_waypoints_constraints_3(self.obj, cont_pts_array, num_cont_pts, scale_factor, des_vel_array)
        return objective
    
# control_points = np.array([[3, 0, 9, 5, 7, 8, 0, 1],
#                             [1, 6, 2, 2, 1, 4, 1, 4]])
# desired_velocities = np.array([[7, 3],
#                                [1, 4]])
# #### desired_velocities = desired_velocities / np.linalg.norm(desired_velocities,2,0)
# scale_factor = 1
# vel_const = WaypointConstraints(2)
# vel_at_waypoints_constraints = vel_const.velocity_at_waypoints_constraints(control_points, scale_factor, desired_velocities)
# print("vel_at_waypoints_constraints: " , vel_at_waypoints_constraints)
