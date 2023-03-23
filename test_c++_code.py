import ctypes 
import pathlib 
import os 
import numpy as np

script_dir = os.path.abspath(os.path.dirname(__file__))
libname_str = os.path.join(script_dir)
libname = pathlib.Path(libname_str)
lib = ctypes.CDLL(libname / "libPathObjectivesAndConstraints.so")

    # WaypointConstraints<2>* WaypointConstraints_2(){return new WaypointConstraints<2>();}
    # float* velocity_at_waypoints_2(WaypointConstraints<2>* obj, float cont_pts[], int num_control_points,
    #         float scale_factor, float desired_velocities[]){return obj->velocity_at_waypoints_constraints(
    #         cont_pts, num_control_points, scale_factor, desired_velocities);}

    # WaypointConstraints<3>* WaypointConstraints_3(){return new WaypointConstraints<3>();}
    # float* velocity_at_waypoints_3(WaypointConstraints<3>* obj, float cont_pts[], int num_control_points,
    #         float scale_factor, float desired_velocities[]){return obj->velocity_at_waypoints_constraints(
    #         cont_pts, num_control_points, scale_factor, desired_velocities);}

class WaypointConstraints(object):

    def __init__(self, dimension):
        ND_POINTER_PY_FLOATS = np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,flags="C")
        ND_POINTER_C_FLOATS = np.ctypeslib.ndpointer(dtype=ctypes.c_float, shape=(dimension,2))
        self._dimension = dimension
        if dimension == 2:
            lib.WaypointConstraints_2.argtypes = [ctypes.c_void_p]
            lib.WaypointConstraints_2.restype = ctypes.c_void_p
            lib.velocity_at_waypoints_constraints_2.argtypes = [ctypes.c_void_p, 
                ND_POINTER_PY_FLOATS, ctypes.c_int, ctypes.c_float, ND_POINTER_PY_FLOATS]
            lib.velocity_at_waypoints_constraints_2.restype = ND_POINTER_C_FLOATS
            self.obj = lib.WaypointConstraints_2(0)
        else: # value == 3
            lib.WaypointConstraints_3.argtypes = [ctypes.c_void_p]
            lib.WaypointConstraints_3.restype = ctypes.c_void_p
            lib.velocity_at_waypoints_constraints_3.argtypes = [ctypes.c_void_p, 
                ND_POINTER_PY_FLOATS, ctypes.c_int, ctypes.c_float, ND_POINTER_PY_FLOATS]
            lib.velocity_at_waypoints_constraints_3.restype = ND_POINTER_C_FLOATS
            self.obj = lib.WaypointConstraints_3(0)

    def velocity_at_waypoints_constraints(self, cont_pts, scale_factor, desired_velocities):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float32')
        des_vel_array = desired_velocities.flatten().astype('float32')
        if self._dimension == 2:
            objective = lib.velocity_at_waypoints_constraints_2(self.obj, cont_pts_array, num_cont_pts, scale_factor, des_vel_array)
        else: # value = 3
            objective = lib.velocity_at_waypoints_constraints_3(self.obj, cont_pts_array, num_cont_pts, scale_factor, des_vel_array)
        return objective
    

# control_points = np.array([[7.91705873, 9.88263331, 0.27303466, 7.50604049, 4.61073475, 5.98801717, 1.52432928, 3.8850049, 1.61195392, 8.22471529],
#                            [5.22947263, 1.33282499, 3.51583204, 8.62435967, 3.03096953, 0.84672315, 0.54028843, 7.24686189, 4.79897482, 5.00498365]])
control_points = np.array([[2, 3, 2, 7, 8, 5, 6, 2],
                           [6, 8, 3, 0, 3, 8, 2, 4],
                           [9, 2, 0, 9, 8, 2, 9, 3]])
desired_velocities = np.array([[1,5],
                               [2,-4]])
desired_velocities = desired_velocities / np.linalg.norm(desired_velocities,2,0)
scale_factor = 1
vel_const = WaypointConstraints(3)
vel_at_waypoints_constraints = vel_const.velocity_at_waypoints_constraints(control_points, scale_factor, desired_velocities)
print("vel_at_waypoints_constraints: " , vel_at_waypoints_constraints)
