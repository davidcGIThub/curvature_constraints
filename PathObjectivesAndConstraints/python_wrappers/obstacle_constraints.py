import ctypes 
import pathlib 
import os 
import numpy as np

script_dir = os.path.abspath(os.path.dirname(__file__))
libname_str = os.path.join(script_dir)
libname = pathlib.Path(libname_str)
lib = ctypes.CDLL(libname / "libPathObjectivesAndConstraints.so")

class ObstacleConstraints(object):

    def __init__(self, dimension):
        ND_POINTER_DOUBLE = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags="C")
        self._dimension = dimension
        if dimension == 2:
            lib.ObstacleConstraints_2.argtypes = [ctypes.c_void_p]
            lib.ObstacleConstraints_2.restype = ctypes.c_void_p
            lib.checkIfObstacleCollides_2.argtypes = [ctypes.c_void_p, ND_POINTER_DOUBLE, ctypes.c_int, 
                ctypes.c_double, ND_POINTER_DOUBLE]
            lib.checkIfObstacleCollides_2.restype = ctypes.c_bool
            lib.getObstacleDistanceToSpline_2.argtypes = [ctypes.c_void_p,ND_POINTER_DOUBLE, ctypes.c_int, 
                ctypes.c_double, ND_POINTER_DOUBLE]
            lib.getObstacleDistanceToSpline_2.restype = ctypes.c_double
            self.obj = lib.ObstacleConstraints_2(0)
        else: # value == 3
            lib.ObstacleConstraints_3.argtypes = [ctypes.c_void_p]
            lib.ObstacleConstraints_3.restype = ctypes.c_void_p
            lib.checkIfObstacleCollides_3.argtypes = [ND_POINTER_DOUBLE, ctypes.c_int, 
                ctypes.c_double, ND_POINTER_DOUBLE]
            lib.checkIfObstacleCollides_3.restype = ctypes.c_bool
            lib.getObstacleDistanceToSpline_3.argtypes = [ND_POINTER_DOUBLE, ctypes.c_int, 
                ctypes.c_double, ND_POINTER_DOUBLE]
            lib.getObstacleDistanceToSpline_3.restype = ctypes.c_double
            self.obj = lib.ObstacleConstraints_3(0)

    def getObstacleDistanceToSpline(self, cont_pts, obstacle_radius, obstacle_center):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float64')
        obstacle_center_array = obstacle_center.flatten().astype('float64')
        if self._dimension == 2:
            distance = lib.getObstacleDistanceToSpline_2(self.obj, cont_pts_array, num_cont_pts, obstacle_radius, obstacle_center_array)
        else: # value = 3
            distance = lib.getObstacleDistanceToSpline_3(self.obj, cont_pts_array, num_cont_pts, obstacle_radius, obstacle_center_array)
        return distance
    
    def checkIfObstacleCollides(self, cont_pts, obstacle_radius, obstacle_center):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float64')
        obstacle_center_array = obstacle_center.flatten().astype('float64')
        if self._dimension == 2:
            collides = lib.checkIfObstacleCollides_2(self.obj, cont_pts_array, num_cont_pts, obstacle_radius, obstacle_center_array)
        else: # value = 3
            collides = lib.checkIfObstacleCollides_3(self.obj, cont_pts_array, num_cont_pts, obstacle_radius, obstacle_center_array)
        return collides

# control_points = np.array([[7.91705873, 9.88263331, 0.27303466, 7.50604049, 4.61073475, 5.98801717, 1.52432928, 3.8850049, 1.61195392, 8.22471529],
#                            [5.22947263, 1.33282499, 3.51583204, 8.62435967, 3.03096953, 0.84672315, 0.54028843, 7.24686189, 4.79897482, 5.00498365]])
# radius = 1
# center = np.array([4.84435679, 6.42836434])
# obst_const = ObstacleConstraints(2)
# distance = obst_const.getObstacleDistanceToSpline(control_points, radius, center)
# print("distance: " , distance)

# collides = obst_const.checkIfObstacleCollides(control_points, radius, center)
# print("collides: " , collides)
