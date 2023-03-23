import ctypes 
import pathlib 
import os 
import numpy as np

script_dir = os.path.abspath(os.path.dirname(__file__))
libname_str = os.path.join(script_dir)
libname = pathlib.Path(libname_str)
lib = ctypes.CDLL(libname / "libPathObjectivesAndConstraints.so")

class CurvatureConstraints(object):

    def __init__(self, dimension):
        ND_POINTER_PY_FLOATS = np.ctypeslib.ndpointer(dtype=np.float32, ndim=1,flags="C")
        self._dimension = dimension
        if dimension == 2:
            lib.CurvatureConstraints_2.argtypes = [ctypes.c_void_p]
            lib.CurvatureConstraints_2.restype = ctypes.c_void_p
            lib.get_spline_curvature_constraint_2.argtypes = [ctypes.c_void_p, 
                ND_POINTER_PY_FLOATS, ctypes.c_int, ctypes.c_float]
            lib.get_spline_curvature_constraint_2.restype = ctypes.c_float
            self.obj = lib.CurvatureConstraints_2(0)
        else: # value == 3
            lib.CurvatureConstraints_3.argtypes = [ctypes.c_void_p]
            lib.CurvatureConstraints_3.restype = ctypes.c_void_p
            lib.get_spline_curvature_constraint_3.argtypes = [ctypes.c_void_p, 
                ND_POINTER_PY_FLOATS, ctypes.c_int, ctypes.c_float]
            lib.get_spline_curvature_constraint_3.restype = ctypes.c_float
            self.obj = lib.CurvatureConstraints_3(0)

    def get_spline_curvature_constraint(self, cont_pts, max_curvature):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float32')
        if self._dimension == 2:
            objective = lib.get_spline_curvature_constraint_2(self.obj, cont_pts_array, num_cont_pts, max_curvature)
        else: # value = 3
            objective = lib.get_spline_curvature_constraint_3(self.obj, cont_pts_array, num_cont_pts, max_curvature)
        return objective
    
# control_points = np.array([[4, 1, 4, 5],
#                            [2, 2, 0, 4],
#                            [7, 0, 1, 7]])
# max_curvature = 2
# curve_const = CurvatureConstraints(3)
# curvature_constraint = curve_const.get_spline_curvature_constraint(control_points, max_curvature)
# print("curvature_constraint: " , curvature_constraint)
