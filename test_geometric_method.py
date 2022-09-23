import numpy as np
from bsplinegenerator.helper_functions import count_number_of_control_points
import time
from max_curvature_evaluators.geometric_max_finder import get_max_curvature
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import calculate_curvature, get_matrix
from bsplinegenerator.bsplines import BsplineEvaluation
from beziercurvegenerator.bezier_curves import BezierCurve
import matplotlib.pyplot as plt

order = 2
control_points = np.random.randint(10, size=(2,order+1)) # random
control_points = np.array([[5, 6, 4],[3, 6, 0]])
# control_points = np.array([[1, 1, 7],[1, 1, 7]])
print("control_points: " , control_points)

bspline = BsplineEvaluation(control_points,order,0,1)
bspline_data, time_data_bspline = bspline.get_spline_data(1000)
bspline_curvature_data, time_data = bspline.get_spline_curvature_data(10000)
F = np.array([[1,1,0],[0,2,0],[0,1,1]])/2
bezier_control_points = np.dot(F,control_points.T).T
print("bezier_control_points: " , bezier_control_points)
max_curvature = get_max_curvature(control_points)
print("max_curvature: ", max_curvature)
bezier_curve = BezierCurve(bezier_control_points,0,1)
curve_data, time_data = bezier_curve.get_curve_data(1000)
bezier_curvature_data, time_data = bezier_curve.get_curvature_data(1000)
print("bezier_max_curvature: " , np.max(bezier_curvature_data))
print("bspline_max_curvature: " , np.max(bspline_curvature_data))
# plt.plot(curve_data[0,:],curve_data[1,:])
# plt.plot(bspline_data[0,:],bspline_data[1,:])
# plt.show()
# bspline.plot_spline_vs_time(1000)
bspline.plot_curvature(1000)
bezier_curve.plot_curvature(1000)
bspline.plot_derivative_magnitude(1000,1)
bspline.plot_derivative_magnitude(1000,2)
# bezier_curve.plot_curvature(1000)
# bezier_curve.plot_curvature(1000)
