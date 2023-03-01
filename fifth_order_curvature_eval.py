import numpy as np
from bsplinegenerator.bsplines import BsplineEvaluation
from max_curvature_evaluators.control_point_method import get_control_point_curvature_bound
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import get_matrix

method = "max_numerator_over_min_denominator"
order = 5
# control_points = np.array([[2,4,5.7,8,8.2,9],[6,4,8,6,9,12],[2,3.4,5,7.5,6.9,9]])
# control_points = np.array([[2,4,5.7,8,8.2,9],[6,4,8,6,9,12]])
control_points = np.array([[1.01,2,3,4.05,2,0.98],[1,2.02,2.99,4,2.03,1]])
# control_points = np.array([[2,2,2,2,2,2],[3,3,3,3,3,3]])
# control_points = np.array([[6,2,2,2,2,2],[8,3,3,3,3,3]])
# control_points = np.array([[6,6,6,2,2,2],[8,8,8,3,3,3]])
# control_points = np.array([[1,2,3,4,2,1],[1,2,3,4,2,1]])
# control_points = np.array([[0,0,0,0,0,0],[1,2,3,4,2,1]])
M = get_matrix(order)
curvature_bound = get_control_point_curvature_bound(control_points, order)

bspline = BsplineEvaluation(control_points, order,0,1)
curvature_data,time_data = bspline.get_spline_curvature_data(1000000)
max_curvature = np.max(curvature_data)

print("curvature_bound: " , curvature_bound)
print("max curvature: " , max_curvature)
# bspline.plot_curvature(1000000)
# bspline.plot_spline(100000)
# bspline.plot_derivative_vs_time(10000,1)
# bspline.plot_derivative_vs_time(10000,2)
