import numpy as np

from max_curvature_evaluators.max_numerator_over_min_denominator import find_curvature_using_max_numerator_over_min_denominator
from bsplinegenerator.bsplines import BsplineEvaluation
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import get_matrix
import matplotlib.pyplot as plt


order = 3
num_control_points = order+1
# dimension = 3
# control_points = np.random.randint(10, size=(dimension,order+1)) # random
control_points = np.array([[1, 3, 2, 0],[1, 3, 2, 0]])
bspline = BsplineEvaluation(control_points,order,0,1)
M = get_matrix(order)

curvature_data, time_data = bspline.get_spline_curvature_data(10000000)
true_max_curvature = np.max(curvature_data)
max_curvature = find_curvature_using_max_numerator_over_min_denominator(control_points, order, M)

print("control_points: " , control_points)
print("max_curvature: ", max_curvature)
print("true_max_curvature: ", true_max_curvature)
bspline.plot_spline(10000)
bspline.plot_curvature(100000)