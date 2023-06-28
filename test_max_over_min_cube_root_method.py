import numpy as np

from max_curvature_evaluators.max_numerator_over_min_denominator import find_curvature_using_max_numerator_over_min_denominator
from bsplinegenerator.bsplines import BsplineEvaluation
from max_curvature_evaluators.helper_files.helper_curvature_evaluations import get_matrix
import matplotlib.pyplot as plt


order = 3
num_control_points = order+1
# dimension = 3
# control_points = np.random.randint(10, size=(dimension,order+1)) # random
control_points = np.array([[1.69984, 3.13518, 3.75944, 3.64158],
                            [4.23739,  3.8813, 4.23739, 5.36086]])
bspline = BsplineEvaluation(control_points,order,0,1)
M = get_matrix(order)

curvature_data, time_data = bspline.get_spline_curvature_data(10000000)
true_max_curvature = np.max(curvature_data)
max_curvature = find_curvature_using_max_numerator_over_min_denominator(control_points, order, M)

print("max_curvature: ", max_curvature)
print("true_max_curvature: ", true_max_curvature)
# bspline.plot_spline(10000)
# bspline.plot_curvature(100000)