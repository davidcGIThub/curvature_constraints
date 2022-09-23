import numpy as np

from max_curvature_evaluators.helper_files.helper_curvature_evaluations import get_matrix
from max_curvature_evaluators.max_at_min_velocity_finder import find_curvature_at_min_velocity_magnitude
from bsplinegenerator.bsplines import BsplineEvaluation
import matplotlib.pyplot as plt
import random


order = 3
num_control_points = order+1
dimension = 3
control_points = np.random.randint(10, size=(dimension,num_control_points)) # random
bspline = BsplineEvaluation(control_points,order,0,1)
velocity_data, time_data = bspline.get_derivative_magnitude_data(1000,1)
acceleration_data, time_data = bspline.get_derivative_magnitude_data(1000,2)
curvature_data, time_data = bspline.get_spline_curvature_data(1000)

M = get_matrix(order)

max_curvature = find_curvature_at_min_velocity_magnitude(control_points, order, M)

curvature_data, time_data = bspline.get_spline_curvature_data(1000)

plt.figure()
plt.plot(time_data, curvature_data)
plt.plot(time_data, time_data*0 + max_curvature)
plt.show()