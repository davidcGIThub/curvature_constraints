import numpy as np
import random
import matplotlib.pyplot as plt
from bsplinegenerator.bsplines import BsplineEvaluation

order = 3
control_points = np.random.randint(10, size=(random.randint(2, 3),order+1)) # random
num_data_points = 1000

if len(control_points) == 1:
    control_points = control_points.flatten()

bspline = BsplineEvaluation(control_points,order,0,1)
curvature_data, time_data = bspline.get_spline_curvature_data(num_data_points)
velocity_data, time_data = bspline.get_derivative_magnitude_data(num_data_points,1)
acceleration_data, time_data = bspline.get_derivative_magnitude_data(num_data_points,2)
# bspline.plot_derivative_magnitude(num_data_points,1)
# bspline.plot_derivative_magnitude(num_data_points,2)
plt.plot(time_data, velocity_data , label="vel")
plt.plot(time_data, acceleration_data , label="accel")
plt.plot(time_data, curvature_data , label="curv")
plt.legend()
plt.show()
