
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator_2nd_order import PathGenerator2nd
import time

## try different initial conditions
## move the max velocity around
## try direction constraint
## add function to get path length
max_curvature = 1.5
order = 2
waypoints = np.array([[0,6],[0,0]])
dimension = np.shape(waypoints)[0]
velocities = np.array([[1,1],[0,0]]) # 0
velocities = np.array([[1,0],[0,1]]) # 1
# velocities = np.array([[0,0],[1,-1]]) # 2
# velocities = np.array([[-1,-1],[0,0]]) # 3
# velocities = np.array([[0,0],[1,1]]) # 4
# velocities = np.array([[-1,1],[0,0]]) # 5
velocities = velocities/np.linalg.norm(velocities,2,0) # normalize velocities
initial_control_points = None




path_gen = PathGenerator2nd(dimension)
start_time = time.time()
control_points, scale_factor = path_gen.generate_path(waypoints, velocities, max_curvature)
end_time = time.time()
print("eval time: " , end_time - start_time)
spline_start_time = 0
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
number_data_points = 1000000
spline_data, time_data = bspline.get_spline_data(number_data_points)
spline_at_knot_points, knot_points = bspline.get_spline_at_knot_points()
curvature_data, time_data = bspline.get_spline_curvature_data(number_data_points)
velocity_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,1)
acceleration_magnitude_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,2)

print("max_curvature: " , np.max(curvature_data))
plt.figure()
plt.plot(spline_data[0,:], spline_data[1,:])
plt.tight_layout()
plt.subplots_adjust(wspace=None, hspace=None)
plt.show()

plt.figure()
plt.plot(time_data, curvature_data)
plt.tight_layout()
plt.subplots_adjust(wspace=None, hspace=None)
plt.title("curvature")
plt.show()

