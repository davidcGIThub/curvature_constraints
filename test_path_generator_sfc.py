import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator import PathGenerator
from path_generation.safe_flight_corridor import SFC, plot_2D_sfcs
import time

sfc_1 = SFC(np.array([[4],[6]]), np.array([[3],[5]]), np.eye(2))
sfc_2 = SFC(np.array([[10],[2]]), np.array([[6],[7]]), np.eye(2))
sfc_3 = SFC(np.array([[4],[8]]), np.array([[9],[10]]), np.eye(2))
sfc_4 = SFC(np.array([[16],[4]]), np.array([[15],[12]]), np.eye(2))
# sfcs = (sfc_1, sfc_2)
# sfcs = (sfc_1, sfc_2, sfc_3)
sfcs = (sfc_1, sfc_2, sfc_3, sfc_4)
# sfcs = None
# waypoint_directions = np.array([[0,0],[1,1]]) # 2
waypoint_directions = np.array([[0,1],[1,0]]) # 2
max_curvature = 1
order = 3
# waypoints = np.array([[2,9],[3,7.5]])
waypoints = np.array([[2,20],[3,12]])
# waypoints = np.array([[2,8],[3,13]])
dimension = np.shape(waypoints)[0]


initial_control_points = None

curvature_method = "roots_numerator_and_denominator"
# curvature_method = "constrain_max_acceleration_and_min_velocity"


path_gen = PathGenerator(dimension)
start_time = time.time()
control_points = path_gen.generate_path(waypoints, waypoint_directions, max_curvature, sfcs)
end_time = time.time()
print("computation time:" , end_time - start_time)
spline_start_time = 0
scale_factor = 1
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)

spline_data, time_data = bspline.get_spline_data(1000)
minvo_cps = bspline.get_minvo_control_points()

# print("sfcs: " , sfcs)
plt.figure()
plot_2D_sfcs(sfcs)
plt.plot(spline_data[0,:], spline_data[1,:])
plt.scatter(waypoints[0,:],waypoints[1,:])
# plt.scatter(minvo_cps[0,:],minvo_cps[1,:])
plt.show()
