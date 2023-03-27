import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator_obst import PathGenerator
import time

waypoint_directions = np.array([[0,0],[1,-1],[0,0]]) # 2
max_curvature = 1
order = 3
waypoints = np.array([[0,0],[0,400],[100,150]])
dimension = np.shape(waypoints)[0]
obstacle_radii = np.array([20])
obstacle_centers = np.array([[3,],[200],[125]])
obstacles = (obstacle_centers, obstacle_radii)

initial_control_points = None

curvature_method = "roots_numerator_and_denominator"
# curvature_method = "constrain_max_acceleration_and_min_velocity"


path_gen = PathGenerator(dimension)
start_time = time.time()
control_points = path_gen.generate_path(waypoints, waypoint_directions, max_curvature, obstacles)
end_time = time.time()
print("computation time:" , end_time - start_time)
spline_start_time = 0
scale_factor = 1
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)

# bspline.plot_spline(1000)
# bspline.plot_curvature(1000)
# spline_data, time_data = bspline.get_spline_data(1000)


radius = obstacle_radii[0]
center = obstacle_centers.flatten()
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = radius  * np.outer(np.cos(u), np.sin(v)) + center[0]
y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
z = radius  * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]

# Plot the surface
ax.plot_surface(x, y, z)
# ax.plot(spline_data[0,:],spline_data[1,:],spline_data[2,:])
# Set an equal aspect ratio
ax.set_aspect('auto')

plt.show()
