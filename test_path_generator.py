import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator import PathGenerator
import time

waypoint_directions = np.array([[0,0],[1,-1],[0,0]]) # 2
max_curvature = 1
order = 3
waypoints = np.array([[0,0],[0,400],[100,150]])
dimension = np.shape(waypoints)[0]

initial_control_points = None

curvature_method = "roots_numerator_and_denominator"
# curvature_method = "constrain_max_acceleration_and_min_velocity"


path_gen = PathGenerator(dimension)
start_time = time.time()
control_points = path_gen.generate_path(waypoints, waypoint_directions, max_curvature)
end_time = time.time()
print("computation time:" , end_time - start_time)
spline_start_time = 0
scale_factor = 1
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)

bspline.plot_spline(1000)
bspline.plot_curvature(1000)
