from max_curvature_evaluators.helper_files.mdm_algorithm_adapted import MDM
import numpy as np
import time

num_points = 6
dim = 2
points = np.random.randint(10, size=(dim,num_points)) + 5.0 # random
# points = np.array([[5,5,5,5,5,5],[6, 6, 6, 6, 6, 6]]) * 1.0
points = np.array([[9.0, 6, 12, 14, 13, 7], [5, 13, 11, 11, 8, 8]])
max_iter = 5000
init_approx_index = 1
mdm_object_2 = MDM(points, dim, max_iter, init_approx_index)
start_time = time.time()
min_distance = mdm_object_2.get_min_distance()
end_time = time.time()
# print("solve time none: " , end_time - start_time)
# mdm_object_2.plot_convex_hull_with_solution()
