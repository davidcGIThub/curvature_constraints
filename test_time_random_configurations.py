import numpy as np
import time
from bsplinegenerator.bsplines import BsplineEvaluation

num_iterations = 100
order = 2

total_time = 0
for i in range(num_iterations):
    control_points = np.random.randint(10, size=(2,order+1)) #random
    start_time = time.time()
    max_curvature = find_curvature_bound()
    end_time = time.time()
    total_time += end_time - start_time
    bspline = BsplineEvaluation(control_points,order,0,1)
    curvature_data, time_data = bspline.get_spline_curvature_data(10000)
    true_max = np.max(curvature_data)
    if true_max < 1e-10 and max_curvature < 1e-10:
        accuracy = 1
    elif true_max < 1e-10:
        accuracy = 0
    elif true_max == np.inf and max_curvature == np.inf:
        accuracy = 1
    elif true_max == np.inf:
        accuracy = 0
    else:
        accuracy = max_curvature/true_max
    average_accuracy += accuracy/num_iterations
    print("average_accuracy: " , average_accuracy)
ave_time = total_time/num_iterations
print("ave_accuracy: " , average_accuracy)
print("ave_time: " , ave_time)


def find_curvature_bound(method):
    if method == "maximize_curvature_equation":
        pass
    elif method == "roots_of_curvature_derivative":
        pass
    elif method == "discrete_evaluations":
        pass
    elif method == "control_point_derivatives":
        pass
    elif method == "geometric method":
        pass
    elif method == "":
        pass
    elif method == "angle_constrained_control_points":
        pass