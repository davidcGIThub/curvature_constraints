import numpy as np
from bsplinegenerator.bsplines import BsplineEvaluation
from max_curvature_evaluators.helper_files.max_curvature_test_functions import test_performance_of_max_curvature_function

num_iterations = 100
order = 2

isRestricted = False
# isRestricted = True
# method = "maximize_curvature_equation"
# method = "roots_of_curvature_derivative"
# method = "discrete_evaluations"
# method = "curvature_at_min_velocity"
method = "max_numerator_over_min_denominator"
isRestricted = False
average_accuracy, std_accuracy, average_time, std_time = \
    test_performance_of_max_curvature_function(order,isRestricted,method,num_iterations)
print("average_accuracy: ", average_accuracy)
print("std_accuracy: ", std_accuracy)
print("average_time: ", average_time)
print("std_time: ", std_time)
