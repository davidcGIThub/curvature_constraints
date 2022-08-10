import numpy as np
from max_curvature_evaluators.mdm_algorithm_adapted import MDM
import matplotlib.pyplot as plt
import time

# points =  np.array([[ 1.,  3.,  2.,  0.],[ 2.,  4., -4.,  0.]])
dimension = 2
num_points = np.random.randint(10)
points = np.random.randint(10, size=(dimension,num_points))*1.0-5
print("points: " , points)
start_time = time.time()
mdm = MDM(points,dimension,True)
solution = mdm.get_closest_point()
print("solution: ", solution)
end_time = time.time()
print("time: " , end_time - start_time)
plt.scatter(points[0,:],points[1,:],label="points")
plt.plot([],[])
plt.scatter(0,0)
plt.plot([0,solution[0]],[0,solution[1]])
plt.show()
