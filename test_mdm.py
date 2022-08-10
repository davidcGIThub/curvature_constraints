import numpy as np
from max_curvature_evaluators.mdm_algorithm_adapted import MDM
import matplotlib.pyplot as plt
import time

points = np.array([[1,4,6,3.5,3.6],[3,6,4,4.2,4.1]])
# points =  np.array([[2., 4., 5., 9., 4.],[1., 0., 3., 8., 4.]])
points = np.array([[ 2., -5., -5.,  1.,  0.],[-3.,  2.,  2.,  3.,  1.]])
dimension = 2
num_points = np.random.randint(10)
points = np.random.randint(10, size=(dimension,num_points))*1.0-5
points = np.array([[ 4.,  4., 4.],[-3., -3., -3.]])
print("points: " , points)
start_time = time.time()
mdm = MDM(points,dimension,True)
solution = mdm.solve()
print("solution: ", solution)
end_time = time.time()
print("time: " , end_time - start_time)
plt.scatter(points[0,:],points[1,:],label="points")
plt.plot([],[])
plt.scatter(0,0)
plt.plot([0,solution[0]],[0,solution[1]])
plt.show()
