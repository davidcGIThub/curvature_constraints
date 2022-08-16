from time import time
from tracemalloc import start
import numpy as np
import matplotlib.pyplot as plt
from bsplinegenerator.bsplines import BsplineEvaluation
from scipy.optimize import minimize, Bounds, LinearConstraint
from max_curvature_evaluators.mdm_algorithm_adapted import MDM
import time
## Control Points ###

order = 3
dimension = 2
num_control_points = order+1
control_points = np.random.randint(100, size=(dimension,num_control_points)) # random
print("control_points: " , control_points)

if len(control_points) == 1:
    control_points = control_points.flatten()

### Parameters
# order = 3
start_time = 0
scale_factor = 1
clamped = False
number_data_points = 1000

### Create B-Spline Object ###
bspline = BsplineEvaluation(control_points, order, start_time, scale_factor, clamped)

####  Evaluate B-Spline Data ###
start_time = time.time()
spline_data, time_data = bspline.get_spline_data(number_data_points)
true_min = np.min(np.linalg.norm(spline_data,2,0))
end_time = time.time()
true_time = end_time - start_time

#### Dr. Beard Method ####
start_time = time.time()
one = np.ones(order+1)
cp_inv = np.linalg.pinv(np.dot(control_points.T,control_points))
numerator = np.dot(control_points,np.dot(cp_inv,one))
denominator = np.dot(one,np.dot(cp_inv,one))
coeficients = np.dot(cp_inv,one) / denominator
beard_bound = np.linalg.norm(numerator/denominator)
end_time = time.time()
beard_time = end_time - start_time

### sqp bound ###
start_time = time.time()
def objective_function(coeficients):
    objective = np.dot(np.dot(coeficients.T,np.dot(control_points.T,control_points)),coeficients)
    return objective
def objective_derivative(coeficients):
    jacobian = 2*np.dot(coeficients.T,np.dot(control_points.T,control_points))
    return jacobian
constraint_matrix = np.ones(num_control_points)
constraint = LinearConstraint(constraint_matrix, lb=1,ub=1)
bounds = Bounds(lb=0.0, ub =1.0)
result = minimize(objective_function, x0=np.ones(num_control_points)/num_control_points, \
            method='SLSQP', bounds=bounds, jac=objective_derivative, constraints=(constraint))
sqp_min_bound = np.sqrt(result.fun)
end_time = time.time()
sqp_time = end_time - start_time

### mdm bound ###
start_time = time.time()
mdm = MDM(np.transpose(control_points), dimension, True)
mdm_min = np.linalg.norm(mdm.get_closest_point())
end_time = time.time()
mdm_time = end_time - start_time

### conservative bounds ####
start_time = time.time()
max_x = np.max(control_points[0,:])
max_y = np.max(control_points[1,:])
min_x = np.min(control_points[0,:])
min_y = np.min(control_points[1,:])
if max_x > 0 and min_x < 0:
    x_dist = 0
else:
    x_dist = np.min((np.abs(max_x),np.abs(min_x)))
if max_y > 0 and min_y < 0:
    y_dist = 0
else:
    y_dist = np.min((np.abs(max_y),np.abs(min_y)))
conservative_bound = np.sqrt(x_dist**2 + y_dist**2)
end_time = time.time()
conservative_time = end_time - start_time

print("mdm_time: " , mdm_time)
print("sqp_min_time: " , sqp_time)
print("conservative time: " , conservative_time)
print("beard_time: " , beard_time)
print("true_time: " , true_time)
print(" ")
print("mdm_min: " , mdm_min)
print("sqp_min_bound: " , sqp_min_bound)
print("conservative bound: " , conservative_bound)
print("beard_bound: " , beard_bound)
print("true_min: " , true_min)


plt.plot(time_data, np.linalg.norm(spline_data,2,0),label = "norm_of_spline")
plt.scatter(np.linspace(0,time_data[-1],order+1),np.linalg.norm(control_points,2,0),label="control_point_norm")
plt.plot(time_data,time_data*0+true_min, label="true min")
plt.plot(time_data,time_data*0+beard_bound, label="beard bound")
plt.plot(time_data,time_data*0+sqp_min_bound, label="sqp min bound" )
plt.plot(time_data,time_data*0+conservative_bound, label="conservative bound")

plt.legend()
plt.show()

# bspline.plot_spline(number_data_points)




