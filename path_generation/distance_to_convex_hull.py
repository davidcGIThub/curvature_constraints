import matplotlib.pyplot as plt
from scipy.optimize import minimize, NonlinearConstraint, Bounds
import numpy as np
import time


def min_norm(point_set):
    # create initial conditions
    # dimension = np.shape(point_set)[0]
    num_points = np.shape(point_set)[1]
    optimization_variables = np.zeros(num_points)
    optimization_variables[0] = 1
    # define constraints and objective function and constraints
    def objectie_func(variables):
        return np.sum(np.dot(point_set,variables.T)**2)
    def constraint_func(variables):
        return np.sum(np.dot(np.eye(num_points),variables.T))
    constraint = NonlinearConstraint(constraint_func, lb= 1, ub=1)
    variable_bounds = Bounds(lb=0, ub = 1)
    # minimize_options = {'disp': True}#, 'maxiter': self.maxiter, 'ftol': tol}
    # perform optimization
    result = minimize(
        objectie_func,
        x0=optimization_variables,
        method='SLSQP', 
        bounds=variable_bounds,
        constraints=(constraint))
    # retrieve data
    optimized_coeficients = np.array(result.x)
    closest_point = np.dot(point_set,optimized_coeficients.T)
    inHull = (np.sum(closest_point) < 0.0000001)
    return optimized_coeficients, closest_point, inHull

n_points = 10
n_dim = 2
# point_set = np.random.rand(n_dim,n_points)*10 - 2
# point_set = np.array([[8, 5, 13, 10,  8, 14],
#                    [13, 9, 10, 10, 12, 10]])
point_set = np.array([[9,  6, 12, 14, 13,  7, 10,  8, 14,  5],
                      [5, 13, 11, 11,  8,  8,  6, 10,  8,  9],
                      [6,  5,  4,  3,  8,  5, 13, 10,  8, 14]])

# x = np.random.rand(n_dim)
start_time = time.time()
coeficients, closest_point, in_hull = min_norm(point_set)
eval_time = time.time() - start_time


print("coeficients: " , coeficients)
print("closest_point: " , closest_point)
print("dist: " , np.linalg.norm(closest_point,2))
print("in_hull: " , in_hull)
print("eval_time: ", eval_time)
# plt.scatter(x[0],x[1])
plt.scatter(point_set[0,:], point_set[1,:])
plt.show()




# def in_hull(points, x):
#     n_points = len(points)
#     c = np.zeros(n_points)
#     A = np.r_[points.T,np.ones((1,n_points))]
#     b = np.r_[x, np.ones(1)]
#     lp = linprog(c, A_eq=A, b_eq=b)
#     return lp.success, lp.fun
# n_points = 10
# n_dim = 2
# Z = np.random.rand(n_points,n_dim)
# x = np.random.rand(n_dim)
# success, variables = in_hull(Z, x)
# print("success: " , success)
# print("variables: " , variables)