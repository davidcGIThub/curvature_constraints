import matplotlib.pyplot as plt
from scipy.optimize import minimize, NonlinearConstraint, Bounds
from bsplinegenerator.bsplines import BsplineEvaluation
import numpy as np
import time


def min_norm(point_set):
    num_points = np.shape(point_set)[1]
    optimization_variables = np.zeros(num_points)
    optimization_variables[0] = 1
    def objectie_func(variables):
        return np.sum(np.dot(point_set,variables.T)**2)
    def constraint_func(variables):
        return np.sum(np.dot(np.eye(num_points),variables.T))
    constraint = NonlinearConstraint(constraint_func, lb= 1, ub=1)
    variable_bounds = Bounds(lb=0, ub = 1)
    result = minimize(
        objectie_func,
        x0=optimization_variables,
        method='SLSQP', 
        bounds=variable_bounds,
        constraints=(constraint))
    optimized_coeficients = np.array(result.x)
    closest_point = np.dot(point_set,optimized_coeficients.T)
    inHull = (np.sum(closest_point) < 0.0000001)
    return optimized_coeficients, closest_point, inHull

n_points = 10
n_dim = 2
radius = 1
point_set = np.random.rand(n_dim,n_points)*10
point = np.random.rand(n_dim,1)*10 - 2
order = 3

bspline = BsplineEvaluation(point_set, order, 0,1)
minvo_cps = bspline.get_minvo_control_points()

point_set_translated = point_set - point
coeficients, closest_point, in_hull = min_norm(point_set_translated)

print("point_set: " , point_set)
print("point: " , point)
# print("dist to point: " , np.linalg.norm(closest_point,2))
print("dist to object: " , np.linalg.norm(closest_point,2) - radius)
print("radius: " , radius)



circle = plt.Circle((point.item(0), point.item(1)), radius, color='b',fill=False)
fig, ax = plt.subplots() 
ax.add_patch(circle)
ax.scatter(point.item(0), point.item(1))
ax.scatter(point_set[0,:], point_set[1,:], color = "g")
ax.scatter(closest_point.item(0) + point.item(0), closest_point.item(1) + point.item(1), color="r")
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