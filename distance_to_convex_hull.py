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
    inHull = (np.linalg.norm(closest_point,2) < 0.00001)
    return optimized_coeficients, closest_point, inHull

n_points = 10
n_dim = 2
radius = 1
# control_point_set= np.random.rand(n_dim,n_points)*10
control_point_set = np.array([[8.33665694, 8.74600799, 0.15411698, 1.03944841, 2.71745767, 5.10167975, 0.49748729, 9.24934927, 2.09060337, 0.3531816],
                            [1.40865962, 5.09002261, 8.24619956, 0.06025564, 7.38338304, 6.72558684, 7.26305315, 8.93062143, 2.91758113, 6.65393249]])
point = np.random.rand(n_dim,1)*10 - 2
point = np.array([[2.6662532], [6.16021986]])
order = 3

bspline = BsplineEvaluation(control_point_set, order, 0,1)
spline_data, time_data = bspline.get_spline_data(1000)
minvo_cps = bspline.get_minvo_control_points()
distance = 1000000000000
closest_minvo_points = np.array([])
closest_point = []

bspline_point_set_translated = control_point_set - point
coeficients_, closest_point_, in_hull_ = min_norm(bspline_point_set_translated)
object_in_hull = np.linalg.norm(closest_point_,2) < radius
# print("object in hull: " , object_in_hull)
if (not object_in_hull):
    closest_point = closest_point_
    mean_bspline_points = np.mean(control_point_set,1)[:,None]
    pseudo_radius = np.max(np.linalg.norm(mean_bspline_points - control_point_set,2,0))
    # print("#### pseudo_radius: ", pseudo_radius)
    distance = np.linalg.norm(closest_point,2) + pseudo_radius - radius
else:
    for i in range(n_points - order):
        points = minvo_cps[:,i*(order+1):i*(order+1)+order+1]
        point_set_translated = points - point
        mean_points = np.mean(points,1)[:,None]
        mean_to_point = np.linalg.norm(mean_points - point,2)
        pseudo_radius = np.max(np.linalg.norm((mean_points - points),2,0))
        coeficients, current_closest_point, in_hull = min_norm(point_set_translated)
        current_distance = np.linalg.norm(current_closest_point,2)
        # print("in_hull: " , in_hull)
        # print("current distance 1: " , current_distance)
        if in_hull:
            current_distance = current_distance - radius - pseudo_radius + mean_to_point
        else:
            current_distance = current_distance - radius
            # print("current distance 2: " , current_distance)
            # print("radius: " , radius)
            # print("point not in hull")

        if current_distance < distance:
            distance = current_distance
            closest_minvo_points = points
            closest_point = current_closest_point
            print("in hull: " , in_hull)
            print("distance: " , distance)
            print("current_distance: " , np.linalg.norm(current_closest_point,2))
            # print("point_set_translated: " , point_set_translated)
            # print("current_closest_point: " , current_closest_point)
            # print("distance to point: " , np.linalg.norm(current_closest_point,2))

print("control_point_set: " , control_point_set)
print("closest_minvo_points: " , closest_minvo_points)
print("point: " , point)
print("dist to object: " , distance)
# print("radius: " , radius)

circle = plt.Circle((point.item(0), point.item(1)), radius, color='b',fill=False)
fig, ax = plt.subplots() 
ax.add_patch(circle)
ax.scatter(point.item(0), point.item(1))
ax.scatter(control_point_set[0,:], control_point_set[1,:], facecolors='none', edgecolors='g')
if object_in_hull:
    ax.scatter(closest_minvo_points[0,:], closest_minvo_points[1,:], color = "y")
ax.plot(spline_data[0,:], spline_data[1,:],color = 'g')
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