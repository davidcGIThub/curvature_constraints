import numpy as np
import matplotlib.pyplot as plt
from bsplinegenerator.bsplines import BsplineEvaluation
from path_following.spline_path_follower import BsplinePathFollower


# control_points = np.array([[-3,-4,-2,-.5,1,0,2,3.5,3,5,6.5],
#                             [.5,3.5,6,5.5,3.7,2,-1,2,5,5.5,5]]) # 2 D

control_points = np.array([[-3,  -4, -2, -.5, 1  ,   0,  2, 3.5, 3],
                           [.5, 3.5,  6, 5.5, 3.7,   2, -1,   2, 5],
                           [ 1, 3.2,  5,   0, 3.3, 1.5, -1, 2.5, 4]]) # 3D


order = 3
dimension = np.shape(control_points)[0]
bspline = BsplineEvaluation(control_points, order, 0,1)
bspline_data, time_data = bspline.get_spline_data(1000)
point = np.random.randint(6, size=(dimension,1))
desired_airspeed = 2
# point = np.array([[2],[6]])

path_follower = BsplinePathFollower(order, path_gain=1, distance_gain=1)
closest_point, path_vector = path_follower.get_closest_point_and_direction_vector(control_points,point)
desired_direction = path_follower.get_desired_direction_vector(closest_point, point, path_vector, desired_airspeed)

print("path_vector: " , path_vector)
print("desired_direction: " , desired_direction)

if dimension == 3:
    fig = plt.figure()
    ax = plt.axes(projection ='3d')
    ax.plot(bspline_data[0,:],bspline_data[1,:],bspline_data[2,:])
    ax.scatter(point[0],point[1],point[2])
    ax.scatter(closest_point[0],closest_point[1],closest_point[2])
    ax.quiver(closest_point.item(0), closest_point.item(1), closest_point.item(2), \
                path_vector.item(0), path_vector.item(1), path_vector.item(2), length=1, normalize=True, color = "r")
    ax.quiver(point.item(0), point.item(1), point.item(2), \
                desired_direction.item(0), desired_direction.item(1), desired_direction.item(2), length=1, normalize=True, color = "k")
    plt.show()
else:
    plt.figure()
    plt.plot(bspline_data[0,:],bspline_data[1,:])
    plt.scatter(point[0],point[1])
    plt.scatter(closest_point[0],closest_point[1])
    plt.arrow(closest_point.item(0), closest_point.item(1), path_vector.item(0), \
              path_vector.item(1), head_width=0.2, head_length=0.2)
    plt.arrow(point.item(0), point.item(1), desired_direction.item(0), \
              desired_direction.item(1), head_width=0.2, head_length=0.2)
    plt.show()
