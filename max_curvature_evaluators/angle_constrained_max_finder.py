import numpy as np
import matplotlib.pyplot as plt
from bsplinegenerator.bsplines import BsplineEvaluation

order = 3
num_control_points = order+1
isBezier = True

if order == 1:
    angle = np.pi/2
elif order == 2:
    angle = np.pi/2
elif order == 3:
    angle = np.pi/4
elif order == 4:
    angle = np.pi/6
elif order == 5:
    angle = np.pi/8

# angle = np.pi/2
length = 3

def create_random_control_points_greater_than_angles(num_control_points,angle,length):
    control_points = np.zeros((2,num_control_points))
    for i in range(num_control_points):
        if i == 0:
            control_points[:,i][:,None] = np.array([[0],[0]])
        elif i == 1:
            random_vec = np.random.rand(2,1)
            next_vec = length*random_vec/np.linalg.norm(random_vec)
            control_points[:,i][:,None] = control_points[:,i-1][:,None] + next_vec
        else:
            new_angle = angle*2*(0.5-np.random.rand())
            R = np.array([[np.cos(new_angle), -np.sin(new_angle)],[np.sin(new_angle), np.cos(new_angle)]])
            prev_vec = control_points[:,i-1][:,None] - control_points[:,i-2][:,None]
            unit_prev_vec = prev_vec/np.linalg.norm(prev_vec)
            next_vec = length*np.dot(R,unit_prev_vec)#*np.random.rand()
            control_points[:,i][:,None] = control_points[:,i-1][:,None] + next_vec
    return control_points

def get_curvature_bound(control_points,isBezier):
    if order == 2:
        constant = 0.35355817313510857
        if isBezier:
            constant = 0.7071078439714923
    elif order == 3:
        constant = 0.7885805074747374
        if isBezier:
            constant = 1.0928304616410691
    elif order == 4:
        constant = 0.9117477044041745
        if isBezier:
            constant = 1.121134469464793
    elif order == 5:
        constant = 0.9436121693107355
        if isBezier:
            constant = 1.1126375193309979
    angles, lengths = get_angles_and_lengths(control_points)
    possible_bounds = np.sin(np.pi-angles)/(length*constant)
    max_curvature_bound = np.max(possible_bounds)
    return max_curvature_bound


def get_angles_and_lengths(control_points):
    vectors = control_points[:,1:] - control_points[:,0:-1]
    next_vectors = vectors[:,1:]
    prev_vectors = -vectors[:,0:-1]
    norm_vectors = np.linalg.norm(vectors,2,0)
    next_norm_vectors = norm_vectors[1:]
    prev_norm_vectors = norm_vectors[0:-1]
    numerator = np.sum(prev_vectors*next_vectors,0)
    denominator = prev_norm_vectors*next_norm_vectors
    angles = np.arccos(numerator/denominator)
    lengths = norm_vectors
    return angles, lengths

control_points = create_random_control_points_greater_than_angles(num_control_points, angle,length)
max_curvature_bound = get_curvature_bound(control_points,isBezier)

spline = BsplineEvaluation(control_points,order,0,1,isBezier)
spline.plot_spline(1000)

curvature, time_data = spline.get_spline_curvature_data(1000)
max_curvature = np.max(curvature)
# print("constant " , np.sin(angle)/(length)/max_curvature)

# print(np.sin(np.pi/2)/(length)/max_curvature)

plt.plot(time_data, curvature)
plt.plot(time_data, time_data*0 + max_curvature_bound)
plt.show()
