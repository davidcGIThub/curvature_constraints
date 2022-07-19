import numpy as np
from bsplinegenerator import bsplines
import matplotlib.pyplot as plt

def create_control_points_with_saturated_angles(num_control_points,angle,length,fixed_angle,fixed_length,ratio=1):
    control_points = np.zeros((2,num_control_points))
    lengths = np.zeros(num_control_points-1)
    angles = np.zeros(num_control_points-2)
    prev_vec = np.array([])
    for i in range(num_control_points):
        if i == 0:
            control_points[:,i][:,None] = np.array([[0],[0]])
        elif i == 1:
            random_vec = np.random.rand(2,1)
            unit_random_vec = random_vec/np.linalg.norm(random_vec)
            if fixed_length:
                new_length = length
            else:
                new_length = length*np.random.rand()
            next_vec = unit_random_vec*new_length
            control_points[:,i][:,None] = control_points[:,i-1][:,None] + next_vec
            prev_vec = next_vec
            lengths[i-1] = new_length
        else:
            if fixed_angle:
                new_angle = angle
            else:  
                new_angle = (angle - np.random.rand()*angle)*np.random.choice((-1,1))
            unit_prev_vec = prev_vec/np.linalg.norm(prev_vec)
            R = np.array([[np.cos(new_angle), -np.sin(new_angle)],[np.sin(new_angle), np.cos(new_angle)]])
            if fixed_length:
                new_length = length
                if i == num_control_points -1:
                    new_length = length*ratio
            else:
                new_length = length*np.random.rand()
            next_vec = new_length*np.dot(R,unit_prev_vec)
            control_points[:,i][:,None] = control_points[:,i-1][:,None] + next_vec
            prev_vec = next_vec
            lengths[i-1] = new_length
            angles[i-2] = new_angle
    return control_points, angles, lengths

def get_min_leg_times_sin_angle(angles,lengths):
    min_leg = np.inf
    for i in range(len(angles)):
        angle = angles[i]
        leg_1 = lengths[i]/np.sin(angle)
        leg_2 = lengths[i+1]/np.sin(angle)
        min_leg = np.min([min_leg,leg_1,leg_2])
    return min_leg

equalizer = 0.353553921985746 # second order

def get_smallest_ratio(num_cps,_angle,_length):
    min_radius = 0
    min_leg = np.inf
    ratio = 1.0
    while min_radius < min_leg*equalizer:
        control_points, angles, lengths = create_control_points_with_saturated_angles(num_cps,_angle,_length,True,True,ratio)
        min_leg = get_min_leg_times_sin_angle(angles,lengths)
        bspline = bsplines.BsplineEvaluation(control_points,num_cps-1,0,1,False)
        curvature_data, time_data = bspline.get_spline_curvature_data(1000)
        min_radius = 1/np.max(curvature_data)
        ratio = ratio-0.001
    return ratio

def get_leg_radius_data(num_cps):
    first_leg_length = 1
    second_leg_lengths = np.linspace(0.1,1,100)
    radii = np.zeros(100)
    angle = np.pi/2
    for i in range(100):
        ratio = second_leg_lengths[i]/first_leg_length
        control_points, angles, lengths = create_control_points_with_saturated_angles(num_cps,angle,first_leg_length,True,True,ratio)
        bspline = bsplines.BsplineEvaluation(control_points,num_cps-1,0,1,False)
        curvature_data, time_data = bspline.get_spline_curvature_data(1000)
        min_radius = 1/np.max(curvature_data)
        radii[i] = min_radius
        # if 100%(i+1) == 0:
        #     bspline.plot_spline(100)
    return second_leg_lengths, radii



        
num_control_points = 6
angle = np.pi/2
length = 5
fixed_angle = True
fixed_length = True
number_of_data_points = 1000

leg_lengths, radii = get_leg_radius_data(num_control_points)
plt.plot(leg_lengths,radii)
plt.show()

# ratio = get_smallest_ratio(num_control_points, angle, length)
ratio = 1

control_points, angles, lengths = create_control_points_with_saturated_angles(num_control_points,angle,length,fixed_angle,fixed_length,ratio)
min_leg = get_min_leg_times_sin_angle(angles,lengths)
bspline = bsplines.BsplineEvaluation(control_points,num_control_points-1,0,1,False)

curvature_data, time_data = bspline.get_spline_curvature_data(number_of_data_points)
min_radius = 1/np.max(curvature_data)
print("min_radius: " , min_radius)
print("min_leg times equalizer: " , min_leg*equalizer)
print("ratio: " , ratio)
bspline.plot_spline(number_of_data_points)
bspline.plot_curvature(number_of_data_points)

