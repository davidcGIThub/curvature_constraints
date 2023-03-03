from time import time
import numpy as np
import random
from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.bsplines import BsplineEvaluation
from max_curvature_evaluators.helper_files.cross_term_control_points import get_bezier_cross_term_norm_bound, get_minvo_cross_term_norm_bound
from bsplinegenerator.bspline_to_bezier import convert_to_bezier_control_points
from bsplinegenerator.bspline_to_minvo import convert_to_minvo_control_points, bezier_to_minvo_control_points
from max_curvature_evaluators.helper_files.mdm_algorithm_adapted import MDM
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt

### Control Points ###
# control_points = np.array([-3,-4,-2,-.5,1,0,2,3.5,3,6,8]) # 1 D
# control_points = np.array([[-3,-4,-2,-.5,1,0,2,3.5,3,5,6.5],
#                            [.5,3.5,6,5.5,3.7,2,-1,2,5,5.5,5]]) # 2 D
# control_points = np.array([[-3,  -4, -2, -.5, 1  ,   0,  2, 3.5, 3],
#                            [.5, 3.5,  6, 5.5, 3.7,   2, -1,   2, 5],
#                            [ 1, 3.2,  5,   0, 3.3, 1.5, -1, 2.5, 4]]) # 3D

# control_points  = create_random_control_points_greater_than_angles(num_control_points,order,1,dimension)

### Parameters
start_time = 0
scale_factor = 1
clamped = False
num_data_points_per_interval = 1000
order = 5
num_control_points = 11
dimension = 2
# bez
area_list_bez = []
min_vel_list_bez = []
max_vel_list_bez = []
min_accel_list_bez = []
max_accel_list_bez = []
max_cross_term_list_bez = []
# minv
area_list_minv = []
min_vel_list_minv = []
max_vel_list_minv = []
min_accel_list_minv = []
max_accel_list_minv = []
max_cross_term_list_minv = []
# minv
area_list_true = []
min_vel_list_true = []
max_vel_list_true = []
min_accel_list_true = []
max_accel_list_true = []
max_cross_term_list_true = []

num_test_cases = 10000

### Create B-Spline Object ###
for i in range(num_test_cases):
    num_control_points = order + 1
    control_points = np.random.randint(100, size=(dimension,num_control_points)) - 50
    bspline = BsplineEvaluation(control_points, order, start_time, scale_factor, clamped)
    ###### control point area######
    # bezier
    bezier_control_points = convert_to_bezier_control_points(control_points)
    # bezier_control_points = np.random.randint(100, size=(dimension,num_control_points)) - 50 ######
    bez_hull = ConvexHull(bezier_control_points.T)
    bez_area = bez_hull.volume 
    area_list_bez.append(bez_area)
    # minvo
    minvo_control_points = convert_to_minvo_control_points(control_points)
    # minvo_control_points = bezier_to_minvo_control_points(bezier_control_points)
    minvo_hull = ConvexHull(minvo_control_points.T)
    minvo_area = minvo_hull.volume
    area_list_minv.append(minvo_area)
    # truth
    spline_data, time_data = bspline.get_spline_data(1000)
    truth_hull = ConvexHull(spline_data.T)
    truth_area = truth_hull.volume
    area_list_true.append(truth_area)

    ###### velocity ########
    velocity_cps = bspline.get_control_point_derivatives(1)
    # bezier
    bezier_vel_cps = convert_to_bezier_control_points(velocity_cps)
    bez_max_vel_bound = np.max(np.linalg.norm(bezier_vel_cps,2,0))
    max_vel_list_bez.append(bez_max_vel_bound)
    mdm_bez_vel = MDM(bezier_vel_cps,dimension,1000,1)
    bezier_min_vel_bound = mdm_bez_vel.get_min_distance()
    min_vel_list_bez.append(bezier_min_vel_bound)
    # minvo
    minvo_vel_cps = convert_to_minvo_control_points(velocity_cps)
    minvo_max_vel_bound = np.max(np.linalg.norm(minvo_vel_cps,2,0))
    max_vel_list_minv.append(minvo_max_vel_bound)
    mdm_minv_vel = MDM(minvo_vel_cps,dimension,500,1)
    minvo_min_vel_bound = mdm_minv_vel.get_min_distance()
    min_vel_list_minv.append(minvo_min_vel_bound)
    # truth
    velocity_data, time_data = bspline.get_spline_derivative_data(10000,1)
    velocity_mag_data = np.linalg.norm(velocity_data, 2, 0)
    max_vel_list_true.append(np.max(velocity_mag_data))
    min_vel_list_true.append(np.min(velocity_mag_data))

    ###### acceleration ########
    acceleration_cps = bspline.get_control_point_derivatives(2)
    #bezier
    bezier_accel_cps = convert_to_bezier_control_points(acceleration_cps)
    bez_max_accel_bound = np.max(np.linalg.norm(bezier_accel_cps,2,0))
    max_accel_list_bez.append(bez_max_accel_bound)
    mdm_bez_accel = MDM(bezier_accel_cps,dimension,500,1)
    bezier_min_accel_bound = mdm_bez_accel.get_min_distance()
    min_accel_list_bez.append(bezier_min_accel_bound)
    #minvo
    minvo_accel_cps = convert_to_minvo_control_points(acceleration_cps)
    minvo_max_accel_bound = np.max(np.linalg.norm(minvo_accel_cps,2,0))
    max_accel_list_minv.append(minvo_max_accel_bound)
    mdm_minv_accel = MDM(minvo_accel_cps,dimension,500,1)
    minvo_min_accel_bound = mdm_minv_accel.get_min_distance()
    min_accel_list_minv.append(minvo_min_accel_bound)
    # truth
    accel_data, time_data = bspline.get_spline_derivative_data(10000,2)
    accel_mag_data = np.linalg.norm(accel_data, 2, 0)
    max_accel_list_true.append(np.max(accel_mag_data))
    min_accel_list_true.append(np.min(accel_mag_data))

    ###### cross term ########
    #bezier
    bez_cross_term_norm_bound = get_bezier_cross_term_norm_bound(control_points, order,dimension)
    max_cross_term_list_bez.append(bez_cross_term_norm_bound)
    #minvo
    minvo_cross_term_norm_bound = get_minvo_cross_term_norm_bound(control_points, order,dimension)
    max_cross_term_list_minv.append(minvo_cross_term_norm_bound)
    #truth
    cross_term_data = np.abs(np.cross(velocity_data.T, accel_data.T))
    max_cross_term_list_true.append(np.max(cross_term_data))


truth_x_axis_scatter_data = np.zeros(num_test_cases)
bez_x_axis_scatter_data = np.zeros(num_test_cases) + 5
minv_x_axis_scatter_data = np.zeros(num_test_cases) + 10

truth_x_axis_data = np.array([0.5,1.5])
bez_x_axis_data = np.array([1.5,2.5])
minv_x_axis_data = np.array([2.5,3.5])

# Mean values for truth
min_vel_mean_true = np.mean(min_vel_list_true) + np.zeros(2)
max_vel_mean_true = np.mean(max_vel_list_true) + np.zeros(2)
min_accel_mean_true = np.mean(min_accel_list_true) + np.zeros(2)
max_accel_mean_true = np.mean(max_accel_list_true) + np.zeros(2)
area_mean_true = np.mean(area_list_true) + np.zeros(2)
max_cross_term_mean_true = np.mean(max_cross_term_list_true) + np.zeros(2)
min_vel_mean = np.mean(min_vel_list_true)
# Mean values for bezier
min_vel_mean_bez = np.mean(min_vel_list_bez) + np.zeros(2)
max_vel_mean_bez = np.mean(max_vel_list_bez) + np.zeros(2)
min_accel_mean_bez = np.mean(min_accel_list_bez) + np.zeros(2)
max_accel_mean_bez = np.mean(max_accel_list_bez) + np.zeros(2)
area_mean_bez = np.mean(area_list_bez) + np.zeros(2)
max_cross_term_mean_bez = np.mean(max_cross_term_list_bez) + np.zeros(2)
min_vel_mean = np.mean(min_vel_list_bez)
# Mean values for minvo
min_vel_mean_minv = np.mean(min_vel_list_minv) + np.zeros(2)
max_vel_mean_minv = np.mean(max_vel_list_minv) + np.zeros(2)
min_accel_mean_minv = np.mean(min_accel_list_minv) + np.zeros(2)
max_accel_mean_minv = np.mean(max_accel_list_minv) + np.zeros(2)
area_mean_minv = np.mean(area_list_minv) + np.zeros(2)
max_cross_term_mean_minv = np.mean(max_cross_term_list_minv) + np.zeros(2)
min_vel_mean = np.mean(min_vel_list_minv)

#median values
def set_plot_box_colors(plot_box):
    colors = ['tab:green', 'tab:red','tab:blue']
    for patch, color in zip(plot_box['boxes'], colors):
        patch.set_facecolor(color)
    for median in plot_box['medians']:
        median.set(color ='black', linewidth = 1)

fig, ax = plt.subplots(2,3)
fig.suptitle(str(order) + " Order Control Point Comparison")
ax[0, 0].set_title('Area of Convex Hull')
ax[0, 0].plot(truth_x_axis_data, area_mean_true, color = "tab:green", label="truth")
ax[0, 0].plot(bez_x_axis_data, area_mean_bez, color = "tab:red", label="bez")
ax[0, 0].plot(minv_x_axis_data, area_mean_minv, color = "tab:blue", label="minvo")
area_data = [area_list_true, area_list_bez, area_list_minv]
pb_area = ax[0, 0].boxplot(area_data,patch_artist = True, notch=True,showfliers=False)
set_plot_box_colors(pb_area)
ax[0,0].legend()

ax[0,1].set_title('Velocity CPs Max Norm')
ax[0, 1].plot(truth_x_axis_data, max_vel_mean_true, color = "tab:green", label="truth")
ax[0, 1].plot(bez_x_axis_data, max_vel_mean_bez, color = "tab:red", label="bez")
ax[0, 1].plot(minv_x_axis_data, max_vel_mean_minv, color = "tab:blue", label="minvo")
max_velocity_data = [max_vel_list_true, max_vel_list_bez, max_vel_list_minv]
pb_max_vel = ax[0, 1].boxplot(max_velocity_data,patch_artist = True, notch=True,showfliers=False)
set_plot_box_colors(pb_max_vel)

ax[0,2].set_title('Acceleration CPs Max Norm')
ax[0, 2].plot(truth_x_axis_data, max_accel_mean_true, color = "tab:green", label="truth")
ax[0, 2].plot(bez_x_axis_data, max_accel_mean_bez, color = "tab:red", label="bez")
ax[0, 2].plot(minv_x_axis_data, max_accel_mean_minv, color = "tab:blue", label="minvo")
max_acceleration_data = [max_accel_list_true, max_accel_list_bez, max_accel_list_minv]
pb_max_accel = ax[0, 2].boxplot(max_acceleration_data,patch_artist = True, notch=True,showfliers=False)
set_plot_box_colors(pb_max_accel)

ax[1,0].set_title('Cross Term CPs Max Norm')
ax[1, 0].plot(truth_x_axis_data, max_cross_term_mean_true, color = "tab:green", label="truth")
ax[1, 0].plot(bez_x_axis_data, max_cross_term_mean_bez, color = "tab:red", label="bez")
ax[1, 0].plot(minv_x_axis_data, max_cross_term_mean_minv, color = "tab:blue", label="minvo")
max_cross_term_data = [max_cross_term_list_true, max_cross_term_list_bez, max_cross_term_list_minv]
pb_cross_term = ax[1, 0].boxplot(max_cross_term_data,patch_artist = True, notch=True,showfliers=False)
set_plot_box_colors(pb_cross_term)

ax[1,1].set_title('Velocity CPs Min Norm')
ax[1, 1].plot(truth_x_axis_data, min_vel_mean_true, color = "tab:green", label="truth")
ax[1, 1].plot(bez_x_axis_data, min_vel_mean_bez, color = "tab:red", label="bez")
ax[1, 1].plot(minv_x_axis_data, min_vel_mean_minv, color = "tab:blue", label="minvo")
min_velocity_data = [min_vel_list_true, min_vel_list_bez, min_vel_list_minv]
pb_min_vel = ax[1, 1].boxplot(min_velocity_data,patch_artist = True, notch=True,showfliers=False)
set_plot_box_colors(pb_min_vel)

ax[1,2].set_title('Acceleration CPs Min Norm')
ax[1, 2].plot(truth_x_axis_data, min_accel_mean_true, color = "tab:green", label="truth")
ax[1, 2].plot(bez_x_axis_data, min_accel_mean_bez, color = "tab:red", label="bez")
ax[1, 2].plot(minv_x_axis_data, min_accel_mean_minv, color = "tab:blue", label="minvo")
min_accel_data = [min_accel_list_true, min_accel_list_bez, min_accel_list_minv]
pb_min_accel = ax[1, 2].boxplot(min_accel_data,patch_artist = True, notch=True,showfliers=False)
set_plot_box_colors(pb_min_accel)
plt.show()



