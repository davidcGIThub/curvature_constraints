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


### Parameters
num_test_cases = 10000
start_time = 0
scale_factor = 1
clamped = False
num_data_points_per_interval = 1000
order = 5
num_control_points = 11
dimension = 2
# bez
area_list_bez = []
max_norm_list_bez = []
min_norm_list_bez = []
bez_area_win = 0
bez_min_norm_win = 0
bez_max_norm_win = 0
# minv
area_list_minv = []
max_norm_list_minv = []
min_norm_list_minv = []
minv_area_win = 0
minv_min_norm_win = 0
minv_max_norm_win = 0
# minv
area_list_true = []
max_norm_list_true = []
min_norm_list_true = []

### Create B-Spline Object ###
for i in range(num_test_cases):
    num_control_points = order + 1
    control_points = np.random.randint(100, size=(dimension,num_control_points)) - 50
    bspline = BsplineEvaluation(control_points, order, start_time, scale_factor, clamped)
    ###### control point area######
    # bezier
    bezier_control_points = convert_to_bezier_control_points(control_points)
    bez_hull = ConvexHull(bezier_control_points.T)
    bez_area = bez_hull.volume 
    area_list_bez.append(bez_area)
    # minvo
    minvo_control_points = convert_to_minvo_control_points(control_points)
    minvo_hull = ConvexHull(minvo_control_points.T)
    minvo_area = minvo_hull.volume
    area_list_minv.append(minvo_area)
    # truth
    spline_data, time_data = bspline.get_spline_data(1000)
    truth_hull = ConvexHull(spline_data.T)
    truth_area = truth_hull.volume
    area_list_true.append(truth_area)
    if minvo_area < bez_area:
        minv_area_win += 1
    elif minvo_area > bez_area:
        bez_area_win += 1

    ###### norm ########
    # bezier
    norm_bez = np.linalg.norm(bezier_control_points,2,0)
    max_norm_list_bez.append(np.max(norm_bez))
    mdm_bez = MDM(bezier_control_points,dimension,1000,1)
    bezier_min_bound = mdm_bez.get_min_distance()
    min_norm_list_bez.append(bezier_min_bound)
    # minvo
    norm_minv = np.linalg.norm(minvo_control_points,2,0)
    max_norm_list_minv.append(np.max(norm_minv))
    mdm_minv = MDM(minvo_control_points,dimension,1000,1)
    minvo_min_bound = mdm_minv.get_min_distance()
    min_norm_list_minv.append(minvo_min_bound)
    # truth
    norm_true = np.linalg.norm(spline_data,2,0)
    max_norm_list_true.append(np.max(norm_true))
    min_norm_list_true.append(np.min(norm_true))
    if np.max(norm_minv) < np.max(norm_bez):
        minv_max_norm_win += 1
    elif np.max(norm_minv) > np.max(norm_bez):
        bez_max_norm_win += 1
    if minvo_min_bound > bezier_min_bound:
        minv_min_norm_win += 1
    elif minvo_min_bound < bezier_min_bound:
        bez_min_norm_win += 1
    
print("Area wins: Minvo = " , minv_area_win, ",  Bezier = " , bez_area_win)
print("Min Norm wins: Minvo = " , minv_min_norm_win, ",  Bezier = " , bez_min_norm_win)
print("Max Norm wins: Minvo = " , minv_max_norm_win, ",  Bezier = " , bez_max_norm_win)

truth_x_axis_data = np.array([0.5,1.5])
bez_x_axis_data = np.array([1.5,2.5])
minv_x_axis_data = np.array([2.5,3.5])

# Mean values for truth
area_mean_true = np.mean(area_list_true) + np.zeros(2)
max_norm_mean_true = np.mean(max_norm_list_true) + np.zeros(2)
min_norm_mean_true = np.mean(min_norm_list_true) + np.zeros(2)

# Mean values for bezier
area_mean_bez = np.mean(area_list_bez) + np.zeros(2)
max_norm_mean_bez = np.mean(max_norm_list_bez) + np.zeros(2)
min_norm_mean_bez = np.mean(min_norm_list_bez) + np.zeros(2)

# Mean values for minvo
area_mean_minv = np.mean(area_list_minv) + np.zeros(2)
max_norm_mean_minv = np.mean(max_norm_list_true) + np.zeros(2)
min_norm_mean_minv = np.mean(min_norm_list_true) + np.zeros(2)

#median values
def set_plot_box_colors(plot_box):
    colors = ['tab:green', 'tab:red','tab:blue']
    for patch, color in zip(plot_box['boxes'], colors):
        patch.set_facecolor(color)
    for median in plot_box['medians']:
        median.set(color ='black', linewidth = 1)

fig, ax = plt.subplots(1,3)
fig.suptitle(str(order) + " Order Control Point Comparison")
ax[0].set_title('Area of Convex Hull')
ax[0].plot(truth_x_axis_data, area_mean_true, color = "tab:green", label="truth")
ax[0].plot(bez_x_axis_data, area_mean_bez, color = "tab:red", label="bez")
ax[0].plot(minv_x_axis_data, area_mean_minv, color = "tab:blue", label="minvo")
area_data = [area_list_true, area_list_bez, area_list_minv]
pb_area = ax[0].boxplot(area_data,patch_artist = True, notch=True,showfliers=False)
set_plot_box_colors(pb_area)
ax[0].legend()

ax[1].set_title('Max Norm of Control Points')
ax[1].plot(truth_x_axis_data, max_norm_mean_true, color = "tab:green", label="truth")
ax[1].plot(bez_x_axis_data, max_norm_mean_bez, color = "tab:red", label="bez")
ax[1].plot(minv_x_axis_data, max_norm_mean_minv, color = "tab:blue", label="minvo")
max_norm_data = [max_norm_list_true, max_norm_list_bez, max_norm_list_minv]
pb_max_norm = ax[1].boxplot(max_norm_data,patch_artist = True, notch=True,showfliers=False)
set_plot_box_colors(pb_max_norm)

ax[2].set_title('Min Norm of Control Points')
ax[2].plot(truth_x_axis_data, min_norm_mean_true, color = "tab:green", label="truth")
ax[2].plot(bez_x_axis_data, min_norm_mean_bez, color = "tab:red", label="bez")
ax[2].plot(minv_x_axis_data, min_norm_mean_minv, color = "tab:blue", label="minvo")
min_norm_data = [min_norm_list_true, min_norm_list_bez, min_norm_list_minv]
pb_min_norm = ax[2].boxplot(min_norm_data,patch_artist = True, notch=True,showfliers=False)
set_plot_box_colors(pb_min_norm)

plt.show()


fig2, ax = plt.subplots(1,3)
fig2.suptitle(str(order) + " Order Control Point Comparison")
ax[0].set_title('Test Cases with Smaller Area')
ax[0].bar([0,1],[bez_area_win,minv_area_win], color =['tab:red', 'tab:blue'])

ax[1].set_title('Test Cases with Smaller Max Norm')
ax[1].bar([0,1],[bez_max_norm_win,minv_max_norm_win], color =['tab:red', 'tab:blue'])

ax[2].set_title('Test Cases with Larger Min Norm')
ax[2].bar([0,1],[bez_min_norm_win, minv_min_norm_win], color =['tab:red', 'tab:blue'])
plt.show()






