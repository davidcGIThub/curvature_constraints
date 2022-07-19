from time import time
import numpy as np
import random
import matplotlib.pyplot as plt
from bsplinegenerator.bsplines import BsplineEvaluation

control_points = np.random.randint(10, size=(random.randint(2, 3),11)) # random

if len(control_points) == 1:
    control_points = control_points.flatten()

### Parameters
order = 3
start_time = 0
scale_factor = 1
derivative_order = 1
clamped = False
num_data_points = 1000

### Create B-Spline Object ###
bspline = BsplineEvaluation(control_points, order, start_time, scale_factor, clamped)
vel_mag, time_data = bspline.get_derivative_magnitude_data(num_data_points,1)
vel, time_data = bspline.get_spline_derivative_data(num_data_points,1)
accel, time_data = bspline.get_spline_derivative_data(num_data_points,2)

time_step = time_data[1] - time_data[0]
start_time = time_data[0] + time_step/2
end_time = time_data[-1] - time_step/2
plt.plot(time_data,2*np.dot(np.transpose(vel),accel).diagonal())
plt.plot(time_data,vel_mag)
vel_mag_deriv = (vel_mag[1:]**2 - vel_mag[0:-1]**2) / (time_data[1:] - time_data[0:-1])
plt.plot(np.linspace(start_time,end_time,num_data_points-1), vel_mag_deriv)
plt.show()