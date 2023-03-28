import numpy as np
import matplotlib.pyplot as plt
from path_generation.safe_flight_corridor import SFC_2D, plot_2D_sfc, plot_2D_sfcs, get2DRotationAndTranslationFromPoints
import time

point_1 = np.array([[3],[4]])
point_2 = np.array([[7],[10]])
point_3 = np.array([[14],[7]])
point_4 = np.array([[20],[14]])
R1, T1, min_len_1 = get2DRotationAndTranslationFromPoints(point_1, point_2)
R2, T2, min_len_2 = get2DRotationAndTranslationFromPoints(point_2, point_3)
R3, T3, min_len_3 = get2DRotationAndTranslationFromPoints(point_3, point_4)

sfc_1 = SFC_2D(np.array([[min_len_1+3],[2]]), T1, R1)
sfc_2 = SFC_2D(np.array([[min_len_2 + 2],[3]]), T2, R2)
sfc_3 = SFC_2D(np.array([[min_len_3],[2]]), T3, R3)
# sfcs = (sfc_1, sfc_2)
sfcs = (sfc_1, sfc_2, sfc_3)
# sfcs = (sfc_1, sfc_2)





initial_control_points = None

curvature_method = "roots_numerator_and_denominator"
# curvature_method = "constrain_max_acceleration_and_min_velocity"




# print("sfcs: " , sfcs)
plt.figure()
plot_2D_sfcs(sfcs)
# plot_2D_sfc(sfc_1)
# plt.scatter(minvo_cps[0,:],minvo_cps[1,:])
plt.show()
