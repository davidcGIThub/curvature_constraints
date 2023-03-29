import numpy as np
import matplotlib.pyplot as plt
from path_generation.safe_flight_corridor import SFC_2D, plot_2D_sfc, plot_2D_sfcs, get2DRotationAndTranslationFromPoints
import time

def get_composite_sfc_rotation_matrix(dimension, intervals_per_corridor, sfcs, num_minvo_cont_pts):
    order = 3
    num_corridors = len(intervals_per_corridor)
    M_len = num_minvo_cont_pts*dimension
    M_rot = np.zeros((M_len, M_len))
    num_cont_pts_per_interval = order + 1
    interval_count = 0
    dim_step = num_minvo_cont_pts
    for corridor_index in range(num_corridors):
        rotation = sfcs[corridor_index].rotation
        num_intervals = intervals_per_corridor[corridor_index]
        for interval_index in range(num_intervals):
            for cont_pt_index in range(num_cont_pts_per_interval):
                index = interval_count*num_cont_pts_per_interval+cont_pt_index
                M_rot[index, index] = rotation[0,0]
                M_rot[index, index + dim_step] = rotation[0,1]
                M_rot[index + dim_step, index] = rotation[1,0]
                M_rot[index + dim_step, index + dim_step] = rotation[1,1]
                if dimension == 3:
                    M_rot[2*dim_step + index, index] = rotation[2,0]
                    M_rot[2*dim_step + index, index + dim_step] = rotation[2,1]
                    M_rot[2*dim_step + index, index + 2*dim_step] = rotation[2,2]
                    M_rot[dim_step + index, index + 2*dim_step] = rotation[1,2]
                    M_rot[index, index + 2*dim_step] = rotation[0,2]
            interval_count += 1
    return M_rot

def get_intervals_per_corridor(num_corridors):
    if num_corridors == 1:
        # 5 intervals
        return (5)
    elif num_corridors == 2:
        # 6 intervals
        return (3,3)
    elif num_corridors == 3:
        # 8 intervals
        return (3,2,3)
    elif num_corridors == 4:
        # 10 intervals
        return (3,2,2,3)
    elif num_corridors == 5:
        # 12 intervals
        return (3,2,2,2,3)

point_1 = np.array([[3],[4]])
point_2 = np.array([[7],[10]])
point_3 = np.array([[14],[7]])
point_4 = np.array([[20],[14]])
dimension = np.shape(point_1)[0]
R1, T1, min_len_1 = get2DRotationAndTranslationFromPoints(point_1, point_2)
R2, T2, min_len_2 = get2DRotationAndTranslationFromPoints(point_2, point_3)
R3, T3, min_len_3 = get2DRotationAndTranslationFromPoints(point_3, point_4)
sfc_1 = SFC_2D(np.array([[min_len_1+3],[2]]), T1, R1)
sfc_2 = SFC_2D(np.array([[min_len_2 + 2],[3]]), T2, R1)
sfc_3 = SFC_2D(np.array([[min_len_3],[2]]), T3, R1)
# sfc_1 = SFC_2D(np.array([[min_len_1+3],[2]]), T1, np.eye(2))
# sfc_2 = SFC_2D(np.array([[min_len_2 + 2],[3]]), T2, np.eye(2))
# sfc_3 = SFC_2D(np.array([[min_len_3],[2]]), T3, np.eye(2))
# sfcs = (sfc_1, sfc_2)
sfcs = (sfc_1, sfc_2, sfc_3)
# sfcs = (sfc_1, sfc_2)

order = 3
intervals_per_corridor = get_intervals_per_corridor(len(sfcs))
num_intervals = np.sum(intervals_per_corridor)
num_minvo_cps = num_intervals*(order+1)
minvo_cps = np.random.randint(10, size=(dimension,num_minvo_cps)) # random
M_rot = get_composite_sfc_rotation_matrix(dimension, intervals_per_corridor, sfcs, num_minvo_cps)


rotated_minvo_cps = np.reshape(M_rot @ minvo_cps.flatten()[:,None],(2,num_minvo_cps))


# print("sfcs: " , sfcs)
plt.figure()
# plot_2D_sfcs(sfcs)
# plot_2D_sfc(sfc_1)
plt.scatter(minvo_cps[0,:],minvo_cps[1,:])
plt.plot(minvo_cps[0,:],minvo_cps[1,:])
plt.scatter(rotated_minvo_cps[0,:],rotated_minvo_cps[1,:])
plt.plot(rotated_minvo_cps[0,:],rotated_minvo_cps[1,:])

plt.show()
