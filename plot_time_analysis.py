import matplotlib.pyplot as plt

discrete_times = [0.03034057498, 0.03108223224, 0.02852445841, 0.02919289088]
maximize_times = [0.00306290054321289, 0.00259621548652648, 0.003908126116, 0.00438338232]
roots_times = [0.00100869607925415, 0.001381050348, 0.001485831261, 0.001637473583]
min_vel_times = [0.0001661572456, 0.0001793963909]
max_over_min_times = [0.000142635107, 0.0002872812748]
control_point_times = [1.63E-04, 0.0002940416336, 0.0004574177265, 0.0007455370426]
geometric_times = [0.0000478022098541259]

discrete_std = [0.0009367860834, 0.0008166451794, 0.0005866512847, 0.001100102859]
maximize_std = [0.0006508972231, 0.001337830823, 0.001523611294, 0.0020393431]
roots_std = [0.0001813139089, 0.0003825943991, 0.0003895071886, 0.0003591138549]
min_vel_std = [0.0000362807957618437, 0.0000263813995834594]
max_over_min_std = [0.0000341948849480549, 0.000030222755724404]
control_point_std = [0.0000317371148401742, 0.0003108215061, 0.001070696344, 0.001643314309]
geometric_std = [0.000011865477316534]

orders_5 = [2,3,4,5]
orders_3 = [2,3]
orders_2 = [2]
["b", "g", "r", "c", "m", "y", "orange",]
"discrete_evaluations", "maximize_curvature_equation", "roots_of_curvature_derivative", \
#     "curvature_at_min_velocity", "max_numerator_over_min_denominator", "control_point_derivatives", "geometric"
plt.figure()
plt.errorbar(orders_5,discrete_times,discrete_std, fmt='-o',capsize=5,color="b",label="discrete_evaluations")
plt.errorbar(orders_5,maximize_times,maximize_std, fmt='-o',capsize=5,color="g",label="maximize_curvature_equation")
plt.errorbar(orders_5,roots_times,roots_std, fmt='-o',capsize=5,color="r",label="roots_of_curvature_derivative")
plt.errorbar(orders_5,control_point_times,control_point_std, fmt='-o',capsize=5,color="y",label="control_point_derivatives")
plt.errorbar(orders_3,min_vel_times,min_vel_std, fmt='-o',capsize=5,color="c",label="curvature_at_min_velocity")
plt.errorbar(orders_3,max_over_min_times,max_over_min_std, fmt='-o',capsize=5,color="m",label="max_numerator_over_min_denominator")
plt.errorbar(orders_2,geometric_times,geometric_std, fmt='-o',capsize=5,color="orange",label="geometric")
plt.xlabel("Polynomial Order")
plt.ylabel("Computationa Time")
plt.title("Time Analysis of Curvature Extrema Evaluation Methods")
plt.yscale("log") 
plt.legend()
plt.show()