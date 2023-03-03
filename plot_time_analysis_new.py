import matplotlib.pyplot as plt

discrete_times = 0.027242960691452028
maximize_times = 0.0033083910942077635
roots_times = 0.001954437971115112
roots_num_den = 0.000274564266204834
control_point_times = 0.000381046941171658

discrete_std = 0.0010185855394797968
maximize_std = 0.0015308729398982664
roots_std = 0.0009874361002188293
roots_num_den_std = 4.6722790318489266e-05
control_point_std = 0.000381046941171658



["b", "g", "r", "c", "m", "y", "orange",]
"discrete_evaluations", "maximize_curvature_equation", "roots_of_curvature_derivative", \
#     "curvature_at_min_velocity", "max_numerator_over_min_denominator", "control_point_derivatives", "geometric"
plt.figure()
plt.errorbar(1,discrete_times,discrete_std, fmt='-o',capsize=7,color="b",label="1000 discrete_evaluations")
plt.errorbar(2,maximize_times,maximize_std, fmt='-o',capsize=7,color="g",label="maximize_curvature_equation")
plt.errorbar(3,roots_times,roots_std, fmt='-o',capsize=7,color="r",label="roots_of_curvature_derivative")
plt.errorbar(4,roots_num_den,roots_num_den_std, fmt='-o',capsize=7,color="c",label="roots of numerator & denominator")
plt.errorbar(5,control_point_times,control_point_std, fmt='-o',capsize=7,color="m",label="control_point_derivatives")


plt.xlabel("Methods of Evaluating Curvature Extrema")
plt.ylabel("Computationa Time")
plt.title("Time Analysis of Curvature Extrema Evaluation Methods")
plt.yscale("log") 
plt.legend()
plt.show()

# method:  discrete_evaluations
# mean:  -2.1345412016627677e-05
# std:  0.0002103983172481965
# median  -7.046656567857838e-16
# high:  2.1220936054942903e-06
# low:  -0.005839707964350349
# average_time:  0.027242960691452028
# std_time:  0.0010185855394797968
 
# method:  maximize_curvature_equation
# mean:  -0.0009557577980586171
# std:  0.0302149316959877
# median  1.4050463566761877e-16
# high:  3.0470391258164547e-05
# low:  -0.9559579322870172
# average_time:  0.0033083910942077635
# std_time:  0.0015308729398982664
 
# method:  roots_of_curvature_derivative
# mean:  2.0325636754746163e-07
# std:  1.4222470474086938e-06
# median  1.1122185806827343e-11
# high:  3.0470391256640454e-05
# low:  -3.5527136788004883e-15
# average_time:  0.001954437971115112
# std_time:  0.0009874361002188293
 
# method:  roots_of_curvature_numerator_and_denominator
# mean:  0.3310935497004202
# std:  0.7072445362689805
# median  0.11768245399634926
# high:  9.802214439277256
# low:  -3.5527136788004883e-15
# average_time:  0.000274564266204834
# std_time:  4.6722790318489266e-05
 
# method:  control_point_derivatives
# mean:  6.655675895774785
# std:  55.940451705203074
# median  0.4799766425962193
# high:  1451.3182422970435
# low:  -0.9271293969022043
# average_time:  0.00020877599716186522
# std_time:  0.000381046941171658