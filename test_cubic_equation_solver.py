import numpy as np
from max_curvature_evaluators.helper_files import cube_root_solver
import time
# c = 2
# d = 3.2
# print("linear problem 1: " , cube_root_solver.__linear_solver(c,d))
# c = 0
# d = 3.2
# print("linear problem 2: " , cube_root_solver.__linear_solver(c,d))
# b = 3
# c = 6
# d = 3
# print("Quadratic problem 1: " , cube_root_solver.__quadratic_solver(b,c,d))
# b = 3
# c = 6
# d = 8
# print("Quadratic problem 2: " , cube_root_solver.__quadratic_solver(b,c,d))
# b = 3
# c = 6
# d = 2
# print("Quadratic problem 3: " , cube_root_solver.__quadratic_solver(b,c,d))
# a = -243
# b = 432
# c = -192
# print("Quadratic problem ?: " , cube_root_solver.__quadratic_solver(a,b,c))
a = 4
b = 3
c = -6
d = -2
start_time = time.time()
print("Cubic problem 1: " , cube_root_solver.__cubic_solver(a,b,c,d))
print("elapsed_time: ", time.time() - start_time)
# a = 4
# b = 3
# c = 6
# d = 2
# print("Cubic problem 2: " , cube_root_solver.__cubic_solver(a,b,c,d))
# a = 4
# b = 4
# c = 1
# d = 0
# print("Cubic problem 3: " , cube_root_solver.__cubic_solver(a,b,c,d))
# a = 3
# b = 6
# c = 4
# d = 0.8888888888888888
# print("Cubic problem 4: " , cube_root_solver.__cubic_solver(a,b,c,d))
# d =  18*a*b*c*d - 4*(b**3)*d + (b**2)*(c**2) - 4*a*(c**3) - 27*(a**2)*(d**2)
# p = b*b - 3*a*c
# print("d: " , d)
# print("p: " , p)


