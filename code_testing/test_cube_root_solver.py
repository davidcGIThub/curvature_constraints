import numpy as np
from max_curvature_evaluators.helper_evaluations.cube_root_solver import solver
import random

a = random.randint(-10,10)
b = random.randint(-10,10)
c = random.randint(-10,10)
d = random.randint(-10,10)
solver(a,b,c,d)