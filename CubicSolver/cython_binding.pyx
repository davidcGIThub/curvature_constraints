
cdef extern from "CubicEquationSolver.hpp":
    array<double,3> solve_equation(double a_term, double b_term, double c_term, double d_term)