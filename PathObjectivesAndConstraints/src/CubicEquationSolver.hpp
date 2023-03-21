#ifndef CUBICEQUATIONSOLVER_HPP
#define CUBICEQUATIONSOLVER_HPP
#include <array>
#include <list>
#include <limits>

namespace CubicEquationSolver
{
    std::array<float,3> solve_equation(float &a_term, float &b_term, float &c_term, float &d_term);

    std::array<float,3> solve_cubic_equation(float &a_term, float &b_term, float &c_term, float &d_term);

    std::array<float,2> solve_quadratic_equation(float &b_term, float &c_term, float &d_term);

    float solve_linear_equation(float &c_term, float &d_term);

    static const float DOUBLEMAX = std::numeric_limits<float>::max();
}

#endif