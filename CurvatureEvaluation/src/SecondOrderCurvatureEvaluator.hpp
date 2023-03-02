#ifndef SECONDORDERCURVATUREEVALUATOR_HPP
#define SECONDORDERCURVATUREEVALUATOR_HPP
#include <array>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "BsplineToBezier.hpp"

namespace SecondOrderCurvatureEvaluator
{
    double evaluate_2D_spline_curvature(std::array<std::array<double,3>, 2> control_points);
    double evaluate_3D_spline_curvature(std::array<std::array<double,3>, 3> control_points);
    double calculate_max_curvature(double &A, double &dot_start, double &dot_end, \
                                   double &norm_start, double &norm_end, double &norm_mid);
}

#endif