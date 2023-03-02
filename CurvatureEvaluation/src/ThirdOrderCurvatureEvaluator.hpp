#ifndef THIRDORDERCURVATUREEVALUATOR_HPP
#define THIRDORDERCURVATUREEVALUATOR_HPP
#include <array>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "CubicEquationSolver.hpp"

namespace ThirdOrderCurvatureEvaluator
{
    extern "C"
    double find_spline_curvature_bound(double cont_pts[], int num_control_points);
    extern "C"
    double find_min_velocity_of_spline(double cont_pts[], int num_control_points, double scale_factor);
    double evaluate_interval_curvature(Eigen::Matrix<double,2,4> &control_points);
    std::array<double,2> find_minimum_velocity_and_time(Eigen::Matrix<double,2,4> &control_points, double &scale_factor);
    double find_maximum_acceleration(Eigen::Matrix<double,2,4> &control_points);
    double find_maximum_cross_term(Eigen::Matrix<double,2,4> &control_points, double &scale_factor);
    Eigen::Vector4d get_2D_cross_coefficients(Eigen::Matrix<double,2,4> &control_points);
    Eigen::Vector4d get_3D_cross_coefficients(Eigen::Matrix<double,3,4> &control_points);
    double calculate_velocity_magnitude(double &t, Eigen::Matrix<double,2,4> &control_points, double &scale_factor);
    double calculate_acceleration_magnitude(double &t, Eigen::Matrix<double,2,4> &control_points);
    double calculate_cross_term_magnitude(double &t, Eigen::Matrix<double,2,4> &control_points, double &scale_factor);
    Eigen::Matrix<double, 4,4> get_third_order_M_matrix();
    Eigen::Vector4d get_third_order_T_derivative_vector(double &t, double &scale_factor);
    Eigen::Vector4d get_third_order_T_second_derivative_vector(double &t);
}
#endif