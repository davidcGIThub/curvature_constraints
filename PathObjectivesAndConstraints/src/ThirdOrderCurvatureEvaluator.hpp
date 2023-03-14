#ifndef THIRDORDERCURVATUREEVALUATOR_HPP
#define THIRDORDERCURVATUREEVALUATOR_HPP
#include <array>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "gtest/gtest_prod.h"

template <int D> // D is the dimension of the spline
class ThirdOrderCurvatureEvaluator
{
    public:
        ThirdOrderCurvatureEvaluator();
        double find_spline_curvature_bound(double cont_pts[], int num_control_points);
        double evaluate_interval_curvature(Eigen::Matrix<double,D,4> &control_points);
    private:
        double find_maximum_cross_term(Eigen::Matrix<double,D,4> &control_points, double &scale_factor);
        Eigen::Vector4d get_2D_cross_coefficients(Eigen::Matrix<double,D,4> &control_points);
        Eigen::Vector4d get_3D_cross_coefficients(Eigen::Matrix<double,D,4> &control_points);
        double calculate_cross_term_magnitude(double &t, Eigen::Matrix<double,D,4> &control_points, double &scale_factor);
        Eigen::Matrix<double,D,4> array_section_to_eigen(double cont_pts[], int step, unsigned int index);
    FRIEND_TEST(ThirdOrderCurvatureTest, CrossCoeficients2D);
    FRIEND_TEST(ThirdOrderCurvatureTest, CrossCoeficients3D);
    FRIEND_TEST(ThirdOrderCurvatureTest, CrossTermMagnitude);
    FRIEND_TEST(ThirdOrderCurvatureTest, MaxCrossTerm);
};
#endif

