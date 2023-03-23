#ifndef THIRDORDERCURVATUREBOUNDS_HPP
#define THIRDORDERCURVATUREBOUNDS_HPP
#include <array>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "gtest/gtest_prod.h"
#include "CBindingHelper.hpp"

template <int D> // D is the dimension of the spline
class ThirdOrderCurvatureBounds
{
    public:
        ThirdOrderCurvatureBounds();
        float get_spline_curvature_bound(float cont_pts[], int &num_control_points);
        Eigen::VectorXf get_interval_curvature_bounds(float cont_pts[], int &num_control_points);
        float evaluate_interval_curvature_bound(Eigen::Matrix<float,D,4> &control_points);
    private:
        CBindingHelper<D> cbind_help{};
        float find_maximum_cross_term(Eigen::Matrix<float,D,4> &control_points, float &scale_factor);
        Eigen::Vector4d get_2D_cross_coefficients(Eigen::Matrix<float,D,4> &control_points);
        Eigen::Vector4d get_3D_cross_coefficients(Eigen::Matrix<float,D,4> &control_points);
        float calculate_cross_term_magnitude(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor);
    FRIEND_TEST(ThirdOrderCurvatureTest, CrossCoeficients2D);
    FRIEND_TEST(ThirdOrderCurvatureTest, CrossCoeficients3D);
    FRIEND_TEST(ThirdOrderCurvatureTest, CrossTermMagnitude);
    FRIEND_TEST(ThirdOrderCurvatureTest, MaxCrossTerm);
};

#endif

