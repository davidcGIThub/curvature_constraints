#ifndef THIRDORDERCURVATUREEVALUATOR_HPP
#define THIRDORDERCURVATUREEVALUATOR_HPP
#include <array>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "CubicEquationSolver.hpp"
#include "gtest/gtest_prod.h"

template <int D> // D is the dimension of the spline
class ThirdOrderCurvatureEvaluator
{
    public:
        ThirdOrderCurvatureEvaluator();
        double find_spline_curvature_bound(double cont_pts[], int num_control_points);
        double evaluate_interval_curvature(Eigen::Matrix<double,D,4> &control_points);
    private:
        double find_min_velocity_of_spline(double cont_pts[], int num_control_points, double scale_factor);
        std::array<double,2> find_minimum_velocity_and_time(Eigen::Matrix<double,D,4> &control_points, double &scale_factor);
        double find_maximum_acceleration(Eigen::Matrix<double,D,4> &control_points);
        double find_maximum_cross_term(Eigen::Matrix<double,D,4> &control_points, double &scale_factor);
        Eigen::Vector4d get_2D_cross_coefficients(Eigen::Matrix<double,D,4> &control_points);
        Eigen::Vector4d get_3D_cross_coefficients(Eigen::Matrix<double,D,4> &control_points);
        double calculate_velocity_magnitude(double &t, Eigen::Matrix<double,D,4> &control_points, double &scale_factor);
        double calculate_acceleration_magnitude(double &t, Eigen::Matrix<double,D,4> &control_points);
        double calculate_cross_term_magnitude(double &t, Eigen::Matrix<double,D,4> &control_points, double &scale_factor);
        Eigen::Matrix<double, 4,4> get_third_order_M_matrix();
        Eigen::Vector4d get_third_order_T_derivative_vector(double &t, double &scale_factor);
        Eigen::Vector4d get_third_order_T_second_derivative_vector(double &t);
        Eigen::Matrix<double,D,4> array_section_to_eigen(double cont_pts[], int step, unsigned int index);
    FRIEND_TEST(ThirdOrderCurvatureTest, CrossCoeficients2D);
    FRIEND_TEST(ThirdOrderCurvatureTest, CrossCoeficients3D);
    FRIEND_TEST(ThirdOrderCurvatureTest, MatrixRetrieval);
    FRIEND_TEST(ThirdOrderCurvatureTest, DerivativeVectorRetrieval);
    FRIEND_TEST(ThirdOrderCurvatureTest, SecondDerivativeVectorRetrieval);
    FRIEND_TEST(ThirdOrderCurvatureTest, DerivativeVectorRetrievalError);
    FRIEND_TEST(ThirdOrderCurvatureTest, SecondDerivativeVectorRetrievalError);
    FRIEND_TEST(ThirdOrderCurvatureTest, VelocityMagnitude);
    FRIEND_TEST(ThirdOrderCurvatureTest, AccelerationMagnitude);
    FRIEND_TEST(ThirdOrderCurvatureTest, CrossTermMagnitude);
    FRIEND_TEST(ThirdOrderCurvatureTest, MinVelocity);
    FRIEND_TEST(ThirdOrderCurvatureTest, MaxAcceleration);
    FRIEND_TEST(ThirdOrderCurvatureTest, MaxCrossTerm);
    FRIEND_TEST(ThirdOrderCurvatureTest, MinVelocityOfSpline);
};
#endif

