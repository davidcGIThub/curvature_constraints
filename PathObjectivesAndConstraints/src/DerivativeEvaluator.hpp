#ifndef DERIVATIVEEVALUATOR_HPP
#define DERIVATIVEEVALUATOR_HPP
#include <array>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "CubicEquationSolver.hpp"
#include "gtest/gtest_prod.h"
#include "CBindingHelper.hpp"

template <int D> // D is the dimension of the spline
class DerivativeEvaluator
{
    public:
        DerivativeEvaluator();
        float find_min_velocity_of_spline(float cont_pts[], int num_control_points, float scale_factor);
        float find_max_acceleration_of_spline(float cont_pts[], int num_control_points, float scale_factor);
        std::array<float,2> find_min_velocity_and_time(Eigen::Matrix<float,D,4> &control_points, float &scale_factor);
        std::array<float,2> find_max_acceleration_and_time(Eigen::Matrix<float,D,4> &control_points, float &scale_factor);
        Eigen::Matrix<float,D,1> calculate_velocity_vector(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor);
        float calculate_velocity_magnitude(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor);
        Eigen::Matrix<float,D,1> calculate_acceleration_vector(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor);
        float calculate_acceleration_magnitude(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor);
    private:
        CBindingHelper<D> cbind_help{};
        Eigen::Matrix<float, 4,4> get_third_order_M_matrix();
        Eigen::Vector4f get_third_order_T_derivative_vector(float &t, float &scale_factor);
        Eigen::Vector4f get_third_order_T_second_derivative_vector(float &t, float &scale_factor);
    FRIEND_TEST(DerivativeTest, MatrixRetrieval);
    FRIEND_TEST(DerivativeTest, DerivativeVectorRetrieval);
    FRIEND_TEST(DerivativeTest, SecondDerivativeVectorRetrieval);
    FRIEND_TEST(DerivativeTest, DerivativeVectorRetrievalError);
    FRIEND_TEST(DerivativeTest, SecondDerivativeVectorRetrievalError);
};
#endif

