#include "gtest/gtest.h"
#include "ThirdOrderCurvatureEvaluator.hpp"

TEST(ThirdOrderCurvatureTest, CrossCoeficients2D)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double c_3 = 288.0;
    double c_2 = -1008.0;
    double c_1 = 1168.0;
    double c_0 = -448.0;
    Eigen::Matrix<double, 2, 4> control_points;
    control_points << 9, 4, 8, 6,
                      1, 5, 5, 5;
    Eigen::Vector4d coeficient_array = c_eval.get_2D_cross_coefficients(control_points);
    double true_c_3 = coeficient_array(3);
    double true_c_2 = coeficient_array(2);
    double true_c_1 = coeficient_array(1);
    double true_c_0 = coeficient_array(0);
    EXPECT_EQ(true_c_3, c_3);
    EXPECT_EQ(true_c_2, c_2);
    EXPECT_EQ(true_c_1, c_1);
    EXPECT_EQ(true_c_0, c_0);
}

TEST(ThirdOrderCurvatureTest, CrossCoeficients3D)
{
    ThirdOrderCurvatureEvaluator<3> c_eval{};
    double c_3 = 830.5;
    double c_2 = -354.75;
    double c_1 = 789.25;
    double c_0 = -154.0;
    Eigen::Matrix<double, 3, 4> control_points;
    control_points << 2, 4, 4, 2,
                      1, 2, 8, 6,
                      2, 3, 5, 0;
    Eigen::Vector4d coeficient_array = c_eval.get_3D_cross_coefficients(control_points);
    double true_c_3 = coeficient_array(3);
    double true_c_2 = coeficient_array(2);
    double true_c_1 = coeficient_array(1);
    double true_c_0 = coeficient_array(0);
    EXPECT_EQ(true_c_3, c_3);
    EXPECT_EQ(true_c_2, c_2);
    EXPECT_EQ(true_c_1, c_1);
    EXPECT_EQ(true_c_0, c_0);
}

TEST(ThirdOrderCurvatureTest, MatrixRetrieval)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double, 4,4> true_M;
    true_M << -0.16666667,  0.5, -0.5, 0.16666667,
               0.5,          -1,    0, 0.66666667,
              -0.5,         0.5,  0.5, 0.16666667,
        0.16666667,           0,    0, 0;
    Eigen::Matrix<double, 4,4> M = c_eval.get_third_order_M_matrix();
    float tolerance = 0.00000001;
    EXPECT_TRUE(true_M.isApprox(M,tolerance));
}

TEST(ThirdOrderCurvatureTest, DerivativeVectorRetrieval)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double scale_factor = 1;
    Eigen::Vector4d true_dt_vector;
    true_dt_vector << 0.36757914, 0.70007536, 1, 0;
    double t = 0.35003768191677603;
    Eigen::Vector4d dt_vector = c_eval.get_third_order_T_derivative_vector(t,scale_factor);
    float tolerance = 0.00000001;
    EXPECT_TRUE(true_dt_vector.isApprox(dt_vector, tolerance));
}

TEST(ThirdOrderCurvatureTest, SecondDerivativeVectorRetrieval)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Vector4d true_ddt_vector;
    double t = 0.5199907026925293;
    true_ddt_vector << 3.11994422, 2, 0, 0;
    Eigen::Vector4d ddt_vector = c_eval.get_third_order_T_second_derivative_vector(t);
    float tolerance = 0.00000001;
    EXPECT_TRUE(true_ddt_vector.isApprox(ddt_vector, tolerance));
}

TEST(ThirdOrderCurvatureTest, DerivativeVectorRetrievalError)
{
    // this tests _that_ the expected exception is thrown
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double scale_factor = 1;
    EXPECT_THROW(
    try
    {
        double t = -0.4;
        c_eval.get_third_order_T_derivative_vector(t,scale_factor);
    }
    catch(std::invalid_argument const& e)
    {
        EXPECT_STREQ("t value should be between 0 and 1", e.what());
        throw;
    }, std::invalid_argument);
}

TEST(ThirdOrderCurvatureTest, SecondDerivativeVectorRetrievalError)
{
    // this tests _that_ the expected exception is thrown
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    EXPECT_THROW(
    try
    {
        double t = 1.3;
        c_eval.get_third_order_T_second_derivative_vector(t);
    }
    catch(std::invalid_argument const& e)
    {
        EXPECT_STREQ("t value should be between 0 and 1", e.what());
        throw;
    }, std::invalid_argument);
}

TEST(ThirdOrderCurvatureTest, VelocityMagnitude)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double scale_factor = 1;
    double true_velocity = 2.8230311463199973;
    Eigen::Matrix<double, 2,4> control_points;
    control_points << 4, 7, 4, 4,
                      3, 5, 8, 4;
    double time_t = 0.35768880134338066;
    double velocity = c_eval.calculate_velocity_magnitude(time_t, control_points,scale_factor);
    double tolerance = 0.000000001;
    EXPECT_NEAR(true_velocity, velocity,tolerance);
}

TEST(ThirdOrderCurvatureTest, AccelerationMagnitude)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double true_acceleration = 2.2568129685516602;
    Eigen::Matrix<double, 2,4> control_points;
    control_points << 8, 3, 6, 4,
                      6, 5, 7, 3;
    double time_t = 0.4696969696969697;
    double acceleration = c_eval.calculate_acceleration_magnitude(time_t, control_points);
    EXPECT_EQ(true_acceleration, acceleration);
}

TEST(ThirdOrderCurvatureTest, CrossTermMagnitude)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double scale_factor = 1;
    double true_cross_term = 6.632352941176468;
    Eigen::Matrix<double, 2,4> control_points;
    control_points << 5, 6, 8, 3,
                      6, 9, 8, 1;
    double time_t = 0.14705882352941177;
    double cross_term = c_eval.calculate_cross_term_magnitude(time_t, control_points,scale_factor);
    EXPECT_EQ(true_cross_term, true_cross_term);
}

TEST(ThirdOrderCurvatureTest, MinVelocity)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double scale_factor = 1;
    double true_min_velocity = 0.12523779030065027;
    Eigen::Matrix<double, 2,4> control_points;
    control_points << 0.1,  3.2,  2.01, 1.1,
                      0.01, 2.98, 2.1,  0.99; 
    std::array<double,2> min_velocity_and_time = c_eval.find_minimum_velocity_and_time(control_points,scale_factor);
    double min_velocity = min_velocity_and_time[0];
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_min_velocity, min_velocity,tolerance);
}

TEST(ThirdOrderCurvatureTest, MaxAcceleration)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double true_max_acceleration = 9.055385138137417;
    Eigen::Matrix<double, 2,4> control_points;
    control_points << 7, 4, 2, 0,
                      3, 1, 8, 7;
    double max_acceleration = c_eval.find_maximum_acceleration(control_points);
    EXPECT_EQ(true_max_acceleration, max_acceleration);
}

TEST(ThirdOrderCurvatureTest, MaxCrossTerm)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double scale_factor = 1;
    double true_max_cross_term = 18.0;
    Eigen::Matrix<double, 2, 4> control_points;
    control_points << 0, 1, 4, 4,
                      6, 1, 4, 3;
    double max_cross_term = c_eval.find_maximum_cross_term(control_points,scale_factor);
    EXPECT_EQ(true_max_cross_term, max_cross_term);
}

TEST(ThirdOrderCurvatureTest, RandomNormalSpline)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double,2,4> control_points;
    control_points << 0, 6, 2, 1,
                      4, 3, 5, 2;
    double true_curvature_bound = 18.447094346279968;
    double curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound,tolerance);
}

TEST(ThirdOrderCurvatureTest, MaxCurvatureAtEndpoint)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double,2,4> control_points;
    control_points << 6, 1, 7, 1,
                      0, 1, 8, 9;
    double true_curvature_bound = 0.75;
    double curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, NoMovement)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double,2,4> control_points;
    control_points << 2,2,2,2,
                      3,3,3,3;
    double true_curvature_bound = 0;
    double curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, SameTwoInitialCPs)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double,2,4> control_points;
    control_points << 2, 2, 6, 8,
                      3, 3, 0, 5;
    double true_curvature_bound = 1.3193938001976513;
    double curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, SameThreeInitialCPs)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double,2,4> control_points;
    control_points << 2, 2, 2, 8,
                      3, 3, 3, 5;
    double true_curvature_bound = 0;
    double curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, SameTwoEndAndInitialCPs)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double,2,4> control_points;
    control_points << 2, 2, 8, 8,
                      3, 3, 5, 5;
    double true_curvature_bound = 0;
    double curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, InfiniteCurvature)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double,2,4> control_points;
    control_points << 1, 3, 2, 0,
                      1, 3, 2, 0;
    double true_curvature_bound = std::numeric_limits<double>::max();
    double curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, VeryLargeCurvature)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double,2,4> control_points;
    control_points << 0.1, 3.2, 2.01, 1.1,
                      0.01, 2.98, 2.1, 0.99;
    double true_curvature_bound = 367.5125371522953;
    double curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, StraightLine)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<double,2,4> control_points;
    control_points << 1,4,6,9,
                      1,4,6,9;
    double true_curvature_bound = 0;
    double curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, SplineCurvatureBound)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double true_max_curvature = 1.6000000063248174;
    int num_control_points = 8;
    double control_points[] = {-0.96055144, -0.04452866, 1.13866608, 2.49667343, 
        3.80204877, 4.76675134, 5.11662433,  4.76675134, -0.13352168, 0.06676084,
        -0.13352168, -0.58215741, -0.93804997, -0.84450962, -0.10254957,  1.25470791};
    double max_curvature = c_eval.find_spline_curvature_bound(control_points, num_control_points);
    double tolerance = 0.0000001;
    EXPECT_NEAR(true_max_curvature, max_curvature, tolerance);
}

TEST(ThirdOrderCurvatureTest, SplineCurvatureBound2)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double true_max_curvature = 2.028118191541948;
    int num_control_points = 8;
    double control_points[] = {-0.89402549, -0.05285741,  1.10545513, 2.47300498,
        3.79358126,  4.76115495, 5.11942253,  4.76115495, -0.10547684,  0.05273842,
        -0.10547684, -0.47275804, -0.79306865, -0.76080139, -0.11946946,  1.23867924};
    double max_curvature = c_eval.find_spline_curvature_bound(control_points, num_control_points);
    double tolerance = 0.0000001;
    EXPECT_NEAR(true_max_curvature, max_curvature, tolerance);
}

TEST(ThirdOrderCurvatureTest, MinVelocityOfSpline)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    double scale_factor = 1;
    double true_min_velocity = 0.7068769432442363;
    int num_control_points = 8;
    double control_points[] = {   -0.89402549, -0.05285741,  1.10545513,  2.47300498,
         3.79358126,  4.76115495,  5.11942253,  4.76115495, -0.10547684,  0.05273842,
        -0.10547684, -0.47275804, -0.79306865, -0.76080139, -0.11946946,  1.23867924};
    double min_velocity = c_eval.find_min_velocity_of_spline(control_points, num_control_points,scale_factor);
    double tolerance = 0.0000001;
    EXPECT_NEAR(true_min_velocity, min_velocity, tolerance);
}
