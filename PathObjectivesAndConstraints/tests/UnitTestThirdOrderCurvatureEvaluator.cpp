#include "gtest/gtest.h"
#include "ThirdOrderCurvatureEvaluator.hpp"

TEST(ThirdOrderCurvatureTest, CrossCoeficients2D)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    float c_3 = 288.0;
    float c_2 = -1008.0;
    float c_1 = 1168.0;
    float c_0 = -448.0;
    Eigen::Matrix<float, 2, 4> control_points;
    control_points << 9, 4, 8, 6,
                      1, 5, 5, 5;
    Eigen::Vector4d coeficient_array = c_eval.get_2D_cross_coefficients(control_points);
    float true_c_3 = coeficient_array(3);
    float true_c_2 = coeficient_array(2);
    float true_c_1 = coeficient_array(1);
    float true_c_0 = coeficient_array(0);
    EXPECT_EQ(true_c_3, c_3);
    EXPECT_EQ(true_c_2, c_2);
    EXPECT_EQ(true_c_1, c_1);
    EXPECT_EQ(true_c_0, c_0);
}

TEST(ThirdOrderCurvatureTest, CrossCoeficients3D)
{
    ThirdOrderCurvatureEvaluator<3> c_eval{};
    float c_3 = 830.5;
    float c_2 = -354.75;
    float c_1 = 789.25;
    float c_0 = -154.0;
    Eigen::Matrix<float, 3, 4> control_points;
    control_points << 2, 4, 4, 2,
                      1, 2, 8, 6,
                      2, 3, 5, 0;
    Eigen::Vector4d coeficient_array = c_eval.get_3D_cross_coefficients(control_points);
    float true_c_3 = coeficient_array(3);
    float true_c_2 = coeficient_array(2);
    float true_c_1 = coeficient_array(1);
    float true_c_0 = coeficient_array(0);
    EXPECT_EQ(true_c_3, c_3);
    EXPECT_EQ(true_c_2, c_2);
    EXPECT_EQ(true_c_1, c_1);
    EXPECT_EQ(true_c_0, c_0);
}

TEST(ThirdOrderCurvatureTest, CrossTermMagnitude)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    float scale_factor = 1;
    float true_cross_term = 6.632352941176468;
    Eigen::Matrix<float, 2,4> control_points;
    control_points << 5, 6, 8, 3,
                      6, 9, 8, 1;
    float time_t = 0.14705882352941177;
    float cross_term = c_eval.calculate_cross_term_magnitude(time_t, control_points,scale_factor);
    EXPECT_EQ(true_cross_term, true_cross_term);
}

TEST(ThirdOrderCurvatureTest, MaxCrossTerm)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    float scale_factor = 1;
    float true_max_cross_term = 18.0;
    Eigen::Matrix<float, 2, 4> control_points;
    control_points << 0, 1, 4, 4,
                      6, 1, 4, 3;
    float max_cross_term = c_eval.find_maximum_cross_term(control_points,scale_factor);
    EXPECT_EQ(true_max_cross_term, max_cross_term);
}

TEST(ThirdOrderCurvatureTest, RandomNormalSpline)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<float,2,4> control_points;
    control_points << 0, 6, 2, 1,
                      4, 3, 5, 2;
    float true_curvature_bound = 18.447094346279968;
    float curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    float tolerance = 0.00001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound,tolerance);
}

TEST(ThirdOrderCurvatureTest, MaxCurvatureAtEndpoint)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<float,2,4> control_points;
    control_points << 6, 1, 7, 1,
                      0, 1, 8, 9;
    float true_curvature_bound = 0.75;
    float curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    float tolerance = 0.0000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, NoMovement)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<float,2,4> control_points;
    control_points << 2,2,2,2,
                      3,3,3,3;
    float true_curvature_bound = 0;
    float curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    float tolerance = 0.0000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, SameTwoInitialCPs)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<float,2,4> control_points;
    control_points << 2, 2, 6, 8,
                      3, 3, 0, 5;
    float true_curvature_bound = 1.3193938001976513;
    float curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    float tolerance = 0.000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, SameThreeInitialCPs)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<float,2,4> control_points;
    control_points << 2, 2, 2, 8,
                      3, 3, 3, 5;
    float true_curvature_bound = 0;
    float curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    float tolerance = 0.00001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, SameTwoEndAndInitialCPs)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<float,2,4> control_points;
    control_points << 2, 2, 8, 8,
                      3, 3, 5, 5;
    float true_curvature_bound = 0;
    float curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    float tolerance = 0.00001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, InfiniteCurvature)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<float,2,4> control_points;
    control_points << 1, 3, 2, 0,
                      1, 3, 2, 0;
    float true_curvature_bound = std::numeric_limits<float>::max();
    float curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    float tolerance = 0.000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, VeryLargeCurvature)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<float,2,4> control_points;
    control_points << 0.1, 3.2, 2.01, 1.1,
                      0.01, 2.98, 2.1, 0.99;
    float true_curvature_bound = 367.5125371522953;
    float curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    float tolerance = 0.001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, StraightLine)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    Eigen::Matrix<float,2,4> control_points;
    control_points << 1,4,6,9,
                      1,4,6,9;
    float true_curvature_bound = 0;
    float curvature_bound = c_eval.evaluate_interval_curvature(control_points);
    float tolerance = 0.00001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound, tolerance);
}

TEST(ThirdOrderCurvatureTest, SplineCurvatureBound)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    float true_max_curvature = 1.6000000063248174;
    int num_control_points = 8;
    float control_points[] = {-0.96055144, -0.04452866, 1.13866608, 2.49667343, 
        3.80204877, 4.76675134, 5.11662433,  4.76675134, -0.13352168, 0.06676084,
        -0.13352168, -0.58215741, -0.93804997, -0.84450962, -0.10254957,  1.25470791};
    float max_curvature = c_eval.find_spline_curvature_bound(control_points, num_control_points);
    float tolerance = 0.000001;
    EXPECT_NEAR(true_max_curvature, max_curvature, tolerance);
}

TEST(ThirdOrderCurvatureTest, SplineCurvatureBound2)
{
    ThirdOrderCurvatureEvaluator<2> c_eval{};
    float true_max_curvature = 2.028118191541948;
    int num_control_points = 8;
    float control_points[] = {-0.89402549, -0.05285741,  1.10545513, 2.47300498,
        3.79358126,  4.76115495, 5.11942253,  4.76115495, -0.10547684,  0.05273842,
        -0.10547684, -0.47275804, -0.79306865, -0.76080139, -0.11946946,  1.23867924};
    float max_curvature = c_eval.find_spline_curvature_bound(control_points, num_control_points);
    float tolerance = 0.000001;
    EXPECT_NEAR(true_max_curvature, max_curvature, tolerance);
}

TEST(ThirdOrderCurvatureTest, SplineCurvatureBound3)
{
    ThirdOrderCurvatureEvaluator<3> c_eval{};
    float true_max_curvature = 2.4603294220963416;
    int num_control_points = 4;
    float control_points[] = {4, 1, 4, 5,
                                2, 2, 0, 4,
                                7, 0, 1, 7};
    float max_curvature = c_eval.find_spline_curvature_bound(control_points, num_control_points);
    float tolerance = 0.000001;
    EXPECT_NEAR(true_max_curvature, max_curvature, tolerance);
}
