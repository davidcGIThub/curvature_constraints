#include "FifthOrderCurvatureEvaluator.hpp"
#include "gtest/gtest.h"

TEST(FifthOrderCurvature, CurvatureBoundRandomControlPoints)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 2, 4, 5.7, 8, 8.2, 9,
                      6, 4, 8  , 6, 9  , 12;
    double true_curvature_bound = 1.2400376904192356;
    double curvature_bound = FifthOrderCurvatureEvaluator::evaluate_interval_curvature_bound(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound,curvature_bound,tolerance);
}

TEST(FifthOrderCurvature, CurvatureBoundStraightLine)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 1, 2, 3, 6, 9, 10,
                      1, 2, 3, 6, 9, 10;
    double true_curvature_bound = 0;
    double curvature_bound = FifthOrderCurvatureEvaluator::evaluate_interval_curvature_bound(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound,curvature_bound,tolerance);
}

TEST(FifthOrderCurvature, CurvatureBoundCusp)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 1, 2, 3, 4, 3, 1,
                      1, 2, 3, 4, 3, 1;
    double true_curvature_bound = std::numeric_limits<double>::max();
    double curvature_bound = FifthOrderCurvatureEvaluator::evaluate_interval_curvature_bound(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_curvature_bound,curvature_bound,tolerance);
}

TEST(FifthOrderCurvature, CurvatureBoundAlmostCuspTrue)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 1.01,2,3,4.05,2,0.98,
                      1,2.02,2.99,4,2.03,1;
    double true_curvature_bound = 1943826.5789539549;
    double curvature_bound = FifthOrderCurvatureEvaluator::evaluate_interval_curvature_bound(control_points);
    double tolerance = 0.000000001;
    EXPECT_NEAR(true_curvature_bound,curvature_bound,tolerance);
}

TEST(FifthOrderCurvature, CurvatureBoundZeroVelocityStart)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 1,1,1,1,1,8,
                      3,3,3,3,3,4;
    double true_curvature_bound = 0;
    double curvature_bound = FifthOrderCurvatureEvaluator::evaluate_interval_curvature_bound(control_points);
    double tolerance = 0.000000001;
    EXPECT_NEAR(true_curvature_bound,curvature_bound,tolerance);
}

TEST(FifthOrderCurvature, CurvatureBoundZeroVelocity)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 1,1,1,1,1,1,
                      3,3,3,3,3,3;
    double true_curvature_bound = 0;
    double curvature_bound = FifthOrderCurvatureEvaluator::evaluate_interval_curvature_bound(control_points);
    double tolerance = 0.000000001;
    EXPECT_NEAR(true_curvature_bound,curvature_bound,tolerance);
}

TEST(FifthOrderCurvature, IsCuspTrue)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 1,2,3,4,2,1,
                      1,2,3,4,2,1;
    bool is_cusp{true};
    EXPECT_EQ(FifthOrderCurvatureEvaluator::check_for_cusp(control_points),is_cusp);
}

TEST(FifthOrderCurvature, IsCuspFalse)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << -1,0,0,0,0,0,
                       1,2,3,4,5,6;
    bool is_cusp{false};
    EXPECT_EQ(FifthOrderCurvatureEvaluator::check_for_cusp(control_points),is_cusp);
}

TEST(FifthOrderCurvature, CrossTermBound)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 2, 4, 5.7, 8, 8.2, 9,
                      6, 4, 8  , 6, 9  , 12; 
    double true_cross_term_bound = 4.187500000000001;
    double cross_term_bound = FifthOrderCurvatureEvaluator::get_cross_term_max_bound(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_cross_term_bound, cross_term_bound, tolerance);
}

TEST(FifthOrderCurvature, AccelerationBound)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 2, 4, 5.7, 8, 8.2, 9,
                      6, 4, 8  , 6, 9  , 12;
    double true_acceleration_bound = 2.6238224872205897;
    double acceleration_bound = FifthOrderCurvatureEvaluator::get_acceleration_max_bound(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_acceleration_bound, acceleration_bound, tolerance);
}

TEST(FifthOrderCurvature, VelocityBound)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 2, 4, 5.7, 8, 8.2, 9,
                      6, 4, 8  , 6, 9  , 12;
    double true_min_velocity = 1.4546207623632545;
    double velocity_bound = FifthOrderCurvatureEvaluator::get_velocity_min_bound(control_points);
    double tolerance = 0.00000000001;
    EXPECT_NEAR(true_min_velocity, velocity_bound, tolerance);
}

TEST(FifthOrderCurvature, FifthOrderCrossTermCPs)
{
    Eigen::Matrix<double, 2, 6> control_points;
    control_points << 2, 4, 5.7, 8, 8.2, 9,
                      6, 4, 8  , 6, 9  , 12;
    Eigen::Matrix<double, 2, 7> true_cross_term_control_points;
    true_cross_term_control_points << -4.17083333, -4.1875, -2.67833333, -0.47083333, 1.645, 3.19166667, 3.81666667,
                                        0        ,       0,           0,           0,     0,          0,          0;
    Eigen::Matrix<double, 2, 7> cross_term_control_points = 
        FifthOrderCurvatureEvaluator::get_fifth_order_cross_term_control_points(control_points);
    double tolerance = 0.00001; 
    EXPECT_TRUE(cross_term_control_points.isApprox(true_cross_term_control_points,tolerance));
}
    


