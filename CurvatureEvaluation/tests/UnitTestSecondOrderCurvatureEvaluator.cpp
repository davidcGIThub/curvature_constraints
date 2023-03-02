#include "SecondOrderCurvatureEvaluator.hpp"
#include "gtest/gtest.h"
#include <iostream>

TEST(SecondOrderCurvatureTests, MethodOne)
{   
    std::array<std::array<double,3>,2> control_points = {{{5, 3, 0},{9, 8, 5}}};
    double true_curvature_bound = 0.2683281572999747;
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_2D_spline_curvature(control_points);
    EXPECT_EQ(true_curvature_bound, curvature_bound);
}

TEST(SecondOrderCurvatureTests, MethodTwo)
{
    std::array<std::array<double,3>,2> control_points = {{{5, 4, 5},{0, 3, 5}}};
    double true_curvature_bound = 0.44721359549995787;
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_2D_spline_curvature(control_points);
    EXPECT_EQ(true_curvature_bound, curvature_bound);
}

TEST(SecondOrderCurvatureTests, MethodThree)
{
    std::array<std::array<double,3>,2> control_points = {{{2, 6, 4},{3, 6, 9}}};
    double true_curvature_bound = 0.6666666666666666;
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_2D_spline_curvature(control_points);
    EXPECT_EQ(true_curvature_bound, curvature_bound);
}

TEST(SecondOrderCurvatureTests, NoMovement)
{
    std::array<std::array<double,3>,2> control_points = {{{3,3,3},{5,5,5}}};
    double true_curvature_bound = 0;
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_2D_spline_curvature(control_points);
    EXPECT_EQ(true_curvature_bound, curvature_bound);
}

TEST(SecondOrderCurvatureTests, SameInitialCPs)
{
    std::array<std::array<double,3>,2> control_points = {{{1, 1, 7},{1, 1, 7}}};
    double true_curvature_bound = 0;
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_2D_spline_curvature(control_points);
    EXPECT_EQ(true_curvature_bound, curvature_bound);
}

TEST(SecondOrderCurvatureTests, SameEndCps)
{
    std::array<std::array<double,3>,2> control_points = {{{1, 7, 7},{1, 7, 7}}};
    double true_curvature_bound = 0;
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_2D_spline_curvature(control_points);
    EXPECT_EQ(true_curvature_bound, curvature_bound);
}

TEST(SecondOrderCurvatureTests, InfiniteCurvature)
{
    std::array<std::array<double,3>,2> control_points = {{{5, 6, 4},{3, 6, 0}}};
    double true_curvature_bound = std::numeric_limits<double>::max();
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_2D_spline_curvature(control_points);
    EXPECT_EQ(true_curvature_bound, curvature_bound);
}

TEST(SecondOrderCurvatureTests, MethodOne3D)
{   
    std::array<std::array<double,3>,3> control_points = {{{6, 1, 4},{8, 9, 8},{2, 2, 5}}};
    double true_curvature_bound = 2.838959066508998;
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_3D_spline_curvature(control_points);
    double tolerance = 0.000000001;
    EXPECT_NEAR(true_curvature_bound, curvature_bound,tolerance);
}

TEST(SecondOrderCurvatureTests, MethodTwo3D)
{
    std::array<std::array<double,3>,3> control_points = {{{3, 7, 7},{0, 7, 8},{0, 5, 4}}};
    double true_curvature_bound = 4.690415759823429;
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_3D_spline_curvature(control_points);
    EXPECT_EQ(true_curvature_bound, curvature_bound);
}

TEST(SecondOrderCurvatureTests, MethodThree3D)
{
    std::array<std::array<double,3>,3> control_points = {{{9, 6, 0},{5, 5, 4},{8, 8, 2}}};
    double true_curvature_bound = 0.6758625033664688;
    double curvature_bound = SecondOrderCurvatureEvaluator::evaluate_3D_spline_curvature(control_points);
    EXPECT_EQ(true_curvature_bound, curvature_bound);
}