#include "gtest/gtest.h"
#include "BsplineToBezier.hpp"

TEST(BsplineToBezierTests, 2DFirstOrderConversion)
{
    Eigen::Matrix<double, 2,2> bspline_control_points;
    bspline_control_points << 0, 2,
                              2,-1;
    Eigen::Matrix<double, 2,2> true_bezier_control_points;
    true_bezier_control_points << 0, 2,
                             2,-1;
    Eigen::Matrix<double, 2,2> bezier_control_points;
    BsplineToBezier<2> converter{};
    bezier_control_points << converter.convert_1st_order_spline(bspline_control_points);
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points));
}

TEST(BsplineToBezierTests, 2DSecondOrderConversion)
{
    Eigen::Matrix<double, 2, 3> bspline_control_points;
    bspline_control_points << 0,   2,   3.5,
                              2,  -1,   2;
    Eigen::Matrix<double, 2,3> true_bezier_control_points;
    true_bezier_control_points << 1,    2,    2.75,
                                0.5,   -1,    0.5;
    Eigen::Matrix<double, 2,3> bezier_control_points;
    BsplineToBezier<2> converter{};
    bezier_control_points << converter.convert_2nd_order_spline(bspline_control_points);
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points));
}

TEST(BsplineToBezierTests, 2DThirdOrderConversion)
{
    Eigen::Matrix<double, 2, 4> bspline_control_points;
    bspline_control_points << 0,  2, 3.5, 3,
                              2, -1,   2, 5;
    Eigen::Matrix<double, 2, 4> true_bezier_control_points;
    true_bezier_control_points << 1.91666667, 2.5,        3,         3.16666667,
                                  0.        , 0  ,        1,         2;
    Eigen::Matrix<double, 2,4> bezier_control_points;
    BsplineToBezier<2> converter{};
    bezier_control_points << converter.convert_3rd_order_spline(bspline_control_points);
    float tolerance = 0.00000001;
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points,tolerance));
}

TEST(BsplineToBezierTests, 2DFourthOrderConversion)
{
    Eigen::Matrix<double, 2, 5> bspline_control_points;
    bspline_control_points << 0,  2, 3.5, 3, 5,
                              2, -1,   2, 5, 5.5;
    Eigen::Matrix<double, 2, 5> true_bezier_control_points;
    true_bezier_control_points <<  2.64583333, 2.95833333, 3.16666667, 3.20833333, 3.27083333,
                                   0.75,       1.25,       2,          2.75,       3.39583333;
    Eigen::Matrix<double, 2, 5> bezier_control_points;
    BsplineToBezier<2> converter{};
    bezier_control_points << converter.convert_4th_order_spline(bspline_control_points);
    float tolerance = 0.00000001;
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points,tolerance));
}

TEST(BsplineToBezierTests, 2DFifthOrderConversion)
{
    Eigen::Matrix<double, 2, 6> bspline_control_points;
    bspline_control_points << 0,  2, 3.5, 3,   5,  6.5,
                              2, -1,   2, 5, 5.5,  5;
    Eigen::Matrix<double, 2, 6> true_bezier_control_points;
    true_bezier_control_points << 3.05,       3.175,      3.25,       3.3,        3.4,        3.5625,
                                  2.02916667, 2.55833333, 3.11666667, 3.63333333, 4.06666667, 4.40833333;
    Eigen::Matrix<double, 2, 6> bezier_control_points;
    BsplineToBezier<2> converter{};
    bezier_control_points << converter.convert_5th_order_spline(bspline_control_points);
    float tolerance = 0.00000001;
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points,tolerance));
}

TEST(BsplineToBezierTests, 3DFirstOrderConversion)
{
    Eigen::Matrix<double, 3,2> bspline_control_points;
    bspline_control_points << 0,   2,
                              2,  -1,
                              1, 3.2;
    Eigen::Matrix<double, 3,2> true_bezier_control_points;
    true_bezier_control_points << 0,   2,
                                  2,  -1,
                                  1,   3.2;
    Eigen::Matrix<double, 3,2> bezier_control_points;
    BsplineToBezier<3> converter{};
    bezier_control_points << converter.convert_1st_order_spline(bspline_control_points);
    float tolerance = 0.00000001;
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points,tolerance));
}

TEST(BsplineToBezierTests, 3DSecondOrderConversion)
{
    Eigen::Matrix<double, 3, 3> bspline_control_points;
    bspline_control_points << 0,   2,   3.5,
                              2,  -1,   2,
                              1, 3.2,   5;
    Eigen::Matrix<double, 3,3> true_bezier_control_points;
    true_bezier_control_points << 1,    2,    2.75,
                                0.5,   -1,     0.5,
                                2.1,  3.2,     4.1;
    Eigen::Matrix<double, 3,3> bezier_control_points;
    BsplineToBezier<3> converter{};
    bezier_control_points << converter.convert_2nd_order_spline(bspline_control_points);
    float tolerance = 0.00000001;
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points,tolerance));
}

TEST(BsplineToBezierTests, 3DThirdOrderConversion)
{
    Eigen::Matrix<double, 3, 4> bspline_control_points;
    bspline_control_points << 0,   2, 3.5, 3,
                              2,  -1,   2, 5,
                              1, 3.2,   5, 0;
    Eigen::Matrix<double, 3, 4> true_bezier_control_points;
    true_bezier_control_points << 1.91666667, 2.5,        3,         3.16666667,
                                  0.        , 0  ,        1,         2,
                                  3.13333333, 3.8,      4.4,        3.86666667;
    Eigen::Matrix<double, 3,4> bezier_control_points;
    BsplineToBezier<3> converter{};
    bezier_control_points << converter.convert_3rd_order_spline(bspline_control_points);
    float tolerance = 0.00000001;
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points,tolerance));
}

TEST(BsplineToBezierTests, 3DFourthOrderConversion)
{
    Eigen::Matrix<double, 3, 5> bspline_control_points;
    bspline_control_points << 0,  2 , 3.5, 3, 5  ,
                              2, -1 ,   2, 5, 5.5,
                              1, 3.2,   5, 0, 3.3;
    Eigen::Matrix<double, 3, 5> true_bezier_control_points;
    true_bezier_control_points <<  2.64583333, 2.95833333, 3.16666667, 3.20833333, 3.27083333,
                                   0.75,       1.25,       2,          2.75,       3.39583333,
                                   3.8,        3.98333333, 3.86666667, 3.18333333, 2.5625;
    Eigen::Matrix<double, 3, 5> bezier_control_points;
    BsplineToBezier<3> converter{};
    bezier_control_points << converter.convert_4th_order_spline(bspline_control_points);
    float tolerance = 0.0000001;
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points,tolerance));
}

TEST(BsplineToBezierTests, 3DFifthOrderConversion)
{
    Eigen::Matrix<double, 3, 6> bspline_control_points;
    bspline_control_points << 0,  2 , 3.5, 3,   5,  6.5,
                              2, -1 ,   2, 5, 5.5,  5,
                              1, 3.2,   5, 0, 3.3,  1.5;
    Eigen::Matrix<double, 3, 6> true_bezier_control_points;
    true_bezier_control_points << 3.05,       3.175,      3.25,       3.3,        3.4,        3.5625,
                                  2.02916667, 2.55833333, 3.11666667, 3.63333333, 4.06666667, 4.40833333,
                                  3.47916667, 3.23166667, 2.82333333, 2.32666667, 1.99333333, 1.8375;
    Eigen::Matrix<double, 3, 6> bezier_control_points;
    BsplineToBezier<3> converter{};
    bezier_control_points << converter.convert_5th_order_spline(bspline_control_points);
    float tolerance = 0.00000001;
    EXPECT_TRUE(bezier_control_points.isApprox(true_bezier_control_points,tolerance));
}

