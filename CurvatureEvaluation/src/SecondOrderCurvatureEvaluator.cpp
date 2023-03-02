#include "SecondOrderCurvatureEvaluator.hpp"
#include <iostream>

namespace SecondOrderCurvatureEvaluator
{
    double evaluate_2D_spline_curvature(std::array<std::array<double,3>,2> points)
    {
        Eigen::Matrix<double, 2, 3> bspline_control_points;
        bspline_control_points << points[0][0], points[0][1], points[0][2],
                                  points[1][0], points[1][1], points[1][2];
        Eigen::Matrix<double, 2, 3> bezier_control_points;
        bezier_control_points << BsplineToBezier::convert_2D_second_order_spline(bspline_control_points);
        Eigen::Vector2d p_0 = bezier_control_points.block(0,0,2,1);
        Eigen::Vector2d p_1 = bezier_control_points.block(0,1,2,1);
        Eigen::Vector2d p_2 = bezier_control_points.block(0,2,2,1);
        Eigen::Vector2d mid = (p_0+p_2)/2.0;
        Eigen::Vector2d leg_start = p_1 - p_0;
        Eigen::Vector2d leg_middle = p_1 - mid;
        Eigen::Vector2d leg_end = p_1 - p_2;
        double cross_product_norm = abs(leg_start(0)*leg_end(1) - leg_start(1)*leg_end(0));
        double A = cross_product_norm/2.0;
        double norm_start = leg_start.norm();
        double norm_end = leg_end.norm();
        double norm_mid = leg_middle.norm();
        double dot_start = leg_start.dot(leg_middle);
        double dot_end = leg_middle.dot(leg_end);
        double max_curvature = calculate_max_curvature(A, dot_start, dot_end, norm_start, norm_end, norm_mid);
        return max_curvature;
    }

    double evaluate_3D_spline_curvature(std::array<std::array<double,3>,3> points)
    {
        Eigen::Matrix<double, 3, 3> bspline_control_points;
        bspline_control_points << points[0][0], points[0][1], points[0][2],
                                  points[1][0], points[1][1], points[1][2],
                                  points[2][0], points[2][1], points[2][2];
        Eigen::Matrix<double, 3, 3> bezier_control_points;
        bezier_control_points << BsplineToBezier::convert_3D_second_order_spline(bspline_control_points);
        Eigen::Vector3d p_0 = bezier_control_points.block(0,0,3,1);
        Eigen::Vector3d p_1 = bezier_control_points.block(0,1,3,1);
        Eigen::Vector3d p_2 = bezier_control_points.block(0,2,3,1);
        Eigen::Vector3d mid = (p_0+p_2)/2.0;
        Eigen::Vector3d leg_start = p_1 - p_0;
        Eigen::Vector3d leg_middle = p_1 - mid; 
        Eigen::Vector3d leg_end = p_1 - p_2;
        double cross_product_norm = (leg_start.cross(leg_end)).norm();
        double A = cross_product_norm/2.0;
        double norm_start = leg_start.norm();
        double norm_end = leg_end.norm();
        double norm_mid = leg_middle.norm();
        double dot_start = leg_start.dot(leg_middle);
        double dot_end = leg_middle.dot(leg_end);
        double max_curvature = calculate_max_curvature(A, dot_start, dot_end, norm_start, norm_end, norm_mid);
        return max_curvature;
    }

    double calculate_max_curvature(double &A, double &dot_start, double &dot_end, \
                                    double &norm_start, double &norm_end, double &norm_mid)
    {
        double max_curvature;
        if (dot_start <= 0)
        {
            if (norm_start == 0)
                max_curvature = 0;
            else
                max_curvature = A/pow(norm_start,3);
        }
        else if (dot_end <= 0)
        {
            if (norm_end == 0)
                max_curvature = 0;
            else
                max_curvature = A/pow(norm_end,3);
        }
        else
        {
            if (A == 0)
                max_curvature = std::numeric_limits<double>::max();
            else
                max_curvature = pow(norm_mid,3)/(A*A);
        }
        return max_curvature;
    }
}