#ifndef BSPLINETOBEZIER_HPP
#define BSPLINETOBEZIER_HPP
#include <array>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

namespace BsplineToBezier
{
    Eigen::Matrix<double, 2,2> convert_2D_first_order_spline(Eigen::Matrix<double, 2,2> bspline_control_points);
    Eigen::Matrix<double, 2,3> convert_2D_second_order_spline(Eigen::Matrix<double, 2,3> bspline_control_points);
    Eigen::Matrix<double, 2,4> convert_2D_third_order_spline(Eigen::Matrix<double, 2,4> bspline_control_points);
    Eigen::Matrix<double, 2,5> convert_2D_fourth_order_spline(Eigen::Matrix<double, 2,5> bspline_control_points);
    Eigen::Matrix<double, 2,6> convert_2D_fifth_order_spline(Eigen::Matrix<double, 2,6> bspline_control_points);
    Eigen::Matrix<double, 3,2> convert_3D_first_order_spline(Eigen::Matrix<double, 3,2> bspline_control_points);
    Eigen::Matrix<double, 3,3> convert_3D_second_order_spline(Eigen::Matrix<double, 3,3> bspline_control_points);
    Eigen::Matrix<double, 3,4> convert_3D_third_order_spline(Eigen::Matrix<double, 3,4> bspline_control_points);
    Eigen::Matrix<double, 3,5> convert_3D_fourth_order_spline(Eigen::Matrix<double, 3,5> bspline_control_points);
    Eigen::Matrix<double, 3,6> convert_3D_fifth_order_spline(Eigen::Matrix<double, 3,6> bspline_control_points);
    Eigen::Matrix<double, 2,2> get_first_order_bspline_to_bezier_conversion_matrix();
    Eigen::Matrix<double,3,3> get_second_order_bspline_to_bezier_conversion_matrix();
    Eigen::Matrix<double,4,4> get_third_order_bspline_to_bezier_conversion_matrix();
    Eigen::Matrix<double,5,5> get_fourth_order_bspline_to_bezier_conversion_matrix();
    Eigen::Matrix<double,6,6> get_fifth_order_bspline_to_bezier_conversion_matrix();
}

#endif