#include "BsplineToBezier.hpp"
#include <iostream>


template<int D>
BsplineToBezier<D>::BsplineToBezier()
{

}

template<int D>
Eigen::Matrix<double, D,2> BsplineToBezier<D>::convert_1st_order_spline(Eigen::Matrix<double, D,2> bspline_control_points)
{
    Eigen::Matrix<double, 2,2> conversion_matrix = get_first_order_bspline_to_bezier_conversion_matrix();
    Eigen::Matrix<double, D,2> bezier_control_points = bspline_control_points*conversion_matrix;
    return bezier_control_points;
}

template<int D>
Eigen::Matrix<double, D,3> BsplineToBezier<D>::convert_2nd_order_spline(Eigen::Matrix<double, D,3> bspline_control_points)
{
    Eigen::Matrix<double, 3,3> conversion_matrix = get_second_order_bspline_to_bezier_conversion_matrix();
    Eigen::Matrix<double, D,3> bezier_control_points = bspline_control_points*conversion_matrix;
    return bezier_control_points;
}

template<int D>
Eigen::Matrix<double, D,4> BsplineToBezier<D>::convert_3rd_order_spline(Eigen::Matrix<double, D,4> bspline_control_points)
{
    Eigen::Matrix<double, 4,4> conversion_matrix = get_third_order_bspline_to_bezier_conversion_matrix();
    Eigen::Matrix<double, D,4> bezier_control_points = bspline_control_points*conversion_matrix;
    return bezier_control_points;
}

template<int D>
Eigen::Matrix<double, D,5> BsplineToBezier<D>::convert_4th_order_spline(Eigen::Matrix<double, D,5> bspline_control_points)
{
    Eigen::Matrix<double, 5,5> conversion_matrix = get_fourth_order_bspline_to_bezier_conversion_matrix();
    Eigen::Matrix<double, D,5> bezier_control_points = bspline_control_points*conversion_matrix;
    return bezier_control_points;
}

template<int D>
Eigen::Matrix<double, D,6> BsplineToBezier<D>::convert_5th_order_spline(Eigen::Matrix<double, D,6> bspline_control_points)
{
    Eigen::Matrix<double, 6,6> conversion_matrix = get_fifth_order_bspline_to_bezier_conversion_matrix();
    Eigen::Matrix<double, D,6> bezier_control_points = bspline_control_points*conversion_matrix;
    return bezier_control_points;
}

template<int D>
Eigen::Matrix<double, 2,2> BsplineToBezier<D>::get_first_order_bspline_to_bezier_conversion_matrix()
{
    Eigen::Matrix<double,2,2> conversion_matrix;
    conversion_matrix << 1, 0,
                            0, 1;
    return conversion_matrix;
}

template<int D>
Eigen::Matrix<double,3,3> BsplineToBezier<D>::get_second_order_bspline_to_bezier_conversion_matrix()
{
    Eigen::Matrix<double,3,3> conversion_matrix;
    conversion_matrix << 0.5, 0   , 0,
                            0.5, 1   , 0.5,
                            0  , 0   , 0.5;
    return conversion_matrix;
}

template<int D>
Eigen::Matrix<double,4,4> BsplineToBezier<D>::get_third_order_bspline_to_bezier_conversion_matrix()
{
    Eigen::Matrix<double,4,4> conversion_matrix;
    conversion_matrix << 1/6.0 , 0    , 0    , 0    ,
                            4/6.0 , 4/6.0, 2/6.0, 1/6.0,
                            1/6.0 , 2/6.0, 4/6.0, 4/6.0,
                            0     , 0    , 0    , 1/6.0;
    return conversion_matrix;
}

template<int D>
Eigen::Matrix<double,5,5> BsplineToBezier<D>::get_fourth_order_bspline_to_bezier_conversion_matrix()
{
    Eigen::Matrix<double,5,5> conversion_matrix;
    conversion_matrix << 1/24.0 , 0      , 0      , 0      , 0      ,
                            11/24.0, 8/24.0 , 4/24.0 , 2/24.0 , 1/24.0 ,
                            11/24.0, 14/24.0, 16/24.0, 14/24.0, 11/24.0,
                            1/24.0 , 2/24.0 , 4/24.0 , 8/24.0 , 11/24.0,
                            0      , 0      , 0      , 0      , 1/24.0 ;
    return conversion_matrix;
}

template<int D>
Eigen::Matrix<double,6,6> BsplineToBezier<D>::get_fifth_order_bspline_to_bezier_conversion_matrix()
{
    Eigen::Matrix<double,6,6> conversion_matrix;
    conversion_matrix << 1/120.0 , 0     , 0     , 0     , 0     , 0     ,
                            26/120.0, 16/120.0, 8/120.0 , 4/120.0 , 2/120.0 , 1/120.0 ,
                            66/120.0, 66/120.0, 60/120.0, 48/120.0, 36/120.0, 26/120.0,
                            26/120.0, 36/120.0, 48/120.0, 60/120.0, 66/120.0, 66/120.0,
                            1/120.0 , 2/120.0 , 4/120.0 , 8/120.0 , 16/120.0, 26/120.0,
                            0     , 0     , 0     , 0     , 0     , 1/120.0     ;
    return conversion_matrix;
}

//Explicit template instantiations
template class BsplineToBezier<2>;
template class BsplineToBezier<3>;