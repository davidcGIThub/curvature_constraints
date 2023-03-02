#ifndef FIFTHORDERCURVATUREEVALUATOR_HPP
#define FIFTHORDERCURVATUREEVALUATOR_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "MDMAlgorithmClass.hpp"
#include "BsplineToBezier.hpp"

namespace FifthOrderCurvatureEvaluator
{
    double evaluate_interval_curvature_bound(Eigen::Matrix<double, 2, 6> &control_points);
    double get_cross_term_max_bound(Eigen::Matrix<double, 2, 6> &control_points);
    double get_acceleration_max_bound(Eigen::Matrix<double, 2, 6> &control_points);
    double get_velocity_min_bound(Eigen::Matrix<double, 2, 6> &control_points);
    Eigen::Matrix<double, 2,7> get_fifth_order_cross_term_control_points(Eigen::Matrix<double,2,6> &control_points);
    bool check_for_cusp(Eigen::Matrix<double, 2, 6> &control_points);
}



#endif