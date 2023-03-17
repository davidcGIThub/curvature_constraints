#include "ConvexHullCollisionChecker.hpp"

template<int D>
ConvexHullCollisionChecker<D>::ConvexHullCollisionChecker(){}

template<int D>
double ConvexHullCollisionChecker<D>::getDistanceToObstacle(Eigen::Matrix<double,D,1> obstacle_center, 
                                                double radius, Eigen::MatrixXd points, int num_points)
{
    double distance_to_point = getDistanceToPoint(obstacle_center, points, num_points);
    double distance_to_obstacle = distance_to_point - radius;
    return distance_to_obstacle;
}

template<int D>
double getDistanceToPoint(Eigen::Matrix<double,D,1> obstacle_center, 
                          Eigen::MatrixXd points, int num_points)
{
    Eigen::MatrixXd translated_points = points - obstacle_center;
    double distance_to_point = mdm.min_norm(translated_points, num_points);
    return distance_to_point;
}