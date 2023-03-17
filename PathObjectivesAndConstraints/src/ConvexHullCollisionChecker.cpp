#include "ConvexHullCollisionChecker.hpp"
#include <iostream>

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
double ConvexHullCollisionChecker<D>::getDistanceToPoint(Eigen::Matrix<double,D,1> obstacle_center, 
                          Eigen::MatrixXd points, int num_points)
{
    Eigen::MatrixXd translated_points = points.colwise() - obstacle_center;
    int max_iter = 5000;
    unsigned int init_index = 1;
    double tolerance = 0.000001;
    double distance_to_point = mdm.min_norm(translated_points, num_points, max_iter, 
                                            init_index, tolerance);
    return distance_to_point;
}

// explicit instantiation
template class ConvexHullCollisionChecker<2>;
template class ConvexHullCollisionChecker<3>;