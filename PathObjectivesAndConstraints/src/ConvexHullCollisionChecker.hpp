#ifndef CONVEXHULLCOLLISIONCHECKER_HPP
#define CONVEXHULLCOLLISIONCHECKER_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "MDMAlgorithmClass.hpp"

template <int D>
class ConvexHullCollisionChecker
{
    public:
        ConvexHullCollisionChecker();
        double getDistanceToObstacle(Eigen::Matrix<double,D,1> obstacle_center, 
                                    double radius, Eigen::MatrixXd points, int num_points);
    private:

        double getDistanceToPoint(Eigen::Matrix<double,D,1> obstacle_center, 
                                   Eigen::MatrixXd points, int num_points);
        MDMAlgorithmClass<D> mdm{};
};

#endif