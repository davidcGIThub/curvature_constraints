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
        float getDistanceToSphere(Eigen::Matrix<float,D,1> &sphere_center, 
                                    float &sphere_radius, Eigen::MatrixXf &spline_points, int &num_spline_points);
        float getConservativeDistanceToSphere(Eigen::Matrix<float,D,1> &sphere_center, 
                                    float &sphere_radius, Eigen::MatrixXf &points, int &num_points);
        bool checkIfCollides(Eigen::Matrix<float,D,1> &sphere_center, 
                                    float &sphere_radius, Eigen::MatrixXf &points, int &num_points);
    private:
        float getDistanceToPoint(Eigen::Matrix<float,D,1> &point, 
                                   Eigen::MatrixXf &spline_points, int &num_spline_points);
        MDMAlgorithmClass<D> mdm{};
};
#endif