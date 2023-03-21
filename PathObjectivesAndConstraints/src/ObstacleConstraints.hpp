#ifndef OBSTACLECONSTRAINTS_HPP
#define OBSTACLECONSTRAINTS_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "ConvexHullCollisionChecker.hpp"
#include "CBindingHelper.hpp"
#include "BsplineToMinvo.hpp"

template <int D>
class ObstacleConstraints
{
    public:
        ObstacleConstraints();
        float* getObstacleDistancesToSpline(float cont_pts[], int num_control_points,
                                          float obstacle_radii[], float obstacle_centers[],
                                          unsigned int num_obstacles);
        float getObstacleDistanceToSpline(float cont_pts[], int num_control_points,
                                float obstacle_radius, float obstacle_center[]);
    private:
        ConvexHullCollisionChecker<D> collision_checker{};
        BsplineToMinvo<D> cp_converter{};
        CBindingHelper<D> helper{};
        float getDistanceToClosestInterval(Eigen::MatrixXf &control_points, int &num_control_points,
                                float &obstacle_radius, Eigen::Matrix<float,D,1> &obstacle_center);
        Eigen::Matrix<float,D,1> getObstacleCenterFromArray(float obstacle_centers[], unsigned int &obstacle_num,
                                                            unsigned int &num_obstacles);
};

#endif