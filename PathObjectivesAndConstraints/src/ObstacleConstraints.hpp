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
        bool* checkIfObstaclesCollide(float cont_pts[], int num_control_points,
                                          float obstacle_radii[], float obstacle_centers[],
                                          unsigned int num_obstacles);
        bool checkIfObstacleCollides(float cont_pts[], int num_control_points,
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

extern "C"
{
    ObstacleConstraints<2>* ObstacleConstraints_2(){return new ObstacleConstraints<2>();}
    bool checkIfObstacleCollides_2(ObstacleConstraints<2>* ex, float cont_pts[], int num_control_points,
                                float obstacle_radius, float obstacle_center[]){return ex->checkIfObstacleCollides(cont_pts,
                                num_control_points,obstacle_radius,obstacle_center);}
    float getObstacleDistanceToSpline_2(ObstacleConstraints<2>* ex, float cont_pts[], int num_control_points,
                                float obstacle_radius, float obstacle_center[]){return ex->getObstacleDistanceToSpline(cont_pts,
                                num_control_points,obstacle_radius,obstacle_center);}
    ObstacleConstraints<3>* ObstacleConstraints_3(){return new ObstacleConstraints<3>();}
    bool checkIfObstacleCollides_3(ObstacleConstraints<3>* ex, float cont_pts[], int num_control_points,
                                float obstacle_radius, float obstacle_center[]){return ex->checkIfObstacleCollides(cont_pts,
                                num_control_points,obstacle_radius,obstacle_center);}
    float getObstacleDistanceToSpline_3(ObstacleConstraints<3>* ex, float cont_pts[], int num_control_points,
                                float obstacle_radius, float obstacle_center[]){return ex->getObstacleDistanceToSpline(cont_pts,
                                num_control_points,obstacle_radius,obstacle_center);}
}

#endif