#ifndef OBSTACLECONSTRAINTS_HPP
#define OBSTACLECONSTRAINTS_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

template <int D>
class ObstacleConstraints
{
    public:
        ObstacleConstraints();
        double* getObstacleDistancesToSpline(double cont_pts[], int num_control_points,
                                          double obstacle_radii[], double obstacle_centers[],
                                          int num_obstacles);
        double getObstacleDistanceToSpline(double cont_pts[], int num_control_points,
                                double obstacle_radius, double obstacle_center)
    private:
        ConvexHullCollisionChecker<D> collision_checker{};
        double getDistanceToClosestInterval(Eigen::MatrixXd control_points, int num_control_points,
                                double  obstacle_radius, Eigen::Matrix<double,D,1> obstacle_center);
        Eigen::Matrix<double,D,1> getObstacleCenterFromArray(double obstacle_centers[], int obstacle_num,
                                                            int num_obstacles);
};

#endif