#include "ObstacleConstraints.hpp"
#include <iostream>

template <int D>
ObstacleConstraints<D>::ObstacleConstraints()
{

}

template <int D>
float* ObstacleConstraints<D>::getObstacleDistancesToSpline(float cont_pts[], int num_control_points,
                                float obstacle_radii[], float obstacle_centers[], unsigned int num_obstacles)
{
    const int order = 3;
    float distance_to_obstacle{0};
    float* distances_array = new float[num_obstacles];
    Eigen::MatrixXf control_points(D,num_control_points);
    control_points = helper.array_to_eigen(cont_pts, num_control_points);
    for (unsigned int i = 0; i < num_obstacles; i++)
    {
        float obstacle_radius = obstacle_radii[i];
        Eigen::Matrix<float,D,1> obstacle_center = getObstacleCenterFromArray(
                obstacle_centers, i, num_obstacles);
        float distance_to_obstacle = collision_checker.getConservativeDistanceToSphere(obstacle_center, obstacle_radius,
                control_points, num_control_points); 
        
        if (distance_to_obstacle < 0)
        {
            distance_to_obstacle = getDistanceToClosestInterval(control_points, num_control_points,
                                obstacle_radius, obstacle_center);
        }
        distances_array[i] = -distance_to_obstacle;
    }
    return distances_array;
}

template <int D>
float ObstacleConstraints<D>::getObstacleDistanceToSpline(float cont_pts[], int num_control_points,
                                float obstacle_radius, float obstacle_center[])
{
    const int order = 3;
    Eigen::MatrixXf control_points(D,num_control_points);
    control_points = helper.array_to_eigen(cont_pts, num_control_points);
    unsigned int obstacle_index = 0;
    unsigned int num_obstacles = 1;
    Eigen::Matrix<float,D,1> obstacle_center_ = getObstacleCenterFromArray(obstacle_center, obstacle_index,
                                                num_obstacles);
    float distance_to_obstacle = collision_checker.getConservativeDistanceToSphere(obstacle_center_, obstacle_radius,
                control_points, num_control_points); 
    if (distance_to_obstacle < 0)
    {
        distance_to_obstacle = getDistanceToClosestInterval(control_points, num_control_points,
                                obstacle_radius, obstacle_center_);
    }
    return -distance_to_obstacle;
}

template <int D>
float ObstacleConstraints<D>::getDistanceToClosestInterval(Eigen::MatrixXf &control_points, int &num_control_points,
                                float &obstacle_radius, Eigen::Matrix<float,D,1> &obstacle_center)
{ 
    const int order = 3;
    const int NUM_POINTS = order+1;
    int num_points = NUM_POINTS;
    float min_distance{std::numeric_limits<float>::max()};
    float distance = 0;

    for (unsigned int i = 0; i < num_control_points-order-1; i++)
    {
        Eigen::Matrix<float, D,NUM_POINTS> bspline_points = control_points.block(0,i,D,order+1);
        Eigen::MatrixXf minvo_points(D,NUM_POINTS);
        minvo_points = cp_converter.convert_3rd_order_spline(bspline_points);
        distance = collision_checker.getDistanceToSphere(obstacle_center, obstacle_radius,
                minvo_points, num_points);
        if (distance < min_distance)
        {
            min_distance = distance;
        }
    }
    return min_distance;
}

template <int D>
Eigen::Matrix<float,D,1> ObstacleConstraints<D>::getObstacleCenterFromArray(float obstacle_centers[], unsigned int &obstacle_num,
                                                            unsigned int &num_obstacles)
{
    Eigen::Matrix<float,D,1> obstacle_center;
    if (D == 2)
    {
        obstacle_center << obstacle_centers[obstacle_num], obstacle_centers[obstacle_num+num_obstacles];
    }
    else //D == 3
    {
        obstacle_center << obstacle_centers[obstacle_num], obstacle_centers[obstacle_num+num_obstacles], 
            obstacle_centers[obstacle_num+num_obstacles*2];
    }
    return obstacle_center;
}

// explicit instantiation
template class ObstacleConstraints<2>;
template class ObstacleConstraints<3>;
