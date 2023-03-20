#include "ObstacleConstraints.hpp"
#include "ConvexHullCollisionChecker.hpp"

template <int D>
ObstacleConstraints<D>::ObstacleConstraints()
{

}

template <int D>
double* ObstacleConstraints<D>::getObstacleDistancesToSpline(double cont_pts[], int num_control_points,
                                double obstacle_radii[], double obstacle_centers[], int num_obstacles)
{
    const int order = 3;
    double distance_to_obstacle{0};
    double* distances_array = new double[num_obstacles];
    Eigen::MatrixXd control_points(D,num_control_points) = array_to_eigen(cont_pts, num_control_points);
    for (unsigned int i = 0; i < num_obstacles; i++)
    {
        double obstacle_radius = obstacle_radii[i];
        Eigen::Matrix<double,D,1> obstacle_center = getObstacleCenterFromArray(
                obstacle_centers, i, num_obstacles);
        double distance_to_obstacle = collision_checker.getConservativeDistanceToSphere(obstacle_center, obstacle_radius,
                control_points, num_control_points); 
        
        if (distance_to_obstacle < 0)
        {
            distance_to_obstacle = getDistanceToClosestInterval(control_points, num_control_points,
                                obstacle_radius, obstacle_center);
        }
        distances_array[i] = distance_to_obstacle;
    }
    return distances_array;
}

template <int D>
double ObstacleConstraints<D>::getObstacleDistanceToSpline(double cont_pts[], int num_control_points,
                                double obstacle_radius, double obstacle_center)
{
    const int order = 3;
    Eigen::MatrixXd control_points(D,num_control_points) = array_to_eigen(cont_pts, num_control_points);
    double distance_to_obstacle = collision_checker.getConservativeDistanceToSphere(obstacle_center, obstacle_radius,
                control_points, num_control_points); 
    if (distance_to_obstacle < 0)
    {
        distance_to_obstacle = getDistanceToClosestInterval(control_points, num_control_points,
                                obstacle_radius, obstacle_center);
    }
    return distance_to_obstacle
}

template <int D>
double ObstacleConstraints<D>::getDistanceToClosestInterval(Eigen::MatrixXd control_points, int num_control_points,
                                double obstacle_radius, Eigen::Matrix<double,D,1> obstacle_center)
{
    const int order = 3;
    int num_points = order+1;
    double min_distance{std::numeric_limits<double>::max()};
    for (unsigned int i = 0; i < num_control_points-3; i++)
    {
        Eigen::Matrix<double, D,order> points = control_points.block(0,i,D,i+order+1);
        double distance = collision_checker.getDistanceToSphere(obstacle_center, obstacle_radius,
                points, num_points);
        if (distance < min_distance)
        {
            min_distance = distance;
        }
    }
    return min_distance;
}

template <int D>
Eigen::Matrix<double,D,1> ObstacleConstraints<D>::getObstacleCenterFromArray(double obstacle_centers[], int obstacle_num,
                                                            int num_obstacles)
{
    Eigen::Matrix<double,D,1> obstacle_center;
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
