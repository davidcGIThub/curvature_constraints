#include "gtest/gtest.h"
#include "RectPrismCollisionChecker.hpp"

TEST(RectPrismCollisionTests, AboveNotColliding)
{
    const int D{2};
    RectPrismCollisionChecker<D> checker{};
    Eigen::Matrix<double,D,1> obstacle_center;
    obstacle_center << 7, 9;
    double obstacle_radius = 2; 
    const int num_points = 5;
    Eigen::Matrix<double, D, num_points> points; 
    points << 1, 3,  -1,  2, 4.8,
              2, 5, 6.9, -4, 6;
    double distance = checker.getDistanceToObstacle(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    double true_distance = 1.0413812651491097;
    double tolerance = 0.0000000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}

TEST(RectPrismCollisionTests, BelowNotColliding)
{
    const int D{2};
    RectPrismCollisionChecker<D> checker{};
    Eigen::Matrix<double,D,1> obstacle_center;
    obstacle_center << -5, -5.1;
    double obstacle_radius = 4; 
    const int num_points = 5;
    Eigen::Matrix<double, D, num_points> points; 
    points << 1, 3,  -1,  2, 4.8,
              2, 5, 6.9, -4, 6;
    double distance = checker.getDistanceToObstacle(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    double true_distance = 0.14849370253830774;
    double tolerance = 0.0000000001;
    EXPECT_NEAR(true_distance, distance, tolerance);

}

TEST(RectPrismCollisionTests, EdgeColliding)
{
    const int D{2};
    RectPrismCollisionChecker<D> checker{};
    Eigen::Matrix<double,D,1> obstacle_center;
    obstacle_center << -2, 0;
    double obstacle_radius = 4; 
    const int num_points = 5;
    Eigen::Matrix<double, D, num_points> points; 
    points << 1.8, 3,  -1,  2.3, 4.8,
              2, 4, 6.9, -4, 5.8;
    double distance = checker.getDistanceToObstacle(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    double true_distance = -3;
    double tolerance = 0.0000000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}

TEST(RectPrismCollisionTests, CenterInsideOfObstacle)
{
    const int D{2};
    RectPrismCollisionChecker<D> checker{};
    Eigen::Matrix<double,D,1> obstacle_center;
    obstacle_center << 0, 0;
    double obstacle_radius = 4; 
    const int num_points = 5;
    Eigen::Matrix<double, D, num_points> points; 
    points << 1.8, 3,  -1,  2.3, 4.8,
              2, 4, 6.9, -4, 5.8;
    double distance = checker.getDistanceToObstacle(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    double true_distance = -4;
    double tolerance = 0.0000000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}

TEST(RectPrismCollisionTests, Colliding3D)
{
    const int D{3};
    RectPrismCollisionChecker<D> checker{};
    Eigen::Matrix<double,D,1> obstacle_center;
    obstacle_center << -1.2, -5, -3;
    double obstacle_radius = 5; 
    const int num_points = 5;
    Eigen::Matrix<double, D, num_points> points; 
    points << 1.8, 3,  -1,  2.3, 4.8,
              2  , 4, 6.9, -4  , 5.8,
              -1 , 5,   2,    6,   7;
    double distance = checker.getDistanceToObstacle(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    double true_distance = -2.7550055679356351;
    double tolerance = 0.0000000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
    
}

TEST(RectPrismCollisionTests, NotColliding3D)
{
    const int D{3};
    RectPrismCollisionChecker<D> checker{};
    Eigen::Matrix<double,D,1> obstacle_center;
    obstacle_center << -8, -5, 20;
    double obstacle_radius = 5; 
    const int num_points = 5;
    Eigen::Matrix<double, D, num_points> points; 
    points << 1.8, 3,  -1,  2.3, 4.8,
              2  , 4, 6.9, -4  , 5.8,
              -1 , 5,   2,    6,   7;
    double distance = checker.getDistanceToObstacle(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    double true_distance = 9.798648586948742;
    double tolerance = 0.0000000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}