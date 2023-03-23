#include "gtest/gtest.h"
#include "ConvexHullCollisionChecker.hpp"

TEST(ConvexHullCollisionTests, NotColliding)
{
    const int D{2};
    ConvexHullCollisionChecker<D> checker{};
    Eigen::Matrix<float,D,1> sphere_center;
    sphere_center << 17, 15;
    float sphere_radius = 3; 
    int num_spline_points = 6;
    Eigen::MatrixXf spline_points(D, num_spline_points); 
    spline_points << 8,  5, 13, 10,  8, 14,
                     8,  9, 10, 10, 12, 10;
    float distance = checker.getDistanceToSphere(sphere_center, sphere_radius, 
                                                    spline_points, num_spline_points);
    float true_distance = 2.8309518976313477;
    float tolerance = 0.00001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}

TEST(ConvexHullCollisionTests, RadiusColliding)
{
    const int D{2};
    ConvexHullCollisionChecker<D> checker{};
    Eigen::Matrix<float,D,1> obstacle_center;
    obstacle_center << 4, 7;
    float obstacle_radius = 5; 
    int num_points = 6;
    Eigen::MatrixXf points(D, num_points);
    points << 8,  4, 13, 10,  8, 14,
              8, 10, 10, 10, 12, 10;
    float distance = checker.getDistanceToSphere(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    float true_distance = -2.316718414998114;
    float tolerance = 0.0000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}

TEST(ConvexHullCollisionTests, CenterColliding)
{
    const int D{2};
    ConvexHullCollisionChecker<D> checker{};
    Eigen::Matrix<float,D,1> obstacle_center;
    obstacle_center << 9, 10.4;
    float obstacle_radius = 3; 
    int num_points = 6;
    Eigen::MatrixXf points(D, num_points); 
    points << 8,  4, 13, 10,  8, 14, 8, 10, 10, 10, 12, 10;
    float distance = checker.getDistanceToSphere(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    float true_distance = -7.859687576256714;
    float tolerance = 0.0000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}

TEST(ConvexHullCollisionTests, Colliding3D)
{
    const int D{3};
    ConvexHullCollisionChecker<D> checker{};
    Eigen::Matrix<float,D,1> obstacle_center;
    obstacle_center << 4, 7, 1.5;
    float obstacle_radius = 4; 
    int num_points = 6;
    Eigen::MatrixXf points(D, num_points); 
    points << 8,   4,    13,  10,     8,  14,
              8,  10,    10,  10,    12,  10,
              2,   7,   3.7,   4,   7.9,  3.4;
    float distance = checker.getDistanceToSphere(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    float true_distance = -0.21699707756002873;
    float tolerance = 0.0000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}

TEST(ConvexHullCollisionTests, NotColliding3D)
{
    const int D{3};
    ConvexHullCollisionChecker<D> checker{};
    Eigen::Matrix<float,D,1> obstacle_center;
    obstacle_center << 2, 4, 0;
    float obstacle_radius = 1; 
    int num_points = 6;
    Eigen::MatrixXf points(D, num_points); 
    points << 8,   4,    13,  10,     8,  14,
              8,  10,    10,  10,    12,  10,
              2,   7,   3.7,   4,   7.9,  3.4;
    float distance = checker.getDistanceToSphere(obstacle_center, obstacle_radius, 
                                                    points, num_points);
    float true_distance = 6.429670248464965;
    float tolerance = 0.0000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}

TEST(ConvexHullCollisionTests, ConservativeNotColliding3D)
{
    const int D{3};
    ConvexHullCollisionChecker<D> checker{};
    Eigen::Matrix<float,D,1> obstacle_center;
    obstacle_center << 2, 4, 0;
    float obstacle_radius = 1; 
    int num_points = 6;
    Eigen::MatrixXf points(D, num_points);  
    points << 8,   4,    13,  10,     8,  14,
              8,  10,    10,  10,    12,  10, 
              2,   7,   3.7,   4,   7.9,  3.4;
    float distance = checker.getConservativeDistanceToSphere(obstacle_center,
                 obstacle_radius, points, num_points);
    float true_distance = 12.404153024776605;
    float tolerance = 0.00001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}

TEST(ConvexHullCollisionTests, ConservativeColliding3D)
{
    const int D{3};
    ConvexHullCollisionChecker<D> checker{};
    Eigen::Matrix<float,D,1> obstacle_center;
    obstacle_center << 9, 10, 5;
    float obstacle_radius = 1; 
    int num_points = 6;
    Eigen::MatrixXf points(D, num_points);  
    points << 8,   4,    13,  10,     8,  14,
              8,  10,    10,  10,    12,  10, 
              2,   7,   3.7,   4,   7.9,  3.4;
    float distance = checker.getConservativeDistanceToSphere(obstacle_center,
                 obstacle_radius, points, num_points);
    float true_distance = -1;
    float tolerance = 0.0000001;
    EXPECT_NEAR(true_distance, distance, tolerance);
}


TEST(ConvexHullCollisionTests, CheckIfCollidesDoes)
{
    const int D{2};
    ConvexHullCollisionChecker<D> checker{};
    Eigen::Matrix<float,D,1> obstacle_center;
    obstacle_center << 7, 7;
    float obstacle_radius = 5; 
    int num_points = 6;
    Eigen::MatrixXf points(D, num_points);  
    points << 8,   4,    13,  10,     8,  14,
              8,  10,    10,  10,    12,  10;
    bool collides = checker.checkIfCollides(obstacle_center,
                 obstacle_radius, points, num_points);
    float true_collides = true;
    EXPECT_EQ(true_collides, collides);
}

TEST(ConvexHullCollisionTests, CheckIfCollidesDoesnt)
{
    const int D{2};
    ConvexHullCollisionChecker<D> checker{};
    Eigen::Matrix<float,D,1> obstacle_center;
    obstacle_center << 4, 3;
    float obstacle_radius = 1; 
    int num_points = 6;
    Eigen::MatrixXf points(D, num_points);  
    points << 8,   4,    13,  10,     8,  14,
              8,  10,    10,  10,    12,  10;
    bool collides = checker.checkIfCollides(obstacle_center,
                 obstacle_radius, points, num_points);
    float true_collides = false;
    EXPECT_EQ(true_collides, collides);
}