#include "MDMAlgorithm.hpp"
#include "gtest/gtest.h"
#include <iostream>
#include <array>
#include <chrono>
#include <ctime> 

TEST(MDMAlgorithmTests, MDMRandomCase1)
{
    Eigen::Matrix<double, 2, 6> points;
    points << 8, 5, 13, 10,  8, 14,
             13, 9, 10, 10, 12, 10; 
    int max_iter = 5000;
    unsigned int init_index = 1;
    double true_min_norm = 10.295630140987;
    double tolerance = 0.000001;
    // auto start = std::chrono::system_clock::now();
    double min_norm = MDM_Algorithm::mdm_algorithm(points, max_iter, init_index, tolerance);
    // auto end = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end-start;
    // std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;
    double error_tolerance = 0.00000001;
    EXPECT_NEAR(min_norm, true_min_norm, tolerance);
}

TEST(MDMAlgorithmTests, MDMRandomCase2)
{
    Eigen::Matrix<double, 2, 6> points;
    points << 10, 8, 14, 9, 14, 5,
              6, 10, 14, 14, 8, 9; 
    int max_iter = 5000;
    unsigned int init_index = 1;
    double true_min_norm = 10.289915108550531;
    double tolerance = 0.000001;
    double min_norm = MDM_Algorithm::mdm_algorithm(points, max_iter, init_index, tolerance);
    double error_tolerance = 0.00000001;
    EXPECT_NEAR(true_min_norm, min_norm, error_tolerance);
}

TEST(MDMAlgorithmTests, MDMStraightLine)
{
    Eigen::Matrix<double, 2, 6> points;
    points << 1,2,3,4,5,6,
              1,2,3,4,5,6; 
    int max_iter = 5000;
    unsigned int init_index = 1;
    double true_min_norm = 1.4142135623730951;
    double tolerance = 0.000001;
    double min_norm = MDM_Algorithm::mdm_algorithm(points, max_iter, init_index, tolerance);
    double error_tolerance = 0.00000001;
    EXPECT_NEAR(true_min_norm, min_norm, error_tolerance);
}

TEST(MDMAlgorithmTests, MDMStraightLine2)
{
    Eigen::Matrix<double, 2, 6> points;
    points << 6,5,4,3,2,1,
              1,2,3,4,5,6; 
    int max_iter = 5000;
    unsigned int init_index = 1;
    double true_min_norm = 4.949747468305833;
    double tolerance = 0.000001;
    double min_norm = MDM_Algorithm::mdm_algorithm(points, max_iter, init_index, tolerance);
    double error_tolerance = 0.00000001;
    EXPECT_NEAR(true_min_norm, min_norm, error_tolerance);
}

TEST(MDMAlgorithmTests, MDMSamePoints)
{
    Eigen::Matrix<double, 2, 6> points;
    points << 5, 5, 5, 5, 5, 5,
              6, 6, 6, 6, 6, 6;
    int max_iter = 5000;
    double tolerance = 0.000001;
    unsigned int init_index = 1;
    double true_min_norm = 7.810249675906654;
    double min_norm = MDM_Algorithm::mdm_algorithm(points, max_iter, init_index, tolerance);
    double error_tolerance = 0.00000001;
    EXPECT_NEAR(true_min_norm, min_norm, error_tolerance);
}

TEST(MDMAlgorithmTests, FixingCode)
{
    Eigen::Matrix<double, 2, 6> points;
    points << 9, 6, 12, 14, 13, 7,
              5, 13, 11, 11, 8, 8;
    int max_iter = 5000;
    double tolerance = 0.000001;
    unsigned int init_index = 1;
    double true_min_norm = 10.261953630166738;
    double min_norm = MDM_Algorithm::mdm_algorithm(points, max_iter, init_index, tolerance);
    double error_tolerance = 0.00000001;
    EXPECT_NEAR(true_min_norm, min_norm, error_tolerance);
}



