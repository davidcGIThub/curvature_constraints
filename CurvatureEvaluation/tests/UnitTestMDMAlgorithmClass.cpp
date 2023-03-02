#include "MDMAlgorithmClass.hpp"
#include "gtest/gtest.h"
#include <iostream>
#include <array>
#include <chrono>
#include <ctime>  

TEST(MDMAlgorithmTests, MDMClassRandomCase1)
{
    Eigen::Matrix<double, 2, 6> points;
    points << 8, 5, 13, 10,  8, 14,
             13, 9, 10, 10, 12, 10;
    int max_iter = 5000;
    unsigned int init_index = 1;
    double true_min_norm = 10.295630140987;
    double tolerance = 0.000001;
    // auto start = std::chrono::system_clock::now();
    MDMAlgorithmClass<6,2> mdm{};
    double min_norm = mdm.min_norm(points, max_iter, init_index, tolerance);
    // auto end = std::chrono::system_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end-start;
    // std::cout << "elapsed time class: " << elapsed_seconds.count() << "s" << std::endl;
    double error_tolerance = 0.00000001;
    EXPECT_NEAR(min_norm, true_min_norm, tolerance);
}
