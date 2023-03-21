#ifndef CBINDINGHELPER_HPP
#define CBINDINGHELPER_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "gtest/gtest_prod.h"

template <int D>
class CBindingHelper
{
    public:
        CBindingHelper();
        Eigen::Matrix<float,D,4> array_section_to_eigen(float cont_pts[], int &num_cps, unsigned int &index);
        Eigen::MatrixXf array_to_eigen(float cont_pts[], int &num_cps);
};

#endif