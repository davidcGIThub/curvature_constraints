#include "MDMAlgorithm.hpp"
#include <vector>
#include <iostream>

namespace MDM_Algorithm
{

    double mdm_algorithm(Eigen::Matrix<double, 2, 6> points, 
        int max_iterations, unsigned int initial_index, double tolerance)
    {
        int iterations{0};
        double delta_p{1.0};
        int num_points{6};
        double p_vector[num_points] = {0};
        std::vector<int> supp_vector;
        int index{0};
        Eigen::Matrix<double,2,1> current_vector = points.block(0,index,2,1);
        supp_vector.push_back(index);
        p_vector[index] = 1;
        while( delta_p > 0.000001 && iterations < max_iterations && supp_vector.size() > 0)
        {
            Eigen::MatrixXd mat;
            Eigen::VectorXd mult;
            mat.resize(supp_vector.size(),2);
            for(unsigned int i = 0; i<supp_vector.size(); i++)
            {
                for(unsigned int j = 0; j<2; j++)
                {
                    mat(i,j) = points(j,supp_vector[i]);
                }
            }
            mult = mat*current_vector;
            int max_index = find_max_index(mult,supp_vector.size());
            max_index = supp_vector[max_index];
            mult = current_vector.transpose()*points;
            int min_index = find_min_index(mult,6);
            Eigen::Matrix<double,2,1> min_point = points.block(0,min_index,2,1);
            Eigen::Matrix<double,2,1> max_point = points.block(0,max_index,2,1);
            Eigen::Matrix<double,2,1> diff = max_point - min_point;
            delta_p = diff.dot(current_vector);
            if (delta_p > tolerance)
            {
                double diff_norm = diff.norm();
                double t_param = delta_p/(p_vector[max_index] * diff_norm * diff_norm);
                if (t_param >= 1)
                {
                    t_param = 1.0;
                }
                current_vector = current_vector - (t_param * p_vector[max_index] * diff);
                supp_vector.clear();
                double temp1 = t_param*p_vector[max_index];
                double temp2 = 1 - t_param;
                p_vector[min_index] += temp1;
                p_vector[max_index] *= temp2;
                for(int i = 0; i< num_points; i++)
                {
                    if(p_vector[i] > tolerance)
                    {
                        supp_vector.push_back(i);
                    }
                }
            iterations += 1;
            }
        }
        double min_norm = current_vector.norm();
        return min_norm;
    }

    int find_max_index(Eigen::VectorXd matrix, int length)
    {
        double max_size = std::numeric_limits<double>::min();
        int max_index = 0;
        for(unsigned int i = 0; i<length; i++)
        {
            if (matrix(i) > max_size)
            {
                max_size = matrix(i);
                max_index = i;
            }
        }
        return max_index;
    }

    int find_min_index(Eigen::VectorXd matrix, int length)
    {
        double min_size = std::numeric_limits<double>::max();
        int min_index = 0;
        for(unsigned int i = 0; i<length; i++)
        {
            if (matrix(i) < min_size)
            {
                min_size = matrix(i);
                min_index = i;
            }
        }
        return min_index;
    }

    template <int N>
    int get_the_max(std::array<int, N> arr)
    {
        return arr.size();
    }
}
