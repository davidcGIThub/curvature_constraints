#include "MDMAlgorithmClass.hpp"
#include <vector>
#include <iostream>



template <int D>
MDMAlgorithmClass<D>::MDMAlgorithmClass(){}

template <int D>
float MDMAlgorithmClass<D>::min_norm(Eigen::MatrixXf &points, int &num_points,
            int &max_iterations, unsigned int &initial_index, float &tolerance)
{
    int iterations{0};
    float delta_p{1.0};
    float p_vector[num_points] = {0};
    std::vector<int> supp_vector;
    int index{0};
    Eigen::Matrix<float,D,1> current_vector = points.block(0,index,D,1);
    supp_vector.push_back(index);
    p_vector[index] = 1;
    while( delta_p > 0.000001 && iterations < max_iterations && supp_vector.size() > 0)
    {
        Eigen::MatrixXf mat;
        Eigen::VectorXf mult;
        mat.resize(supp_vector.size(),D);
        for(unsigned int i = 0; i<supp_vector.size(); i++)
        {
            for(unsigned int j = 0; j<D; j++)
            {
                mat(i,j) = points(j,supp_vector[i]);
            }
        }
        mult = mat*current_vector;
        int length = supp_vector.size();
        int max_index = find_max_index(mult,length);
        max_index = supp_vector[max_index];
        mult = current_vector.transpose()*points;
        int min_index = find_min_index(mult,num_points);
        Eigen::Matrix<float,D,1> min_point = points.block(0,min_index,D,1);
        Eigen::Matrix<float,D,1> max_point = points.block(0,max_index,D,1);
        Eigen::Matrix<float,D,1> diff = max_point - min_point;
        delta_p = diff.dot(current_vector);
        if (delta_p > tolerance)
        {
            float diff_norm = diff.norm();
            float t_param = delta_p/(p_vector[max_index] * diff_norm * diff_norm);
            if (t_param >= 1)
            {
                t_param = 1.0;
            }
            current_vector = current_vector - (t_param * p_vector[max_index] * diff);
            supp_vector.clear();
            float temp1 = t_param*p_vector[max_index];
            float temp2 = 1 - t_param;
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
    float min_norm = current_vector.norm();
    return min_norm;
}

template <int D>
int MDMAlgorithmClass<D>::find_max_index(Eigen::VectorXf &matrix, int &length)
{
    float max_size = std::numeric_limits<float>::min();
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

template <int D>
int MDMAlgorithmClass<D>::find_min_index(Eigen::VectorXf &matrix, int &length)
{
    float min_size = std::numeric_limits<float>::max();
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

template class MDMAlgorithmClass<2>;
template class MDMAlgorithmClass<3>;
