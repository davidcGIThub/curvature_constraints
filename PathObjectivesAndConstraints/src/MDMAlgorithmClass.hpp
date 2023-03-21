#ifndef MDM_ALOGIRTHM_CLASS_HPP
#define MDM_ALGORITHM_CLASS_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>


template <int D>
class MDMAlgorithmClass
{
    public:
        MDMAlgorithmClass();
        float min_norm(Eigen::MatrixXf &points, int &num_points,
            int &max_iterations, unsigned int &initial_index,float &tolerance);
    private:
        int find_max_index(Eigen::VectorXf &matrix, int &length);
        int find_min_index(Eigen::VectorXf &matrix, int &length);
};
#endif