#ifndef MDM_ALOGIRTHM_CLASS_HPP
#define MDM_ALGORITHM_CLASS_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>


template <int N, int D>
class MDMAlgorithmClass
{
    public:
        MDMAlgorithmClass();
        double min_norm(Eigen::Matrix<double, D, N> points, 
            int max_iterations, unsigned int initial_index,double tolerance);
    private:
        int get_the_max(std::array<int, N> arr);
        int find_max_index(Eigen::VectorXd matrix, int length);
        int find_min_index(Eigen::VectorXd matrix, int length);
};
#endif