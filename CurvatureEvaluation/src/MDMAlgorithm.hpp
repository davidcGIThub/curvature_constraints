#ifndef MDM_ALOGIRTHM_JPP
#define MDM_ALGORITHM_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <array>

namespace MDM_Algorithm
{
    template <int N>
    int get_the_max(std::array<int, N> arr);
    double mdm_algorithm(Eigen::Matrix<double, 2, 6> points, 
        int max_iterations, unsigned int initial_index,double tolerance);
    int find_max_index(Eigen::VectorXd matrix, int length);
    int find_min_index(Eigen::VectorXd matrix, int length);
}
#endif

// template <typename T> T myMax(T x, T y)
// {
//     return (x > y) ? x : y;
// }
  
// int main()
// {
//     cout << myMax<int>(3, 7) << endl; // Call myMax for int
//     cout << myMax<double>(3.0, 7.0)
//          << endl; // call myMax for double
//     cout << myMax<char>('g', 'e')
//          << endl; // call myMax for char
  
//     return 0;
// }