#ifndef OBJECTIVEFUNCTIONS_HPP
#define OBJECTIVEFUNCTIONS_HPP
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include "CBindingHelper.hpp"

template <int D>
class ObjectiveFunctions
{
    public:
        ObjectiveFunctions();
        float minimize_acceleration_and_time(float cont_pts[], int num_control_points, float scale_factor);
        float minimize_distance_and_time(float cont_pts[], int num_control_points, float scale_factor);
    private:
        CBindingHelper<D> cbind_help{};
        float minimize_acceleration(float cont_pts[], int &num_control_points);
        float minimize_distance(float cont_pts[], int &num_control_points);
};


#endif