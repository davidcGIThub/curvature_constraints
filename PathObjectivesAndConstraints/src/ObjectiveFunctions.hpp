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
        double minimize_acceleration_and_time(double cont_pts[], int num_control_points, double scale_factor);
        double minimize_distance_and_time(double cont_pts[], int num_control_points, double scale_factor);
    private:
        CBindingHelper<D> cbind_help{};
        double minimize_acceleration(double cont_pts[], int &num_control_points);
        double minimize_distance(double cont_pts[], int &num_control_points);
};


#endif