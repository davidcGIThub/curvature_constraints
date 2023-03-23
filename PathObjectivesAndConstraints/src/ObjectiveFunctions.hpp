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

extern "C"
{
    ObjectiveFunctions<2>* ObjectiveFunctions_2(){return new ObjectiveFunctions<2>();}
    float minimize_acceleration_and_time_2(ObjectiveFunctions<2>* obj, float cont_pts[], 
        int num_control_points, float scale_factor){return obj->minimize_acceleration_and_time(
            cont_pts, num_control_points, scale_factor);}

    ObjectiveFunctions<3>* ObjectiveFunctions_3(){return new ObjectiveFunctions<3>();}
    float minimize_acceleration_and_time_3(ObjectiveFunctions<3>* obj, float cont_pts[], 
        int num_control_points, float scale_factor){return obj->minimize_acceleration_and_time(
            cont_pts, num_control_points, scale_factor);}
}

#endif