#ifndef WAYPOINTCONSTRAINTS_HPP
#define WAYPOINTCONSTRAINTS_HPP
#include "DerivativeEvaluator.hpp"
#include "CBindingHelper.hpp"

template<int D>
class WaypointConstraints
{
    public:
        WaypointConstraints();
        double* velocity_at_waypoints_constraints(double cont_pts[], int num_control_points,
                                    double scale_factor, double desired_velocities[]);
};

extern "C"
{
    WaypointConstraints<2>* WaypointConstraints_2(){return new WaypointConstraints<2>();}
    double* velocity_at_waypoints_constraints_2(WaypointConstraints<2>* obj, double cont_pts[], int num_control_points,
            double scale_factor, double desired_velocities[]){return obj->velocity_at_waypoints_constraints(
            cont_pts, num_control_points, scale_factor, desired_velocities);}

    WaypointConstraints<3>* WaypointConstraints_3(){return new WaypointConstraints<3>();}
    double* velocity_at_waypoints_constraints_3(WaypointConstraints<3>* obj, double cont_pts[], int num_control_points,
            double scale_factor, double desired_velocities[]){return obj->velocity_at_waypoints_constraints(
            cont_pts, num_control_points, scale_factor, desired_velocities);}
}

#endif