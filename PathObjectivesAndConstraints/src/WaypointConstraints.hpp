#ifndef WAYPOINTCONSTRAINTS_HPP
#define WAYPOINTCONSTRAINTS_HPP
#include "DerivativeEvaluator.hpp"
#include "CBindingHelper.hpp"

template<int D>
class WaypointConstraints
{
    public:
        WaypointConstraints();
        double* velocity_at_waypoints(double cont_pts[], int num_control_points,
                                    double scale_factor, double desired_velocities[]);
};

#endif