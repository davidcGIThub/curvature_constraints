#ifndef WAYPOINTCONSTRAINTS_HPP
#define WAYPOINTCONSTRAINTS_HPP
#include "DerivativeEvaluator.hpp"
#include "CBindingHelper.hpp"

template<int D>
class WaypointConstraints
{
    public:
        WaypointConstraints();
        float* velocity_at_waypoints(float cont_pts[], int num_control_points,
                                    float scale_factor, float desired_velocities[]);
};

#endif