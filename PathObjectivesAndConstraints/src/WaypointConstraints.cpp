#include "WaypointConstraints.hpp"
#include <iostream>

template<int D>
WaypointConstraints<D>::WaypointConstraints()
{

}

template<int D>
double* WaypointConstraints<D>::velocity_at_waypoints_constraints(double cont_pts[], int num_cps, double scale_factor, double desired_velocities[])
{
    int order = 3;
    double* constraints = new double[D*2];
    DerivativeEvaluator<D> d_dt_eval{};
    CBindingHelper<D> helper{};
    double t0 = 0;
    double tf = scale_factor;
    unsigned int start_index = 0;
    unsigned int end_index = num_cps - order -1;
    Eigen::Matrix<double,D,4> initial_control_points = helper.array_section_to_eigen(cont_pts, num_cps, start_index);
    Eigen::Matrix<double,D,4> final_control_points = helper.array_section_to_eigen(cont_pts, num_cps, end_index);
    Eigen::Matrix<double,D,1> velocity_vector_0 = d_dt_eval.calculate_velocity_vector(t0, initial_control_points, scale_factor);
    Eigen::Matrix<double,D,1> velocity_vector_f = d_dt_eval.calculate_velocity_vector(tf, final_control_points, scale_factor);
    for(unsigned int i = 0; i < D; i++)
    {
        constraints[2*i] = velocity_vector_0.coeff(i) - desired_velocities[2*i];
        constraints[2*i+1] = velocity_vector_f.coeff(i) - desired_velocities[2*i+1];
    }
    return constraints;
}


//explicit instantiation
template class WaypointConstraints<2>;
template class WaypointConstraints<3>;
