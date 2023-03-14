#include "WaypointConstraints.hpp"
#include <iostream>

template<int D>
WaypointConstraints<D>::WaypointConstraints()
{

}

template<int D>
double* WaypointConstraints<D>::velocity_at_waypoints(double cont_pts[], int num_cps, double scale_factor, double desired_velocities[])
{
    int order = 3;
    double* constraints = new double[D*2];
    DerivativeEvaluator<D> d_dt_eval{};
    CBindingHelper<D> helper{};
    double t0 = 0;
    double tf = scale_factor;
    Eigen::Matrix<double,D,4> initial_control_points = helper.array_section_to_eigen(cont_pts, num_cps, 0);
    Eigen::Matrix<double,D,4> final_control_points = helper.array_section_to_eigen(cont_pts, num_cps, num_cps - 3 -1);
    Eigen::Matrix<double,D,1> velocity_vector_0 = d_dt_eval.calculate_velocity_vector(t0, initial_control_points, scale_factor);
    Eigen::Matrix<double,D,1> velocity_vector_f = d_dt_eval.calculate_velocity_vector(tf, final_control_points, scale_factor);
    for(unsigned int i = 0; i < D; i++)
    {
        constraints[2*i] = velocity_vector_0.coeff(i) - desired_velocities[2*i];
        constraints[2*i+1] = velocity_vector_f.coeff(i) - desired_velocities[2*i+1];
    }
    return constraints;
}

// int* arr2_return()
// {
//     int* information = new int[9];
//     for(int k=0; k<9; k++)
//     {
//         information[k] = k;
//     }
//     return information;
// }

    // def __create_waypoint_velocity_constraint(self, velocities):
    //     def velocity_constraint_function(variables):
    //         constraints = np.zeros(self._dimension*2)
    //         control_points = np.reshape(variables[0:self._num_control_points*self._dimension], \
    //         (self._dimension,self._num_control_points))
    //         segement_1_control_points = control_points[:,0:self._order+1]
    //         segement_2_control_points = control_points[:,self._num_control_points-self._order-1:]
    //         scale_factor = variables[-1]
    //         T_0 = get_T_derivative_vector(self._order,0,0,1,scale_factor)
    //         T_f = get_T_derivative_vector(self._order,scale_factor,0,1,scale_factor)
    //         start_velocity = np.dot(segement_1_control_points,np.dot(self._M,T_0)).flatten()
    //         end_velocity = np.dot(segement_2_control_points,np.dot(self._M,T_f)).flatten()
    //         desired_start_velocity = velocities[:,0]
    //         desired_end_velocity = velocities[:,1]
    //         constraints[0:self._dimension] = start_velocity - desired_start_velocity
    //         constraints[self._dimension:] = end_velocity - desired_end_velocity
    //         return constraints
    //     lower_bound = 0
    //     upper_bound = 0
    //     velocity_vector_constraint = NonlinearConstraint(velocity_constraint_function, lb= lower_bound, ub=upper_bound)
    //     return velocity_vector_constraint

//explicit instantiation
template class WaypointConstraints<2>;
template class WaypointConstraints<3>;
