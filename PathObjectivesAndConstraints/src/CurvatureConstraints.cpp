#include "CurvatureConstraints.hpp"

template <int D>
CurvatureConstraints<D>::CurvatureConstraints()
{

}

template <int D>
float CurvatureConstraints<D>::get_spline_curvature_constraint(float cont_pts[], 
    int num_cont_pts, float max_curvature)
{
    float constraint;
    constraint = curv_bound.get_spline_curvature_bound(cont_pts, num_cont_pts) - max_curvature;
    return constraint;
}

template <int D>
float* CurvatureConstraints<D>::get_interval_curvature_constraints(float cont_pts[], int num_cont_pts, 
    float max_curvature)
{
    int order = 3;
    int num_intervals = num_cont_pts - order;
    float* constraints = new float[num_intervals];
    Eigen::VectorXf curvature_bounds(num_intervals);
    curvature_bounds = curv_bound.get_interval_curvature_bounds(cont_pts, num_cont_pts);
    for(unsigned int i = 0; i < num_intervals; i++)
    {
        constraints[i] = curvature_bounds.coeff(i) - max_curvature;
    }
    return constraints;
}

// explicit instantiation
template class CurvatureConstraints<2>;
template class CurvatureConstraints<3>;