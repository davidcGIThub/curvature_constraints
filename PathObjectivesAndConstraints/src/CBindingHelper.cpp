#include "CBindingHelper.hpp"

template<int D>
CBindingHelper<D>::CBindingHelper()
{
    
}

template<int D>
Eigen::Matrix<double,D,4> CBindingHelper<D>::array_section_to_eigen(double cont_pts[], int num_cps, unsigned int index)
{
    Eigen::Matrix<double,D,4> interval_control_points;
    if (D == 2)
    {
        interval_control_points << cont_pts[index], cont_pts[index+1], cont_pts[index+2], cont_pts[index+3],
            cont_pts[num_cps+index], cont_pts[num_cps+index+1], cont_pts[num_cps+index+2], cont_pts[num_cps+index+3];
    }
    else
    {
        interval_control_points << cont_pts[index], cont_pts[index+1], cont_pts[index+2], cont_pts[index+3],
            cont_pts[num_cps+index], cont_pts[num_cps+index+1], cont_pts[num_cps+index+2], cont_pts[num_cps+index+3],
            cont_pts[2*num_cps+index], cont_pts[2*num_cps+index+1], cont_pts[2*num_cps+index+2], cont_pts[2*num_cps+index+3];
    }
    return interval_control_points;
}

// explicit instantiations
template class CBindingHelper<2>;
template class CBindingHelper<3>;