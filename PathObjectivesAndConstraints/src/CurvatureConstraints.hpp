#ifndef CURVATURECONSTRAINTS_HPP
#define CURVATURECONSTRAINTS_HPP
#include "ThirdOrderCurvatureBounds.hpp"
#include "CBindingHelper.hpp"

template <int D>
class CurvatureConstraints
{
    public:
        CurvatureConstraints();
        float get_spline_curvature_constraint(float cont_pts[], int num_cont_pts, float max_curvature);
        float* get_interval_curvature_constraints(float cont_pts[], int num_cont_pts, float max_curvature);

    private:
        ThirdOrderCurvatureBounds<D> curv_bound{};
        CBindingHelper<D> cbind_help{};

};


extern "C"
{
    CurvatureConstraints<2>* CurvatureConstraints_2(){return new CurvatureConstraints<2>();}
    float get_spline_curvature_constraint_2(CurvatureConstraints<2>* obj, float cont_pts[], 
        int num_cont_pts, float max_curvature){return obj->get_spline_curvature_constraint(
        cont_pts, num_cont_pts, max_curvature);}

    CurvatureConstraints<3>* CurvatureConstraints_3(){return new CurvatureConstraints<3>();}
    float get_spline_curvature_constraint_3(CurvatureConstraints<3>* obj, float cont_pts[],
        int num_cont_pts, float max_curvature){return obj->get_spline_curvature_constraint(
        cont_pts, num_cont_pts, max_curvature);}
}

#endif