#include "ThirdOrderCurvatureEvaluator.hpp"
#include "CubicEquationSolver.hpp"
#include "DerivativeEvaluator.hpp"
#include <iostream>
#include <stdexcept>

template <int D>
ThirdOrderCurvatureEvaluator<D>::ThirdOrderCurvatureEvaluator()
{

}

template <int D>
double ThirdOrderCurvatureEvaluator<D>::find_spline_curvature_bound(double cont_pts[], int num_control_points)
{
    double max_curvature{0};
    double curvature;
    for (unsigned int i = 0; i < num_control_points-3; i++)
    {
        Eigen::Matrix<double,D,4> interval_control_points = array_section_to_eigen(cont_pts, num_control_points, i);
        curvature = evaluate_interval_curvature(interval_control_points);
        if (curvature > max_curvature)
        {
            max_curvature = curvature;
        }
    }
    return max_curvature;
}

template<int D>
Eigen::Matrix<double,D,4> ThirdOrderCurvatureEvaluator<D>::array_section_to_eigen(double cont_pts[], int step, unsigned int index)
{
    Eigen::Matrix<double,D,4> interval_control_points;
    if (D == 2)
    {
        interval_control_points << cont_pts[index], cont_pts[index+1], cont_pts[index+2], cont_pts[index+3],
            cont_pts[step+index], cont_pts[step+index+1], cont_pts[step+index+2], cont_pts[step+index+3];
    }
    else
    {
        interval_control_points << cont_pts[index], cont_pts[index+1], cont_pts[index+2], cont_pts[index+3],
            cont_pts[step+index], cont_pts[step+index+1], cont_pts[step+index+2], cont_pts[step+index+3],
            cont_pts[2*step+index], cont_pts[2*step+index+1], cont_pts[2*step+index+2], cont_pts[2*step+index+3];
    }
    return interval_control_points;
}

template <int D>
double ThirdOrderCurvatureEvaluator<D>::evaluate_interval_curvature(Eigen::Matrix<double,D,4> &control_points)
{
    double scale_factor = 1;
    DerivativeEvaluator<D> d_dt_eval{};
    std::array<double,2> min_velocity_and_time = d_dt_eval.find_min_velocity_and_time(control_points,scale_factor);
    double min_velocity = min_velocity_and_time[0];
    double time_at_min_velocity = min_velocity_and_time[1];
    double max_cross_term = find_maximum_cross_term(control_points,scale_factor);
    std::array<double,2> max_acceleration_and_time = d_dt_eval.find_max_acceleration_and_time(control_points, scale_factor);
    double max_acceleration = max_acceleration_and_time[0];
    double curvature_bound;
    if (min_velocity <= 1.0e-15)
    {
        double acceleration_at_min_vel = 
            d_dt_eval.calculate_acceleration_magnitude(time_at_min_velocity,
                control_points, scale_factor);
        if(acceleration_at_min_vel == 0)
        {
            curvature_bound = 0;
        }
        else
        {
            curvature_bound = std::numeric_limits<double>::max();
        }
    }
    else
    {
        curvature_bound = max_acceleration/(min_velocity*min_velocity);
        double curvature = max_cross_term/(min_velocity*min_velocity*min_velocity);
        if (curvature < curvature_bound)
        {
            curvature_bound = curvature;
        }
    }
    return curvature_bound;
}

template <int D>
double ThirdOrderCurvatureEvaluator<D>::find_maximum_cross_term(Eigen::Matrix<double,D,4> &control_points, double &scale_factor)
{
    Eigen::Vector4d coeficients;
    if (D == 2) {coeficients = get_2D_cross_coefficients(control_points);}
    else {coeficients = get_3D_cross_coefficients(control_points);}
    std::array<double,3> roots = CubicEquationSolver::solve_equation(coeficients(0),
        coeficients(1), coeficients(2), coeficients(3));
    double t0 = 0;
    double tf = 1.0;
    double max_cross_term = calculate_cross_term_magnitude(t0,control_points,scale_factor);
    double cross_term_at_tf = calculate_cross_term_magnitude(tf,control_points,scale_factor);
    if (cross_term_at_tf > max_cross_term)
    {
        max_cross_term = cross_term_at_tf;
    }
    for(int index = 0; index < 3; index++)
    {
        double root = roots[index];
        if(root > 0 && root < 1.0)
        {
            double cross_term =  calculate_cross_term_magnitude(root, control_points,scale_factor);
            if (cross_term > max_cross_term)
            {
                max_cross_term = cross_term;
            }
        }
    }
    return max_cross_term;
}

template <int D>
Eigen::Vector4d ThirdOrderCurvatureEvaluator<D>::get_2D_cross_coefficients(Eigen::Matrix<double,D,4> &control_points)
{
    double p0x = control_points(0,0);
    double p0y = control_points(1,0);
    double p1x = control_points(0,1);
    double p1y = control_points(1,1);
    double p2x = control_points(0,2);
    double p2y = control_points(1,2);
    double p3x = control_points(0,3);
    double p3y = control_points(1,3);
    double c_3 = ((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) 
        *((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - 
        (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2);
    double c_2 = - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y))*
        ((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - 
        (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2) - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*
        (p0y - 3*p1y + 3*p2y - p3y))*((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)* 
        (p0y - 3*p1y + 3*p2y - p3y));
    double c_1 = ((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y - 2*p1y + p2y)*
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) + 
        pow((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y),2);
    double c_0 = -((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y/2 - p2y/2)* 
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y));
    Eigen::Vector4d coeficients;
    coeficients << c_0, c_1, c_2, c_3;
    return coeficients;
}

template <int D>
Eigen::Vector4d ThirdOrderCurvatureEvaluator<D>::get_3D_cross_coefficients(Eigen::Matrix<double,D,4> &control_points)
{
    double p0x = control_points(0,0);
    double p0y = control_points(1,0);
    double p0z = control_points(2,0);
    double p1x = control_points(0,1);
    double p1y = control_points(1,1);
    double p1z = control_points(2,1);
    double p2x = control_points(0,2);
    double p2y = control_points(1,2);
    double p2z = control_points(2,2);
    double p3x = control_points(0,3);
    double p3y = control_points(1,3);
    double p3z = control_points(2,3);
    double c_3 = ((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y))
        *((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - 
        (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2) + ((p0z - 2*p1z + p2z)*
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0z - 3*p1z + 3*p2z - p3z))* 
        ((p0x*p1z)/2 - (p1x*p0z)/2 - p0x*p2z + p2x*p0z + (p0x*p3z)/2 + (3*p1x*p2z)/2 - (3*p2x*p1z)/2 - 
        (p3x*p0z)/2 - p1x*p3z + p3x*p1z + (p2x*p3z)/2 - (p3x*p2z)/2) + ((p0z - 2*p1z + p2z)* 
        (p0y - 3*p1y + 3*p2y - p3y) - (p0y - 2*p1y + p2y)*(p0z - 3*p1z + 3*p2z - p3z))*((p0y*p1z)/2 - 
        (p1y*p0z)/2 - p0y*p2z + p2y*p0z + (p0y*p3z)/2 + (3*p1y*p2z)/2 - (3*p2y*p1z)/2 - (p3y*p0z)/2 - 
        p1y*p3z + p3y*p1z + (p2y*p3z)/2 - (p3y*p2z)/2);
    double c_2 = - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y))*
        ((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - 
        (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2) - ((p0z/2 - p2z/2)*
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0z - 3*p1z + 3*p2z - p3z))*((p0x*p1z)/2 - 
        (p1x*p0z)/2 - p0x*p2z + p2x*p0z + (p0x*p3z)/2 + (3*p1x*p2z)/2 - (3*p2x*p1z)/2 - (p3x*p0z)/2 - 
        p1x*p3z + p3x*p1z + (p2x*p3z)/2 - (p3x*p2z)/2) - ((p0z/2 - p2z/2)*(p0y - 3*p1y + 3*p2y - p3y) - 
        (p0y/2 - p2y/2)*(p0z - 3*p1z + 3*p2z - p3z))*((p0y*p1z)/2 - (p1y*p0z)/2 - p0y*p2z + p2y*p0z + 
        (p0y*p3z)/2 + (3*p1y*p2z)/2 - (3*p2y*p1z)/2 - (p3y*p0z)/2 - p1y*p3z + p3y*p1z + (p2y*p3z)/2 - 
        (p3y*p2z)/2) - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*
        (p0y - 3*p1y + 3*p2y - p3y))*((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)* 
        (p0y - 3*p1y + 3*p2y - p3y)) - ((p0z/2 - p2z/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*
        (p0z - 3*p1z + 3*p2z - p3z))*((p0z - 2*p1z + p2z)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*
        (p0z - 3*p1z + 3*p2z - p3z)) - ((p0z/2 - p2z/2)*(p0y - 3*p1y + 3*p2y - p3y) - (p0y/2 - p2y/2)*
        (p0z - 3*p1z + 3*p2z - p3z))*((p0z - 2*p1z + p2z)*(p0y - 3*p1y + 3*p2y - p3y) - (p0y - 2*p1y + p2y)
        *(p0z - 3*p1z + 3*p2z - p3z));
    double c_1 = ((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y - 2*p1y + p2y)*
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) + ((p0z/2 - p2z/2)*
        (p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0z - 2*p1z + p2z))*((p0z - 2*p1z + p2z)* 
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0z - 3*p1z + 3*p2z - p3z)) + ((p0z/2 - p2z/2)*
        (p0y - 2*p1y + p2y) - (p0y/2 - p2y/2)*(p0z - 2*p1z + p2z))*((p0z - 2*p1z + p2z)* 
        (p0y - 3*p1y + 3*p2y - p3y) - (p0y - 2*p1y + p2y)*(p0z - 3*p1z + 3*p2z - p3z)) + 
        pow((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y),2) + 
        pow((p0z/2 - p2z/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0z - 3*p1z + 3*p2z - p3z),2) + 
        pow((p0z/2 - p2z/2)*(p0y - 3*p1y + 3*p2y - p3y) - (p0y/2 - p2y/2)*(p0z - 3*p1z + 3*p2z - p3z),2);
    double c_0 = -((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y/2 - p2y/2)* 
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y)) - ((p0z/2 - p2z/2)* 
        (p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0z - 2*p1z + p2z))*((p0z/2 - p2z/2)* 
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0z - 3*p1z + 3*p2z - p3z)) - 
        ((p0z/2 - p2z/2)*(p0y - 2*p1y + p2y) - (p0y/2 - p2y/2)*(p0z - 2*p1z + p2z))*((p0z/2 - p2z/2)* 
        (p0y - 3*p1y + 3*p2y - p3y) - (p0y/2 - p2y/2)*(p0z - 3*p1z + 3*p2z - p3z));
    Eigen::Vector4d coeficients;
    coeficients << c_0, c_1, c_2, c_3;
    return coeficients;
}

template <int D>
double ThirdOrderCurvatureEvaluator<D>::calculate_cross_term_magnitude(double &t, Eigen::Matrix<double,D,4> &control_points, double &scale_factor)
{
    DerivativeEvaluator<D> d_dt_eval{};
    Eigen::Matrix<double,D,1> velocity_vector = d_dt_eval.calculate_velocity_vector(t, control_points, scale_factor);
    Eigen::Matrix<double,D,1> acceleration_vector = d_dt_eval.calculate_acceleration_vector(t, control_points, scale_factor);
    double cross_term_magnitude;
    if (D == 2)
    { 
        cross_term_magnitude = abs(velocity_vector(0)*acceleration_vector(1) - velocity_vector(1)*acceleration_vector(0));
    }
    else
    {
        float x = velocity_vector(1)*acceleration_vector(2) - velocity_vector(2)*acceleration_vector(1);
        float y = velocity_vector(2)*acceleration_vector(0) - velocity_vector(0)*acceleration_vector(2);
        float z = velocity_vector(0)*acceleration_vector(1) - velocity_vector(1)*acceleration_vector(0);
        cross_term_magnitude = sqrt(x*x + y*y + z*z);
    }
    return cross_term_magnitude;
}

//Explicit template instantiations
template class ThirdOrderCurvatureEvaluator<2>;
template class ThirdOrderCurvatureEvaluator<3>;