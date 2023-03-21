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
float ThirdOrderCurvatureEvaluator<D>::find_spline_curvature_bound(float cont_pts[], int num_control_points)
{
    float max_curvature{0};
    float curvature;
    for (unsigned int i = 0; i < num_control_points-3; i++)
    {
        Eigen::Matrix<float,D,4> interval_control_points = cbind_help.array_section_to_eigen(cont_pts, num_control_points, i);
        curvature = evaluate_interval_curvature(interval_control_points);
        if (curvature > max_curvature)
        {
            max_curvature = curvature;
        }
    }
    return max_curvature;
}

template <int D>
float ThirdOrderCurvatureEvaluator<D>::evaluate_interval_curvature(Eigen::Matrix<float,D,4> &control_points)
{
    float scale_factor = 1;
    DerivativeEvaluator<D> d_dt_eval{};
    std::array<float,2> min_velocity_and_time = d_dt_eval.find_min_velocity_and_time(control_points,scale_factor);
    float min_velocity = min_velocity_and_time[0];
    float time_at_min_velocity = min_velocity_and_time[1];
    float max_cross_term = find_maximum_cross_term(control_points,scale_factor);
    std::array<float,2> max_acceleration_and_time = d_dt_eval.find_max_acceleration_and_time(control_points, scale_factor);
    float max_acceleration = max_acceleration_and_time[0];
    float curvature_bound;
    if (min_velocity <= 1.0e-6)
    {
        float acceleration_at_min_vel = 
            d_dt_eval.calculate_acceleration_magnitude(time_at_min_velocity,
                control_points, scale_factor);
        if(acceleration_at_min_vel == 0)
        {
            curvature_bound = 0;
        }
        else
        {
            curvature_bound = std::numeric_limits<float>::max();
        }
    }
    else
    {
        curvature_bound = max_acceleration/(min_velocity*min_velocity);
        float curvature = max_cross_term/(min_velocity*min_velocity*min_velocity);
        if (curvature < curvature_bound)
        {
            curvature_bound = curvature;
        }
    }
    return curvature_bound;
}

template <int D>
float ThirdOrderCurvatureEvaluator<D>::find_maximum_cross_term(Eigen::Matrix<float,D,4> &control_points, float &scale_factor)
{
    Eigen::Vector4d coeficients;
    if (D == 2) {coeficients = get_2D_cross_coefficients(control_points);}
    else {coeficients = get_3D_cross_coefficients(control_points);}
    float a_term = coeficients(0);
    float b_term = coeficients(1);
    float c_term = coeficients(2);
    float d_term = coeficients(3);
    std::array<float,3> roots = CubicEquationSolver::solve_equation(a_term,
        b_term, c_term, d_term);
    float t0 = 0;
    float tf = 1.0;
    float max_cross_term = calculate_cross_term_magnitude(t0,control_points,scale_factor);
    float cross_term_at_tf = calculate_cross_term_magnitude(tf,control_points,scale_factor);
    if (cross_term_at_tf > max_cross_term)
    {
        max_cross_term = cross_term_at_tf;
    }
    for(int index = 0; index < 3; index++)
    {
        float root = roots[index];
        if(root > 0 && root < 1.0)
        {
            float cross_term =  calculate_cross_term_magnitude(root, control_points,scale_factor);
            if (cross_term > max_cross_term)
            {
                max_cross_term = cross_term;
            }
        }
    }
    return max_cross_term;
}

template <int D>
Eigen::Vector4d ThirdOrderCurvatureEvaluator<D>::get_2D_cross_coefficients(Eigen::Matrix<float,D,4> &control_points)
{
    float p0x = control_points(0,0);
    float p0y = control_points(1,0);
    float p1x = control_points(0,1);
    float p1y = control_points(1,1);
    float p2x = control_points(0,2);
    float p2y = control_points(1,2);
    float p3x = control_points(0,3);
    float p3y = control_points(1,3);
    float c_3 = ((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) 
        *((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - 
        (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2);
    float c_2 = - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y))*
        ((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - 
        (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2) - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*
        (p0y - 3*p1y + 3*p2y - p3y))*((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)* 
        (p0y - 3*p1y + 3*p2y - p3y));
    float c_1 = ((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y - 2*p1y + p2y)*
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) + 
        pow((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y),2);
    float c_0 = -((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y/2 - p2y/2)* 
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y));
    Eigen::Vector4d coeficients;
    coeficients << c_0, c_1, c_2, c_3;
    return coeficients;
}

template <int D>
Eigen::Vector4d ThirdOrderCurvatureEvaluator<D>::get_3D_cross_coefficients(Eigen::Matrix<float,D,4> &control_points)
{
    float p0x = control_points(0,0);
    float p0y = control_points(1,0);
    float p0z = control_points(2,0);
    float p1x = control_points(0,1);
    float p1y = control_points(1,1);
    float p1z = control_points(2,1);
    float p2x = control_points(0,2);
    float p2y = control_points(1,2);
    float p2z = control_points(2,2);
    float p3x = control_points(0,3);
    float p3y = control_points(1,3);
    float p3z = control_points(2,3);
    float c_3 = ((p0y - 2*p1y + p2y)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y))
        *((p0x*p1y)/2 - (p1x*p0y)/2 - p0x*p2y + p2x*p0y + (p0x*p3y)/2 + (3*p1x*p2y)/2 - (3*p2x*p1y)/2 - 
        (p3x*p0y)/2 - p1x*p3y + p3x*p1y + (p2x*p3y)/2 - (p3x*p2y)/2) + ((p0z - 2*p1z + p2z)*
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0z - 3*p1z + 3*p2z - p3z))* 
        ((p0x*p1z)/2 - (p1x*p0z)/2 - p0x*p2z + p2x*p0z + (p0x*p3z)/2 + (3*p1x*p2z)/2 - (3*p2x*p1z)/2 - 
        (p3x*p0z)/2 - p1x*p3z + p3x*p1z + (p2x*p3z)/2 - (p3x*p2z)/2) + ((p0z - 2*p1z + p2z)* 
        (p0y - 3*p1y + 3*p2y - p3y) - (p0y - 2*p1y + p2y)*(p0z - 3*p1z + 3*p2z - p3z))*((p0y*p1z)/2 - 
        (p1y*p0z)/2 - p0y*p2z + p2y*p0z + (p0y*p3z)/2 + (3*p1y*p2z)/2 - (3*p2y*p1z)/2 - (p3y*p0z)/2 - 
        p1y*p3z + p3y*p1z + (p2y*p3z)/2 - (p3y*p2z)/2);
    float c_2 = - ((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y))*
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
    float c_1 = ((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y - 2*p1y + p2y)*
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0y - 3*p1y + 3*p2y - p3y)) + ((p0z/2 - p2z/2)*
        (p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0z - 2*p1z + p2z))*((p0z - 2*p1z + p2z)* 
        (p0x - 3*p1x + 3*p2x - p3x) - (p0x - 2*p1x + p2x)*(p0z - 3*p1z + 3*p2z - p3z)) + ((p0z/2 - p2z/2)*
        (p0y - 2*p1y + p2y) - (p0y/2 - p2y/2)*(p0z - 2*p1z + p2z))*((p0z - 2*p1z + p2z)* 
        (p0y - 3*p1y + 3*p2y - p3y) - (p0y - 2*p1y + p2y)*(p0z - 3*p1z + 3*p2z - p3z)) + 
        pow((p0y/2 - p2y/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0y - 3*p1y + 3*p2y - p3y),2) + 
        pow((p0z/2 - p2z/2)*(p0x - 3*p1x + 3*p2x - p3x) - (p0x/2 - p2x/2)*(p0z - 3*p1z + 3*p2z - p3z),2) + 
        pow((p0z/2 - p2z/2)*(p0y - 3*p1y + 3*p2y - p3y) - (p0y/2 - p2y/2)*(p0z - 3*p1z + 3*p2z - p3z),2);
    float c_0 = -((p0y/2 - p2y/2)*(p0x - 2*p1x + p2x) - (p0x/2 - p2x/2)*(p0y - 2*p1y + p2y))*((p0y/2 - p2y/2)* 
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
float ThirdOrderCurvatureEvaluator<D>::calculate_cross_term_magnitude(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor)
{
    DerivativeEvaluator<D> d_dt_eval{};
    Eigen::Matrix<float,D,1> velocity_vector = d_dt_eval.calculate_velocity_vector(t, control_points, scale_factor);
    Eigen::Matrix<float,D,1> acceleration_vector = d_dt_eval.calculate_acceleration_vector(t, control_points, scale_factor);
    float cross_term_magnitude;
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