#include "DerivativeEvaluator.hpp"
#include "CubicEquationSolver.hpp"
#include <iostream>
#include <stdexcept>

template <int D>
DerivativeEvaluator<D>::DerivativeEvaluator()
{

}

template <int D>
float DerivativeEvaluator<D>::find_min_velocity_of_spline(float cont_pts[], int num_control_points, float scale_factor)
{   
    float min_velocity = std::numeric_limits<float>::max();
    float velocity;
    int step = num_control_points;
    for (unsigned int i = 0; i < num_control_points-3; i++)
    {
        Eigen::Matrix<float,D,4> interval_control_points = cbind_help.array_section_to_eigen(cont_pts, num_control_points, i);
        std::array<float,2> vel_and_time = find_min_velocity_and_time(interval_control_points, scale_factor);
        velocity = vel_and_time[0];
        if (velocity < min_velocity)
        {
            min_velocity = velocity;
        }
    }
    return min_velocity;
}

template <int D>
float DerivativeEvaluator<D>::find_max_acceleration_of_spline(float cont_pts[], int num_control_points, float scale_factor)
{   
    float max_acceleration = std::numeric_limits<float>::min();
    float acceleration;
    int step = num_control_points;
    for (unsigned int i = 0; i < num_control_points-3; i++)
    {
        Eigen::Matrix<float,D,4> interval_control_points = cbind_help.array_section_to_eigen(cont_pts, num_control_points, i);
        std::array<float,2> accel_and_time = find_max_acceleration_and_time(interval_control_points, scale_factor);
        acceleration = accel_and_time[0];
        if (acceleration > max_acceleration)
        {
            max_acceleration = acceleration;
        }
    }
    return max_acceleration;
}

template <int D>
std::array<float,2> DerivativeEvaluator<D>::find_min_velocity_and_time(Eigen::Matrix<float,D,4> &control_points, float &scale_factor)
{
    Eigen::Matrix<float, 4,4> M = get_third_order_M_matrix();
    Eigen::Matrix<float, 4,4> J = M.transpose()*control_points.transpose()*control_points*M;
    float A = 36*J(0,0);
    float B = 12*J(0,1) + 24*J(1,0);
    float C = 8*J(1,1) + 12*J(2,0);
    float D_ = 4*J(2,1);
    std::array<float,3> roots = CubicEquationSolver::solve_equation(A, B, C, D_);
    float t0 = 0;
    float tf = 1.0;
    float time_at_min = t0;
    float min_velocity = calculate_velocity_magnitude(t0,control_points,scale_factor);
    float velocity_at_tf = calculate_velocity_magnitude(tf,control_points,scale_factor);
    if (velocity_at_tf < min_velocity)
    {
        min_velocity = velocity_at_tf;
        float time_at_min = tf;
    }
    for(int index = 0; index < 3; index++)
    {
        float root = roots[index];
        if(root > 0 && root < 1.0)
        {
            float velocity =  calculate_velocity_magnitude(root, control_points,scale_factor);
            if (velocity < min_velocity)
            {
                min_velocity = velocity;
                float time_at_min = root;
            }
        }
    }
    std::array<float,2> min_velocity_and_time = {min_velocity, time_at_min};
    return min_velocity_and_time;
}

template <int D>
std::array<float,2> DerivativeEvaluator<D>::find_max_acceleration_and_time(Eigen::Matrix<float,D,4> &control_points, float &scale_factor)
{
    float t0 = 0;
    float tf = 1.0;
    float time_at_max = 0;
    float max_acceleration = calculate_acceleration_magnitude(t0,control_points,scale_factor);
    float acceleration_at_tf = calculate_acceleration_magnitude(tf, control_points,scale_factor);
    if (acceleration_at_tf > max_acceleration)
    {
        max_acceleration = acceleration_at_tf;
        time_at_max = tf;
    }
    std::array<float,2> max_acceleration_and_time = {max_acceleration, time_at_max};
    return max_acceleration_and_time;
}

template <int D>
Eigen::Matrix<float,D,1> DerivativeEvaluator<D>::calculate_velocity_vector(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor)
{
    Eigen::Vector4f dT = get_third_order_T_derivative_vector(t, scale_factor);
    Eigen::Matrix<float,4,4> M = get_third_order_M_matrix();
    Eigen::Matrix<float,D,1> velocity_vector = control_points*M*dT;
    return velocity_vector;
}

template <int D>
float DerivativeEvaluator<D>::calculate_velocity_magnitude(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor)
{
    Eigen::Matrix<float,D,1> velocity_vector = calculate_velocity_vector(t, control_points, scale_factor);
    float velocity_magnitude = velocity_vector.norm();
    return velocity_magnitude;
}

template <int D>
Eigen::Matrix<float,D,1> DerivativeEvaluator<D>::calculate_acceleration_vector(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor)
{
    Eigen::Vector4f ddT = get_third_order_T_second_derivative_vector(t, scale_factor);
    Eigen::Matrix<float, 4,4> M = get_third_order_M_matrix();
    Eigen::Matrix<float,D,1> acceleration_vector = control_points*M*ddT;
    return acceleration_vector;
}

template <int D>
float DerivativeEvaluator<D>::calculate_acceleration_magnitude(float &t, Eigen::Matrix<float,D,4> &control_points, float &scale_factor)
{
    Eigen::Matrix<float,D,1> acceleration_vector = calculate_acceleration_vector(t, control_points, scale_factor);
    float acceleration_magnitude = acceleration_vector.norm();
    return acceleration_magnitude;
}

template <int D>
Eigen::Matrix<float, 4,4> DerivativeEvaluator<D>::get_third_order_M_matrix()
{
    Eigen::Matrix<float, 4,4> M;
    M <<  -1/6.0 ,   1/2.0 , -1/2.0 , 1/6.0,
            1/2.0 ,      -1 ,  0     , 2/3.0,
            -1/2.0 ,   1/2.0 ,  1/2.0 , 1/6.0,
            1/6.0 ,  0      ,  0     , 0;
    return M;
}

template <int D>
Eigen::Vector4f DerivativeEvaluator<D>::get_third_order_T_derivative_vector(float &t, float &scale_factor)
{
    float alpha = scale_factor;
    Eigen::Vector4f t_vector;
    if (scale_factor == 1.0)
    {
        if (t < 0 || t > 1)
        {
            throw std::invalid_argument("t value should be between 0 and 1");
        }
        else
        {
            t_vector << 3*t*t , 2*t , 1, 0;
        }
    }
    else
    {
        t_vector << 3*t*t/(alpha*alpha*alpha) , 2*t/(alpha*alpha) , 1/alpha , 0;
    }
    return t_vector;
}

template <int D>
Eigen::Vector4f DerivativeEvaluator<D>::get_third_order_T_second_derivative_vector(float &t,  float &scale_factor)
{
    float alpha = scale_factor;
    Eigen::Vector4f t_vector;
    if (scale_factor == 1)
    {
        if (t < 0 || t > 1)
        {
            throw std::invalid_argument("t value should be between 0 and 1");
        }
        else
        {
            t_vector << 6*t , 2 , 0 , 0;
        }
    }
    else
    {
        t_vector << 6*t/(alpha*alpha*alpha ) , 2/(alpha*alpha) , 0 , 0;
    }
    return t_vector;
}

//Explicit template instantiations
template class DerivativeEvaluator<2>;
template class DerivativeEvaluator<3>;