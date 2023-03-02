#include "CurvatureBounds3rdOrder.hpp"
#include <stdexcept>

template <int D>
CurvatureBounds3rdOrder<D>::CurvatureBounds3rdOrder(){}

template <int D>
double CurvatureBounds3rdOrder<D>::find_spline_curvature_bound(double cont_pts[], int &num_control_points)
{   
    double max_curvature{0};
    double curvature;
    for (int i = 0; i < num_control_points-3; i++)
    {
        int step = num_control_points;
        Eigen::Matrix<double,D,4> interval_control_points;
        if (D == 2)
        {
            interval_control_points << cont_pts[i], cont_pts[i+1], cont_pts[i+2], cont_pts[i+3],
                                    cont_pts[step+i], cont_pts[step+i+1], cont_pts[step+i+2], cont_pts[step+i+3];
        }
        else // D == 3
        {
            interval_control_points << cont_pts[i], cont_pts[i+1], cont_pts[i+2], cont_pts[i+3],
                                    cont_pts[step+i], cont_pts[step+i+1], cont_pts[step+i+2], cont_pts[step+i+3];
                                    cont_pts[2*step+i], cont_pts[2*step+i+1], cont_pts[2*step+i+2], cont_pts[2*step+i+3];
        }
        curvature = evaluate_interval_curvature(interval_control_points);
        if (curvature > max_curvature)
        {
            max_curvature = curvature;
        }
    }
    return max_curvature;
}

template <int D>
double CurvatureBounds3rdOrder<D>::find_min_velocity_of_spline(double cont_pts[], int &num_control_points, double &scale_factor)
{   
    double min_velocity = std::numeric_limits<double>::max();
    double velocity;
    for (int i = 0; i < num_control_points-3; i++)
    {
        int step = num_control_points;
        Eigen::Matrix<double,2,4> interval_control_points; 
        if (D == 2)
        {
            interval_control_points << cont_pts[i], cont_pts[i+1], cont_pts[i+2], cont_pts[i+3],
                                    cont_pts[step+i], cont_pts[step+i+1], cont_pts[step+i+2], cont_pts[step+i+3];
        }
        else // D == 3
        {
            interval_control_points << cont_pts[i], cont_pts[i+1], cont_pts[i+2], cont_pts[i+3],
                                    cont_pts[step+i], cont_pts[step+i+1], cont_pts[step+i+2], cont_pts[step+i+3];
                                    cont_pts[2*step+i], cont_pts[2*step+i+1], cont_pts[2*step+i+2], cont_pts[2*step+i+3];
        }
        std::array<double,2> vel_and_time = ThirdOrderCurvatureEvaluator::find_minimum_velocity_and_time(interval_control_points, scale_factor);
        velocity = vel_and_time[0];
        if (velocity < min_velocity)
        {
            min_velocity = velocity;
        }
    }
    return min_velocity;
}

template <int D>
double CurvatureBounds3rdOrder<D>::evaluate_interval_curvature(Eigen::Matrix<double,D,4> &control_points)
{
    double scale_factor = 1;
    std::array<double,2> min_velocity_and_time = ThirdOrderCurvatureEvaluator::find_minimum_velocity_and_time(control_points,scale_factor);
    double min_velocity = min_velocity_and_time[0];
    double time_at_min_velocity = min_velocity_and_time[1];
    double max_cross_term = ThirdOrderCurvatureEvaluator::find_maximum_cross_term(control_points,scale_factor);
    double max_acceleration = ThirdOrderCurvatureEvaluator::find_maximum_acceleration(control_points);
    double curvature_bound;
    if (min_velocity <= 1.0e-15)
    {
        double acceleration_at_min_vel = 
            ThirdOrderCurvatureEvaluator::calculate_acceleration_magnitude(time_at_min_velocity,
                control_points);
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
std::array<double,2> CurvatureBounds3rdOrder<D>::find_minimum_velocity_and_time(Eigen::Matrix<double,D,4> &control_points, double &scale_factor)
{
    Eigen::Matrix<double, 4,4> M = get_third_order_M_matrix();
    Eigen::Matrix<double, 4,4> J = M.transpose()*control_points.transpose()*
        control_points*M;
    double A = 36*J(0,0);
    double B = 12*J(0,1) + 24*J(1,0);
    double C = 8*J(1,1) + 12*J(2,0);
    double D = 4*J(2,1);
    std::array<double,3> roots = CubicEquationSolver::solve_equation(A, B, C, D);
    double t0 = 0;
    double tf = 1.0;
    double time_at_min = t0;
    double min_velocity = ThirdOrderCurvatureEvaluator::calculate_velocity_magnitude(t0,control_points,scale_factor);
    double velocity_at_tf = ThirdOrderCurvatureEvaluator::calculate_velocity_magnitude(tf,control_points,scale_factor);
    if (velocity_at_tf < min_velocity)
    {
        min_velocity = velocity_at_tf;
        double time_at_min = tf;
    }
    for(int index = 0; index < 3; index++)
    {
        double root = roots[index];
        if(root > 0 && root < 1.0)
        {
            double velocity =  ThirdOrderCurvatureEvaluator::calculate_velocity_magnitude(root, control_points,scale_factor);
            if (velocity < min_velocity)
            {
                min_velocity = velocity;
                double time_at_min = root;
            }
        }
    }
    std::array<double,2> min_velocity_and_time = {min_velocity, time_at_min};
    return min_velocity_and_time;
}

template <int D>
double CurvatureBounds3rdOrder<D>::find_maximum_acceleration(Eigen::Matrix<double,D,4> &control_points)
{
    double t0 = 0;
    double tf = 1.0;
    double max_acceleration = ThirdOrderCurvatureEvaluator::calculate_acceleration_magnitude(t0,control_points);
    double acceleration_at_tf = ThirdOrderCurvatureEvaluator::calculate_acceleration_magnitude(tf, control_points);
    if (acceleration_at_tf > max_acceleration)
    {
        max_acceleration = acceleration_at_tf;
    }
    return max_acceleration;
}

template <int D>
double CurvatureBounds3rdOrder<D>::find_maximum_cross_term(Eigen::Matrix<double,D,4> &control_points, double &scale_factor)
{
    Eigen::Vector4d coeficients = get_2D_cross_coefficients(control_points);
    std::array<double,3> roots = CubicEquationSolver::solve_equation(coeficients(0),
        coeficients(1), coeficients(2), coeficients(3));
    double t0 = 0;
    double tf = 1.0;
    double max_cross_term = ThirdOrderCurvatureEvaluator::calculate_cross_term_magnitude(t0,control_points,scale_factor);
    double cross_term_at_tf = ThirdOrderCurvatureEvaluator::calculate_cross_term_magnitude(tf,control_points,scale_factor);
    if (cross_term_at_tf > max_cross_term)
    {
        max_cross_term = cross_term_at_tf;
    }
    for(int index = 0; index < 3; index++)
    {
        double root = roots[index];
        if(root > 0 && root < 1.0)
        {
            double cross_term =  ThirdOrderCurvatureEvaluator::calculate_cross_term_magnitude(root, control_points,scale_factor);
            if (cross_term > max_cross_term)
            {
                max_cross_term = cross_term;
            }
        }
    }
    return max_cross_term;
}

template <int D>
Eigen::Vector4d CurvatureBounds3rdOrder<D>::get_2D_cross_coefficients(Eigen::Matrix<double,2,4> &control_points)
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
Eigen::Vector4d CurvatureBounds3rdOrder<D>::get_3D_cross_coefficients(Eigen::Matrix<double,3,4> &control_points)
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
double CurvatureBounds3rdOrder<D>::calculate_velocity_magnitude(double &t, Eigen::Matrix<double,D,4> &control_points, double &scale_factor)
{
    Eigen::Vector4d dT = ThirdOrderCurvatureEvaluator::get_third_order_T_derivative_vector(t, scale_factor);
    Eigen::Matrix<double, 4,4> M = get_third_order_M_matrix();
    Eigen::Vector2d velocity_vector = control_points*M*dT;
    double velocity_magnitude = velocity_vector.norm();
    return velocity_magnitude;
}

template <int D>
double CurvatureBounds3rdOrder<D>::calculate_acceleration_magnitude(double &t, Eigen::Matrix<double,D,4> &control_points)
{
    Eigen::Vector4d ddT = ThirdOrderCurvatureEvaluator::get_third_order_T_second_derivative_vector(t);
    Eigen::Matrix<double, 4,4> M = get_third_order_M_matrix();
    Eigen::Vector2d acceleration_vector = control_points*M*ddT;
    double acceleration_magnitude = acceleration_vector.norm();
    return acceleration_magnitude;
}

template <int D>
double CurvatureBounds3rdOrder<D>::calculate_cross_term_magnitude(double &t, Eigen::Matrix<double,D,4> &control_points, double &scale_factor)
{
    Eigen::Matrix<double, 4,4> M = get_third_order_M_matrix();
    Eigen::Vector4d ddT = ThirdOrderCurvatureEvaluator::get_third_order_T_second_derivative_vector(t);
    Eigen::Vector2d acceleration_vector = control_points*M*ddT;
    Eigen::Vector4d dT = ThirdOrderCurvatureEvaluator::get_third_order_T_derivative_vector(t,scale_factor);
    Eigen::Vector2d velocity_vector = control_points*M*dT;
    double cross_term_magnitude = velocity_vector(0)*acceleration_vector(1) 
        - velocity_vector(1)*acceleration_vector(0);
    return abs(cross_term_magnitude);
}

template <int D>
Eigen::Matrix<double, 4,4> CurvatureBounds3rdOrder<D>::get_third_order_M_matrix()
{
    Eigen::Matrix<double, 4,4> M;
    M <<  -1/6.0 ,   1/2.0 , -1/2.0 , 1/6.0,
            1/2.0 ,      -1 ,  0     , 2/3.0,
            -1/2.0 ,   1/2.0 ,  1/2.0 , 1/6.0,
            1/6.0 ,  0      ,  0     , 0;
    return M;
}

template <int D>
Eigen::Vector4d CurvatureBounds3rdOrder<D>::get_third_order_T_derivative_vector(double &t, double &scale_factor)
{
    Eigen::Vector4d t_vector;
    if (t < 0 || t > 1)
    {
        throw std::invalid_argument("t value should be between 0 and 1");
    }
    else
    {
        t_vector << 3*t*t , 2*t , 1 , 0;
    }
    return t_vector;
}

template <int D>
Eigen::Vector4d CurvatureBounds3rdOrder<D>::get_third_order_T_second_derivative_vector(double &t)
{
    Eigen::Vector4d t_vector;
    if (t < 0 || t > 1)
    {
        throw std::invalid_argument("t value should be between 0 and 1");
    }
    else
    {
        t_vector << 6*t , 2 , 0 , 0;
    }
    return t_vector;
}