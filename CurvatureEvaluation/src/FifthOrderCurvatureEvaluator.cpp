#include "FifthOrderCurvatureEvaluator.hpp"
#include <iostream>

namespace FifthOrderCurvatureEvaluator
{
    double evaluate_interval_curvature_bound(Eigen::Matrix<double, 2, 6> &control_points)
    {
        double max_cross_term = get_cross_term_max_bound(control_points);
        double max_acceleration = get_acceleration_max_bound(control_points);
        double min_velocity = get_velocity_min_bound(control_points);
        double max_curvature{0};
        if (min_velocity <= 1e-10)
        {
            if(check_for_cusp(control_points))
            {
                max_curvature = std::numeric_limits<double>::max();
            }
            else
            {
                max_curvature = 0;
            }
        }
        else
        {
            max_curvature = max_acceleration/(min_velocity*min_velocity);
            double cross_method_bound = max_cross_term/(min_velocity*min_velocity*min_velocity);
            if (cross_method_bound < max_curvature)
            {
                max_curvature = cross_method_bound;
            }
        }
        return max_curvature;
    }

    double get_cross_term_max_bound(Eigen::Matrix<double, 2, 6> &control_points)
    {
        double scale_factor = 1.0;
        Eigen::Matrix<double, 2, 7> cross_term_control_points =
            get_fifth_order_cross_term_control_points(control_points);
        double max_cross_term = cross_term_control_points.colwise().norm().maxCoeff();
        return max_cross_term;
    }
    
    double get_acceleration_max_bound(Eigen::Matrix<double, 2, 6> &control_points)
    {
        double scale_factor = 1.0;
        Eigen::Matrix<double, 2, 5> velocity_control_points =
            (control_points.block(0,1,2,5) - control_points.block(0,0,2,5))/scale_factor;
        Eigen::Matrix<double, 2, 4> acceleration_control_points =
            (velocity_control_points.block(0,1,2,4) - velocity_control_points.block(0,0,2,4))
            /scale_factor;
        Eigen::Matrix<double, 2, 4> bezier_acceleration_control_points = 
            BsplineToBezier::convert_2D_third_order_spline(acceleration_control_points);
        double max_acceleration = bezier_acceleration_control_points.colwise().norm().maxCoeff();
        return max_acceleration;

    }

    double get_velocity_min_bound(Eigen::Matrix<double, 2, 6> &control_points)
    {
        double scale_factor = 1.0;
        Eigen::Matrix<double, 2, 5> velocity_control_points =
            (control_points.block(0,1,2,5) - control_points.block(0,0,2,5))/scale_factor;
        Eigen::Matrix<double, 2, 5> bezier_velocity_control_points = 
            BsplineToBezier::convert_2D_fourth_order_spline(velocity_control_points);
        MDMAlgorithmClass<5,2> mdm{};
        int max_iter = 500;
        unsigned int init_index = 1;
        double tolerance = 0.000001;
        double min_velocity= mdm.min_norm(bezier_velocity_control_points, max_iter, init_index, tolerance);
        return min_velocity;
    }

    Eigen::Matrix<double, 2,7> get_fifth_order_cross_term_control_points(Eigen::Matrix<double,2,6> &control_points)
    {
        Eigen::Matrix<double, 2, 1> P_0 = control_points.block(0,0,2,1);
        Eigen::Matrix<double, 2, 1> P_1 = control_points.block(0,1,2,1);
        Eigen::Matrix<double, 2, 1> P_2 = control_points.block(0,2,2,1);
        Eigen::Matrix<double, 2, 1> P_3 = control_points.block(0,3,2,1);
        Eigen::Matrix<double, 2, 1> P_4 = control_points.block(0,4,2,1);
        Eigen::Matrix<double, 2, 1> P_5 = control_points.block(0,5,2,1);
        Eigen::Matrix<double, 2, 1> p1;
        Eigen::Matrix<double, 2, 1> p2;
        Eigen::Matrix<double, 2, 1> p3;
        Eigen::Matrix<double, 2, 1> p4;
        Eigen::Matrix<double, 2, 1> p5;
        p1 = P_0/4 - P_1/2 + P_3/2 - P_4/4;
        p2 = P_0/24 - (5*P_1)/24 + (5*P_2)/12 - (5*P_3)/12 + (5*P_4)/24 - P_5/24;
        p3 = P_0/6 - (2*P_1)/3 + P_2 - (2*P_3)/3  + P_4/6;
        p4 = P_0/6  + P_1/3      - P_2        + P_3/3      + P_4/6;
        p5 = P_0/24 + (5*P_1)/12 - (5*P_3)/12 - P_4/24;
        Eigen::Matrix<double, 2,2> Y1;
        Eigen::Matrix<double, 2,2> Y2;
        Eigen::Matrix<double, 2,2> Y3;
        Eigen::Matrix<double, 2,2> Y4;
        Y1 << 1,0,0,0;
        Y4 << 1,0,0,0;
        Y2 << 0,1,0,0;
        Y3 << 0,1,0,0;
        Eigen::Matrix<double, 2,1> c_0;
        Eigen::Matrix<double, 2,1> c_1;
        Eigen::Matrix<double, 2,1> c_2;
        Eigen::Matrix<double, 2,1> c_3;
        Eigen::Matrix<double, 2,1> c_4;
        Eigen::Matrix<double, 2,1> c_5;
        Eigen::Matrix<double, 2,1> c_6;
        c_6 = (Y1*(3*p3)).cwiseProduct(Y2*p2) - (Y2*(3*p3)).cwiseProduct(Y1*p2) + \
            (Y3*p3).cwiseProduct(Y4*4*p2) - (Y4*p3).cwiseProduct(Y3*4*p2);
        c_5 = (Y1*p1).cwiseProduct(Y2*4*p2) - (Y2*p1).cwiseProduct(Y1*4*p2) + \
            (Y4*p3).cwiseProduct(Y3*3*p3) - (Y3*p3).cwiseProduct(Y4*3*p3) + \
            (Y3*2*p1).cwiseProduct(Y4*p2) - (Y4*2*p1).cwiseProduct(Y3*p2);
        c_4 = (Y2*p4).cwiseProduct(Y1*4*p2) - (Y1*p4).cwiseProduct(Y2*4*p2) - \
            (Y3*p4).cwiseProduct(Y4*p2)   + (Y4*p4).cwiseProduct(Y3*p2) + \
            (Y3*p1).cwiseProduct(Y4*3*p3) - (Y4*p1).cwiseProduct(Y3*3*p3) - \
            (Y3*2*p1).cwiseProduct(Y4*p3) + (Y4*2*p1).cwiseProduct(Y3*p3);
        c_3 = (Y4*p4)  .cwiseProduct(Y3*3*p3) - (Y4*3*p3).cwiseProduct(Y3*p4) - \
            (Y3*p5)  .cwiseProduct(Y4*4*p2) + (Y4*p5)  .cwiseProduct(Y3*4*p2) - \
            (Y4*2*p1).cwiseProduct(Y3*p1)   + (Y3*2*p1).cwiseProduct(Y4*p1) + \
            (Y4*p3)  .cwiseProduct(Y3*p4)   - (Y3*p3)  .cwiseProduct(Y4*p4);
        c_2 = (Y1*2*p1).cwiseProduct(Y2*p4) - (Y2*2*p1).cwiseProduct(Y1*p4) + \
            (Y3*p1).cwiseProduct(Y4*p4) - (Y4*p1).cwiseProduct(Y3*p4) + \
            (Y3*p5).cwiseProduct(Y4*3*p3) - (Y4*p5).cwiseProduct(Y3*3*p3);
        c_1 = (Y2*2*p1).cwiseProduct(Y1*p5) - (Y1*2*p1).cwiseProduct(Y2*p5);
        c_0 = (Y3*p5).cwiseProduct(Y4*p4) - (Y4*p5).cwiseProduct(Y3*p4);
        Eigen::Matrix<double, 2,1> a_0;
        Eigen::Matrix<double, 2,1> a_1;
        Eigen::Matrix<double, 2,1> a_2;
        Eigen::Matrix<double, 2,1> a_3;
        Eigen::Matrix<double, 2,1> a_4;
        Eigen::Matrix<double, 2,1> a_5;
        Eigen::Matrix<double, 2,1> a_6;
        a_0 = c_0;
        a_1 = (c_1 + 6*a_0)/6;
        a_2 = (c_2 - 15*a_0 + 30*a_1)/15;
        a_3 = (c_3 + 60*a_2 + 20*a_0 - 60*a_1)/20;
        a_4 = (c_4 - 15*a_0 + 60*a_1 - 90*a_2 + 60*a_3)/15;
        a_5 = (c_5 - 30*a_1 + 6*a_0  + 60*a_2 - 60*a_3 + 30*a_4)/6;
        a_6 =  c_6 - a_0    + 6*a_1  - 15*a_2 + 20*a_3 - 15*a_4 + 6*a_5;
        Eigen::Matrix<double, 2, 7> cross_term_control_points;
        cross_term_control_points << a_0, a_1, a_2, a_3, a_4, a_5, a_6;
        return cross_term_control_points;
    }

    bool check_for_cusp(Eigen::Matrix<double, 2, 6> &control_points)
    {
        Eigen::Matrix<double, 2, 5> difference =
                (control_points.block(0,1,2,5) - control_points.block(0,0,2,5));
        bool is_cusp{false};
        for(int i = 0; i < 2; i++)
        {
            if ((difference.block(0,0,i,5).array() > 0.0).any() && (difference.block(0,0,i,5).array() < 0.0).any())
            {
                is_cusp = true;
                break;
            }
        }
        return is_cusp;
    }
}