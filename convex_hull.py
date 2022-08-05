import numpy as np 
import matplotlib.pyplot as plt

# Mitchell–Dem'yanov–Malozemov Algorithm

points = np.array([[2,5,6],[5,3,5]])
# points = np.array([[-1,4,7],[-3,6,6]])
points = np.random.randint(10, size=(2,5))
print("points: " , points)

def find_minimum_norm_of_convex_hull(points):
    improvement = 100
    num_points = np.shape(points)[1]
    alpha = np.ones(num_points)/num_points
    u = np.dot(points,alpha)
    u_prev = u
    # print("u: ", u)
    while improvement > 0.01:
        plt.plot(np.array([0,u[0]]),[0,u[1]])
        plt.scatter(points[0,:],points[1,:],label="points")
        plt.show()
        points_unit = points / np.linalg.norm(points,2,0)
        dot_product_of_points = np.dot(np.transpose(points),u)
        # print("dot_product_of_points: " , np.shape(dot_product_of_points))
        point_most_parallel_index = np.argmax(dot_product_of_points)
        # print("point_most_parallel_index: " , point_most_parallel_index)
        point_most_perpindicular_index = np.argmin(dot_product_of_points)
        # print("point_most_perpindicular_index: " , point_most_perpindicular_index) 
        point_tilde = points[:,point_most_parallel_index]
        print("point_tilde: " , point_tilde)
        point_bar = points[:,point_most_perpindicular_index]
        print("point_bar: " , point_bar)
        alpha_tilde = alpha[point_most_parallel_index]
        # print("alpha_tilde: " , alpha_tilde)
        term = (point_bar-point_tilde)*alpha_tilde
        # print("term: " , term)
        s_extrema = -np.dot(u,term)/np.dot(u,u)
        phi_values = np.zeros(100)
        s_values = np.linspace(0,1,100)
        for i in range(100):
            s = s_values[i]
            phi_values[i] = np.dot(u+s*term,u+s*term)
        plt.plot(s_values,phi_values)
        plt.show()
        # print("s_extrema: " , s_extrema)
        start_s_value = np.dot(u,u)
        # print("start_s_value: " , start_s_value)
        end_s_value = np.dot(u+term,u+term)
        # print("end_s_value: " , end_s_value)
        mid_s_value = np.dot(u+0.5*term,u+0.5*term)
        # print("mid_s_value: " , mid_s_value)
        s_extrema_value = np.dot(u+s_extrema*term,u+s_extrema*term)
        # print("s_extrema_value: " , s_extrema_value)
        min_index = np.argmin((start_s_value,end_s_value,s_extrema_value))
        # print("min_index: " , min_index)
        s_k  = np.array([0,1,s_extrema])[min_index]
        # print("s_k: " , s_k)
        u_prev = u
        # print("u_prev: " , u_prev)
        u = u_prev + s_k*term
        # print("u: " , u)
        A = np.concatenate((np.ones((1,num_points)),points),0)
        # print("A: " , A)
        b = np.concatenate((np.array([1]),u),0)
        # print("b: " , b)
        alpha = np.dot(np.linalg.pinv(A),b)
        # print("alpha: " , alpha)
        improvement = np.dot(u_prev,u_prev) - np.dot(u,u)
        print("### improvement #### : " , improvement)
    return u

closest_point = find_minimum_norm_of_convex_hull(points)
plt.plot(np.array([0,closest_point[0]]),[0,closest_point[1]])
plt.scatter(points[0,:],points[1,:],label="points")
plt.show()

