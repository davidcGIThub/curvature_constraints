import numpy as np
from scipy.spatial import ConvexHull
from bsplinegenerator.helper_functions import count_number_of_control_points
from scipy.spatial import ConvexHull, convex_hull_plot_2d
import matplotlib.pyplot as plt


class MDM(object):
    __class__ = 'MDM'
    __doc__ = """
        This is an implementation the accelerated Mitchell-Demyanov-Malozemov 
        method for finding nearest to coordinates beginning point. Also plots
        convex-hull and optimal solution in 2- and 3-dimensional cases.
        
        #############################################################
        ######################   ARGUMENTS    ######################
        #############################################################
        _dim: int 
            Dimensionality of the problem, usually dim = 2 or 3 for plotting purposes
            
        _points: np.ndarray of shape (dim, number_of_points)
            Points which form convex hull
            
        is_accelerated: boolean
            Flag for using accelerated either not-accelerated MDM-method
            
        supp_vector: list of floats
            Support vector (i.e. {i \in 0 : dim - 1 | p[i] > 0} )
            
        max_iter: int (default value is 500)
            Maximum iterations we can evaluate on

        init_approx_index: int (default is 1)
            Index of element of convex hull, which forms first approximation 
            of v vector. It can be changed for lowering iterations sake
            For first approximation we'll just take point from a board of 
            hull - cause it's easily reduced
        #############################################################
        #############################################################
        #############################################################
    """

    def __init__(self, points, dim, max_iter = 500, init_approx_index = 1):
        self._dim = dim
        self._num_points = count_number_of_control_points(points)
        self._points = np.transpose(points).copy()
        self._A_matrix = points.copy()
        self.max_iter = max_iter
        self.init_approx_index = init_approx_index

    def get_min_distance(self):
        closest_point = self.get_closest_point()
        min_distance = np.linalg.norm(closest_point)
        # print("")
        # print("get min distance()")
        # print("closest_point: " , closest_point)
        # print("min_distance: " , min_distance)
        # print("")
        return min_distance

    def get_closest_point(self):
        # print("")
        if self._num_points == 1:
            closest_point = self._points[0,:]
            return closest_point
        elif self._num_points == 2:
            closest_point = self.solve_two_points()
            return closest_point
        else:
            return self.solve_three_or_more_points()

    def solve_two_points(self):
        point1 = self._points[0,:]
        point2 = self._points[1,:]
        if np.array_equal(point1,point2):
            closest_point = point1
            return closest_point
        point_difference = point1 - point2
        alpha = np.dot(point1,point_difference)/np.dot(point_difference,point_difference)
        if alpha > 1:
            closest_point = point2
        elif alpha < 0:
            closest_point = point1
        else:
            closest_point = alpha*point2 + (1-alpha)*point1
        return closest_point

    def solve_three_or_more_points(self):
        # print("")
        # print("Solve for three points or more()")
        # print("")
        # print("points: " , self._points)
        # print("max_iterations: " , self.max_iter)
        # print("initial index: " , self.init_approx_index)
        iterations = 0
        # print("iterations : ", iterations)
        delta_p = 1
        # print("delta_p : ", delta_p)
        p_vector = [0 for i in range(0, len(self._points))]
        # print("p_vector: ", p_vector)
        supp_vector = []
        # print("supp_vector : ", supp_vector)
        index = 0
        # print("index : ", index)
        vector_current = self._points[index,:].copy()
        # print("vector_current: ", vector_current)
        supp_vector.append(index)
        # print("supp_vector : ", supp_vector)
        p_vector[index] = 1
        # print("p_vector : ", p_vector)
        #then we need to find vect_{k+1} iteratively
        while delta_p > 0.000001 and iterations < self.max_iter and len(supp_vector) != 0:
            mat = self._points[supp_vector]
            # print("mat: " , mat)
            mult = np.dot(mat, vector_current)
            # print("mult: " , mult)
            ind_max = np.argmax(mult)           # finding max for indices in supp_vector
            # print("ind_max : ", ind_max)
            ind_max = supp_vector[ind_max]      # finding max general in our mult product
            # print("ind_max : ", ind_max)
            mult = np.matmul(vector_current, self._A_matrix)
            # print("mult : ", mult)
            ind_min = np.argmin(mult)   # i''_k
            # print("ind_min : ", ind_min)
            diff = self._points[ind_max] - self._points[ind_min]
            # print("diff : ", diff)
            delta_p = np.dot(diff, vector_current)
            # print("delta_p : ", delta_p)
            if delta_p > 0.000001:                  #if not bigger, then we've found a solution
                t_param = delta_p /(p_vector[ind_max] * (np.linalg.norm(diff)) ** 2)  # recounting all variables
                # print("t_param : ", t_param)
                if t_param >= 1:
                    t_param = 1
                # print("vector_current: " , vector_current)
                # print("t_param : ", t_param)
                # print("p_vector[ind_max] : ", p_vector[ind_max])
                # print("diff: " , diff)
                # print("t_param * p_vector[ind_max] : ", t_param * p_vector[ind_max])
                # print("t_param * p_vector[ind_max] * diff: " , t_param * p_vector[ind_max] * diff)
                vector_current -= t_param * p_vector[ind_max] * diff
                # print("vector_current : ", vector_current)
                supp_vector = [] # recounting
                # print("supp_vector : ", supp_vector)
                temp1 = t_param * p_vector[ind_max]
                # print("temp1 : ", temp1)
                temp2 = (1 - t_param)
                # print("temp2 : ", temp2)
                p_vector[ind_min] += temp1
                # print("p_vector : ", p_vector)
                p_vector[ind_max] *= temp2
                # print("p_vector : ",p_vector )
                for i in range(len(p_vector)):
                    if p_vector[i] > 0.0000001:
                        supp_vector.append(i)
                        # print("supp_vector : ", supp_vector)
            iterations += 1
            # print("iterations : ", iterations)
            closest_point = vector_current
            # print("closest_point : ", closest_point)
        # print("")
        return closest_point

    def plot_convex_hull_with_solution(self):
        hull = ConvexHull(self._points)
        plt.plot(self._points[:,0], self._points[:,1], 'o')
        for simplex in hull.simplices:
            plt.plot(self._points[simplex, 0], self._points[simplex, 1], 'k-')
        closest_point = self.get_closest_point()
        plt.scatter(closest_point[0],closest_point[1],color="tab:red")
        min_distance = np.linalg.norm(closest_point)
        # print("min distance: " , min_distance)
        plt.show()

        