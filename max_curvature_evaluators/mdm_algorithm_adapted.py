import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

class MDM(object):
    __class__ = 'MDM'
    __doc__ = """
        This is an implementation the accelerated Mitchell-Demyanov-Malozemov method for 
        finding nearest to coordinates beginning point.
        Also plots convex-hull and optimal solution in 2- and 3-dimensional cases.
        
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
            Index of element of convex hull, which forms first approximation of v vector.
            It can be changed for lowering iterations sake
            For first approximation we'll just take point from a board of hull - cause it's easy reduced
        #############################################################
        #############################################################
        #############################################################
    """

    def __init__(self, points, dim, accel, max_iter = 500, init_approx_index = 1):
        self._dim = dim
        self._points = points.copy()
        self._hull = ConvexHull(points)
        self._A_matrix = points.copy().transpose()
        self.is_accelerated = accel
        self.iterations = None
        self.delta_p = None
        self.p_vector = None
        self.vector_current = None
        self.supp_vector = None
        self.max_iter = max_iter
        self.init_approx_index = init_approx_index


    def solve(self):
        V = 0
        iterations = 0
        delta_p = 1
        p_vector = [0 for i in range(0, len(self._points))]
        supp_vector = []
        t_param_vector =[]
        MIN_set = []
        MAX_set = []
        diff_vector = []                            # for cycles finding
        P_vectors = []                              # matrix of p_vectors
        V_vectors = []
        cycle_constructed = False
        cycle_is_constructing = False
        special_upd_done = False                    #whether special update Wn = W + lambdaV is done
        cycle_current_size = 0                      #we will search actual size of cycle
        vector_current = self._points[self._hull.vertices[self.init_approx_index]].copy()    #need copy() there for non-changing _points
        supp_vector.append(self._hull.vertices[self.init_approx_index])  # approximation => get vect_0
        p_vector[self._hull.vertices[self.init_approx_index]] = 1                         # working right concat
        #then we need to find vect_{k+1} iteratively
        while delta_p > 0.000001 and iterations < self.max_iter and len(supp_vector) != 0:
            if self.is_accelerated is True and cycle_constructed is True and special_upd_done is False:
                for i in range(cycle_size):  #constructing V as linear combination of D's that we used previously
                    V += -1 * t_param_vector[cycle_start + i - 1] * diff_vector[cycle_start + i -1]
                p_vector = P_vectors[cycle_start - 1]                   #returning to value where cycle had begun
                vector_current = V_vectors[cycle_start - 1]       #returning
                supp_vector = []                        #returning
                for i in range(len(p_vector)):          #returning
                    if p_vector[i] > 0.0000001:
                        supp_vector.append(i)

                lambda_t = -np.dot(vector_current, V) / np.linalg.norm(V) ** 2
                for i in range(cycle_size):
                    if t_param_vector[i] > 0:
                        if lambda_t > (1 - p_vector[MIN_set[i]]) / t_param_vector[i]:
                            lambda_t = (1 - p_vector[MIN_set[i]]) / t_param_vector[i]
                    elif t_param_vector[i] < 0:
                        if lambda_t > -p_vector[MAX_set[i]] / t_param_vector[i]:
                            lambda_t = -p_vector[MAX_set[i]] / t_param_vector[i]
                vector_current += lambda_t * V
                for i in range(cycle_size):
                    p_vector[MAX_set[i]] -= lambda_t * t_param_vector[i]
                    p_vector[MIN_set[i]] += lambda_t * t_param_vector[i]
                special_upd_done = True     #once it's done we're forgiving about that
            mult = np.dot(self._points[supp_vector], vector_current)
            ind_max = np.argmax(mult)           # finding max for indices in supp_vector
            ind_max = supp_vector[ind_max]      # finding max general in our mult product
            mult = np.matmul(vector_current, self._A_matrix)
            ind_min = np.argmin(mult)                                                 # i''_k
            if self.is_accelerated is True and cycle_constructed is False:
                MIN_set.append(ind_min)
                MAX_set.append(ind_max)
            diff = self._points[ind_max] - self._points[ind_min]
            delta_p = np.dot(diff, vector_current)
            if delta_p > 0.000001:                  #if not bigger, then we've found a solution
                t_param = delta_p /(p_vector[ind_max] * (np.linalg.norm(diff)) ** 2)  # recounting all variables
                if t_param >= 1:
                    t_param = 1
                ###################
                #Block for cycle finding. We must remember that MDM can improve wrong choices for updating vectors
                #It could be labeled as "false cycle construction"
                #for one cycles
                #####################
                if self.is_accelerated is True:          #if using accelerated MDM-method
                    if iterations > 0 and cycle_is_constructing is False:  #constructing cycle(active finding cycle, i mean, active-active)
                        contains = np.where(np.all(diff_vector == diff, axis = 1))[0]    #finds if diff_vector contains diff
                        if len(contains) != 0:      #found first element of cycle
                            cycle_is_constructing = True        #cycle is constructing now
                            cycle_start = contains[len(contains)-1]                     #index of first element of cycle; not changing
                            cycle_size = iterations - cycle_start         #not changing
                            cycle_current_size += 1         #this var for checking if all variables actually are cycle
                    elif cycle_is_constructing is True and cycle_constructed is False:
                        if cycle_current_size < cycle_size:
                            indices_match = np.where(np.all(diff_vector == diff, axis=1))[0]
                            if len(indices_match) > 0 and  indices_match[len(indices_match) - 1] == (cycle_start + cycle_current_size):
                                cycle_current_size += 1     #checking if cycle exists
                                if cycle_current_size == cycle_size:
                                    cycle_constructed = True
                            else:
                                cycle_is_constructing = False
                                cycle_current_size = 0
                    P_vectors.append(p_vector.copy())
                    V_vectors.append(vector_current.copy())
                    t_param_vector.append(t_param)
                    diff_vector.append(diff)
                vector_current -= t_param * p_vector[ind_max] * diff
                supp_vector = []                #recounting
                temp1 = t_param * p_vector[ind_max]
                temp2 = (1 - t_param)
                p_vector[ind_min] += temp1
                p_vector[ind_max] *= temp2
                for i in range(len(p_vector)):
                    if p_vector[i] > 0.0000001:
                        supp_vector.append(i)
            iterations += 1
        return vector_current
