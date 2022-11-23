import numpy as np
import time

def solver(a,b,c,d):
    if a == 0:
        if b == 0:
            if c == 0:
                roots = []
            else:
                roots = __linear_solver(c,d)
        else:
            roots = __quadratic_solver(b,c,d)
    else:
        roots = __cubic_solver(a,b,c,d)
    return roots

def __linear_solver(c,d):
    if c == 0:
        roots = []
    else:
        root = -d/c
        roots = [root]
    return roots

def __quadratic_solver(b,c,d):
    if c*c - 4*b*d == 0:
        root = -c/(2*b)
        roots = [root]
    elif c*c - 4*b*d < 0:
        roots = []
    else:
        root1 = (-c + np.sqrt(c*c - 4*b*d))/(2*b)
        root2 = (-c - np.sqrt(c*c - 4*b*d))/(2*b)
        roots = [root1, root2]
    return roots

def __cubic_solver(A,B,C,D):
    d = 18*A*B*C*D - 4*(B**3)*D + (B**2)*(C**2) - 4*A*(C**3) - 27*(A**2)*(D**2)
    if d > 0:
        Q_ = (3*(C/A) - (B/A)**2)/9
        R_ = (9*(B/A)*(C/A) - 27*(D/A) -2*(B/A)**3)/54
        term = R_/np.sqrt((-Q_)**3)
        theta_ = np.arccos(term)
        t1 = 2*np.sqrt(-Q_)*np.cos(theta_/3) - (B/A)/3
        t2 = 2*np.sqrt(-Q_)*np.cos((theta_+2*np.pi)/3) - (B/A)/3
        t3 = 2*np.sqrt(-Q_)*np.cos((theta_+4*np.pi)/3) - (B/A)/3
        roots = [t1,t2,t3]
    elif d < 0:
        P = B**2 - 3*A*C
        Q = 9*A*B*C - 2*(B**3) - 27*(A**2)*D
        cube_root = lambda x: x**(1./3.) if x >= 0 else -(-x)**(1./3.)
        term_1 = Q/2 + np.sqrt((Q**2)/4 - P**3)
        term_2 = Q/2 - np.sqrt((Q**2)/4 - P**3)
        N = cube_root(term_1) + cube_root(term_2)
        t1 = -B/(3*A) + N/(3*A)
        roots = [t1]
    else:
        P = B**2 - 3*A*C
        if P == 0:
            t1 = -B/(3*A)
            roots = [t1]
        else:
            t1 = (9*A*D - B*C)/(2*P)
            t2 = (4*A*B*C - 9*A**2*D - B**3)/(A*P)
            roots = [t1,t2]
    return roots