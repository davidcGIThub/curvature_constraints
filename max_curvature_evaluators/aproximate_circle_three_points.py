# https://www.geeksforgeeks.org/equation-of-circle-when-three-points-on-the-circle-are-given/

# Python3 implementation of the approach
from math import sqrt
import matplotlib.pyplot as plt
 
# Function to find the circle on
# which the given three points lie
def findCircle(x1, y1, x2, y2, x3, y3) :
    x12 = x1 - x2
    x13 = x1 - x3

    y12 = y1 - y2
    y13 = y1 - y3
 
    y31 = y3 - y1
    y21 = y2 - y1
 
    x31 = x3 - x1
    x21 = x2 - x1
 
    # x1^2 - x3^2
    sx13 = pow(x1, 2) - pow(x3, 2)
 
    # y1^2 - y3^2
    sy13 = pow(y1, 2) - pow(y3, 2)

    sx21 = pow(x2, 2) - pow(x1, 2)
    sy21 = pow(y2, 2) - pow(y1, 2)
 
    f = (((sx13) * (x12) + (sy13) *
          (x12) + (sx21) * (x13) +
          (sy21) * (x13)) // (2 *
          ((y31) * (x12) - (y21) * (x13))))
             
    g = (((sx13) * (y12) + (sy13) * (y12) +
          (sx21) * (y13) + (sy21) * (y13)) //
          (2 * ((x31) * (y12) - (x21) * (y13))))
 
    c = (-pow(x1, 2) - pow(y1, 2) -
         2 * g * x1 - 2 * f * y1)
 
    # eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0
    # where centre is (h = -g, k = -f) and
    # radius r as r^2 = h^2 + k^2 - c
    h = -g
    k = -f
    sqr_of_r = h * h + k * k - c
 
    # r is the radius
    r = round(sqrt(sqr_of_r), 5)
 
    print("Centre = (", h, ", ", k, ")")
    print("Radius = ", r)
    return r, [h,k]
 
# Driver code
if __name__ == "__main__" :
     
    x1 = 1 ; y1 = 1
    x2 = 2 ; y2 = 4
    x3 = 5 ; y3 = 3
    radius, center = findCircle(x1, y1, x2, y2, x3, y3)
    circle=plt.Circle(center,radius,color = 'b',fill=False)
    fig,ax = plt.subplots()
    ax.scatter([x1,x2,x3],[y1,y2,y3])
    ax.add_patch(circle)
    plt.show()
 
# This code is contributed by Ryuga