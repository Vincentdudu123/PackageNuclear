"""
This file contains some functions in python which I think it will be useful in the future
"""
import numpy as np

def integral_midpoint_1D(f,ax,bx,nr_x):
    '''
    f:     function to be integrated, one variable
    ax,bx: lower and upper bound of the variable x
    nr_x:  number of intervals in x
    '''
    dx = (bx-ax)/nr_x
    midpoints = np.arange(nr_x)*dx+0.5*dx+ax
    integral = 0
    for xx in midpoints:
        integral += f(xx)
    return integral*dx


def integral_midpoint_2D(f,ax,bx,ay,by,nr_x,nr_y):
    '''
    f:           function to be integrated, contains two variables
    ax,bx ay,by: lower and upper bounds of x and y
    nr_x,nr_y:   number of intervals in x & y
    '''
    g = lambda y: integral_midpoint_1D(lambda x: f(x,y),ax,bx,nr_x)
    return integral_midpoint_1D(g,ay,by,nr_y)

def integral_midpoint_3D(f,ax,bx,ay,by,az,bz,nr_x,nr_y,nr_z):

    '''
    f:                 function to be integrated, contains three variables
    ax,bx ay,by az,bz: lower and upper bounds of x and y and z
    nr_x,nr_y,nr_z:    number of intervals in x & y & z
    '''   

    g = lambda x,y: integral_midpoint_1D(lambda z: f(x,y,z), az, bz, nr_z)
    return integral_midpoint_2D(g,ax,bx,ay,by,nr_x,nr_y)
    

if __name__=="__main__":
 
    import matplotlib.pyplot as plt

    
#    sin3 = lambda x,y,z: np.sin(x)*np.sin(y)*np.sin(z)
#    print(integral_midpoint_3D(sin3,0,np.pi,0,np.pi,0,np.pi,100,100,100))

#    A = np.array([(3.0,2,1,-2),(-1,4,5,4),(2,-8,10,3),(-2,-8,10,0.1)])
#    answer = np.ones(4)
#    b = np.dot(A,answer)
#    print(b)
#    row_order = LU_factor(A)

#    print(LU_solve(A,b,row_order))
    

 