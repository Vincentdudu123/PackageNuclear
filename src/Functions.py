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
    
def LU_factor(A):
    """
    factor L dot U = A
    """
    [Nrow, Ncol] = A.shape
    assert Nrow == Ncol
    N=Nrow
    for column in range(0,N):
        for row in range(column+1,N):
            mod_row = A[row]
            factor = mod_row[column]/A[column,column]
            mod_row=mod_row-factor*A[column,:]

            mod_row[column]=factor
            mod_row = mod_row[column:N]
            A[row,column:N]=mod_row
    return

def LU_solve(A,b):
    """LU_factor first and then solve"""
    [Nrow, Ncol] = A.shape
    assert Nrow == Ncol    
    N=Nrow
    LU_factor(A)

    x=np.zeros(N)
    # temporary vector for L^-1 b
    y=np.zeros(N)

    for row in range(N):
        RHS=b[row]
        for column in range(0,row):
            RHS-= y[column]*A[row,column]
        y[row]=RHS

    for row in range(N-1,-1,-1):
        RHS = y[row]
        for column in range(row+1,N):
            RHS -= x[column]*A[row,column]
        x[row] = RHS/A[row,row]

    return x


if __name__=="__main__":
 
    import matplotlib.pyplot as plt

    
    sin3 = lambda x,y,z: np.sin(x)*np.sin(y)*np.sin(z)
    print(integral_midpoint_3D(sin3,0,np.pi,0,np.pi,0,np.pi,100,100,100))

    A = np.array([(3.0,2,1),(-1,4,5),(2,-8,10)])
    print(A)
    b = np.array([6,8,4])
    print(LU_solve(A,b))
 