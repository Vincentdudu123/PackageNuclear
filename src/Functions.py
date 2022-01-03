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
    
def LU_factor(A,LOUD=True):
    """
    factor L dot U = A, modifies matrix A to be composed of L and U in which the diagonals of L are 1.0
    return row order
    """
    [Nrow, Ncol] = A.shape
    assert Nrow == Ncol
    N=Nrow
    # create scale factors
    s = np.zeros(N)
    count = 0
    row_order = np.arange(N)
    for row in A:
        s[count] = np.max(np.fabs(row))
        count+=1
    if LOUD:
    #    print("s= ",s)
        print("Original Matrix is \n", A)
    
    for column in range(0,N):
        #swap rows if needed
        largest_pos = np.argmax(np.fabs(A[column:N,column])/s[column]) + column
        if largest_pos != column:
            if LOUD:
                print("Swapping row", column, "with row", largest_pos)
            swap_rows(A,column,largest_pos)
            #changes to RHS
            tmp=row_order[column]
            row_order[column]=row_order[largest_pos]
            row_order[largest_pos]=tmp
            #re/order s
            tmp=s[column]
            s[column]=s[largest_pos]
            s[largest_pos]=tmp
        
        for row in range(column+1,N):
            mod_row = A[row]
            factor = mod_row[column]/A[column,column]
            mod_row=mod_row-factor*A[column,:]
            # put the factor in the correct place in the modified row we need
            mod_row[column]=factor
            # only take the part of the modified row we need
            mod_row = mod_row[column:N]
            A[row,column:N]=mod_row
    if LOUD:
        print("after swapping and factorizing\n",A)
    return row_order

def LU_solve(A,b,row_order):
    """take the pivoting order"""
    [Nrow, Ncol] = A.shape
    assert Nrow == Ncol    
    assert b.size++Ncol
    assert row_order.max()==Ncol-1
    N=Nrow
    # reorder the equations
    tmp = b.copy()

    for row in range(N):
        b[row_order[row]]=tmp[row]

    x=np.zeros(N)
    # temporary vector for L^-1 b
    y=np.zeros(N)
    #forward solve
    for row in range(N):
        RHS=b[row]
        for column in range(0,row):
            RHS -= y[column]*A[row,column]
        y[row]=RHS

    for row in range(N-1,-1,-1):
        RHS = y[row]
        for column in range(row+1,N):
            RHS -= x[column]*A[row,column]
        x[row] = RHS/A[row,row]

    return x

def swap_rows(A,a,b):
    """Rows two rows in a matrix, switch row a with row b, a and b both index of rows"""
    assert (a>=0) and (b>=0)
    N = A.shape[0] # number of rows
    assert (a<N) and (b<N) 
    temp = A[a,:].copy()
    A[a,:] = A[b,:].copy()
    A[b,:] = temp.copy() 

def Jacobi_fast(A,b,x0=np.array([]),tol=1.0e-6,max_iter=1000,LOUD=False):
    """Solve a Linear system using Jacobi method based on initial guess x0"""
    [Nrow, Ncol]=A.shape()
    assert Nrow==Ncol
    N = Nrow
    converged=False
    iteration=1
    if x0.size==0:
        x0=np.random.rand(N)
    x=x0.copy()
    while not converged:
        x0=x.copy()
        x=(b-np.dot(A,x0)+A.diagonal()*x0)/A.diagonal()
        relative_change=np.linalg.norm(x-x0)/np.linalg.norm(x)
        if LOUD:
            print("iteration",iteration, ": Relative Change =", relative_change)
        if (relative_change<tol) or (iteration >=max_iter):
            converged=True
        iteration+=1
    return x

if __name__=="__main__":
 
    import matplotlib.pyplot as plt

    
#    sin3 = lambda x,y,z: np.sin(x)*np.sin(y)*np.sin(z)
#    print(integral_midpoint_3D(sin3,0,np.pi,0,np.pi,0,np.pi,100,100,100))

    A = np.array([(3.0,2,1,-2),(-1,4,5,4),(2,-8,10,3),(-2,-8,10,0.1)])
    answer = np.ones(4)
    b = np.dot(A,answer)
#    print(b)
    row_order = LU_factor(A)

    print(LU_solve(A,b,row_order))
 