import numpy as np

from scipy.linalg import lu, inv

def gausselim(A,B):
    """
    Solve Ax = B using Gaussian elimination and LU decomposition.
    A = LU   decompose A into lower and upper triangular matrices
    LUx = B  substitute into original equation for A
    Let y = Ux and solve:
    Ly = B --> y = (L^-1)B  solve for y using "forward" substitution
    Ux = y --> x = (U^-1)y  solve for x using "backward" substitution
    :param A: coefficients in Ax = B
    :type A: numpy.ndarray of size (m, n)
    :param B: dependent variable in Ax = B
    :type B: numpy.ndarray of size (m, 1)
    """
    # LU decomposition with pivot
    pl, u = lu(A, permute_l=True)
    # forward substitution to solve for Ly = B
    y = np.zeros(B.size)
    for m, b in enumerate(B.flatten()):
        y[m] = b
        # skip for loop if m == 0
        if m:
            for n in range(m):
                y[m] -= y[n] * pl[m,n]
        y[m] /= pl[m, m]

    # backward substitution to solve for y = Ux
    x = np.zeros(B.size)
    lastidx = B.size - 1  # last index
    for midx in range(B.size):
        m = B.size - 1 - midx  # backwards index
        x[m] = y[m]
        if midx:
            for nidx in range(midx):
                n = B.size - 1  - nidx
                x[m] -= x[n] * u[m,n]
        x[m] /= u[m, m]
    return x



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
    assert b.size==Ncol
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
    [Nrow, Ncol]=A.shape
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