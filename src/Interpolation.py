
def lagrange_interp(a,f,x):
    """Compute a lagrange interpolant
    a: array of n points
    f: array of the value of f(a) at the n points"""
    answer = 0 
    assert a.size==f.size
    n=a.size
    for i in range(n):
        product = 1
        for j in range(n):
            if i!=j:
                product *= (x-a[j])/(a[i]-a[j])
        answer += product*f[i]
    return answer

def spline_interp(x,f,option=3):
    """Interpolation using CubicSpline
    option1: natural spline condition where the second derivative of the start and end points equals zero
    option2: clamped spline conditions where the first derivative of the start and end points equals zero
    option3: not a knot spline conditions where the third derivative of the start and end points equals to the same value as to the second and second to the last points"""
    assert x.size==f.size

    from scipy.interpolate import CubicSpline
    if option==1:
        return CubicSpline(x,f,bc_type='natural')
    elif option == 2:
        return CubicSpline(x,f,bc_type='clamped')
    else:
        return CubicSpline(x,f)


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

    import numpy as np

    a = np.linspace(0.2,np.pi*1.8,5)
    data = np.zeros((5,2))
    data[:,0]=a
    data[:,1]=np.sinh(a)

    splineFunc= spline_interp(data[:,0],data[:,1])

    points=200
    X=np.linspace(0,np.pi*2,points)

    plt.plot(X,splineFunc(X),linestyle='--')
    plt.plot(X,np.sinh(X),linestyle='-')