# find the root of a given function and given interval
# close root: ridder method
# open root: inexact_newton
import numpy as np

from Matrix import *

def ridder(f,a,b,epsilon=1.0e-6,iter_max=200):
    """find the root of function f in interval [a,b]
    """

    assert b>a
    fa=f(a)
    fb=f(b)
    try:
        assert fa*fb<0
    except:
        print("The function has the same sign at the border of the given interval")
        quit()
    delta = b-a
    iterations=0
    residual = 1.0
    while np.fabs(residual)>epsilon:
        c = 0.5*(a+b)
        d = 0.0
        fc = f(c)
        if (fa-fb>0):
            d=c+(c-a)*fc/np.sqrt(fc**2-fa*fb)
        else:
            d=c-(c-a)*fc/np.sqrt(fc**2-fa*fb)
        fd=f(d)
        #
        if fa*fd<0:
            b=d
            fb=fd
        elif fb*fd<0:
            a=d
            fa=fd
        residual = fd
        iterations+=1
        if iterations==iter_max:
            print("maximum iteration number %s reached" % str(iter_max))
            print("the closest prediction is %s at point %s" % (str(residual), str(d)))
            break
    print("It took ", iterations, " interations")
    return d

def inexact_newton(f,x0=0,delta=1.0e-7,epsilon=1.0e-6, iter_max=200, LOUD=False):
    """find the open root of a arbitrary function f,
    x0 initial guess
    delta for estimating the slope
    """
    x=x0
    iterations=0
    fx=f(x)
    while (np.fabs(fx)>epsilon):
        fxdelta=f(x+delta)
        slope=(fxdelta-fx)/delta
        x=x-fx/slope
        fx=f(x)
        iterations+=1
        if iterations==iter_max:
            print("maximum iteration number %s reached" % str(iter_max))
            print("the closest prediction is %s at point %s" % (str(fx), str(x)))
    print("It took ", iterations, " Iterations")
    return x

def newton_system(f, x0, delta=1.0e-7,epsilon=1.0e-6, iter_max=200):
    """find the open root of a system f using jacobian matrix J(xi)(xi+1 - xi)=-F(xi)
    given initial guess of x as x0 (np.array)
    J(xi)= matrix of partial derivative with respect to xi 
    """
    
    def jacobian(f,x,delta):
        N=x0.size
        J=np.zeros((N,N))
        fx=f(x)
        x_plusdelta=x.copy()
        idelta=1/delta
        #print(x_plusdelta)
        for i in range(N):
            x_plusdelta[i] += delta
            col = (f(x_plusdelta)-fx)/delta
            x_plusdelta[i]=x[i]
            J[:,i]=col
        return J
    x=x0
    iterations=0
    fx=f(x)
    while (np.linalg.norm(fx)>epsilon):
        J=jacobian(f,x,delta)

        RHS = -fx

#        row_order = LU_factor(J,LOUD=False)
        delta_x = np.linalg.lstsq(J,RHS,rcond=0)[0]
        #print("the solved deltax is:",delta_x)
        x += delta_x

        fx=f(x)
        iterations += 1

        if iterations==iter_max:
            break
    print('It took ', iterations, ' Iterations')
    return x
            

if __name__=="__main__":
#    f = lambda x: 3*x**3+2*x**2-5*x-20
#
#    print(ridder(f,-5,2))
#    print(inexact_newton(f,1.5))

    def Reactor(x, D=9.21, nuSigf =0.157, Siga=0.1532):
        answer=np.zeros((3))
        answer[0] = (np.pi/(x[0]+2*D))**2 + \
                    (np.pi/(x[1]+2*D))**2 + \
                    (np.pi/(x[2]+2*D))**2 - (nuSigf-Siga)/D
        answer[1] = 2*(x[0]*x[1]+x[0]*x[2]+x[1]*x[2])-1.2e6
        answer[2] = x[0]-x[1]
        return answer

    x0=np.array([7000.,7000.,100.])
    x=newton_system(Reactor, x0, epsilon=1e-8, delta=1.0e-10)