"""
This file contains few methods doing linear regression using numpy (np.linalg.lstsq(A,points)[0] which returns a set of parameters of a0 a1 a2.... 
number of ax depends on the number of attributes in A),
Matrix.py 
"""
import numpy as np

def linearRegress(data,points):
    if type(data)==dict:
        print("type of data is Dictionary")
        for i in range(len(list(data.values()))):
            assert len(points)==len(list(data.values())[i])
        size=len(list(data.values())[0])
        A=np.ones(size)
        keywords=list(data.keys())
        for i in keywords:
            A=np.vstack([A,data[i]])
        A=A.transpose()
        print(A)
        solution = np.linalg.lstsq(A,points,rcond=0)[0]
        print("The solution of this function is: ", solution[0])
        for i in range(len(keywords)):
            print(" and ", keywords[i], solution[i+1])


    elif type(data)==list:
        print("type of data is a list")
        for i in range(len(data.values())):
            assert len(points)==len(data[i])
        size=len(data[0])
        A=np.ones(size)
        for i in range(len(data)):
            A=np.vstack([A,data[i]])
        A=A.transpose()
        solution = np.linalg.lstsq(A,points,rcond=0)[0]
        print("The solution of this function is: ", solution[0])
        for i in range(len(data)):
            print(" and ", solution[i+1])

    else:
        print("Non-supported data struction, must be either dictionary or list")
        assert "1"=="0"
    #mean absolute error
    yhat=np.dot(A,solution)
    errorAbs=np.mean(np.fabs(points-yhat))
    return solution,errorAbs        

if __name__=="__main__":
    A={"a":[486,714,628,581,523,587,602,558,564,537,299,379,541],"b":[1,1,2,1,0,2,2,5,1,3,2,1,0],"c":[2,2,1,3,2,1,1,3,4,1,0,0,2]}
    a,b=linearRegress(A,[52,65,42,42,45,41,41,56,57,51,10,21,52])
    print(b)
#
#    b=[1,1,2,1,0,2,2,5,1,3,2,1,0]
#    c=[2,2,1,3,2,1,1,3,4,1,0,0,2]
#    bc=[]
#    for i in range(len(b)):
#        bc.append(b[i]*c[i])
#    A={"a":[486,714,628,581,523,587,602,558,564,537,299,379,541],"bc":bc}
#    a,b=linearRegress(A,[52,65,42,42,45,41,41,56,57,51,10,21,52])
#    print(b)
#


    


