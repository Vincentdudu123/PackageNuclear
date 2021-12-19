"""
This file contains some functions in python which I think it will be useful in the future
"""


class PKModel(object):
    def __init__(self,n=1.0,beta=[0.000215,0.001424,0.001274,0.002568,0.000748,0.000273],llambda=[0.0124,0.0305,0.111,0.301,1.14,3.01],promptGeT=2.0e-5):
        self.n=n
        self.beta=beta
        self.llambda=llambda
        self.promptGeT=promptGeT
        self.precCon=[]

        self.sumbeta=0.0
        for i in self.beta:
            self.sumbeta+=i

        self.check()
        self.calcpreCon()

    
    def check(self):
        if len(self.beta)!=len(self.llambda):
            print("The lenght of beta and llambda are not equal.")
            

    def calcpreCon(self):
        print(len(self.beta))
        for i in range(len(self.beta)):
            self.precCon.append(self.beta[i]*self.n/self.promptGeT/self.llambda[i])
    
    def nextTimeStep(self,rho,dt=0.0001):

    
        if dt>0.0001:
            a = int(dt/0.0001)
            b = dt%0.0001

            dt = 0.0001
        else:
            a = 1
            b = -1.0

            dt=dt

        for i in range(a):

            self.dn=rho/self.promptGeT*self.n*dt

            self.dCon=[]

            for i in range(len(self.beta)):
                self.dCon.append((self.beta[i]/self.promptGeT*self.n-self.llambda[i]*self.precCon[i])*dt)
                self.dn-=self.dCon[-1]
            #print(self.dn)
            #print(self.precCon)
            self.synchronize()
        if b>0:
            self.dn=rho/self.promptGeT*self.n*b

            self.dCon=[]

            for i in range(len(self.beta)):
                self.dCon.append((self.beta[i]/self.promptGeT*self.n-self.llambda[i]*self.precCon[i])*b)
                self.dn-=self.dCon[-1]
            #print(self.dn)
            #print(self.precCon)
            self.synchronize()            
        
    
    def synchronize(self):
        self.n+=self.dn
        for i in range(len(self.beta)):
            self.precCon[i]+=self.dCon[i]


if __name__=="__main__":
 
    import matplotlib.pyplot as plt

    reactor = PKModel(n=10000)
    #print(reactor.precCon)
    n=[]
    time=[]
    n.append(reactor.n)
    time.append(0)

    dt=0.011
    for i in range(100):
        time.append(dt*(i+1))
        reactor.nextTimeStep(-5.0e-3,dt)
        n.append(reactor.n)
    
    reactor2 = PKModel(n=10000)
    n2=[10000]
    time2=[0]
    dt=0.0001
    for i in range(100000):
        time2.append(dt*(i+1))
        reactor2.nextTimeStep(-5.0e-3,dt)
        n2.append(reactor2.n)


    fig = plt.figure()

    ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])

    ax.plot(time, n)
    ax.plot(time2, n2)

    reactor3 = PKModel(n=10000)
    time=[0,0.1,1,5,10]
    dt=[]
    n3 = [10000]
    for i in range(len(time)-1):
        dt=time[i+1]-time[i]
        reactor3.nextTimeStep(-5.0e-3,dt)
        n3.append(reactor3.n)
    ax.plot(time, n3,"bo")   
    plt.show()
    