"""
This file contains some functions in python which I think it will be useful in the future
"""

class Nuclide(object):

    def __init__(self,name="a",llambda=0.0,Concentration=0.0):
        self.Concentration=Concentration
        self.name=name
        self.parent={}
        self.daughter={}
        self.activity=0.0
        self.Con=[]
        self.llambda=llambda

    def resetNumber(self):
        for i in self.daughter:
            i["number"]=0.0

    def addParent(self,parentNuclide):
        if parentNuclide.name not in self.parent.keys():
            self.parent[parentNuclide.name]=parentNuclide

    def addDaughter(self,daughterNuclide,yyield):
        if daughterNuclide.name not in self.daughter.keys():
            self.daughter[daughterNuclide.name]={ "nuclide": daughterNuclide, "yyield": yyield, "number": 0.0}
            daughterNuclide.addParent(self)
        else:
            self.daughter[daughterNuclide.name]["yyield"] = yyield

    def decay(self,dt):
        self.activity=0.0
        for i in self.daughter.keys():
            self.daughter[i]["number"]=self.daughter[i]["yyield"]*dt*self.Concentration*self.llambda
        self.activity=self.llambda*self.Concentration
    
    def synchronize(self):
        for i in self.daughter.keys():
            self.daughter[i]["nuclide"].Concentration+=self.daughter[i]["number"]
            self.Concentration-=self.daughter[i]["number"]
            self.daughter[i]["number"]=0.0
        self.prepare()

    def prepare(self):
        self.Con.append(self.Concentration)

if __name__=="__main__":
 
    import matplotlib.pyplot as plt

 

    U235=Nuclide("U235",0.02)
    U238=Nuclide("U238",0.01)
    U232=Nuclide("U232",0.03)
    U233=Nuclide("U233")

    U235.addDaughter(U232,1)
    U238.addDaughter(U232,1.0)
    U232.addDaughter(U233,1)

    U235.Concentration=100.0
    U238.Concentration=50.0


    Nuclides=[U235,U232,U233,U238]
    for j in Nuclides:
        j.prepare()
    Time=[0.0]
    dt=0.01

    for i in range(10000):
        Time.append(dt*(i+1))
        for j in Nuclides:
            j.decay(dt)
        for j in Nuclides:
            j.synchronize()


    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])

    ax.plot(Time, U235.Con,label="U235")
    ax.plot(Time, U232.Con,label="U232")
    ax.plot(Time, U238.Con,label="U238")
    ax.plot(Time, U233.Con,label="U233")
    ax.legend()
    plt.show()
 