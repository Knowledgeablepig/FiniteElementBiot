import numpy as np
import matplotlib as mpl
import math as math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


#Even more crucial I should have [0,L],[0,T], p and h values here.
startingSpace = 0
startingTime = 0
finalTime = 2
finalSpace = 10
h = .1
p = .001
amountX = (int) (finalSpace/h)
amountT = (int) (finalTime/p)
porosity = 1

#Alright I should have functions and ICs here

#Sources
#Due to indexing F[0] and S[0] are both F(t_1,x) and S(t_1,x), applegies
F = np.zeros(shape = (amountT,amountX))
S = np.zeros(shape = (amountT,amountX))
for i in range(0,amountT):
    for j in range(0,amountX):
        tNow = startingTime + p*i
        xNow = startingSpace + h*(j)
        #Adjust function here
        #CURRENT FUNCTION x^3t^3
        F[i][j] = math.sin(math.pow(xNow,3)+math.pow(tNow,3))

for i in range(0,amountT):
    for j in range(0,amountX):
        tNow = startingTime + p*i
        xNow = startingSpace + h*(j)
        #Adjust function here
        #CURRENT FUNCTION x^2t^2
        S[i][j] = math.cos(math.pow(xNow,2)+math.pow(tNow,2))

#Initial Conditions
uNaught = np.zeros(amountX+1)
#Boundary Condition
uNaught[0] = 0
uNaught[amountX] = 0
for j in range(1,amountX):
    #RemainingICTerms
    xNow = startingSpace + h*(j)
    uNaught[j] = math.sin(math.pi*xNow*.5)

# pNaught = np.zeros(amountX+1)
# #Boundary Condition
#     pNaught[0] = 0
# for j in range(1,amountX)
#     #RemainingICTerms
#     pNow = startingSpace + h*(j)
#     pNaught[i] = math.pow(xNow,2)

#SolutionInitializations
uSol = np.zeros(shape = (amountT+1,amountX+1))
pSol = np.zeros(shape = (amountT+1,amountX+1))

for j in range(1,amountX):
    uSol[0][j] = uNaught[j]


#FindingaSolution

def pUpdate(timeStep):
    #Initializations in case of Boundary Conditions
    pSol[timeStep][0] = 0
    #Updatestep
    for j in range(1,amountX):
        pSol[timeStep][j] = pSol[timeStep][j-1] + h*F[timeStep][j] + (uSol[timeStep][j+1] - 2*uSol[timeStep][j] + uSol[timeStep][j-1])/h
    pSol[timeStep][amountX] = pSol[timeStep][amountX-1]    
    return

def uUpdate(timeStep):
    #Initializations in case of Boundary Conditions
    uSol[timeStep][0] = 0
    for j in range(0,amountX-1):
       uSol[timeStep][j+1] = uSol[timeStep-1][j+1] + uSol[timeStep][j] - uSol[timeStep-1][j] + p*S[timeStep-1][j] + porosity*p*(pSol[timeStep-1][j+1] - 2*pSol[timeStep-1][j] + pSol[timeStep-1][j-1])/h
    uSol[timeStep][amountX] = uSol[timeStep][amountX-1]
    return

def FullRun():
    print("eat several cats")
    pUpdate(0)
    for i in range(1,amountT):
        uUpdate(i)
        pUpdate(i)
    return

def GraphP(timeStep):

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.set_xlim(startingSpace,finalSpace)
    # ax.set_ylim(startingTime,finalTime)
    # ax.set_xlabel('x')
    # ax.set_ylabel('time')
    # plt.show()

    fig, ax = plt.subplots()
    xArray = np.zeros(amountX+1)
    plotArray = np.zeros(amountX+1)
    for j in range(0,amountX+1):
        plotArray[j] = pSol[timeStep][j]
        xArray[j] = h*j
    line1, = ax.plot(xArray,plotArray,label = "Discrete Solution", lw = 4, ms = 20)
    plt.show()
    return

def GraphU(timeStep):

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    # ax.set_xlim(startingSpace,finalSpace)
    # ax.set_ylim(startingTime,finalTime)
    # ax.set_xlabel('x')
    # ax.set_ylabel('time')
    # plt.show()

    fig2, ax2 = plt.subplots()
    xArray = np.zeros(amountX+1)
    plotArray = np.zeros(amountX+1)
    for j in range(0,amountX+1):
        plotArray[j] = uSol[timeStep][j]
        xArray[j] = h*j
    line1, = ax2.plot(xArray,plotArray,label = "Discrete Solution", lw = 4, ms = 20)
    plt.show()
    return    


FullRun()
GraphU(amountT-1)