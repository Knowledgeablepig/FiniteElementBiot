import numpy as np
import matplotlib as mpl
import math as math
import matplotlib.pyplot as plt

#Even more crucial I should have [0,L],[0,T], p and h values here.
startingSpace = 0
startingTime = 0
finalTime = 2
finalSpace = 10
h = .01
p = .001
amountX = (int) (finalSpace/h)
amountT = (int) (finalTime/p)
porosity = 1
R = 100 #This is the connecting resistor value
C = .002 #This is the connecting capacitor value
a = .05

#Alright I should have functions and ICs here

#Sources
#Due to indexing F[0][j], S[0][j], s[0] are F(p,j*h), S(p,j*h) and s(p)
F = np.zeros(shape = (amountT,amountX))
Ftilde = np.zeros(shape = (amountT,amountX))
S = np.zeros(shape = (amountT,amountX))
s = np.zeros(amountT)
for i in range(0,amountT):
    for j in range(0,amountX):
        tNow = startingTime + p*(i+1)
        xNow = startingSpace + h*(j+1)
        #Adjust function here
        #CURRENT FUNCTION sin(x+t^3)
        F[i][j] = math.sin(xNow +math.pow(tNow,3))

for i in range(0,amountT):
    for j in range(0,amountX):
        tNow = startingTime + p*(i+1)
        xNow = startingSpace + h*(j+1)
        #Adjust function here
        #CURRENT FUNCTION x^3t^3
        Ftilde[i][j] = math.sin(math.pow(xNow,3)+math.pow(tNow,3))*3*math.pow(tNow,2)

for i in range(0,amountT):
    for j in range(0,amountX):
        tNow = startingTime + p*(i+1)
        xNow = startingSpace + h*(j+1)
        #Adjust function here
        #CURRENT FUNCTION cos(x^2+t^2)
        S[i][j] = math.cos(math.pow(xNow,2)+math.pow(tNow,2))

for i in range(0,amountT):
  tNow = startingTime + p*(i+1)
  #Adjust function Here
  #CURRENT FUNCTION t
  s[i] = tNow

#Initial Conditions
uNaught = np.zeros(amountX+1)
piNaught = .23

#Boundary Condition
uNaught[0] = 0
uNaught[amountX] = 0
for j in range(1,amountX):
    #RemainingICTerms
    xNow = startingSpace + h*(j)
    uNaught[j] = math.sin(math.pi*xNow*.5)

#SolutionInitializations
uSol = np.zeros(shape = (amountT+1,amountX+1))
pSol = np.zeros(shape = (amountT+1,amountX+1))
piSol = np.zeros(amountT+1)


for j in range(1,amountX):
    uSol[0][j] = uNaught[j]


#BackwardEuler Matrix Intialization
matrixA = np.zeros(shape = (amountX,amountX))
matrixA[0][0] = 1 + 2*porosity*p/(h*h)
matrixA[0][1] = -porosity*p/(h*h)
matrixA[amountX-1][amountX-2] = (1 - h/R)*(p/(R*C))
matrixA[amountX-1][amountX-1] = (1 - (a + 1/(R*C))*p - p*h/(R*R*C))
for i in range(1,amountX-1):
    matrixA[i][i-1] =  -porosity*p/(h*h)
    matrixA[i][i] = 1 + 2*porosity*p/(h*h)
    matrixA[i][i+1] =  -porosity*p/(h*h)


#FindingaSolution

def pAndpiUpdate(timeStep):
    #Define the variable vectors in Ax + f = y
    xVector = np.zeros(amountX)
    yVector = np.zeros(amountX)
    fVector = np.zeros(amountX)

    #Define y
    for j in range(0,amountX-1):
        yVector[j] = pSol[timeStep][j+1]
    yVector[finalSpace-1] = piSol[timeStep]

    #Define f
    for j in range(0,amountX-1):
        fVector[j] = Ftilde[timeStep+1][j] + S[timeStep+1][j]
    fVector[finalSpace-1] = s[timeStep+1]
    
    #Turning f + Ax = y into Ax = y-f
    for j in range(0,amountX):
        yVector[j] = yVector[j] - fVector[j]

    #Solve the linear system and label solutions correctly
    xVector = np.linalg.solve(matrixA,yVector)

    if(timeStep == 4):
        for j in range(0,amountX):
            print(xVector[j])

    for j in range(0,amountX-1):
        pSol[timeStep+1][j+1] = xVector[j]
    piSol[timeStep+1] = xVector[amountX-1]
    return


#Initial Condition is in terms of u, solves for p(0,x)
def pIntialize(timeStep):
    #Initializations in case of Boundary Conditions
    pSol[timeStep][0] = 0
    #Updatestep
    for j in range(1,amountX):
        pSol[timeStep][j] = pSol[timeStep][j-1] + h*F[timeStep][j] + (uSol[timeStep][j+1] - 2*uSol[timeStep][j] + uSol[timeStep][j-1])/h
        pSol[timeStep][amountX] = (piSol[timeStep]*(h/R) + pSol[timeStep][amountX-1])*(1-h/R)
    return

#p and pi Update function does not solve for u, so calling this returns the corresponding displacement
def uSolver(timeStep):
    for j in range(0,amountX-1):
        uSol[timeStep][j+1] = -1*h*(F[timeStep][j] + pSol[timeStep][j]) + uSol[timeStep][j]
    return

#Run Function
def FullRun():
    print("Pet Several Dogs")
    pIntialize(0)   
    for i in range(0,amountT-1):
        pAndpiUpdate(i)
        uSolver(i)
    return

#Graph pressure function
def GraphP(timeStep):

    fig, ax = plt.subplots()
    xArray = np.zeros(amountX+1)
    plotArray = np.zeros(amountX+1)
    for j in range(0,amountX+1):
        plotArray[j] = pSol[timeStep][j]
        xArray[j] = h*j
    line1, = ax.plot(xArray,plotArray,label = "Discrete Solution", lw = 4, ms = 20)
    plt.show()
    return

#Graph displacement function
def GraphU(timeStep):

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
GraphP(amountT-1)