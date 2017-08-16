# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 11:06:38 2017

@author: jazicoan
"""

"""
Note:
- MW = mainWing
- SW = servoWing
- CL = lift coefficient
- CD = drag coefficient
- uppercase variable names are often abbreviation of the function used to calculate it
"""

import numpy as np
from drawingFunctions import drawWingSailAngle, drawWingSailIRT, drawEquilibriumAngles, drawArrow, drawEvolutionMWAngle, drawHull
from integrationSchemes import eulerScheme, rk2Scheme,rk2Step
from angleFunctions import wrapTo2pi, degTorad, radTodeg, listRadTodeg, listDegTorad, listWrapTo2pi, apparentWindCoordinatesToMWCoordinates, MWCoordinatesToBoatCoordinates
import matplotlib.pyplot as plt

def reynolds(Speed,Lcarac,nu):
    return(Speed*Lcarac/nu)





X = 2

"""Constants"""
global trueWindAngle
trueWindAngle              = float(input('Enter true wind angle angle in degrees:'))
trueWindAngle              = wrapTo2pi(degTorad(trueWindAngle))
trueWindSpeed              = 7
rho                        = 1
nuAir                      = 15.6*(10**(-6))

#Main wing
MWChord                    = 0.74
MWSpan                     = 2.785
aspectRatioMW              = MWChord/MWSpan
MWArea                     = 0.1184*MWSpan
constPartWindForceMW       = (1/2)*rho*trueWindSpeed**2*MWArea
MWThickness                = 0.16
MWDesignedLiftCoefficient  = 2*np.pi*aspectRatioMW/(aspectRatioMW+2)
print(MWDesignedLiftCoefficient)
dontKnowHowToName          = (aspectRatioMW+1)/(aspectRatioMW+2)
constPartWindForceMW       = (1/2)*rho*trueWindSpeed**2*MWArea
distanceMWB                = 0.11
distanceMWA                = -0.10
MWMass                     = 11


#Servo wing
SWChord                   = 0.3
SWSpan                    = 0.8
aspectRatioSW             = SWChord/SWSpan
SWThickness               = 0.055
SWArea                    = SWThickness*SWSpan
distanceTail              = -0.812 #meter
maxTailAngle              = 26
constPartWindForceSW      = (1/2)*rho*trueWindSpeed**2*SWArea
SWDesignedLiftCoefficient = 2*np.pi*aspectRatioSW/(aspectRatioSW+2)
SWMass                    = 2



ReServo                   = reynolds(trueWindSpeed, SWChord, nuAir)
thickOnChordServo         = SWThickness/SWChord
CoeffDragLamServo         = 1.328/np.sqrt(ReServo)
CoeffDragServo            = 2*CoeffDragLamServo*(1+thickOnChordServo)+thickOnChordServo**2

ReMW                      = reynolds(trueWindSpeed, MWChord, nuAir)
thickOnChordMW            = MWThickness/MWChord
CoeffDragLamMW            = 1.328/np.sqrt(ReMW)
CoeffDragMW               = 2*CoeffDragLamMW*(1+thickOnChordMW)+thickOnChordMW**2



momentOfInertiaStructureOnMast  = SWMass*distanceTail**2+MWMass*0.2**2
#momentOfInertiaStructureOnMast  = (np.pi/2)*1.9*(Rext-Rint)**4*MWSpan + (np.pi/2)*1.7*(Rext-Rint)**4*SWSpan
#momentOfInertiaStructureOnMast  = (np.pi/2)*1.9*(35.5-33.5)**4*MWSpan + (np.pi/2)*1.7*(12.3-10.3)**4*SWSpan


#print(momentOfInertiaStructureOnMast)


MWAngle               = float(input('Enter main wing angle in degrees:'))
MWAngle               = wrapTo2pi(degTorad(MWAngle))
dotMWAngle            = wrapTo2pi(degTorad(0))
SWAngle               = float(input('Enter tail angle in degrees (range =-+26):'))
SWAngle               = wrapTo2pi(degTorad(SWAngle))






####################################################################################################
#                                             Physics                                              #
####################################################################################################

def angleOfAttackMW(MWAngle,trueWindAngle):
    #print(MWAngle,trueWindAngle)
    return(wrapTo2pi(trueWindAngle-MWAngle-np.pi))




def angleOfLiftForceMW(state,trueWindAngle):
    alpha = angleOfAttackMW(state[0][0],trueWindAngle)
    return(wrapTo2pi(trueWindAngle-np.pi/2))



def aerodynamicForcesCFD(alpha,tailAngle):
    """ The force are returned in the wind coordinate system we will use symetry for negative angles"""
    liftForceMW = 200.162*abs(wrapTo2pi(alpha)) # experimental formula

    liftForceSW = 17.677*abs(wrapTo2pi(alpha)) +12.67*abs(wrapTo2pi(tailAngle)) # experimental formula 

    dragForceMW = -101.086*abs(wrapTo2pi(alpha))

    dragForceSW = -10.18*abs(wrapTo2pi(alpha)) - 4.54*abs(tailAngle)
    if alpha > 0:
        liftForceMW = -liftForceMW
    if wrapTo2pi(alpha+tailAngle) > 0:
        liftForceSW = -liftForceSW
    return (liftForceMW, liftForceSW, dragForceMW, dragForceSW)


def aerodynamicForcesTheory(alpha,tailAngle):
    liftForceMW = constPartWindForceMW*MWDesignedLiftCoefficient*abs(wrapTo2pi(alpha))*5.91
    #print (constPartWindForceMW*MWDesignedLiftCoefficient*5.91)

    liftForceSW = constPartWindForceSW*SWDesignedLiftCoefficient*abs(wrapTo2pi(dontKnowHowToName*alpha-tailAngle))*5.91
    #print(5.91*constpartwindforcesw*swdesignedliftcoefficient, "wr", wrapto2pi(dontknowhowtoname*alpha))
    dragForceMW = -constPartWindForceMW*MWDesignedLiftCoefficient*abs(wrapTo2pi(alpha))**2*5.91/2

    dragForceSW = -constPartWindForceSW*SWDesignedLiftCoefficient*abs(wrapTo2pi(dontKnowHowToName*alpha-tailAngle))**2*5.91/2
    
    if alpha > 0:
        #print(alpha)
        liftForceMW= -liftForceMW
    if wrapTo2pi(alpha-tailAngle) > 0:
        #print(alpha+tailAngle)
        liftForceSW = -liftForceSW
    
    return (liftForceMW, liftForceSW, dragForceMW, dragForceSW)

def windCoordinatesToMWCoordinates(liftForce,dragForce,alpha):
    y = liftForce*np.cos(alpha)-dragForce*np.sin(alpha)
    x = liftForce*np.sin(alpha)+dragForce*np.cos(alpha)
    return (x,y)
    

def equationθpp(state, truewindAngle,positionAerodynamicCenter ='A', aerodynamicForces = 'T' ):
    alpha      = angleOfAttackMW(state[0][0],truewindAngle)
    #print(alpha)
    tailAngle  = state[2][0]
    #print (tailAngle)
    if aerodynamicForces == 'T':
        liftForceMW,liftForceSW,dragForceMW,dragForceSW = aerodynamicForcesTheory(alpha,tailAngle)
    else:
        liftForceMW,liftForceSW,dragForceMW,dragForceSW = aerodynamicForcesCFD(alpha,tailAngle)

    #print('lMW',liftForceMW, 'dMW',dragForceMW)
    #print('lSW',liftForceSW, 'dSW',dragForceSW)
    
    xMW_TailRelatedPart,yMW_TailRelatedPart = windCoordinatesToMWCoordinates(liftForceSW,dragForceSW,alpha)
    tailRelatedMoment   = yMW_TailRelatedPart*distanceTail
    #print('moment tail',tailRelatedMoment) 
    tailRelatedPart     = tailRelatedMoment/momentOfInertiaStructureOnMast
    #print('tail part',tailRelatedPart) 
    if positionAerodynamicCenter == 'O':
        return(tailRelatedPart)
    else: 
        yMW_MWRelatedPart = liftForceMW*np.cos(alpha)-dragForceMW*np.sin(alpha)
        #xMW_MWRelatedPart = liftForceMW*np.sin(alpha)+dragForceMW*np.cos(alpha)

        if positionAerodynamicCenter =='A':
            MWRelatedMoment = yMW_MWRelatedPart*distanceMWA
        else :
            MWRelatedMoment = yMW_MWRelatedPart*distanceMWB
        #print('moment MW',MWRelatedMoment)
        MWRelatedPart     = MWRelatedMoment/momentOfInertiaStructureOnMast
        #print('part MW',MWRelatedPart)
        #print('tailRelatedPart: ',tailRelatedPart,'MWRelatedPart: ', MWRelatedPart)
        return(tailRelatedPart+MWRelatedPart)
    


global wingState
wingState = np.array([MWAngle,dotMWAngle,SWAngle])
wingState = wingState.reshape((3,1))
#print(wingState)

###########################State equation of the angle of the wingsail""############################

def evolution(X,t,trueWindAngle=0,positionAerodynamicCenter='A', AeroDynamicForces = 'T'):
    Xdot      = [0]*3
    rest      = equationθpp(X, trueWindAngle,positionAerodynamicCenter, AeroDynamicForces)
    Xdot[0]   = wrapTo2pi(X[1][0])
    Xdot[1]   = rest
    Xdot      = np.array([[wrapTo2pi(Xdot[0])],[wrapTo2pi(Xdot[1])],[Xdot[2]]])
    Xdot      = Xdot.reshape((3,1))
#   print('Xdot: ', Xdot)
    return(Xdot)



####################################################################################################
#                          Parameters to determin equilibrium position                             #
####################################################################################################


def equilibriumAngleForDifferentStartAngles(stop = 200,trueWindAngle = 0, positionAerodynamicCenter = 'A', aerodynamicForces = 'T'):
    equilibriumAnglesFDSA = [] #FDSA = ForDifferentStartAngles
    for i in range(-30,31):
        print(i)
        startTime     = 0
        stopTime      = stop
        stepSize      = 0.01
        state         = [0]*3
        state[0]      = wrapTo2pi(degTorad(i))
        state         = np.array([state])
        state         = state.reshape((3,1))
        times, states = rk2Scheme(state, startTime, stopTime, stepSize, evolution, trueWindAngle, positionAerodynamicCenter,aerodynamicForces)
        states        = states[0][-3:]
        equilibriumAnglesFDSA.append(np.mean(states))
    return(equilibriumAnglesFDSA)


def evolutionMWAngleForDifferentStartAngles(wide = (-30,-19),stop = 200,trueWindAngle = 0, positionAerodynamicCenter = 'A', aerodynamicForces = 'T'):
    evolutionAnglesFDSA = [] #FDSA = ForDifferentStartAngles
    times =[]
    for i in range(wide[0],wide[1]):
        print(i)
        startTime     = 0
        stopTime      = stop
        stepSize      = 0.01
        state         = [0]*3
        state[0]      = wrapTo2pi(degTorad(i))
        state         = np.array([state])
        state         = state.reshape((3,1))
        times, states = rk2Scheme(state, startTime, stopTime, stepSize, evolution, trueWindAngle, positionAerodynamicCenter,aerodynamicForces)
        evolutionAnglesFDSA.append(states)
    return(evolutionAnglesFDSA,times, wide)



""" Results will mean nothing until the problem of how to implement angles is fixed """
def equilibriumAngleForDifferentTailAngles(stop = 200, trueWindAngle = 0, positionAerodynamicCenter = 'A', aerodynamicForces = 'T'):
    equilibriumAnglesFDTA = [] #FDTA = ForDifferentTailAngles
    for i in range (-maxTailAngle,maxTailAngle+1):
        state        = [0]*3
        state[2]     = wrapTo2pi(degTorad(i))
        state        = np.array([state])
        state        = state.reshape((3,1))
        startTime    = 0
        stopTime     = stop
        stepSize     = 0.01
        times,states = rk2Scheme(state, startTime, stopTime, stepSize, evolution, trueWindAngle, positionAerodynamicCenter,aerodynamicForces)
        states       = states[0][-3:]
#       print(states)
        equilibriumAnglesFDSWA.append(np.mean(states))
    return(equilibriumAnglesFDSWA)

# Needs to be fixed
def evolutionMWAngleForDifferentTailAngles(wide = (-26,-15), stop = 300,trueWindAngle = 0, positionAerodynamicCenter = 'A', aerodynamicForces = 'T'):
    evolutionAnglesFDTA = [] #FDSA = ForDifferentStartAngles
    times = []
    for i in range(wide[0],wide[1]):
        print(i)
        startTime     = 0
        stopTime      = stop
        stepSize      = 0.01
        state         = [0]*3
        state[2]      = wrapTo2pi(degTorad(i))
        state         = np.array([state])
        state         = state.reshape((3,1))
        #print(state)
        times, states = rk2Scheme(state, startTime, stopTime, stepSize, evolution, trueWindAngle, positionAerodynamicCenter,aerodynamicForces)
        evolutionAnglesFDTA.append(states)
    return(evolutionAnglesFDTA,times,wide)

####################################################################################################
#                       Determin the Lift in function of the angle of attack                       #
####################################################################################################

def lift_dragInfunctionOfAlpha():
    theLifts  = []
    theAngles = []
    theDrags  = []
    for alpha in range (-26,27):
        liftForceMW,useless1, dragForceMW, useless2 = aerodynamicForcesTheory(degTorad(alpha),0)
        
        theLifts.append((liftForceMW))
        theDrags.append((dragForceMW))
        theAngles.append(alpha)
    plt.figure()
    plt.plot(theAngles, theLifts, label = 'Lift')
    plt.plot(theAngles, theDrags, label = 'Drag')
    plt.legend()
    plt.title('Lift and drag of the main wing in function of the angle of attack')
    return(theAngles, theLifts, theDrags)



   
    


 
    
