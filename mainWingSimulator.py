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

"""
UNCOMMENT MAIN BLOC AT THE END OF THE FILE FOR TESTS
"""



X = 2

"""Constants"""
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
    if alpha < 0:
        liftForceMW = -liftForceMW
    if wrapTo2pi(alpha+tailAngle) < 0:
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
    if wrapTo2pi(alpha+tailAngle) > 0:
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



if __name__ == '__main__':

    """
    dt                          = 0.01
    t                           = 0
    fig                         = plt.figure()
    positionOfAerodynamicCenter = (str(input('Enter position of aerodynamic center B/O/A:'))).upper()
    aerodynamicForces           = (str(input('Enter type of aerodynamic forces   (theory)T/CFD:'))).upper()
    ax                          = fig.add_subplot(111, aspect ='equal')
    while t+dt <=1700:
        t                                 = t+dt
        plt.cla()
        ax.set_xlim(-15,15)
        ax.set_ylim(-15,15)
##      print(t)
##      wingState = wingState + dt*evolution(wingState,1)
##      wingState[0],wingState[1] = wrapTo2pi(wingState[0]),wrapTo2pi(wingState[1]
        wingState                         = wingState + dt*rk2Step(wingState,t,dt,evolution,trueWindAngle,positionOfAerodynamicCenter)
        wingState[0][0]                   = wrapTo2pi(wingState[0][0])
        wingState[1][0]                   = wrapTo2pi(wingState[1][0])
##      print('main wing angle: ', wingState[0][0])
##      print('true wind angle: ', trueWindAngle)
        alpha                             = angleOfAttackMW(wingState[0],trueWindAngle)
        liftMW, useless1,dragMW,useless3  = aerodynamicForcesTheory(alpha, SWAngle)
        angleTotalForceMW                 = float(np.arctan2(liftMW,dragMW))
        moduleTotalForceMW                = float(np.sqrt(liftMW**2+dragMW**2))
        #print(radTodeg(angleLift), lift)
        drawArrow(0,0,angleTotalForceMW,moduleTotalForceMW,'k')
        drawWingSailIRT(wingState[0][0],wingState[2][0],trueWindAngle,trueWindSpeed)
        drawHull(0)
        plt.pause(0.0001)
    plt.show()
    """


    # ===== Equilibrium position in function of the position of the aerodynamic center (Theory forces) ===== #
    
    """
    #plt.figure()
    
    timesRK2,statesRK2O = rk2Scheme(wingState,0,300, 0.01, evolution, trueWindAngle,'O')
    anglesInDegreesRK2O = listRadTodeg(statesRK2O[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2O,'on rot axis','RK2')
    """
    
    """
    timesRK2,statesRK2B = rk2Scheme(wingState,0,300, 0.01, evolution, trueWindAngle,'B')
    anglesInDegreesRK2B = listRadTodeg(statesRK2B[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2B,'before rot axis','RK2')
    """

    
    """
    timesEuler,statesEulerA = eulerScheme(wingState,0,300,0.01,evolution,trueWindAngle)
    anglesInDegreesEulerA = listRadTodeg(statesEulerA[0,:])
    drawWingSailAngle(timesEuler,anglesInDegreesEulerA,'euler','fds')
    """
    """
    timesRK2,statesRK2A = rk2Scheme(wingState,0,300, 0.01, evolution, trueWindAngle)
    anglesInDegreesRK2A = listRadTodeg(statesRK2A[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2A,'after rot axis','RK2 (theoretical forces)')
    
    plt.legend()
    plt.show()
    """
    #plt.savefig('Simulation_pics/Comparison position aerodynamic center for theoretical aerodynamic forces.png')
    
    

    # ==== Equilibrium position in function of the position of the aerodynamic center (Experimental forces) ===== #
    
    """
    plt.figure()
    timesRK2,statesRK2O = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'O','CFD')  
    anglesInDegreesRK20 = listRadTodeg(statesRK2O[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK20,'on rot axis','RK2')
    """
    
    """
    timesRK2,statesRK2B = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'B','CFD')
    anglesInDegreesRK2B = listRadTodeg(statesRK2B[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2B,'before rot axis','RK2')
    """
    
    """
    timesRK2,statesRK2A = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'A','CFD')
    anglesInDegreesRK2A = listRadTodeg(statesRK2A[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2A,'after rot axis','RK2 (experimental forces)')
    plt.legend()
    plt.show()
    """
    
    """
    #plt.savefig('Simulation_pics/Comparison position aerodynamic center for experimental aerodynamic forces2.png')
    """


    
    # ===== Equilibrium position in function of the forces (aerodynamic center behind the mast) ===== #

    """
    plt.figure()
    timesRK2,statesRK2A = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'A','CFD')
    drawWingSailAngle(timesRK2,statesRK2A[0,:],'experimental forces','RK2')
    timesRK2,statesRK2A = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle)
    drawWingSailAngle(timesRK2,statesRK2A[0,:],'theoretical forces','RK2')
     plt.legend()
    plt.savefig('Simulation_pics/Comparison type of forces aerodynamic center behind the mast.png')
    """



    # ===== Evolution of the MW angle for different Start angles ===== #
    """
    wide = (-9,1)
    evolutionAnglesFDSA,times, wide = evolutionMWAngleForDifferentStartAngles(wide,200,trueWindAngle,'A', 'T') 
    drawEvolutionMWAngle(evolutionAnglesFDSA,times,wide)
    plt.show()
    #plt.savefig('Simulation_pics/evolution of main wing angles for different start angles theoretical forces_aerodynamic center behind the mast')
    """
    
    """
    wide = (11,21)
    evolutionAnglesFDSA, times, wide = evolutionMWAngleForDifferentStartAngles(wide,200,trueWindAngle,'A', 'CFD') 
    drawEvolutionMWAngle(evolutionAnglesFDSA,times,wide)
    plt.show()
    """
    
    #plt.savefig('Simulation_pics/evolution of main wing angles for different start angles experimental forces_aerodynamic center behind the mast')


    
    # ===== Equilibrium position for different Start angles ===== #

    # Warning!! run the previous section before in order to put a reasonnable stop point
    
    """
    equilibriumAnglesFDSA  = equilibriumAngleForDifferentStartAngles(200,trueWindAngle,'A','T')
    drawEquilibriumAngles( equilibriumAnglesFDSA)
    plt.savefig('Simulation_pics/equilibrium angles for different start angles theoretical forces_aerodynamic center behind the mast.png')
    """
    
    """
    equilibriumAnglesFDSA  = equilibriumAngleForDifferentStartAngles(200,trueWindAngle,'A','CFD')
    drawEquilibriumAngles( equilibriumAnglesFDSA)
    plt.savefig('Simulation_pics/equilibrium angles for different start angles experimental forces_aerodynamic center behind the mast.png')
    """


    
    # =====  Evolution of the MW angle for different tail angles ===== #

    """
    wide = (-14,0) 
    evolutionAnglesFDTA,times, wide = evolutionMWAngleForDifferentTailAngles(wide, 200,trueWindAngle,'A', 'T') 
    drawEvolutionMWAngle(evolutionAnglesFDTA,times, wide,'FDTA')
    plt.show()
    """
    
    #plt.savefig('Simulation_pics/evolution of main wing angles for different tail  angles theoretical forces_aerodynamic center behind the mast')
    

    """
    wide = (-14,0)
    evolutionAnglesFDTA,times, wide = evolutionMWAngleForDifferentTailAngles(wide, 200,trueWindAngle,'A', 'CFD') 
    drawEvolutionMWAngle(evolutionAnglesFDTA,times,wide,'FDTA')
    plt.show()
    """
    #plt.savefig('Simulation_pics/evolution of main wing angles for different tail angles experimental forces_aerodynamic center behind the mast')
    

    
    # ===== Equilibrium position for different servo wing angles (start at 0) ===== #

    """
    plt.figure()
    equilibriumAnglesFDTA = equilibriumAngleForDifferentTailAngles(200,trueWindAngle,'A','T')
    drawEquilibriumAngles( equilibriumAnglesFDTA,'FDTA')
    plt.savefig('Simulation_pics/equilibrium angles for different tail angles theoretical forces_aerodynamic center behind the mast.png')
    """

    # ===== Lift force on main wing in function of the angle of attack ===== #
    """
    #In wind coordinate system : 
    angles , lifts, drags = lift_dragInfunctionOfAlpha()
    
    #print(lifts,drags)
    print(angles)
    # in MW coordinate system
    xMW_forces, yMW_forces = [], []
    for i in range (len(angles)):
        Xforce, Yforce = apparentWindCoordinatesToMWCoordinates(drags[i], lifts[i], degTorad(angles[i]))
        xMW_forces.append(Xforce)
        yMW_forces.append(Yforce)
    plt.figure()
    plt.plot(angles,xMW_forces,label='xMW_force')
    plt.plot(angles,yMW_forces,label='yMW_force')
    plt.legend()
    #print(xMW_forces,yMW_forces)
    
    
    # in boat coordinate system
    MWAngle = int(input('Angle main wing: '))
    xBoat_forces, yBoat_forces = [], []
    for i in range (len(xMW_forces)):
        Xforce, Yforce = apparentWindCoordinatesToMWCoordinates(xMW_forces[i], yMW_forces[i], MWAngle)
        xBoat_forces.append(Xforce)
        yBoat_forces.append(Yforce)
    plt.figure()
    plt.plot(angles,xBoat_forces,label='xBoat_force')
    plt.plot(angles,yBoat_forces,label='yBoat_force')
    plt.legend()
    plt.title('forces in Boat coordinate system with a MW angle of ' + str(MWAngle)+'°')

    plt.show()
    """
    # =====  calculate best orientation for main wing ===== #

    #In wind coordinate system : 
    angles , lifts, drags = lift_dragInfunctionOfAlpha()

    #apparent wind angle in boat coordinate system
    phi = int(input('Give angle of the apparent wind in the boat coordinate system: ')) 
    phi = wrapTo2pi(degTorad(phi))

    # calculating forces in boat coordinate system
    xBoat_forces,yBoat_forces = [],[]
    for i in range(len(lifts)):
        xBoat_force = drags[i]*np.cos(phi)-lifts[i]*np.sin(phi)
        yBoat_force = lifts[i]*np.cos(phi)+drags[i]*np.sin(phi)
        xBoat_forces.append(xBoat_force)
        yBoat_forces.append(yBoat_force)

    maxi_xForce = max(xBoat_forces)
    index_maxi  = xBoat_forces.index(maxi_xForce)
    print('maxi_xForce: ', maxi_xForce)
    print('index_maxi: ', index_maxi-26) # best angle of attack in degrees if we want max(xForce)
    
    plt.figure()
    plt.plot(angles,xBoat_forces,label = 'xboat_force')
    plt.plot(angles,yBoat_forces,label = 'yboat_force')
    plt.legend()
    
    # calculate the thrust/drift ratio
    xForcesOnyForces = []
    for i in range (len(xBoat_forces)):
        if yBoat_forces[i] !=0:
            xOny = xBoat_forces[i]/abs(yBoat_forces[i]) # if x < 0 useless because it pushes the boat backwards
            xForcesOnyForces.append(xOny)
        else:
            xForcesOnyForces.append(0)

    maxi_xOny_ratio = max(xForcesOnyForces)
    index_maxi_xOny = xForcesOnyForces.index(maxi_xOny_ratio)
    maxi_xOny       = xBoat_forces[index_maxi_xOny]

    print('maxi_xOny: ',maxi_xOny)
    print('index_maxi_xOny: ', index_maxi_xOny-26) # best angle of attack if we want the greatest thrust/drift ratio
    
    plt.figure ()
    plt.plot(angles, xForcesOnyForces, label = 'xForcesOnyForces')
    plt.legend()

    # return best angle to go where we want
    
    minXforce = 2 # newtons
    MWDeviation = index_maxi_xOny-26 # degrees
    if maxi_xOny < minXforce :
        MWDeviation = index_maxi-26

    UMWAngle = radTodeg(phi) -  MWDeviation #degrees, best angle for the main wing
    Utail    = MWDeviation                  # command in degrees to give to the tail
    print('Control MW angle: ' , UMWAngle)
    print('control tail angle: ', Utail)
    plt.show()
    
    
