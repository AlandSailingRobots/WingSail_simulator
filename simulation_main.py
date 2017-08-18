import numpy as np
from drawingFunctions import drawWingSailAngle, drawWingSailIRT, drawEquilibriumAngles, drawArrow, drawEvolutionMWAngle, drawHull
from angleFunctions import wrapTo2pi, degTorad, radTodeg, listRadTodeg, listDegTorad, listWrapTo2pi, apparentWindCoordinatesToMWCoordinates, MWCoordinatesToBoatCoordinates
import matplotlib.pyplot as plt
from WingSimulator import *

"""
UNCOMMENT MAIN BLOC AT THE END OF THE FILE FOR TESTS
"""


def simulate_IRT():
    """ This function draws the state of the wing sail in real time """
    global wingState  
    global trueWindAngle  
    dt                          = 0.01 # step size for numerical integration
    t                           = 0    # start time 
    fig                         = plt.figure()
    positionOfAerodynamicCenter = (str(input('Enter position of aerodynamic center B/O/A:'))).upper() 
    aerodynamicForces           = (str(input('Enter type of aerodynamic forces   (theory)T/CFD:'))).upper()
    ax                          = fig.add_subplot(111, aspect ='equal')
    while t+dt <=1700:
        t                                 = t+dt
        plt.cla()
        ax.set_xlim(-15,15)
        ax.set_ylim(-15,15)
        wingState                         = wingState + dt*rk2Step(wingState,t,dt,evolution,trueWindAngle,positionOfAerodynamicCenter)
        wingState[0][0]                   = wrapTo2pi(wingState[0][0])
        wingState[1][0]                   = wrapTo2pi(wingState[1][0])
        alpha                             = angleOfAttackMW(wingState[0],trueWindAngle)

        # We need the lift and the drag on the main wing to draw the resulstant of the force 
        liftMW, useless1,dragMW,useless3  = aerodynamicForcesTheory(alpha, SWAngle) 

        # calculating the angle of the resultant
        angleTotalForceMW                 = float(np.arctan2(liftMW,dragMW))
        # calculating the module of the resultant
        moduleTotalForceMW                = float(np.sqrt(liftMW**2+dragMW**2))

        # drawing resultant
        drawArrow(0,0,angleTotalForceMW,moduleTotalForceMW,'k')
        # drawing the wing sail
        drawWingSailIRT(wingState[0][0],wingState[2][0],trueWindAngle,trueWindSpeed)
        # drawing a fake hull for display
        drawHull(0)
        plt.pause(0.0001)
    plt.show()
    return()



# ===== Equilibrium position in function of the position of the aerodynamic center (Theory forces) ===== #


def calculationsTheoryForces():
    """ This functions is here to allow automatic tests of the wing sail with different center 
    of pressure postition with the forces of the model obtained with theory. You obtain a curve with the evolution of the angle of the main wing in the hull 
    coordinate system""" 
    global wingState
    global trueWindAngle 
    plt.figure()

    # center of pressure on the mast
    timesRK2,statesRK2O = rk2Scheme(wingState,0,300, 0.01, evolution, trueWindAngle,'O')
    anglesInDegreesRK2O = listRadTodeg(statesRK2O[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2O,'on rot axis','RK2')
       
    # center of pressure before the mast   
    timesRK2,statesRK2B = rk2Scheme(wingState,0,300, 0.01, evolution, trueWindAngle,'B')
    anglesInDegreesRK2B = listRadTodeg(statesRK2B[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2B,'before rot axis','RK2')

    #  center of pressure after the mast
    timesRK2,statesRK2A = rk2Scheme(wingState,0,300, 0.01, evolution, trueWindAngle)
    anglesInDegreesRK2A = listRadTodeg(statesRK2A[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2A,'after rot axis','RK2 (theoretical forces)')
    
    plt.legend()
    #plt.savefig('Comparison position aerodynamic center for theoretical aerodynamic forces.png')
    plt.show()
    return()

# ==== Equilibrium position in function of the position of the aerodynamic center (Experimental forces) ===== #

def calculationsCFDForces():
    """ This function is here to allow automatic tests of the wing sail with different center 
    of pressure postition with the forces of the model obtained with CFD. You obtain a curve with the evolution of the angle of the main wing in the hull 
    coordinate system""" 
    global wingState
    global trueWindAngle 
    plt.figure()
    
    # center of pressure on the mast
    timesRK2,statesRK2O = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'O','CFD')  
    anglesInDegreesRK20 = listRadTodeg(statesRK2O[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK20,'on rot axis','RK2')
    
    # center of pressure before the mast 
    timesRK2,statesRK2B = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'B','CFD')
    anglesInDegreesRK2B = listRadTodeg(statesRK2B[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2B,'before rot axis','RK2')
    
    #  center of pressure after the mast
    timesRK2,statesRK2A = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'A','CFD')
    anglesInDegreesRK2A = listRadTodeg(statesRK2A[0,:])
    drawWingSailAngle(timesRK2,anglesInDegreesRK2A,'after rot axis','RK2 (experimental forces)')
    
    plt.legend()
    #plt.savefig('Comparison position aerodynamic center for experimental aerodynamic forces.png')
    plt.show()
    return()

# ===== Evolution of the MW angle for different Start angles ===== #

def calculationsEvolutionForDifferentStartAngles(wide = (-10,10), positionOfAerodynamicCenter = 'A',forcesType ='T'):
    """ This function calculates the evolution of the main wing angle in function of the angle you have at the start at the simulation
    mainly used to check stability"""
    global trueWindAngle 
    evolutionAnglesFDSA,times, wide = evolutionMWAngleForDifferentStartAngles(wide,200,trueWindAngle,positionOfAerodynamicCenter, forcesType) 
    drawEvolutionMWAngle(evolutionAnglesFDSA,times,wide)
    plt.show()
    #if forcesType == 'T':
        #plt.savefig('evolution of main wing angles for different start angles theoretical forces_aerodynamic center behind the mast')
    #else:
        #plt.savefig('evolution of main wing angles for different start angles experimental forces_aerodynamic center behind the mast')
    return()


# ===== Equilibrium position for different Start angles ===== #

# Warning!! run the previous section before in order to put a reasonnable stop point

def calculationsEquilibriumForDifferentStartAngles(stopTime=200, positionOfAerodynamicCenter = 'A',forcesType ='T'):   
    """ This function returns the equilibrium position of the main wing angle in function of the angle you have at the start at the simulation
    mainly used to check stability"""
    global trueWindAngle 
    equilibriumAnglesFDSA  = equilibriumAngleForDifferentStartAngles(stopTime,trueWindAngle,positionOfAerodynamicCenter, forcesType)
    drawEquilibriumAngles( equilibriumAnglesFDSA)
    plt.show()
    #if forcesType == 'T':
        #plt.savefig('equilibrium angles for different start angles theoretical forces_aerodynamic center behind the mast.png')
    #else:
        #plt.savefig('equilibrium angles for different start angles experimental forces_aerodynamic center behind the mast.png')
    return()


# =====  Evolution of the MW angle for different tail angles ===== #
def calculationsEvolutionForDifferentTailAngles(wide =(-10,10), positionOfAerodynamicCenter = 'A', forcesType = 'T'):    
    """ This function returns the evolution of the main wing angle in function of the angle you have for the tail wing"""
    global trueWindAngle 
    evolutionAnglesFDTA,times, wide = evolutionMWAngleForDifferentTailAngles(wide, 200,trueWindAngle,positionOfAerodynamicCenter, forcesType) 
    drawEvolutionMWAngle(evolutionAnglesFDTA,times, wide,'FDTA')
    plt.show()
    #if forcesType == 'T':
        #plt.savefig('evolution of main wing angles for different tail  angles theoretical forces_aerodynamic center behind the mast')
    #else:
        #plt.savefig('evolution of main wing angles for different tail angles experimental forces_aerodynamic center behind the mast')
    return()

# ===== Equilibrium position for different servo wing angles (start at 0) ===== #

# Warning!! run the previous section before in order to put a reasonnable stop point

def calculationsEquilibriumForDifferentTailAngles(stopTime=200, positionOfAerodynamicCenter = 'A',forcesType ='T'):   
    """ This function returns the evolution of the main wing angle in function of the angle you have for the tail wing"""
    global trueWindAngle
    equilibriumAnglesFDTA = equilibriumAngleForDifferentTailAngles(stopTime,trueWindAngle,positionOfAerodynamicCenter,forcesType)
    drawEquilibriumAngles( equilibriumAnglesFDTA,'FDTA')
    #if forcesType =='T':
        #plt.savefig('equilibrium angles for different tail angles theoretical forces_aerodynamic center behind the mast.png')
    #else:
        #plt.savefig('equilibrium angles for different tail angles experimental forces_aerodynamic center behind the mast.png')
    return()

    # ===== Lift force on main wing in function of the angle of attack ===== #
def liftAndDragForceOnMWInFunctionOfAlpha():
    """ This function returns two lists with the drag and the lift (and create and display a graph) for different position of the main 
    wing regarding the wind"""
    #In wind coordinate system 
    angles , lifts, drags = lift_dragInfunctionOfAlpha()
    
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

    plt.print('lifts: ',lift, 'drags: ',drags)
    plt.show()
    return(lifts,drags)

   # =====  calculate best orientation for main wing ===== #
def calculateTailOrder():
    """ In this function the first 3 inputs are not useful you can type any number only
    'Give angle of the apparent wind in the boat coordinate system: ' is used """

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
    return(Utail)

if __name__ == '__main__':

    #simulate_IRT()
    #calculationsTheoryForces()
    #calculateTailOrder()