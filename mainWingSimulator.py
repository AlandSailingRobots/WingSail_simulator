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
from angleFunctions import wrapTo2pi, degTorad, radTodeg, listRadTodeg, listDegTorad, listWrapTo2pi
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
dontKnowHowToName          = (aspectRatioMW+1)/(aspectRatioMW+2)
constPartWindForceMW       = (1/2)*rho*trueWindSpeed**2*MWArea
distanceMWB                = 0.11
distanceMWA                = -0.02
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


print(momentOfInertiaStructureOnMast)


angleMaxLiftCoeff     = degTorad(15)
MWAngle               = float(input('Enter main wing angle in degrees:'))
MWAngle               = wrapTo2pi(degTorad(MWAngle))
dotMWAngle            = wrapTo2pi(degTorad(0))
SWAngle               = float(input('Enter tail angle in degrees (range =-+26):'))
SWAngle               = wrapTo2pi(degTorad(SWAngle))






####################################################################################################
#                                             Physics                                              #
####################################################################################################

def angleOfAttackMW(MWAngle,trueWindAngle):
	return(-wrapTo2pi(trueWindAngle-MWAngle-np.pi))




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
		liftForceMW, liftForceSW = -liftForceMW, -liftForceSW
	return (liftForceMW, liftForceSW, dragForceMW, dragForceSW)


def aerodynamicForcesTheory(alpha,tailAngle):
	liftForceMW = constPartWindForceMW*MWDesignedLiftCoefficient*abs(wrapTo2pi(alpha))*5.91
	#print (constPartWindForceMW*MWDesignedLiftCoefficient*5.91)

	liftForceSW = constPartWindForceSW*SWDesignedLiftCoefficient*abs(wrapTo2pi(dontKnowHowToName*alpha-tailAngle))*5.91
	#print(5.91*constpartwindforcesw*swdesignedliftcoefficient, "wr", wrapto2pi(dontknowhowtoname*alpha))
	dragForceMW = 0

	dragForceSW = 0
	if alpha < 0:
		liftForceMW, liftForceSW = -liftForceMW, -liftForceSW
	return (liftForceMW, liftForceSW, dragForceMW, dragForceSW)




def equationθpp(state, truewindAngle,positionAerodynamicCenter ='A', aerodynamicForces = 'T' ):
	alpha      = angleOfAttackMW(state[0][0],truewindAngle)
	tailAngle  = state[2][0]
	#print (tailAngle)
	if aerodynamicForces == 'T':
		liftForceMW,liftForceSW,dragForceMW,dragForceSW = aerodynamicForcesTheory(alpha,tailAngle)
	else:
		liftForceMW,liftForceSW,dragForceMW,dragForceSW = aerodynamicForcesCFD(alpha,tailAngle)

	yMW_TailRelatedPart = liftForceSW*np.cos(alpha)-dragForceSW*np.sin(alpha)
	#xMW_TailRelatedPart = liftForceSW*np.sin(alpha)+dragForceSW*np.cos(alpha)
	tailRelatedMoment   = yMW_TailRelatedPart*distanceTail
	tailRelatedPart     = tailRelatedMoment/momentOfInertiaStructureOnMast

	if positionAerodynamicCenter == 'O':
		return(tailRelatedPart)
	else: 
		yMW_MWRelatedPart = liftForceMW*np.cos(alpha)-dragForceMW*np.sin(alpha)
		#xMW_MWRelatedPart = liftForceMW*np.sin(alpha)+dragForceMW*np.cos(alpha)

		if positionAerodynamicCenter =='A':
		    MWRelatedMoment = yMW_MWRelatedPart*distanceMWA
		else :
			MWRelatedMoment = yMW_MWRelatedPart*distanceMWB

		MWRelatedPart     = MWRelatedMoment/momentOfInertiaStructureOnMast
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
#	print('Xdot: ', Xdot)
	return(Xdot)



####################################################################################################
#                                  Parameters Testing Methods                                      #
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
		print(i)
		state        = [0]*3
		state[2]     = wrapTo2pi(degTorad(i))
		state         = np.array([state])
		state        = state.reshape((3,1))
		startTime    = 0
		stopTime     = stop
		stepSize     = 0.01
		times,states = rk2Scheme(state, startTime, stopTime, stepSize, evolution, trueWindAngle, positionAerodynamicCenter,aerodynamicForces)
		states       = states[0][-3:]
#		print(states)
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
		stepSiz0      = 0.01
		state         = [0]*3
		state[2]      = wrapTo2pi(degTorad(i))
		state         = np.array([state])
		state         = state.reshape((3,1))
		#print(state)
		times, states = rk2Scheme(state, startTime, stopTime, stepSize, evolution, trueWindAngle, positionAerodynamicCenter,aerodynamicForces)
		evolutionAnglesFDTA.append(states)
	return(evolutionAnglesFDTA,times,wide)



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
##		print(t)
##		wingState = wingState + dt*evolution(wingState,1)
##		wingState[0],wingState[1] = wrapTo2pi(wingState[0]),wrapTo2pi(wingState[1]
		wingState                         = wingState + dt*rk2Step(wingState,t,dt,evolution,trueWindAngle,positionOfAerodynamicCenter)
		wingState[0][0]                   = wrapTo2pi(wingState[0][0])
		wingState[1][0]                   = wrapTo2pi(wingState[1][0])
##		print('main wing angle: ', wingState[0][0])
##		print('true wind angle: ', trueWindAngle)
		alpha                             = angleOfAttackMW(wingState[0],trueWindAngle)
		lift, useless1,useless2,useless3  = aerodynamicForcesCFD(alpha, SWAngle)
		angleLift                         = angleOfLiftForceMW(wingState,trueWindAngle)
		print(radTodeg(angleLift), lift)
		drawArrow(0,0,angleLift,lift/5,'k')
		drawWingSailIRT(wingState[0][0],wingState[2][0],trueWindAngle,trueWindSpeed)
		drawHull(0)
		plt.pause(0.0001)
	plt.show()
	"""


	# ===== Equilibrium position in function of the position of the aerodynamic center (Theory forces) ===== #
	
	"""
	plt.figure()

	timesRK2,statesRK2O = rk2Scheme(wingState,0,00, 0.01, evolution, trueWindAngle,'O')
	anglesInDegreesRK2O = listRadTodeg(statesRK2O[0,:])
	drawWingSailAngle(timesRK2,anglesInDegreesRK2O,'on rot axis','RK2')

	timesRK2,statesRK2B = rk2Scheme(wingState,0,300, 0.01, evolution, trueWindAngle,'B')
	anglesInDegreesRK2B = listRadTodeg(statesRK2B[0,:])
	drawWingSailAngle(timesRK2,anglesInDegreesRK2B,'before rot axis','RK2')
	
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
	
	timesRK2,statesRK2B = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'B','CFD')
	anglesInDegreesRK2B = listRadTodeg(statesRK2B[0,:])
	drawWingSailAngle(timesRK2,anglesInDegreesRK2B,'before rot axis','RK2')
	
	timesRK2,statesRK2A = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'A','CFD')
	anglesInDegreesRK2A = listRadTodeg(statesRK2A[0,:])
	drawWingSailAngle(timesRK2,anglesInDegreesRK2A,'after rot axis','RK2 (experimental forces)')
	plt.legend()
	plt.show()
	"""
	#plt.savefig('Simulation_pics/Comparison position aerodynamic center for experimental aerodynamic forces2.png')
	
	
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
	
	wide = (-9,1)
	evolutionAnglesFDSA,times, wide = evolutionMWAngleForDifferentStartAngles(wide,200,trueWindAngle,'A', 'T') 
	drawEvolutionMWAngle(evolutionAnglesFDSA,times,wide)
	plt.show()
	#plt.savefig('Simulation_pics/evolution of main wing angles for different start angles theoretical forces_aerodynamic center behind the mast')
	
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
	evolutionAnglesFDTA,times, wide = evolutionMWAngleForDifferentTailAngles(wide, 200,trueWindAngle,'A', 'T') 
	drawEvolutionMWAngle(evolutionAnglesFDTA,times, wide,'FDTA')
	plt.savefig('Simulation_pics/evolution of main wing angles for different tail  angles theoretical forces_aerodynamic center behind the mast')
	"""

	"""
	evolutionAnglesFDTA,times, wide = evolutionMWAngleForDifferentTailAngles(wide, 200,trueWindAngle,'A', 'CFD') 
	drawEvolutionMWAngle(evolutionAnglesFDTA,times,wide,'FDTA')
	plt.savefig('Simulation_pics/evolution of main wing angles for different tail angles experimental forces_aerodynamic center behind the mast')
	"""


	
	# ===== Equilibrium position for different servo wing angles (start at 0) ===== #
	"""
	plt.figure()
	equilibriumAnglesFDTA = equilibriumAngleForDifferentTailAngles(200,trueWindAngle,'A','T')
	drawEquilibriumAngles( equilibriumAnglesFDTA,'FDTA')
	plt.savefig('Simulation_pics/equilibrium angles for different tail angles theoretical forces_aerodynamic center behind the mast.png')
	"""
