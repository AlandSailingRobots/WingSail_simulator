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
from angleFunctions import wrapTo2pi, degTorad, radTodeg
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
distanceCanard            = 0.700 #meter
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



angleMaxLiftCoeff     = degTorad(15)
MWAngle               = float(input('Enter main wing angle in degrees:'))
MWAngle               = wrapTo2pi(degTorad(MWAngle))
dotMWAngle            = wrapTo2pi(degTorad(0))
tailAngle             = float(input('Enter tail angle in degrees (range =-+26):'))
tailAngle             = wrapTo2pi(tailAngle)





####################################################################################################
#                                             Physics                                              #
####################################################################################################

def angleOfAttackMW(MWAngle,trueWindAngle):
	return(wrapTo2pi(trueWindAngle-MWAngle-np.pi))




def angleOfLiftForceMW(state,trueWindAngle):
	alpha              = angleOfAttackMW(state[0][0],trueWindAngle)
	if alpha<=0:
		return(wrapTo2pi(-alpha+(np.pi/2)))
	else :
		return(wrapTo2pi(-alpha-(np.pi/2)))



def aerodynamicForcesCFD(alpha,tailAngle):
	liftForceMW = 200.162*alpha # experimental formula

	liftForceSW = 17.677*(wrapTo2pi(alpha-tailAngle)) # experimental formula

	dragForceMW = 101.086*alpha-10.596

	dragForceSW = 0 #28.356*(alpha-tailAngle)**2-5.767*(alpha-tailAngle)+0.068

	return (liftForceMW, liftForceSW, dragForceMW, dragForceSW)


def aerodynamicForcesTheory(alpha,tailAngle):
	liftForceMW = constPartWindForceMW*MWDesignedLiftCoefficient*alpha*5.91
	#print (constPartWindForceMW*MWDesignedLiftCoefficient*5.91)


	liftForceSW = constPartWindForceSW*SWDesignedLiftCoefficient*wrapTo2pi(dontKnowHowToName*alpha-tailAngle)*5.91
	#print(5.91*constpartwindforcesw*swdesignedliftcoefficient, "wr", wrapto2pi(dontknowhowtoname*alpha))
	dragforcemw = 0

	dragforcesw = 0

	return (liftforcemw, liftforcesw, dragforcemw, dragforcesw)






def equationθpp(state, truewindangle,positionAerodynamicCenter ='a', aerodynamicForces = 't' ):
	alpha      = angleOfAttackMW(state[0][0],truewindangle)
	tailangle  = state[2][0]
	if aerodynamicForces == 't':
		liftforcemw,liftForceSW,dragForceMW,dragForceSW = aerodynamicForcesTheory(alpha,tailAngle)  
	else:
		liftForceMW,liftForceSW,dragForceMW,dragForceSW = aerodynamicForcesCFD(alpha,tailAngle)   
	tailRelatedPartLift     = distanceTail*liftForceSW
#	tailRelatedPartDrag     = 0
#	MWRelatedPartDrag = 0
#	Uncomment the following line to add drag to the similation
	tailRelatedPartDrag     = dragForceSW
	tailRelatedPart         = (tailRelatedPartLift*np.cos(alpha)-tailRelatedPartDrag*np.sin(alpha))/momentOfInertiaStructureOnMast
	if positionAerodynamicCenter == 'O':
		return(tailRelatedPart)
	elif positionAerodynamicCenter == 'B':
		MWRelatedPartLift   = distanceMWB*liftForceMW
#		MWRelatedPartDrag   = 0
#		Uncomment the following line to add drag to the similation
		MWRelatedPartDrag   = distanceMWB*dragForceMW
		MWRelatedPart       = (MWRelatedPartLift*np.cos(alpha)-MWRelatedPartDrag*np.sin(alpha))/momentOfInertiaStructureOnMast
		return(tailRelatedPart+MWRelatedPart)
	elif positionAerodynamicCenter == 'A':
		MWRelatedPartLift   = distanceMWA*liftForceMW
#		MWRelatedPartDrag   = 0
#		Uncomment the following line to add drag to the similation
		MWRelatedPartDrag   = distanceMWA*dragForceMW
		MWRelatedPart       = (MWRelatedPartLift*np.cos(alpha)-MWRelatedPartDrag*np.sin(alpha))/momentOfInertiaStructureOnMast
		return(tailRelatedPart+MWRelatedPart)







wingState = np.array([MWAngle,dotMWAngle,tailAngle])
wingState = wingState.reshape((3,1))
#print(wingState)

###########################State equation of the angle of the wingsail""############################

def evolution(X,t,trueWindAngle=0,positionAerodynamicCenter='A', AeroDynamicForces = 'T'):
	Xdot      = [0]*3
	rest      = equationθpp(X, trueWindAngle,positionAerodynamicCenter, AeroDynamicForces)
	Xdot[0]   = wrapTo2pi(X[1][0])
	Xdot[1]   = -rest
	Xdot      = np.array([[wrapTo2pi(Xdot[0])],[wrapTo2pi(Xdot[1])],[Xdot[2]]])
	Xdot      = Xdot.reshape((3,1))
#	print('Xdot: ', Xdot)
	return(Xdot)



####################################################################################################
#                                  Parameters Testing Methods                                      #
####################################################################################################


def equilibriumAngleForDifferentStartAngles(wingState, stop=200,trueWindAngle = 0, positionAerodynamicCenter = 'A', aerodynamicForces = 'T'):
	equilibriumAnglesFDSA = [] #FDSA = ForDifferentStartAngles
	for i in range(-30,31):
		print(i)
		startTime     = 0
		stopTime      = stop
		stepSize      = 0.01
		wingState[0]  = wrapTo2pi(degTorad(i))
		wingState[1]  = 0
		times, states = rk2Scheme(wingState, startTime, stopTime, stepSize, evolution, trueWindAngle, positionAerodynamicCenter,aerodynamicForces)
		states        = states[0][-3:]
		equilibriumAnglesFDSA.append(np.mean(states))
	return(equilibriumAnglesFDSA)


def evolutionMWAngleForDifferentStartAngles(wingState, stop=200,trueWindAngle = 0, positionAerodynamicCenter = 'A', aerodynamicForces = 'T'):
	evolutionAnglesFDSA = [] #FDSA = ForDifferentStartAngles
	times =[]
	for i in range(-30,31):
		print(i)
		startTime     = 0
		stopTime      = stop
		stepSize      = 0.01
		wingState[0]  = wrapTo2pi(degTorad(i))
		wingState[1]  = 0
		times, states = rk2Scheme(wingState, startTime, stopTime, stepSize, evolution, trueWindAngle, positionAerodynamicCenter,aerodynamicForces)
		evolutionAnglesFDSA.append(states)
	return(evolutionAnglesFDSA,times)



""" Results will mean nothing until the problem of how to implement angles is fixed """
def equilibriumAngleForDifferentTailAngles(stop = 200, trueWindAngle = 0, positionAerodynamicCenter = 'A', aerodynamicForces = 'T'):
	equilibriumAnglesFDTA = [] #FDTA = ForDifferentTailAngles
	for i in range (-maxTailAngle,maxTailAngle+1):
		print(i)
		tailAngle    = wrapTo2pi(degTorad(i))
		wingState[0] = 0
		wingState[1] = 0
		wingState[2] = tailAngle
		startTime    = 0
		stopTime     = stop
		stepSize     = 0.01
		times,states = rk2Scheme(wingState, startTime, stopTime, stepSize, evolution, trueWindAngle, positionAerodynamicCenter,aerodynamicForces)
		states       = states[0][-3:]
#		print(states)
		equilibriumAnglesFDSWA.append(np.mean(states))
	return(equilibriumAnglesFDSWA)



if __name__ == '__main__':

	"""
	dt                          = 0.01
	t                           = 0
	fig                         = plt.figure()
	positionOfAerodynamicCenter = (str(input('Enter position of aerodynamic center B/O/A:'))).upper()
	aerodynamicForces           = (str(input('Enter type of aerodynamic forces   (theory)T/CFD:'))).upper()
	ax                          = fig.add_subplot(111, aspect ='equal')
	while t+dt <=1700:
		t               = t+dt
		plt.cla()
		ax.set_xlim(-15,15)
		ax.set_ylim(-15,15)
##		print(t)
##		wingState = wingState + dt*evolution(wingState,1)
##		wingState[0],wingState[1] = wrapTo2pi(wingState[0]),wrapTo2pi(wingState[1]
		wingState       = wingState + dt*rk2Step(wingState,t,dt,evolution,trueWindAngle,positionOfAerodynamicCenter)
		wingState[0][0] = wrapTo2pi(wingState[0][0])
		wingState[1][0] = wrapTo2pi(wingState[1][0])
##		print('main wing angle: ', wingState[0][0])
##		print('true wind angle: ', trueWindAngle)
##		angleLift       = angleOfLiftForceMW(wingState,trueWindAngle)
##		lift            = liftForceOnMW(wingState,trueWindAngle)
##		print(angleLift,lift)
##		drawArrow(0, 0, angleLift, lift,'b')
		drawWingSailIRT(wingState[0][0],wingState[2][0],trueWindAngle,trueWindSpeed)
		drawHull(0)
		plt.pause(0.0001)
	plt.show()
	"""



	# ===== Equilibrium position in function of the position of the aerodynamic center (Theory forces) ===== #
	
	"""
	plt.figure()
	timesRK2,statesRK2O = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'O')
	drawWingSailAngle(timesRK2,statesRK2O[0,:],'on rot axis','RK2')

	timesRK2,statesRK2B = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'B')
	drawWingSailAngle(timesRK2,statesRK2B[0,:],'before rot axis','RK2')

	timesRK2,statesRK2A = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle)
	drawWingSailAngle(timesRK2,statesRK2A[0,:],'after rot axis','RK2 (theoretical forces)')
	plt.legend()
	plt.savefig('Simulation_pics/Comparison position aerodynamic center for theoretical aerodynamic forces.png')
	"""
	

	# ==== Equilibrium position in function of the position of the aerodynamic center (Experimental forces) ===== #
	
	
	plt.figure()
	timesRK2,statesRK2O = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'O','CFD')
	drawWingSailAngle(timesRK2,statesRK2O[0,:],'on rot axis','RK2')

	timesRK2,statesRK2B = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'B','CFD')
	drawWingSailAngle(timesRK2,statesRK2B[0,:],'before rot axis','RK2')
	
	timesRK2,statesRK2A = rk2Scheme(wingState,0,200, 0.01, evolution, trueWindAngle,'A','CFD')
	drawWingSailAngle(timesRK2,statesRK2A[0,:],'after rot axis','RK2 (experimental forces)')
	plt.legend()
	plt.savefig('Simulation_pics/Comparison position aerodynamic center for experimental aerodynamic forces2.png')
	
	
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
	evolutionAnglesFDSA,times = evolutionMWAngleForDifferentStartAngles(wingState,200,trueWindAngle,'A', 'T') 
	drawEvolutionMWAngle(evolutionAnglesFDSA,times)
	plt.savefig('Simulation_pics/evolution of main wing angles for different start angles theoretical forces_aerodynamic center behind the mast')
	"""

	"""
	evolutionAnglesFDSA,times = evolutionMWAngleForDifferentStartAngles(wingState,200,trueWindAngle,'A', 'CFD') 
	drawEvolutionMWAngle(evolutionAnglesFDSA,times)
	plt.savefig('Simulation_pics/evolution of main wing angles for different start angles experimental forces_aerodynamic center behind the mast')
	"""
	# ===== Equilibrium position for different Start angles ===== #

	# Warning!! run the previous section before in order to put a reasonnable stop point
	
	"""
	equilibriumAnglesFDSA  = equilibriumAngleForDifferentStartAngles(wingState,200,trueWindAngle,'A','T')
	drawEquilibriumAngles( equilibriumAnglesFDSA)
	plt.savefig('Simulation_pics/equilibrium angles for different start angles theoretical forces_aerodynamic center behind the mast.png')
	"""
	
	"""
	equilibriumAnglesFDSA  = equilibriumAngleForDifferentStartAngles(wingState,200,trueWindAngle,'A','CFD')
	drawEquilibriumAngles( equilibriumAnglesFDSA)
	plt.savefig('Simulation_pics/equilibrium angles for different start angles experimental forces_aerodynamic center behind the mast.png')
	"""

	
	# ===== Equilibrium position for different servo wing angles (start at 0)===== #
	"""
	plt.figure()
	equilibriumAnglesFDTA = equilibriumAngleForDifferentTailAngles(wingState,200,trueWindAngle,'A','T')
	drawEquilibriumAngles( equilibriumAnglesFDTA,'FDTA')
	plt.savefig('Simulation_pics/equilibrium angles for different tail angles theoretical forces_aerodynamic center behind the mast.png')
	"""
	
