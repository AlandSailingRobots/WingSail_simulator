# -*- coding: utf-8 -*-
"""
Created on Tue Jun 13 14:24:28 2017

@author: jazicoan
"""



"""
UNCOMMENT MAIN BLOC AT THE END OF THE FILE FOR TESTS
"""





# test of simulation of a wing sail
import numpy as np
import matplotlib.pyplot as plt


X = 0

"""Constants"""
apparentWindAngle           = np.pi/6
apparentWindSpeed           = 2
servoWingArea               = X #to find
MainWingArea                = X #to find
rho                         = 1
distanceCanard              = X #to find
distanceTail                = X #to find
maxCanardAngle              = 13
epsilon                     = X #to find
constPartWindForceServoWing = (1/2)*rho*apparentWindSpeed**2*servoWingArea
constPartWindForceMainWing  = (1/2)*rho*apparentWindSpeed**2*MainWingArea
#Should Check if useful or not
#servoWingLiftCoefficient    = X #to find
#servoWingDragCoefficient    = X #to find

#if lucky need to complete it with aerodynamics measurements
CZs                         = np.ones((1,360))
CXs                         = np.ones((1,360))

def degTorad(angle):
	return(angle*np.pi/180)
	
def radTodeg(angle):
	return(angle*180/np.pi)


"""Parameters for experiment"""
canardAngle                 = degTorad(5)
tailAngle                   = -2*canardAngle


"""Parameters that will change"""
mainWingAngle               = X
dotMainWingAngle            = X


####################################################################################################
#                          Aerodynamic Coefficients Search Methods                                 #
####################################################################################################


def findCX(alphaInDegree):
	upper = int(radTodeg(alphaInDegree)+1)
	lower = int(alphaInDegree)
	if abs(upper-alphaInDegree)>abs(lower-alphaInDegree):
		approxAngleOfAttack = lower
	else:
		approxAngleOfAttack = upper
	global CXs
	CX = CXs[0][approxAngleOfAttack]
	return(CX)

def findCZ(alphaInDegree):
	upper = int(radTodeg(alphaInDegree)+1)
	lower = int(alphaInDegree)
	if abs(upper-alphaInDegree)>abs(lower-alphaInDegree):
		approxAngleOfAttack = lower
	else:
		approxAngleOfAttack = upper
	global CZs
	CZ = CZs[0][approxAngleOfAttack]
	return(CZ)

####################################################################################################
#                                       Main wing parameters                                       #
####################################################################################################

def angleOfAttackMainWing():
	global mainWingAngle
	global apparentWindAngle
	return(mainWingAngle+apparentWindAngle)

def liftCoefficientMainWing():
	return(findCZ(radTodeg(angleOfAttackMainWing())))

def liftForceOnMainwing():
	global constPartWindForceMainWing
	return(constPartWindForceMainWing*liftCoefficientMainWing())

def dragCoefficientMainWing():
	return(findCX(radTodeg(angleOfAttackMainWing())))

def dragForceOnMainWing():
	global constPartWindForceMainWing
	return(constPartWindForceMainWing*dragCoefficientMainWing())



####################################################################################################
#                                          Tail parameters                                         #
####################################################################################################
def angleOfAttackTail():
	global apparentWindAngle
	global tailAngle
	return(angleOfAttackMainWing()-np.sign(np.cos(apparentWindAngle))*tailAngle)

def dragCoefficientTail():
	return(findCX(radTodeg(angleOfAttackTail())))

def liftCoefficientTail():
	return(findCZ(radTodeg(angleOfAttackTail())))

def orientationTail():
	global tailAngle
	global epsilon
	return(angleOfAttackMainWing()+tailAngle-epsilon/2)

def liftForceOnTail():
	global constPartWindForceServoWing
	return(constPartWindForceServoWing*liftCoefficientTail()*orientationTail())

def momentOfLiftForceOnTail():
	global distanceTail
	return(distanceTail*liftForceOnTail())

#########################################Need Check#################################################
def dragForceOnTail():
	global constPartWindForceServoWing
	return(constPartWindForceServoWing*dragCoefficientTail()*np.sin(tailAngle))

def momentOfDragForceOnTail():
	global distanceTail
	return(distanceTail*dragForceOnTail())
####################################################################################################




####################################################################################################
#                                         Canard parameters                                        #
####################################################################################################

def angleOfAttackCanard():
	global apparentWindAngle
	global canardAngle
	return(angleOfAttackMainWing()+np.sign(np.cos(apparentWindAngle))*canardAngle)


#Should be checked
def dragCoefficientCanard():
	return(findCX(angleOfAttackCanard()))

#Should be checked
def liftCoefficientCanard():
	return(findCZ(angleOfAttackCanard()))


def orientationCanard():
	global canardAngle
	return(angleOfAttackMainWing()+canardAngle)


def liftForceOnCanard():
	global constPartWindForceServoWing
	return(constPartWindForceServoWing*liftCoefficientCanard()*orientationCanard())

def momentOfLiftForceOnCanard():
	global distanceCanard
	return(distanceCanard*liftForceOnCanard())


#########################################Need Check#################################################
def dragForceOnCanard():
	global constPartWindForceServoWing
	return(constPartWindForceServoWing*dragCoefficientCanard()*np.sin(canardAngle))

def momentOfDragForceOnCanard():
	global distanceCanard
	return(distanceCanard*dragForceOnCanard())
####################################################################################################



startingState=np.array([[degTorad(mainWingAngle)],[degTorad(dotMainWingAngle)]])


###########################State equation of the angle of the wingsail""############################
def evolution(X,t):
	Xdot    = [0]*2
	Xdot[0] = X[1]
	Xdot[1] = -momentOfLiftForceOnTail() - momentOfLiftForceOnCanard()
	return(np.array([[Xdot[0]],[Xdot[1]]]))


####################################################################################################
#                                         Integration Schemes                                      #
####################################################################################################
def eulerScheme(state,start,stop,stepSize):
	t      = start
	times  = [t]
	states = state
	while t+stepSize<=stop:
		state = state+stepSize*evolution(state,t)
		times.append(t)
		np.hstack((states,state))
		t = t + stepSize
	return(times,states)

def rk2Scheme(state,start,stop,stepSize):
	t      = start
	times  = [t]
	states = state
	while t+stepSize<=stop:
		state = state+stepSize*evolution(state+stepSize*evolution(state,t),t+stepSize/2)
		times.append(t)
		np.hstack((states,state))
		t = t + stepSize
	return(times,states)

def rk4Scheme(state,start,stop,stepSize):
	t      = start
	times  = [t]
	states = state
	while t+stepSize<=stop:
		k1 = evolution(state,t)
		k2 = evolution(state + k1/2, t+stepSize/2)
		k3 = evolution(state + k2/2, t+stepSize/2)
		k4 = evolution(state + k3/2, t+stepSize)
		state =state + stepSize*(k1/6+k2/3+k3/3+k4/6)
		times.append(t)
		np.hstack((states,state))
		t = t + stepSize
	return(times,states)

####################################################################################################
#                                  Parameters Testing Methods                                      #
####################################################################################################

def equilibriumAngleForDifferentStartAngles():
	equilibriumAnglesFDSA = [] #FDSA = ForDifferentStartAngles
	for i in range(360):
		startingState = np.array([[degTorad(i)],[0]])
		startTime     = 0
		stopTime      = 10
		stepSize      = 0.1
		times, states = eulerScheme(startingState,startTime,stopTime,stepSize)
		states        = states[0][-3:]
		equilibriumAnglesFDSA.append(np.mean(states))
	return(equilibriumAnglesFDSA)

def equilibriumAngleForDiffrentServoWingAngles():
	equilibriumAnglesFDSWA = [] #FDSWA = ForDifferentServoWingAngles
	global canardAngle
	global tailAngle
	for i in range (-maxCanardAngle,maxCanardAngle+1):
		canardAngle  = i
		tailAngle    = -2*canardAngle
		startingState = np.array([[0],[0]])
		startTime     = 0
		stopTime      = 10
		stepSize      = 0.1
		times, states = eulerScheme(startingState,startTime,stopTime,stepSize)
		states        = states[0][-3:]
		equilibriumAnglesFDSWA.append(np.mean(states))
	return(equilibriumAnglesFDSWA)



####################################################################################################
#                                         Drawing Methods                                          #
####################################################################################################
def drawMainWingsailAngle(times,states):
	plt.figure()
	plt.plot(times,states)
	plt.xlabel("time")
	plt.ylabel("angle")
	plt.title("angle in funtcion of time")
	plt.show()
	return()

def drawEquilibriumAngles(equilibriumAngles):
	title = input("Enter title of Figure :")
	if type(title) != str:
		title =str(title)
	plt.figure()
	plt.axis([-1,1,-1,1])
	plt.plot(np.cos(equilibriumAngles),np.sin(equilibriumAngles),'bo')
	plt.plot(0,0,'r+')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.title(title.upper())
	plt.show()
	return()







if __name__ =='__main__':
#	times,states = eulerScheme(startingState,0,10,0.1)
#	equilibriumAnglesFDSWA = equilibriumAngleForDiffrentServoWingAngles()
	drawEquilibriumAngles( [1,2,3,4])
