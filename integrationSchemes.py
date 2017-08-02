# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 11:09:58 2017

@author: jazicoan
"""

"""
UNCOMMENT MAIN BLOC AT THE END OF THE FILE FOR TESTS
"""



import numpy as np
from angleFunctions import wrapTo2pi
import matplotlib.pyplot as plt

def eulerScheme(state, start, stop, stepSize, evolutionFunction,trueWindAngle=0, positionAerodynamicCenter='A', aerodynamicForces = 'T'):
	t      = start
	times  = [t]
	stateList = state
	while t+stepSize<=stop:
		state = state+stepSize*evolutionFunction(state,t,trueWindAngle,positionAerodynamicCenter, aerodynamicForces)
		state[0] = wrapTo2pi(state[0])
		state[1] = wrapTo2pi(state[1])
		t = t + stepSize
		times.append(t)
		stateList = np.concatenate((stateList,state),axis=1)
	return(times,stateList)


def rk2Step(state, t, stepSize, evolutionFunction, trueWindAngle=0,positionAerodynamicCenter='A', aerodynamicForces = 'T'):
	k1 = evolutionFunction(state,t,trueWindAngle,positionAerodynamicCenter, aerodynamicForces)
	k2 = evolutionFunction(state+stepSize*k1,t+stepSize/2,trueWindAngle,positionAerodynamicCenter, aerodynamicForces)
	return(k2)

def rk2Scheme(state, start, stop, stepSize, evolutionFunction,trueWindAngle=0,positionAerodynamicCenter='A',  aerodynamicForces = 'T'):
	t      = start
	times  = [t]
	stateList = state
	while t+stepSize<=stop:
#		state = state+stepSize*evolutionFunction(state+stepSize*evolutionFunction(state,t,trueWindAngle),t+stepSize/2,trueWindAngle)
		state = state+stepSize*rk2Step(state, t, stepSize, evolutionFunction, trueWindAngle,positionAerodynamicCenter, aerodynamicForces)
		state[0] = wrapTo2pi(state[0])
		state[1] = wrapTo2pi(state[1])
		t = t + stepSize
		times.append(t)
		stateList = np.concatenate((stateList,state),axis=1)
	return(times,stateList)



def rk4(y0, t0, t1, h, evolutionFunction, trueWindAngle = 0, positionAerodynamicCenter='A', aerodynamicForces = 'T'):
	vt = [t0] 
	t  = t0
	vy = y0
	y  = y0
	while t <= t1:
		k1 = h * evolutionFunction(y,t,trueWindAngle,positionAerodynamicCenter,aerodynamicForces)
		k2 = h * evolutionFunction(y + 0.5 * k1, t + 0.5 * h, trueWindAngle,positionAerodynamicCenter,aerodynamicForces)
		k3 = h * evolutionFunction(y + 0.5 * k2, t + 0.5 * h, trueWindAngle,positionAerodynamicCenter,aerodynamicForces)
		k4 = h * evolutionFunction(y + k3, t + h, trueWindAngle,positionAerodynamicCenter,aerodynamicForces)
		t  = t+h
		vt.append(t)
		vy = np.concatenate( (vy,y + (k1 + k2 + k2 + k3 + k3 + k4) / 6),axis=1)
	return vt, vy


def evolution(X,t,trueWindAngle=0,positionAerodynamicCenter='A',aerodynamicForces='T'):
	k = 1
	g = 60
	l = 1
	Xdot = [0]*2
	Xdot[0] = X[1]
	Xdot[1] = -(g/l)*np.sin(X[0])-k*X[1]
	Xdot    = np.array([Xdot])
	Xdot    = Xdot.reshape((2,1))

	return Xdot




if __name__=='__main__':
	state = np.array([[0.3],[0]])
	times,stateListE   = rk2Scheme2(state,0,100,0.01,evolution)
	timesRK2,stateListRK2 = rk2Scheme(state,0,100,0.01,evolution)
	plt.figure()
	plt.plot(timesRK2,stateListRK2[0,:],'g', label='Rk2 step 1/2')
	plt.plot(times,stateListE[0,:],'b',label='Rk2 step 1')
	plt.xlabel('time')
	plt.ylabel('theta')
	plt.legend()
	plt.title('Evolution of theta in function of time')
	plt.figure()
	plt.plot(stateListRK2[0,:],stateListRK2[1,:],'g')
	plt.plot(stateListE[0,:],stateListE[1,:],'b')
	plt.xlabel('position')
	plt.ylabel('speed')
	plt.title('phasis diagram')
	plt.show()
