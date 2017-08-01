# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 11:00:56 2017

@author: jazicoan
"""
import numpy as np
import matplotlib.pyplot as plt
from angleFunctions import wrapTo2pi,degTorad,radTodeg, listRadTodeg, listDegTorad, listWrapTo2pi



"""
UNCOMMENT MAIN BLOC AT THE END OF THE FILE FOR TESTS
"""


def drawArrow(x,y,θ,L,color):
	e=0.2
	M1=L*np.array([[0,1,1-e,1,1-e],[0,0,-e,0,e]])
	M=np.append(M1,[[1,1,1,1,1]],axis=0)
	R=np.array([[np.cos(θ),-np.sin(θ),x],[np.sin(θ),np.cos(θ),y],[0,0,1]])
	M=R@M
	plt.plot(M[0,:],M[1,:],color)

def drawWingSailAngle(times,states,label,schemeName):
	plt.plot(times,states,label = label)
	plt.xlabel("time")
	plt.ylabel("angle")
	plt.title("angle in function of time using "+str(schemeName)+" scheme")
	return()

def drawEquilibriumAngles(equilibriumAngles,ident='FDSA'):
	title = input("Enter title of Figure :")
	if type(title) != str:
		title =str(title)
	plt.figure()
	plt.axis([-2,2,-2,2])
	if ident == 'FDTA':
		for i in range(len(equilibriumAngles)):
	#		print(i)
			plt.text(np.cos(equilibriumAngles[i]),np.sin(equilibriumAngles[i]),str(i-26))
	else:
		for i in range(len(equilibriumAngles)):
	#		print(i)
			plt.text(np.cos(equilibriumAngles[i]),np.sin(equilibriumAngles[i]),str(i-30))
#	plt.plot(np.cos(equilibriumAngles),np.sin(equilibriumAngles),'bo')
	plt.plot([-2,2],[0,0],'k')
	plt.plot([0,0],[-2,2],'k')
	x,y = [],[]
	circle = np.linspace(0,7,360)
	for i in circle:
		x.append(np.cos(i))
		y.append(np.sin(i))
	plt.plot(x,y,'k')
	plt.xlabel('x')
	plt.ylabel('y')
	plt.axis('equal')
	plt.title(title.upper())
	return()

def drawEvolutionMWAngle(evolutionAngles,times,wide,ident = 'FDSA'):
	title = input("Enter title of Figure :")
	if type(title) != str:
		title =str(title) 
	plt.figure()
	if ident == 'FDTA':
		for i in range(len(evolutionAngles)):
			
			plt.plot(times,listRadTodeg(evolutionAngles[i][0,:]),label='TA= '+str(i+wide[0]))
	else:
		for i in range(len(evolutionAngles)):
			plt.plot(times,listRadTodeg(evolutionAngles[i][0,:]),label='SA= '+str(i+wide[0]))
	plt.legend(ncol = 4 ,fontsize = 4, )
	plt.axis([0,times[-1]+10,-50,50])
	plt.title(title.upper())
	return()


def drawHull(heading):

	hull = np.array([[ 6,  3, -1, -7, -7, -1,  3,  6],
	                 [ 0,  2,  3,  1, -1, -3, -2,  0],
	                 [ 1,  1,  1,  1,  1,  1,  1,  1]])

	rotationMatrix = np.array([[np.cos(heading), -np.sin(heading), 0],
                               [np.sin(heading),  np.cos(heading), 0],
							   [              0,                0, 1]])

	hull = rotationMatrix.dot(hull)

	plt.plot(hull[0,:],hull[1,:],'k')
	return ()



def drawWingSailIRT(theta, tailAngle, trueWindAngle, trueWindSpeed):
	tailAngle = wrapTo2pi(tailAngle)
	
	mainWing               = np.array([[ 2, 0, -3,  0,  2],
	                                   [ 0, 1,  0, -1,  0],
	                                   [ 1, 1,  1,  1,  1]])

	servoWing              = np.array([[ 0.5,   0, -1.5,    0, 0.5],
	                                   [ 0  , 0.3,    0, -0.3,   0],
	                                   [ 1  ,   1,    1,    1,   1]])


	
	rotationMatrixMainWing = np.array([[ np.cos(theta), -np.sin(theta),0],
	                                   [ np.sin(theta),  np.cos(theta),0],
	                                   [             0,              0,1]])
#	print("canardAngle: ",canardAngle)
#	print("tailAngle: ", tailAngle)
#	canardDeviation        = wrapTo2pi(theta+canardAngle)
	tailDeviation          = wrapTo2pi(theta+tailAngle)
	rotationMatrixTail     = np.array([[ np.cos(tailDeviation), -np.sin(tailDeviation),   -5.5*np.cos(theta)],
	                                   [ np.sin(tailDeviation),  np.cos(tailDeviation),   -5.5*np.sin(theta)],
	                                   [                     0,                      0,                    1]])

#	rotationMatrixCanard   = np.array([[ np.cos(canardDeviation), -np.sin(canardDeviation),    4.8*np.cos(theta)],
#	                                   [ np.sin(canardDeviation),  np.cos(canardDeviation),    4.8*np.sin(theta)],
#	                                   [                       0,                        0,                    1]])

	mainWing               = rotationMatrixMainWing.dot(mainWing)
	tail                   = rotationMatrixTail.dot(servoWing)
#	canard                 = rotationMatrixCanard.dot(servoWing)
	plt.plot(mainWing[0,:],mainWing[1,:],'r')
	plt.plot(tail[0,:],tail[1,:],'g')
#	plt.plot(canard[0,:],canard[1,:],'b')
	plt.plot(0,0,'r+')
#	print('theta: ',theta)
#	print('canardDeviation: ',canardDeviation)
#	print('tailDeviation: ',tailDeviation)
	plt.text(-11,13,'Wind:')
	drawArrow(-9.0,10.0,trueWindAngle,trueWindSpeed,'k')
#	drawArrow(4.8*np.cos(theta),4.8*np.sin(theta),canardDeviation,1,'b')
	drawArrow(-5.5*np.cos(theta),-5.5*np.sin(theta),tailDeviation,1,'g')



#if __name__=='__main__':

####################################################################################################
#                               drawMainWingsailAngle Test Block                                   #
####################################################################################################
#	times = [0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08]
#	states = np.sin(times)
#	drawMainWingsailAngle(times,states,'rk2')
#	plt.show()



####################################################################################################
#                               drawEquilibriumAngles Test Block                                   #
####################################################################################################

#FDSWA mode

#	equilibriumAnglesFDSWA = [-1.46,-1.35,-1.28,-1.15,-0.99, -0.84,-0.72,-0.68,-0.51,-0.42,-0.37,-0.27,-0.13,0,0.13,0.27,0.37,0.42,0.51,0.62,0.73,0.84,0.95,1.08,1.26,1.38,1.45]
#	drawEquilibriumAngles(equilibriumAnglesFDSWA,'FDSWA')
#	plt.show()



#FDSA mode

#	equilibriumAnglesFDSA = [0,0.001,0.45,0.005,0.002,0.41,0.456]
#	drawEquilibriumAngles(equilibriumAnglesFDSA)
#	plt.show()

####################################################################################################
#                                  drawWingSailIRT Test Block                                      #
####################################################################################################

#	theta         = wrapTo2pi(degTorad(30))
#	canardAngle   = wrapTo2pi(degTorad(13))
#	trueWindAngle = wrapTo2pi(degTorad(180))
#	trueWindSpeed = 3
#	fig = plt.figure()
#	ax  = fig.add_subplot(111, aspect ='equal')
#	ax.set_xlim(-10,10)
#	ax.set_ylim(-10,10)
#	drawWingSailIRT(theta, canardAngle, trueWindAngle, trueWindSpeed)
#	plt.show()
	
