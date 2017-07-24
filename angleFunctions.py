# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 11:02:38 2017

@author: jazicoan
"""
import numpy as np

def degTorad(angle):
	return((angle*np.pi/180))

	
def radTodeg(angle):
	return(angle*180/np.pi)

def wrapTo2pi(theta):
	theta = 2.0*np.arctan(np.tan(theta/2.0))
	return(theta)

"""
print(degTorad(5))
print(degTorad(10))
print(degTorad(15))
print(degTorad(20))
print(degTorad(25))
print(degTorad(30))
"""
