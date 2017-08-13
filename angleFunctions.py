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

def listWrapTo2pi(l):
    u =[]
    for i in l:
        u.append(wrapTo2pi(i))
    return(l)

def listRadTodeg(l):
    u=[]
    for i in l:
        u.append(radTodeg(i))
    return(u)

def listDegTorad(l):
    u = []
    for i in l:
        u.append(degTorad(i))
    return(l)


def apparentWindCoordinatesToMWCoordinates(FxW,FyW,alpha):
    Fy = FyW*np.cos(alpha)-FxW*np.sin(alpha)
    Fx = FyW*np.sin(alpha)+FxW*np.cos(alpha)
    return (Fx,Fy)

def MWCoordinatesToBoatCoordinates(FxMW,FyMW,MWAngle):
    Fx = FxMW*np.cos(MWAngle)+FyMW*np.sin(MWAngle)
    Fy = FyMW*np.cos(MWAngle)-FxMW*np.sin(MWAngle)
    return (Fx,Fy)
"""
print(degTorad(5))
print(degTorad(10))
print(degTorad(15))
print(degTorad(20))
print(degTorad(25))
print(degTorad(30))
"""
