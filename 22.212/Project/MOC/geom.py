import matplotlib.pyplot as plt
from math import pi
from colors import *



class xPlane:
    def __init__(self,x,bc):
        self.x  = x; self.bc = bc

class yPlane:
    def __init__(self,y,bc):
        self.y  = y; self.bc = bc
        
class circle:
    def __init__(self,x,y,r,mat,color,idVal):
        self.x = x; self.y = y; self.r = r; self.mat = mat; self.color = color
        self.id = idVal

class box:
    def __init__(self,xL,yL,radii,modMat,fuelMat,idNum,xi,yi,xPlanes,yPlanes):
        cellColors, circleColors = findColors(radii)
        self.id = idNum 
        self.L = xPlanes[xi]
        self.R = xPlanes[xi+1]
        self.D = yPlanes[yi]
        self.U = yPlanes[yi+1]
        self.color = [cellColors[idNum]]
        circles = []
        for i,r in enumerate(radii):
            circles.append(circle(xL,yL,r,fuelMat,circleColors[i],i+1))
        self.mat = modMat; self.C = circles;

def getVolumes(numRegionsPerCell,cellCircles,sideLength):
    volumes = [0.0]*numRegionsPerCell
    
    volumes[0] = sideLength - pi*cellCircles[0].r**2
    for i in range(1,numRegionsPerCell-1):
        volumes[i] = pi*cellCircles[i-1].r**2 - pi*cellCircles[i].r**2
    volumes[numRegionsPerCell-1] = pi*cellCircles[numRegionsPerCell-2].r**2

    return volumes





