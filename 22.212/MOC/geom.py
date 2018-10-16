import matplotlib.pyplot as plt
from math import pi



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
    def __init__(self,x,y,radii,modMat,fuelMat,circleColor,color,idNum,L,R,D,U):
        self.id = idNum 
        self.L = L 
        self.R = R 
        self.U = U
        self.D = D
        self.color = color
        circles = []
        for i,r in enumerate(radii):
            circles.append(circle(x,y,r,fuelMat,circleColor[i],i+1))
        self.mat = modMat; self.C = circles;

def getVolumes(numRegionsPerCell,cellCircles,sideLength):
    volumes = [0.0]*numRegionsPerCell
    if numRegionsPerCell == 6:
        volumes[0] = sideLength - pi*cellCircles[0].r**2
        volumes[1] = pi*cellCircles[0].r**2 - pi*cellCircles[1].r**2
        volumes[2] = pi*cellCircles[1].r**2 - pi*cellCircles[2].r**2
        volumes[3] = pi*cellCircles[2].r**2 - pi*cellCircles[3].r**2
        volumes[4] = pi*cellCircles[3].r**2 - pi*cellCircles[4].r**2
        volumes[5] = pi*cellCircles[4].r**2
    if numRegionsPerCell == 5:
        volumes[0] = sideLength - pi*cellCircles[0].r**2
        volumes[1] = pi*cellCircles[0].r**2 - pi*cellCircles[1].r**2
        volumes[2] = pi*cellCircles[1].r**2 - pi*cellCircles[2].r**2
        volumes[3] = pi*cellCircles[2].r**2 - pi*cellCircles[3].r**2
        volumes[4] = pi*cellCircles[3].r**2
    if numRegionsPerCell == 3:
        volumes[0] = sideLength - pi*cellCircles[0].r**2
        volumes[1] = pi*cellCircles[0].r**2 - pi*cellCircles[1].r**2
        volumes[2] = pi*cellCircles[1].r**2
    if numRegionsPerCell == 2:
        volumes[0] = sideLength - pi*cellCircles[0].r**2
        volumes[1] = pi*cellCircles[0].r**2 
    return volumes





