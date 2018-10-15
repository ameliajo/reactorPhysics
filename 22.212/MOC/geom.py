import matplotlib.pyplot as plt



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
    def __init__(self,x,y,Length,radii,modMat,fuelMat,circleColor,color,idNum,L,R,D,U):
        self.id = idNum 
        self.sideLength = Length
        self.L = L 
        self.R = R 
        self.U = U
        self.D = D
        self.color = color
        circles = []
        for i,r in enumerate(radii):
            circles.append(circle(x,y,r,fuelMat,circleColor[i],i+1))
        self.mat = modMat; self.C = circles;


