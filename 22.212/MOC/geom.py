import matplotlib.pyplot as plt



class xPlane:
    def __init__(self,x,bc):
        self.x  = x; self.bc = bc

class yPlane:
    def __init__(self,y,bc):
        self.y  = y; self.bc = bc
        
class circle:
    def __init__(self,x,y,r,mat,color):
        self.x = x; self.y = y; self.r = r; self.mat = mat; self.color = color

class box:
    def __init__(self,x,y,L,radii,modMat,fuelMat,color):
        self.L = xPlane(x-L/2.0,'ref')
        self.R = xPlane(x+L/2.0,'ref')
        self.U = yPlane(y+L/2.0,'ref')
        self.D = yPlane(y-L/2.0,'ref')
        circles = []
        for i,r in enumerate(radii):
            circles.append(circle(x,y,r,fuelMat,color[i]))
        self.mat = modMat; self.C = circles;



colors = ["crimson","deeppink","fuchsia","violet","darkviolet","rebeccapurple","slateblue","slategrey","dodgerblue","skyblue","cadetblue","c","darkslategray","mediumturquoise","mediumaquamarine","mediumseagreen","g","lightgreen","lawngreen", "yellowgreen", "olive", "lemonchiffon", "darkgoldenrod", "orange", "navajowhite", "darkorange", "peachpuff", "coral", "mistyrose", "maroon", "lightcoral", "lightgray", "grey", "k"]
