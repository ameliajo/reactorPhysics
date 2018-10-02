import matplotlib.pyplot as plt



class xPlane:
    def __init__(self,x,bc):
        self.x  = x; self.bc = bc

class yPlane:
    def __init__(self,y,bc):
        self.y  = y; self.bc = bc
        
class circle:
    def __init__(self,x,y,r,mat):
        self.x = x; self.y = y; self.r = r; self.mat = mat

class box:
    def __init__(self,L,R,U,D,mat,circle):
        self.L = L; self.R = R; self.U = U; self.D = D; 
        self.mat = mat; self.C = circle;
    def plot(self,ax):
        plt.plot([self.L.x,self.L.x,self.R.x,self.R.x,self.L.x],
                 [self.U.y,self.D.y,self.D.y,self.U.y,self.U.y])
        c = plt.Circle((self.C.x, self.C.y), self.C.r)
        ax.add_artist(c)


