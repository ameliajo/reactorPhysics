import matplotlib.pyplot as plt
import random
from geom import *

random.seed(3)
fig, ax = plt.subplots() 

class ray:
    def __init__(self,x,y,sin,cos,length):
        self.x0 = x
        self.y0 = y
        self.sin = sin
        self.cos = cos
        self.l = length

def plotIntersection(r,t):
    x_int = r.x0 + t * r.cos
    y_int = r.y0 + t * r.sin
    plt.plot(x_int,y_int,'bo')



circle1 = circle(0,0,0.5,'U')
box1 = box(xPlane(-1,'vac'),xPlane(1,'vac'),yPlane(1,'vac'),yPlane(-1,'vac'),'H',circle1)

x0 = (box1.R.x-box1.L.x)*(random.random()-0.5)
y0 = (box1.U.y-box1.D.y)*(random.random()-0.5)
cos = 2.0*random.random()-1.0
sin = 2.0*random.random()-1.0
r1 = ray(x0,y0,sin,cos,10) # should be 4 length


x = [r1.x0,r1.x0+r1.cos*r1.l]
y = [r1.y0,r1.y0+r1.sin*r1.l]

plt.plot(r1.x0,r1.y0,'ro')
plt.plot(x,y,'r')




# Does my ray intersect with the left boundary? 
L = box1.L
t = (L.x - r1.x0)/r1.cos
# If t > 0 and real, then it does intersect
if ( t > 0 ):
    plotIntersection(r1,t)

# Does my ray intersect with the top boundary? 
U = box1.U
t = (U.y - r1.y0)/r1.sin
if ( t > 0 ):
    plotIntersection(r1,t)

# Does my ray intersect with the right boundary? 
R = box1.R
t = (R.x - r1.x0)/r1.cos
print(t)
# If t > 0 and real, then it does intersect
if ( t > 0 ):
    plotIntersection(r1,t)

# Does my ray intersect with the top boundary? 
D = box1.D
t = (D.y - r1.y0)/r1.sin
print(t)
if ( t > 0 ):
    plotIntersection(r1,t)



#x0 = r1.x
#y0 = r1.y
#x1 = r1.x + r1.cos*r1.l
#y1 = r1.y + r1.sin*r1.l

#plt.plot(x,y,'r')

#if ( x1 < box1.L.x ):
#    print("goes left")
#    r1.cos = - r1.cos
#if ( x1 > box1.R.x ):
#    print("goes right")
#if ( y1 < box1.D.y ):
#    print("goes down")
#if ( y1 > box1.U.y ):
#    print("goes up")

#x = [box1.L.x,box1.L.x+r1.cos*r1.l]
#y = [box1.U.y,box1.U.y+r1.sin*r1.l]
#plt.plot(x,y,'r')





box1.plot(ax)





plt.xlim(-4,4)
plt.ylim(-4,4)
plt.show()



