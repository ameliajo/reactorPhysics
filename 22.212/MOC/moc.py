import matplotlib.pyplot as plt
import random

class ray:
    def __init__(self,x,y,sin,cos,length):
        ray.x = x
        ray.y = y
        ray.sin = sin
        ray.cos = cos
        ray.l = length

fig, ax = plt.subplots() 
random.seed(3)
xL = -2.0
xR =  2.0
yD = -2.0
yU =  2.0
plt.plot([xL,xL,xR,xR,xL],[yU,yD,yD,yU,yU])
circle1 = plt.Circle((0, 0), 1.0)
ax.add_artist(circle1)


ray_begin_x = (xR-xL)*(random.random()-0.5)
ray_begin_y = (yU-yD)*(random.random()-0.5)
ray_cos = 2.0*random.random()-1.0
ray_sin = 2.0*random.random()-1.0
r1 = ray(ray_begin_x,ray_begin_y,ray_sin,ray_cos,10)

x = [r1.x,r1.x+r1.cos*r1.l]
y = [r1.y,r1.y+r1.sin*r1.l]

x0 = r1.x
y0 = r1.y
x1 = r1.x + r1.cos*r1.l
y1 = r1.y + r1.sin*r1.l

plt.plot(x,y,'r')

if ( x1 < xL ):
    print("goes left")
    r1.cos = - r1.cos
if ( x1 > xR ):
    print("goes right")
if ( y1 < yD ):
    print("goes down")
if ( y1 > yU ):
    print("goes up")

x = [xL,xL+r1.cos*r1.l]
y = [y1,y1+r1.sin*r1.l]
plt.plot(x,y,'r')







plt.xlim(-4,4)
plt.ylim(-4,4)
plt.show()



