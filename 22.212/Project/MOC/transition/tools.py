import numpy as np
import cmath
import re

def checktol(x,y,tol):
    """Check absolute difference between two values and compare to 
    a defined tolerance. Return boolean. """
    err = abs(abs(x-y)/x)
    if err>tol:
        return False
    elif err<tol:
        return True

def normalize(x):
    return x/np.linalg.norm(x)

def crossCircle(r, u, circle):
    c = circle
    sin = u[1]
    cos = u[0]
    x = r[0]
    y = r[1]
            
    A = sin**2 + cos**2
    B = 2*(x*cos - c.x0*cos + y*sin - c.y0*sin)
    C = x**2 - 2*x*c.x0 + c.x0*c.x0 + y**2 - 2*y*c.y0 + c.y0**2 - c.r**2

    intersections = []
    for pm in [1.0,-1.0]:
        if B*B - 4*A*C < 0: continue
        t = (-B + pm*(B*B-4*A*C)**0.5)/(2*A)
        if (abs(t.imag) < 1.0e-20 and t.real > 0.0):
            t = t.real
            xInt, yInt = r[0] + t*u[0], r[1] + t*u[1]
            intersections.append({"x":xInt,"y":yInt,"t":t,"surface":circle})

    return intersections

def crossXPlane(r,u,xPlane):
        cos = u[0]
        sin = u[1]
        t = (xPlane.x0-r[0])/cos
        return [{"x":r[0]+t*cos,"y":r[1]+t*sin,"t":t,"surface":xPlane}] if t > 0 else []


def crossYPlane(r,u,yPlane):
        cos = u[0]
        sin = u[1]
        t = (yPlane.y0-r[1])/sin
        return [{"x":r[0]+t*cos,"y":r[1]+t*sin,"t":t,"surface":yPlane}] if t > 0 else []


def get_trailing_numbers(s, zero=False):
    m = re.search(r'\d+$', s)
    if m:
        return int(m.group())
    elif zero:
        return 0
    else:
        return None

