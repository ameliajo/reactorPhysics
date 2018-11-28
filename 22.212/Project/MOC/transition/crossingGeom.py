
def crossCircle(r, u, circle):
    c = circle
    x = r[0]
    y = r[1]
            
    A = u[1]**2 + u[0]**2
    B = 2*(x*u[0]- c.x0*u[0]+ y*u[1]- c.y0*u[1])
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
        t = (xPlane.x0-r[0])/u[0]
        return [{"x":r[0]+t*u[0],"y":r[1]+t*u[1],"t":t,"surface":xPlane}] if t > 0 else []


def crossYPlane(r,u,yPlane):
        t = (yPlane.y0-r[1])/u[1]
        return [{"x":r[0]+t*u[0],"y":r[1]+t*u[1],"t":t,"surface":yPlane}] if t > 0 else []

