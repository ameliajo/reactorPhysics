
def crossCircle(ray, circle):
    c = circle
    x,   y = ray.r
    ux, uy = ray.u
            
    A = ux**2 + uy**2
    B = 2.0*(x*ux - c.x0*ux + y*uy - c.y0*uy)
    C = x**2 - 2.0*x*c.x0 + c.x0*c.x0 + y**2 - 2.0*y*c.y0 + c.y0**2 - c.r**2

    intersections = []
    for pm in [1.0,-1.0]:
        if B*B - 4*A*C < 0: continue
        t = (-B + pm*(B*B-4*A*C)**0.5)/(2*A)
        if (abs(t.imag) < 1.0e-20 and t.real > 0.0):
            t = t.real
            xInt, yInt = x + t*ux, y + t*uy
            intersections.append({"x":xInt,"y":yInt,"t":t,"surface":circle})

    return intersections

def crossXPlane(ray,xPlane):
    x,   y = ray.r
    ux, uy = ray.u
 
    t = (xPlane.x0-x)/ux
    return [{"x":x+t*ux,"y":y+t*uy,"t":t,"surface":xPlane}] if t > 0 else []


def crossYPlane(ray,yPlane):
    x,   y = ray.r
    ux, uy = ray.u
 
    t = (yPlane.y0-y)/uy
    return [{"x":x+t*ux,"y":y+t*uy,"t":t,"surface":yPlane}] if t > 0 else []

