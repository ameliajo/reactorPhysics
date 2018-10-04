import matplotlib.pyplot as plt

def plotRay(r,color):
    x = [r.x0,r.xMax]
    y = [r.y0,r.yMax]
    plt.plot(r.x0,r.y0,color,marker='o')
    plt.plot(x,y,color)

def plotRaySegment(r,color,firstIntersection):
    x = [r.x0,firstIntersection["x-int"]]
    y = [r.y0,firstIntersection["y-int"]]
    plt.plot(x,y,color)



